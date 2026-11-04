# Author: Jakub Zapletal
# Date: 2026-11-04
# Inspects and benchmarks the national DTM files, derives their valid-data footprint, estimates runtimes
# and disk needs, and optionally converts them to a tiled layout with overviews (whole file or corridor).

import json
import os
import re
import shutil
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from osgeo import gdal, gdal_array

import gis

# --- Run variables (no command-line parameters; DEFENSES_<NAME> environment variables override them) ---
CONFIG_PATH = "config/defenses.toml"
COUNTRY = "XX"
FORCE = False  # recompute the footprint and redo a finished conversion
RUN_BENCHMARK = True

KEYS = [
    "general.work_dir",
    "general.output_dir",
    "general.n_workers",
    "general.gdal_cachemax_mb",
    "general.warp_memory_limit_mb",
    "general.log_level",
    "benchmark.n_windows",
    "benchmark.window_px",
    "benchmark.sequential_rows",
    "benchmark.seed",
    "country.dtm",
    "country.working_resolution_m",
    "country.unit_core_size_m",
    "country.unit_halo_m",
    "country.convert_dtm",
    "country.dtm_conversion",
    "country.dtm_compress",
    "country.footprint_resolution_m",
]
OVERVIEW_FACTORS = [2, 4, 8, 16, 32]

if __name__ == "__main__":
    run = gis.Run.env_overrides(
        CONFIG_PATH=CONFIG_PATH,
        COUNTRY=COUNTRY,
        FORCE=FORCE,
        RUN_BENCHMARK=RUN_BENCHMARK,
    )
    cfg = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], KEYS)
    general, country, bench = cfg["general"], cfg["country"], cfg["benchmark"]
    work_dir = Path(general["work_dir"]) / run["COUNTRY"]
    out_dir = Path(general["output_dir"]) / run["COUNTRY"]
    (work_dir / "dtm").mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    log = gis.Run.logger("prepare_dtm", out_dir / "logs", general["log_level"])
    gis.Config.warn_placeholders(cfg)
    gdal.SetCacheMax(general["gdal_cachemax_mb"] * 2**20)
    resolution = country["working_resolution_m"]
    compress = gis.Raster.best_compression(country["dtm_compress"])
    if compress != country["dtm_compress"]:
        log.warning(
            "Compression %s is not available in this GDAL build, using %s.",
            country["dtm_compress"],
            compress,
        )

    # --- 1. Inspection ---
    infos = [gis.Raster.inspect(path) for path in country["dtm"]]
    gis.Raster.build_vrt(
        country["dtm"], work_dir / "dtm" / "dtm_source.vrt"
    )  # stops on CRS mismatch
    for info in infos:
        log.info(
            "DTM %s: %d x %d px, pixel %.3f x %.3f, %s, nodata %s, compression %s, block %d x %d%s, "
            "overviews %d, file %.2f GB (raw %.2f GB)",
            info["path"],
            info["width"],
            info["height"],
            info["pixel_x"],
            abs(info["pixel_y"]),
            info["data_type"],
            info["nodata"],
            info["compression"],
            info["block_x"],
            info["block_y"],
            " (strip layout)" if info["strip_layout"] else "",
            info["overviews"],
            info["file_size_gb"],
            info["raw_size_gb"],
        )
        if info["strip_layout"] or info["overviews"] == 0:
            log.warning(
                "%s: %s; window reads at working resolution will be slow (see the benchmark and dtm_conversion).",
                Path(info["path"]).name,
                "strip layout" if info["strip_layout"] else "no overviews",
            )
        if info["nodata"] is None:
            log.warning(
                "%s has no nodata value; every pixel is treated as valid.",
                Path(info["path"]).name,
            )
        if resolution < abs(info["pixel_x"]) - 1e-9:
            raise SystemExit(
                f"working_resolution_m {resolution} is finer than the DTM pixel {info['pixel_x']}."
            )
    if len({round(abs(i["pixel_x"]), 6) for i in infos}) > 1:
        log.warning("DTM files differ in pixel size; the VRT uses the finest one.")
    pd.DataFrame(infos).drop(columns=["crs_wkt"]).to_csv(
        out_dir / "dtm_inspection.csv", index=False
    )

    # --- 2. Valid-data footprint ---
    footprint_path = out_dir / "dtm_footprint.gpkg"
    if footprint_path.exists() and not run["FORCE"]:
        footprints = gpd.read_file(footprint_path)
        log.info("Footprint reused from %s (set FORCE to recompute).", footprint_path)
    else:
        frames = []
        for info in infos:
            start = time.perf_counter()
            part = gis.Raster.footprint(info["path"], country["footprint_resolution_m"])
            part["dtm"] = Path(info["path"]).stem
            part["path"] = info["path"]
            frames.append(part)
            log.info(
                "Footprint of %s derived in %.1f s.",
                Path(info["path"]).name,
                time.perf_counter() - start,
            )
        footprints = gpd.GeoDataFrame(
            pd.concat(frames, ignore_index=True), crs=frames[0].crs
        )
        footprints["area_km2"] = footprints.area / 1e6
        gis.Vector.write_gpkg(footprints, footprint_path)
    for info in infos:
        extent_km2 = (
            info["width"]
            * info["height"]
            * abs(info["pixel_x"] * info["pixel_y"])
            / 1e6
        )
        area_km2 = footprints.loc[footprints["path"] == info["path"]].area.sum() / 1e6
        log.info(
            "Valid data of %s: %.1f km2 of %.1f km2 raster extent (%.1f %%).",
            Path(info["path"]).name,
            area_km2,
            extent_km2,
            100 * area_km2 / extent_km2,
        )

    # --- 3. Benchmark ---
    estimates = {"country": run["COUNTRY"], "dtm": []}
    bench_path = out_dir / "dtm_benchmark.csv"
    if run["RUN_BENCHMARK"]:
        rng = np.random.default_rng(bench["seed"])
        records = []
        for info in infos:
            ds = gdal.Open(info["path"])
            band = ds.GetRasterBand(1)
            gt = ds.GetGeoTransform()
            win_x, win_y = min(bench["window_px"], info["width"]), min(
                bench["window_px"], info["height"]
            )
            bytes_window = win_x * win_y * info["bytes_per_pixel"]
            polygon = shapely.union_all(
                footprints.loc[footprints["path"] == info["path"]].geometry.values
            )
            # Window centres drawn inside the valid footprint, so that the windows read real data.
            xs = rng.uniform(
                gt[0] + win_x / 2 * gt[1],
                gt[0] + (info["width"] - win_x / 2) * gt[1],
                200 * bench["n_windows"],
            )
            ys = rng.uniform(
                gt[3] + (info["height"] - win_y / 2) * gt[5],
                gt[3] + win_y / 2 * gt[5],
                200 * bench["n_windows"],
            )
            inside = shapely.contains_xy(polygon, xs, ys)
            centres = list(zip(xs[inside], ys[inside]))[: bench["n_windows"]]
            if len(centres) < bench["n_windows"]:
                log.warning(
                    "Only %d benchmark windows fit inside the footprint of %s.",
                    len(centres),
                    info["path"],
                )
            overview_messages = []
            sample_window = (0, 0)
            for k, (x, y) in enumerate(centres):
                col0 = int(
                    np.clip((x - gt[0]) / gt[1] - win_x / 2, 0, info["width"] - win_x)
                )
                row0 = int(
                    np.clip((y - gt[3]) / gt[5] - win_y / 2, 0, info["height"] - win_y)
                )
                sample_window = sample_window if k else (col0, row0)
                start = time.perf_counter()
                native = band.ReadAsArray(col0, row0, win_x, win_y)
                t_native = time.perf_counter() - start
                valid_share = (
                    float(np.mean(native != info["nodata"]))
                    if info["nodata"] is not None
                    else 1.0
                )
                # Same ground window at the working resolution, snapped outward to the working grid.
                xmin, ymax = gt[0] + col0 * gt[1], gt[3] + row0 * gt[5]
                bounds = (
                    xmin,
                    ymax - np.ceil(win_y * abs(gt[5]) / resolution) * resolution,
                    xmin + np.ceil(win_x * gt[1] / resolution) * resolution,
                    ymax,
                )
                if k == 0:
                    # Untimed single-threaded warp with debug output: the warp kernel reports the source
                    # window it reads; a window smaller than the native one proves that an overview was used.
                    gdal.SetConfigOption("CPL_DEBUG", "ON")
                    gdal.PushErrorHandler(
                        lambda err_class, err_no, msg: overview_messages.append(msg)
                    )
                    try:
                        gis.Raster.read_window(
                            info["path"],
                            bounds,
                            resolution,
                            "average",
                            general["warp_memory_limit_mb"],
                            multithread=False,
                        )
                    finally:
                        gdal.PopErrorHandler()
                        gdal.SetConfigOption("CPL_DEBUG", None)
                start = time.perf_counter()
                warped, _ = gis.Raster.read_window(
                    info["path"],
                    bounds,
                    resolution,
                    "average",
                    general["warp_memory_limit_mb"],
                )
                t_warp = time.perf_counter() - start
                records.append(
                    {
                        "dtm": Path(info["path"]).name,
                        "window": k,
                        "col0": col0,
                        "row0": row0,
                        "valid_share": valid_share,
                        "native_s": t_native,
                        "native_mb_s": bytes_window / 1e6 / t_native,
                        "warp_s": t_warp,
                        "warp_native_equivalent_mb_s": bytes_window / 1e6 / t_warp,
                        "warp_output_cells": warped.size,
                    }
                )
            del native, warped
            source_widths = [
                int(m.group(1))
                for m in map(
                    re.compile(r"Src=\d+,\d+,(\d+)x(\d+)").search, overview_messages
                )
                if m
            ]
            if not source_widths:
                log.warning(
                    "Could not determine from the GDAL debug output whether gdal.Warp used an "
                    "overview for %s.",
                    Path(info["path"]).name,
                )
            elif min(source_widths) < 0.9 * win_x:
                log.info(
                    "gdal.Warp read %s from an overview: source windows of %d px instead of %d "
                    "native px (factor %.0f).",
                    Path(info["path"]).name,
                    min(source_widths),
                    win_x,
                    win_x / min(source_widths),
                )
            else:
                log.info(
                    "gdal.Warp read %s at full resolution (source window %d px, %d overviews): "
                    "every native pixel is decoded.",
                    Path(info["path"]).name,
                    min(source_widths),
                    info["overviews"],
                )

            # Sequential full-width read (pattern of a conversion pass).
            n_rows = min(bench["sequential_rows"], info["height"])
            row0 = int(rng.integers(0, info["height"] - n_rows + 1))
            start = time.perf_counter()
            block = band.ReadAsArray(0, row0, info["width"], n_rows)
            t_seq = time.perf_counter() - start
            seq_mb_s = info["width"] * n_rows * info["bytes_per_pixel"] / 1e6 / t_seq
            del block

            # Compression ratio and speed of the target layout, measured on one native window.
            sample = band.ReadAsArray(sample_window[0], sample_window[1], win_x, win_y)
            is_float = np.issubdtype(sample.dtype, np.floating)
            mem = gdal.GetDriverByName("MEM").Create("", win_x, win_y, 1, band.DataType)
            mem.GetRasterBand(1).WriteArray(sample)
            if info["nodata"] is not None:
                mem.GetRasterBand(1).SetNoDataValue(info["nodata"])
            start = time.perf_counter()
            gdal.Translate(
                "/vsimem/compression_sample.tif",
                mem,
                format="GTiff",
                creationOptions=[
                    "TILED=YES",
                    "BLOCKXSIZE=512",
                    "BLOCKYSIZE=512",
                    f"COMPRESS={compress}",
                    f"PREDICTOR={3 if is_float else 2}",
                    "NUM_THREADS=ALL_CPUS",
                ],
            )
            t_compress = time.perf_counter() - start
            ratio = gdal.VSIStatL("/vsimem/compression_sample.tif").size / bytes_window
            gdal.Unlink("/vsimem/compression_sample.tif")

            table = pd.DataFrame(
                [r for r in records if r["dtm"] == Path(info["path"]).name]
            )
            summary = {
                "path": info["path"],
                "native_s_median": float(table["native_s"].median()),
                "native_s_p90": float(table["native_s"].quantile(0.9)),
                "native_mb_s_median": float(table["native_mb_s"].median()),
                "warp_s_median": float(table["warp_s"].median()),
                "warp_s_p90": float(table["warp_s"].quantile(0.9)),
                "warp_native_equivalent_mb_s_median": float(
                    table["warp_native_equivalent_mb_s"].median()
                ),
                "window_ground_km2": win_x * win_y * abs(gt[1] * gt[5]) / 1e6,
                "sequential_mb_s": seq_mb_s,
                "compression_ratio": ratio,
                "compression_mb_s": bytes_window / 1e6 / t_compress,
                "raw_gb": info["raw_size_gb"],
            }
            estimates["dtm"].append(summary)
            log.info(
                "Benchmark %s (%d windows of %d x %d px = %.2f km2): native read %.2f s/window (p90 %.2f), %.0f MB/s; "
                "warp to %.1f m %.2f s/window (p90 %.2f), %.0f MB/s native equivalent; sequential read %.0f MB/s; "
                "%s compression ratio %.2f at %.0f MB/s.",
                Path(info["path"]).name,
                len(table),
                win_x,
                win_y,
                summary["window_ground_km2"],
                summary["native_s_median"],
                summary["native_s_p90"],
                summary["native_mb_s_median"],
                resolution,
                summary["warp_s_median"],
                summary["warp_s_p90"],
                summary["warp_native_equivalent_mb_s_median"],
                seq_mb_s,
                compress,
                ratio,
                summary["compression_mb_s"],
            )
        pd.DataFrame(records).to_csv(bench_path, index=False)
        log.info(
            "Note: the operating system file cache can make repeated benchmarks faster than a cold read."
        )
    elif (out_dir / "dtm_estimates.json").exists():
        estimates = json.loads((out_dir / "dtm_estimates.json").read_text())

    # --- 4. Estimates from the measurements ---
    units_csv = out_dir / "units.csv"
    corridor_path = out_dir / "corridor.gpkg"
    if estimates["dtm"]:
        bench_km2 = np.mean([d["window_ground_km2"] for d in estimates["dtm"]])
        warp_s = np.mean([d["warp_s_median"] for d in estimates["dtm"]])
        if units_csv.exists():
            units = pd.read_csv(units_csv, dtype={"unit_id": str})
            window_km2 = units["est_window_cells"] * resolution**2 / 1e6
            read_total = float((window_km2 / bench_km2 * warp_s).sum())
            log.info(
                "Terrain read estimate: %d units, %.0f s per unit on average, %s in total on %d workers "
                "(source: median warp s/window x unit window area / benchmark window area; assumes reads scale with "
                "workers, which the disk may not allow).",
                len(units),
                read_total / len(units),
                gis.Run.format_seconds(read_total / general["n_workers"]),
                general["n_workers"],
            )
            estimates["terrain_read_s_total"] = read_total
        else:
            side_km = (country["unit_core_size_m"] + 2 * country["unit_halo_m"]) / 1000
            log.info(
                "Terrain read estimate for one unit window of %.1f x %.1f km: %.0f s (median warp s/window x area "
                "ratio). Run build_units.py for the per-country total.",
                side_km,
                side_km,
                side_km**2 / bench_km2 * warp_s,
            )
        corridor = (
            gpd.read_file(corridor_path).union_all().buffer(country["unit_halo_m"])
            if corridor_path.exists()
            else None
        )
        for summary in estimates["dtm"]:
            raw_gb = summary["raw_gb"]
            disk_a = (
                raw_gb * summary["compression_ratio"] * (1 + 1 / 3)
            )  # AVERAGE overviews 2..32 add about 1/3
            read_s = raw_gb * 1e3 / summary["sequential_mb_s"]
            compress_s = raw_gb * 1e3 / summary["compression_mb_s"]
            summary.update(
                {
                    "option_a_disk_gb": disk_a,
                    "option_a_read_s": read_s,
                    "option_a_compress_s": compress_s,
                }
            )
            log.info(
                "Option A for %s: about %.1f GB on disk (raw %.1f GB x measured ratio %.2f + 1/3 overviews); one "
                "sequential pass %s (read at %.0f MB/s) to %s (plus compression at %.0f MB/s); overview building "
                "not included.",
                Path(summary["path"]).name,
                disk_a,
                raw_gb,
                summary["compression_ratio"],
                gis.Run.format_seconds(read_s),
                summary["sequential_mb_s"],
                gis.Run.format_seconds(read_s + compress_s),
                summary["compression_mb_s"],
            )
            if corridor is not None:
                part = footprints.loc[footprints["path"] == summary["path"]].union_all()
                share = (
                    corridor.intersection(part).area / part.area if part.area else 0.0
                )
                summary.update(
                    {"corridor_share": share, "option_b_disk_gb": disk_a * share}
                )
                log.info(
                    "Option B for %s: corridor plus halo covers %.1f %% of the valid data; about %.1f GB on disk; the source is "
                    "still read once in full (%s), compression of the corridor part only (%s).",
                    Path(summary["path"]).name,
                    100 * share,
                    disk_a * share,
                    gis.Run.format_seconds(read_s),
                    gis.Run.format_seconds(compress_s * share),
                )
            else:
                log.info(
                    "Option B estimate needs corridor.gpkg; run build_units.py first."
                )
        (out_dir / "dtm_estimates.json").write_text(json.dumps(estimates, indent=2))

    # --- 5. Optional conversion ---
    targets = []
    if country["convert_dtm"]:
        option = country["dtm_conversion"]
        if option == "B" and not corridor_path.exists():
            raise SystemExit(
                "dtm_conversion = 'B' needs corridor.gpkg; run build_units.py first."
            )
        corridor = (
            gpd.read_file(corridor_path).union_all().buffer(country["unit_halo_m"])
            if option == "B"
            else None
        )
        for info in infos:
            target = work_dir / "dtm" / f"{Path(info['path']).stem}_{option}.tif"
            marker = target.with_suffix(".done")
            if marker.exists() and not run["FORCE"]:
                log.info(
                    "Conversion of %s already finished: %s",
                    Path(info["path"]).name,
                    target,
                )
                targets.append(str(target))
                continue
            need_gb = next(
                (
                    s.get("option_b_disk_gb" if option == "B" else "option_a_disk_gb")
                    for s in estimates["dtm"]
                    if s["path"] == info["path"]
                ),
                None,
            )
            free_gb = shutil.disk_usage(target.parent).free / 1e9
            if need_gb is not None and need_gb * 1.2 > free_gb:
                raise SystemExit(
                    f"Not enough disk space for {target}: about {need_gb:.1f} GB needed, "
                    f"{free_gb:.1f} GB free."
                )
            src = gdal.Open(info["path"])
            band = src.GetRasterBand(1)
            gt = src.GetGeoTransform()
            numpy_type = gdal_array.GDALTypeCodeToNumericTypeCode(band.DataType)
            is_float = np.issubdtype(numpy_type, np.floating)
            creation = [
                "TILED=YES",
                "BLOCKXSIZE=512",
                "BLOCKYSIZE=512",
                f"COMPRESS={compress}",
                f"PREDICTOR={3 if is_float else 2}",
                "BIGTIFF=YES",
                "NUM_THREADS=ALL_CPUS",
            ]
            start = time.perf_counter()
            last_report = [start]

            def progress(complete, message, data):
                now = time.perf_counter()
                if now - last_report[0] > 60:
                    last_report[0] = now
                    eta = (now - start) / max(complete, 1e-9) * (1 - complete)
                    log.info(
                        "  %s: %.1f %% done, ETA %s",
                        data,
                        100 * complete,
                        gis.Run.format_seconds(eta),
                    )
                return 1

            if option == "A":
                gdal.Translate(
                    str(target),
                    src,
                    format="GTiff",
                    creationOptions=creation,
                    callback=progress,
                    callback_data=target.name,
                )
            else:
                part = corridor.intersection(
                    footprints.loc[footprints["path"] == info["path"]].union_all()
                )
                if part.is_empty:
                    log.warning(
                        "The corridor does not touch %s; no copy written.", info["path"]
                    )
                    continue
                xmin, ymin, xmax, ymax = part.bounds
                col0 = max(int(np.floor((xmin - gt[0]) / gt[1])), 0)
                col1 = min(int(np.ceil((xmax - gt[0]) / gt[1])), src.RasterXSize)
                row0 = max(int(np.floor((ymax - gt[3]) / gt[5])), 0)
                row1 = min(int(np.ceil((ymin - gt[3]) / gt[5])), src.RasterYSize)
                width, height = col1 - col0, row1 - row0
                nodata = band.GetNoDataValue()
                if nodata is None:
                    nodata = (
                        float("nan") if is_float else float(np.iinfo(numpy_type).min)
                    )
                dst = gdal.GetDriverByName("GTiff").Create(
                    str(target), width, height, 1, band.DataType, creation
                )
                dst_gt = (
                    gt[0] + col0 * gt[1],
                    gt[1],
                    0.0,
                    gt[3] + row0 * gt[5],
                    0.0,
                    gt[5],
                )
                dst.SetGeoTransform(dst_gt)
                dst.SetProjection(src.GetProjection())
                dst.GetRasterBand(1).SetNoDataValue(nodata)
                bytes_per_row = width * gdal.GetDataTypeSize(band.DataType) // 8
                rows_per_pass = max(
                    512,
                    general["warp_memory_limit_mb"]
                    * 2**20
                    // bytes_per_row
                    // 512
                    * 512,
                )
                for r in range(0, height, rows_per_pass):
                    n_rows = min(rows_per_pass, height - r)
                    data = band.ReadAsArray(col0, row0 + r, width, n_rows)
                    band_gt = (
                        dst_gt[0],
                        dst_gt[1],
                        0.0,
                        dst_gt[3] + r * dst_gt[5],
                        0.0,
                        dst_gt[5],
                    )
                    band_box = (
                        band_gt[0],
                        band_gt[3] + n_rows * band_gt[5],
                        band_gt[0] + width * band_gt[1],
                        band_gt[3],
                    )
                    piece = shapely.clip_by_rect(part, *band_box)
                    inside = (
                        gis.Raster.rasterize_geometries(
                            [piece],
                            [1],
                            data.shape,
                            band_gt,
                            src.GetProjection(),
                            all_touched=True,
                            dtype=np.uint8,
                        )
                        > 0
                    )
                    data[~inside] = nodata
                    dst.GetRasterBand(1).WriteArray(data, 0, r)
                    progress((r + n_rows) / height, "", target.name)
                dst.FlushCache()
                dst = None
            t_convert = time.perf_counter() - start
            for name, value in [
                ("COMPRESS_OVERVIEW", compress),
                ("PREDICTOR_OVERVIEW", "3" if is_float else "2"),
                ("BIGTIFF_OVERVIEW", "YES"),
                ("GDAL_NUM_THREADS", "ALL_CPUS"),
            ]:
                gdal.SetConfigOption(name, value)
            start = time.perf_counter()
            converted = gdal.Open(str(target), gdal.GA_Update)
            converted.BuildOverviews(
                "AVERAGE",
                OVERVIEW_FACTORS,
                callback=progress,
                callback_data=target.name,
            )
            converted = None
            for name in [
                "COMPRESS_OVERVIEW",
                "PREDICTOR_OVERVIEW",
                "BIGTIFF_OVERVIEW",
                "GDAL_NUM_THREADS",
            ]:
                gdal.SetConfigOption(name, None)
            t_overviews = time.perf_counter() - start
            marker.write_text(
                json.dumps(
                    {
                        "source": info["path"],
                        "option": option,
                        "convert_s": t_convert,
                        "overviews_s": t_overviews,
                    }
                )
            )
            log.info(
                "Converted %s -> %s (%.2f GB) in %s, overviews in %s.",
                Path(info["path"]).name,
                target.name,
                os.path.getsize(target) / 1e9,
                gis.Run.format_seconds(t_convert),
                gis.Run.format_seconds(t_overviews),
            )
            targets.append(str(target))

    # --- 6. DTM used by the processing scripts ---
    vrt = gis.Raster.build_vrt(
        targets if targets else country["dtm"], work_dir / "dtm" / "dtm.vrt"
    )
    log.info(
        "Processing DTM: %s -> %s",
        vrt,
        "converted copies" if targets else "original files",
    )
