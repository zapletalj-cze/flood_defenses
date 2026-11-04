# Author: Jakub Zapletal
# Date: 2026-11-04
# Inspects and benchmarks the national DTM files (and DSM files when configured), derives their valid-data footprint,
# estimates runtimes and disk needs, and optionally converts the DTM to a tiled layout with overviews.

import json
import os
import shutil
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from osgeo import gdal, gdal_array, osr

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
    "country.dsm",
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
    list_info = [gis.Raster.inspect(path) for path in country["dtm"]]
    gis.Raster.build_vrt(
        country["dtm"], work_dir / "dtm" / "dtm_source.vrt"
    )  # stops on CRS mismatch
    for info in list_info:
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
    if len({round(abs(i["pixel_x"]), 6) for i in list_info}) > 1:
        log.warning("DTM files differ in pixel size; the VRT uses the finest one.")
    pd.DataFrame(list_info).drop(columns=["crs_wkt"]).to_csv(
        out_dir / "dtm_inspection.csv", index=False
    )

    # --- 2. Valid-data footprint ---
    fp_path = out_dir / "dtm_footprint.gpkg"
    if fp_path.exists() and not run["FORCE"]:
        gdf_fp = gpd.read_file(fp_path)
        log.info("Footprint reused from %s (set FORCE to recompute).", fp_path)
    else:
        list_fp = []
        for info in list_info:
            start = time.perf_counter()
            gdf_part = gis.Raster.footprint(
                info["path"], country["footprint_resolution_m"]
            )
            gdf_part["dtm"] = Path(info["path"]).stem
            gdf_part["path"] = info["path"]
            list_fp.append(gdf_part)
            log.info(
                "Footprint of %s derived in %.1f s.",
                Path(info["path"]).name,
                time.perf_counter() - start,
            )
        gdf_fp = gpd.GeoDataFrame(
            pd.concat(list_fp, ignore_index=True), crs=list_fp[0].crs
        )
        gdf_fp["area_km2"] = gdf_fp.area / 1e6
        gis.Vector.write_gpkg(gdf_fp, fp_path)
    for info in list_info:
        extent_km2 = (
            info["width"]
            * info["height"]
            * abs(info["pixel_x"] * info["pixel_y"])
            / 1e6
        )
        area_km2 = gdf_fp.loc[gdf_fp["path"] == info["path"]].area.sum() / 1e6
        log.info(
            "Valid data of %s: %.1f km2 of %.1f km2 raster extent (%.1f %%).",
            Path(info["path"]).name,
            area_km2,
            extent_km2,
            100 * area_km2 / extent_km2,
        )

    # --- 2b. DSM (Stage D): same inspection, same CRS as the DTM, footprint, VRT ---
    list_dsm_info = [gis.Raster.inspect(path) for path in country["dsm"]]
    gdf_dsm_fp = gpd.GeoDataFrame(columns=["path", "geometry"], geometry="geometry")
    if list_dsm_info:
        dtm_srs = osr.SpatialReference(wkt=list_info[0]["crs_wkt"])
        for info in list_dsm_info:
            if not osr.SpatialReference(wkt=info["crs_wkt"]).IsSame(dtm_srs):
                raise SystemExit(
                    f"DSM {info['path']} has a different CRS than the DTM; reproject it first."
                )
            log.info(
                "DSM %s: %d x %d px, pixel %.3f x %.3f, %s, nodata %s, compression %s, block %d x %d%s, "
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
        gis.Raster.build_vrt(country["dsm"], work_dir / "dsm" / "dsm.vrt")
        pd.DataFrame(list_dsm_info).drop(columns=["crs_wkt"]).to_csv(
            out_dir / "dsm_inspection.csv", index=False
        )
        dsm_fp_path = out_dir / "dsm_footprint.gpkg"
        if dsm_fp_path.exists() and not run["FORCE"]:
            gdf_dsm_fp = gpd.read_file(dsm_fp_path)
        else:
            list_fp = []
            for info in list_dsm_info:
                gdf_part = gis.Raster.footprint(
                    info["path"], country["footprint_resolution_m"]
                )
                gdf_part["path"] = info["path"]
                list_fp.append(gdf_part)
            gdf_dsm_fp = gpd.GeoDataFrame(
                pd.concat(list_fp, ignore_index=True), crs=list_fp[0].crs
            )
            gis.Vector.write_gpkg(gdf_dsm_fp, dsm_fp_path)

    # --- 3. Benchmark (DTM and, when configured, DSM: same windows, same measures) ---
    estimates = {"country": run["COUNTRY"], "dtm": [], "dsm": []}
    if run["RUN_BENCHMARK"]:
        rng = np.random.default_rng(bench["seed"])
        for kind, list_kinfo, gdf_kfp in [
            ("dtm", list_info, gdf_fp),
            ("dsm", list_dsm_info, gdf_dsm_fp),
        ]:
            list_rec = []
            for info in list_kinfo:
                polygon = shapely.union_all(
                    gdf_kfp.loc[gdf_kfp["path"] == info["path"]].geometry.values
                )
                list_rows, summary, list_sw = gis.Raster.benchmark_reads(
                    info,
                    polygon,
                    bench["n_windows"],
                    bench["window_px"],
                    resolution,
                    general["warp_memory_limit_mb"],
                    bench["sequential_rows"],
                    compress,
                    rng,
                )
                if len(list_rows) < bench["n_windows"]:
                    log.warning(
                        "Only %d benchmark windows fit inside the footprint of %s.",
                        len(list_rows),
                        info["path"],
                    )
                if not list_sw:
                    log.warning(
                        "Could not determine from the GDAL debug output whether gdal.Warp used an "
                        "overview for %s.",
                        Path(info["path"]).name,
                    )
                elif min(list_sw) < 0.9 * summary["window_px"]:
                    log.info(
                        "gdal.Warp read %s from an overview: source windows of %d px instead of %d "
                        "native px (factor %.0f).",
                        Path(info["path"]).name,
                        min(list_sw),
                        summary["window_px"],
                        summary["window_px"] / min(list_sw),
                    )
                else:
                    log.info(
                        "gdal.Warp read %s at full resolution (source window %d px, %d overviews): "
                        "every native pixel is decoded.",
                        Path(info["path"]).name,
                        min(list_sw),
                        info["overviews"],
                    )
                list_rec += list_rows
                estimates[kind].append(summary)
                log.info(
                    "Benchmark %s %s (%d windows of %d x %d px = %.2f km2): native read %.2f s/window (p90 %.2f), "
                    "%.0f MB/s; warp to %.1f m %.2f s/window (p90 %.2f), %.0f MB/s native equivalent; sequential "
                    "read %.0f MB/s; %s compression ratio %.2f at %.0f MB/s.",
                    kind.upper(),
                    Path(info["path"]).name,
                    len(list_rows),
                    summary["window_px"],
                    summary["window_px"],
                    summary["window_ground_km2"],
                    summary["native_s_median"],
                    summary["native_s_p90"],
                    summary["native_mb_s_median"],
                    resolution,
                    summary["warp_s_median"],
                    summary["warp_s_p90"],
                    summary["warp_native_equivalent_mb_s_median"],
                    summary["sequential_mb_s"],
                    compress,
                    summary["compression_ratio"],
                    summary["compression_mb_s"],
                )
            if list_kinfo:
                pd.DataFrame(list_rec).to_csv(
                    out_dir / f"{kind}_benchmark.csv", index=False
                )
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
            df_units = pd.read_csv(units_csv, dtype={"unit_id": str})
            window_km2 = df_units["est_window_cells"] * resolution**2 / 1e6
            read_total = float((window_km2 / bench_km2 * warp_s).sum())
            log.info(
                "Terrain read estimate: %d units, %.0f s per unit on average, %s in total on %d workers "
                "(source: median warp s/window x unit window area / benchmark window area; assumes reads scale with "
                "workers, which the disk may not allow).",
                len(df_units),
                read_total / len(df_units),
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
                part = gdf_fp.loc[gdf_fp["path"] == summary["path"]].union_all()
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
    list_tgt = []
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
        for info in list_info:
            target = work_dir / "dtm" / f"{Path(info['path']).stem}_{option}.tif"
            marker = target.with_suffix(".done")
            if marker.exists() and not run["FORCE"]:
                log.info(
                    "Conversion of %s already finished: %s",
                    Path(info["path"]).name,
                    target,
                )
                list_tgt.append(str(target))
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
            list_co = [
                "TILED=YES",
                "BLOCKXSIZE=512",
                "BLOCKYSIZE=512",
                f"COMPRESS={compress}",
                f"PREDICTOR={3 if is_float else 2}",
                "BIGTIFF=YES",
                "NUM_THREADS=ALL_CPUS",
            ]
            start = time.perf_counter()
            list_last = [start]

            def progress(complete, message, data):
                now = time.perf_counter()
                if now - list_last[0] > 60:
                    list_last[0] = now
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
                    creationOptions=list_co,
                    callback=progress,
                    callback_data=target.name,
                )
            else:
                part = corridor.intersection(
                    gdf_fp.loc[gdf_fp["path"] == info["path"]].union_all()
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
                    str(target), width, height, 1, band.DataType, list_co
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
            list_tgt.append(str(target))

    # --- 6. DTM used by the processing scripts ---
    vrt = gis.Raster.build_vrt(
        list_tgt if list_tgt else country["dtm"], work_dir / "dtm" / "dtm.vrt"
    )
    log.info(
        "Processing DTM: %s -> %s",
        vrt,
        "converted copies" if list_tgt else "original files",
    )
