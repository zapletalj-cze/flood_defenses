# Author: Jakub Zapletal
# Date: 2026-11-04
# Gets GloFAS historical daily discharge for the country (EWDS through cdsapi when installed, or local files) and the
# static upstream area map, and stores the river cells with enough upstream area compactly in one .npz file.

import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from osgeo import gdal
from shapely.geometry import box

import gis

# --- Run variables (no command-line parameters; DEFENSES_<NAME> environment variables override them) ---
CONFIG_PATH = "config/defenses.toml"
COUNTRY = "XX"
FORCE = False  # download years again that already have a file

KEYS = [
    "general.work_dir",
    "general.output_dir",
    "general.log_level",
    "sources.glofas_enabled",
    "sources.glofas_dataset_id",
    "sources.glofas_request",
    "sources.glofas_time_keys",
    "sources.glofas_area_key",
    "sources.glofas_first_year",
    "sources.glofas_last_year",
    "sources.glofas_file_extension",
    "sources.glofas_files",
    "sources.glofas_discharge_variable",
    "sources.glofas_upstream_area_file",
    "sources.glofas_upstream_area_url",
    "sources.glofas_upstream_area_variable",
    "sources.glofas_upstream_area_to_km2",
    "sources.glofas_min_upstream_area_km2",
    "sources.glofas_valid_time_ends_day",
    "sources.request_pause_s",
    "sources.max_retries",
    "evidence.glofas_search_radius_km",
]

if __name__ == "__main__":
    run = gis.Run.env_overrides(CONFIG_PATH=CONFIG_PATH, COUNTRY=COUNTRY, FORCE=FORCE)
    cfg = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], KEYS)
    general, sources = cfg["general"], cfg["sources"]
    work_dir = Path(general["work_dir"]) / run["COUNTRY"] / "glofas"
    out_dir = Path(general["output_dir"]) / run["COUNTRY"]
    log = gis.Run.logger("download_glofas", out_dir / "logs", general["log_level"])
    gis.Config.warn_placeholders(cfg)
    if not sources["glofas_enabled"]:
        log.info("glofas_enabled is false; nothing to do.")
        raise SystemExit(0)
    footprint_path = out_dir / "dtm_footprint.gpkg"
    if not footprint_path.exists():
        raise SystemExit("dtm_footprint.gpkg is missing; run prepare_dtm.py first.")
    footprint = gpd.read_file(footprint_path)
    extent = gpd.GeoSeries([box(*footprint.total_bounds)], crs=footprint.crs)
    west, south, east, north = (
        extent.segmentize(extent.iloc[0].length / 400).to_crs(4326).total_bounds
    )
    # Margin of the search radius, so that river cells next to the DTM footprint are kept.
    margin_lat = cfg["evidence"]["glofas_search_radius_km"] / 111.0
    margin_lon = margin_lat / max(np.cos(np.radians((south + north) / 2)), 0.1)
    west, east, south, north = (
        west - margin_lon,
        east + margin_lon,
        south - margin_lat,
        north + margin_lat,
    )

    # --- Discharge files: local list or download per year ---
    files = sources["glofas_files"]
    if not files:
        if not sources["glofas_request"]:
            raise SystemExit(
                "sources.glofas_request is empty. Paste the request of the EWDS form for "
                f"{sources['glofas_dataset_id']} (button 'Show API request', without time and area) into the "
                "config, or download the files manually and list them in sources.glofas_files."
            )
        try:
            import cdsapi
        except ImportError:
            raise SystemExit(
                "The cdsapi client is not installed (it is not part of the environment and needs approval). "
                "Download the files through the EWDS form and list them in sources.glofas_files."
            )
        client = cdsapi.Client()
        keys = sources["glofas_time_keys"]
        work_dir.mkdir(parents=True, exist_ok=True)
        for year in range(
            sources["glofas_first_year"], sources["glofas_last_year"] + 1
        ):
            target = work_dir / f"discharge_{year}.{sources['glofas_file_extension']}"
            if target.exists() and not run["FORCE"]:
                continue
            request = dict(sources["glofas_request"])
            request[keys["year"]] = [str(year)]
            request[keys["month"]] = [f"{m:02d}" for m in range(1, 13)]
            request[keys["day"]] = [f"{d:02d}" for d in range(1, 32)]
            request[sources["glofas_area_key"]] = [
                float(np.ceil(north * 100) / 100),
                float(np.floor(west * 100) / 100),
                float(np.floor(south * 100) / 100),
                float(np.ceil(east * 100) / 100),
            ]
            log.info("Requesting GloFAS discharge %d.", year)
            client.retrieve(sources["glofas_dataset_id"], request, str(target))
            time.sleep(sources["request_pause_s"])
        files = sorted(
            str(p)
            for p in work_dir.glob(f"discharge_*.{sources['glofas_file_extension']}")
        )

    # --- Upstream area and the river cells kept ---
    upstream_path = sources["glofas_upstream_area_file"]
    if not upstream_path:
        if not sources["glofas_upstream_area_url"]:
            raise SystemExit(
                "Set sources.glofas_upstream_area_file or glofas_upstream_area_url (not verified)."
            )
        upstream_path = gis.Http.download(
            sources["glofas_upstream_area_url"],
            work_dir / Path(sources["glofas_upstream_area_url"]).name,
            sources["max_retries"],
            sources["request_pause_s"],
        )
    variable = sources["glofas_upstream_area_variable"]
    upstream = gdal.Open(
        f'NETCDF:"{upstream_path}":{variable}' if variable else str(upstream_path)
    )
    up_gt = upstream.GetGeoTransform()
    col0 = max(int(np.floor((west - up_gt[0]) / up_gt[1])), 0)
    col1 = min(int(np.ceil((east - up_gt[0]) / up_gt[1])), upstream.RasterXSize)
    row0 = max(int(np.floor((north - up_gt[3]) / up_gt[5])), 0)
    row1 = min(int(np.ceil((south - up_gt[3]) / up_gt[5])), upstream.RasterYSize)
    up_band = upstream.GetRasterBand(1)
    area = up_band.ReadAsArray(col0, row0, col1 - col0, row1 - row0).astype(np.float64)
    if up_band.GetNoDataValue() is not None:
        area[area == up_band.GetNoDataValue()] = np.nan
    area *= sources["glofas_upstream_area_to_km2"]
    rows, cols = np.nonzero(
        np.isfinite(area) & (area >= sources["glofas_min_upstream_area_km2"])
    )
    lon = up_gt[0] + (col0 + cols + 0.5) * up_gt[1]
    lat = up_gt[3] + (row0 + rows + 0.5) * up_gt[5]
    log.info(
        "GloFAS cells in the area: %d with upstream area >= %.0f km2.",
        len(rows),
        sources["glofas_min_upstream_area_km2"],
    )
    if not len(rows):
        raise SystemExit(
            "No GloFAS cell reaches glofas_min_upstream_area_km2 in the area."
        )

    # --- Daily discharge of the kept cells ---
    dates, series = [], []
    for path in files:
        name = (
            f'NETCDF:"{path}":{sources["glofas_discharge_variable"]}'
            if sources["glofas_discharge_variable"]
            else str(path)
        )
        ds = gdal.Open(name)
        gt = ds.GetGeoTransform()
        c = np.floor((lon - gt[0]) / gt[1]).astype(int)
        r = np.floor((lat - gt[3]) / gt[5]).astype(int)
        if (
            (c < 0).any()
            or (r < 0).any()
            or (c >= ds.RasterXSize).any()
            or (r >= ds.RasterYSize).any()
        ):
            raise SystemExit(
                f"{path} does not cover all GloFAS cells of the area; request a larger area."
            )
        x0, y0 = int(c.min()), int(r.min())
        width, height = int(c.max() - x0 + 1), int(r.max() - y0 + 1)
        metadata = ds.GetMetadata()
        for b in range(1, ds.RasterCount + 1):
            band = ds.GetRasterBand(b)
            band_meta = band.GetMetadata()
            if "GRIB_VALID_TIME" in band_meta:
                stamp = pd.Timestamp(
                    int(band_meta["GRIB_VALID_TIME"].split()[0]), unit="s"
                )
            else:
                dimension = next(
                    (
                        k
                        for k in band_meta
                        if k.startswith("NETCDF_DIM_") and "time" in k.lower()
                    ),
                    None,
                )
                if dimension is None:
                    raise SystemExit(f"No time information in band {b} of {path}.")
                units = metadata.get(f"{dimension[len('NETCDF_DIM_'):]}#units", "")
                unit, _, reference = units.partition(" since ")
                step = {
                    "days": "D",
                    "day": "D",
                    "hours": "h",
                    "hour": "h",
                    "seconds": "s",
                    "second": "s",
                }.get(unit.strip())
                if step is None or not reference:
                    raise SystemExit(f"Unsupported time units {units!r} in {path}.")
                stamp = pd.Timestamp(reference.strip()) + pd.to_timedelta(
                    float(band_meta[dimension]), unit=step
                )
            values = band.ReadAsArray(x0, y0, width, height).astype(np.float32)
            if band.GetNoDataValue() is not None:
                values[values == np.float32(band.GetNoDataValue())] = np.nan
            # A 24 h value stamped at the end of its period belongs to the day that ends (config switch).
            if sources["glofas_valid_time_ends_day"]:
                stamp -= pd.Timedelta(seconds=1)
            dates.append(np.datetime64(stamp.date(), "D"))
            series.append(values[r - y0, c - x0])
        log.info("Read %s: %d time steps.", Path(path).name, ds.RasterCount)
    if not dates:
        raise SystemExit("No discharge data read.")
    dates = np.array(dates)
    order = np.argsort(dates, kind="stable")
    dates, unique = np.unique(dates[order], return_index=True)
    discharge = np.vstack(series)[order][unique]
    expected = (dates[-1] - dates[0]).astype(int) + 1
    if len(dates) < expected:
        log.warning(
            "%d days missing between %s and %s.",
            expected - len(dates),
            dates[0],
            dates[-1],
        )
    target = out_dir / "glofas_discharge.npz"
    np.savez_compressed(
        target,
        lon=lon,
        lat=lat,
        upstream_area_km2=area[rows, cols],
        dates=dates.astype("datetime64[D]").astype(np.int64),
        discharge=discharge.astype(np.float32),
    )
    log.info(
        "glofas_discharge.npz: %d cells, %d days (%s to %s), %.1f MB.",
        len(lon),
        len(dates),
        dates[0],
        dates[-1],
        target.stat().st_size / 1e6,
    )
