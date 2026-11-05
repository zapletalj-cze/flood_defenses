# Author: Jakub Zapletal
# Date: 2026-11-04
# Searches Copernicus GFM scenes over the processing units (one STAC search, plain HTTP and JSON), screens every scene
# once over the windows of all units it covers and copies the unit windows of scenes with enough flooded cells.

import datetime
import logging
import os
import re
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import url2pathname

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from osgeo import gdal
from shapely.geometry import shape

import gis

# --- Run variables (no command-line parameters; DEFENSES_<NAME> environment variables override them) ---
CONFIG_PATH = "config/defenses.toml"
COUNTRY = "XX"
FORCE = False  # screen the scenes of the selected units again, also those screened in earlier runs
UNITS = []  # unit ids; empty = all units

KEYS = [
    "general.work_dir",
    "general.output_dir",
    "general.warp_memory_limit_mb",
    "general.log_level",
    "country.working_resolution_m",
    "sources.gfm_stac_url",
    "sources.gfm_collection",
    "sources.gfm_assets",
    "sources.gfm_date_from",
    "sources.gfm_date_to",
    "sources.gfm_page_limit",
    "sources.gfm_token_env",
    "sources.gfm_screen_resolution_m",
    "sources.gfm_download_threads",
    "sources.gfm_nodata_values",
    "sources.request_pause_s",
    "sources.max_retries",
    "evidence.min_valid_pixels",
]
REQUIRED_ASSETS = ["flood", "exclusion", "reference_water"]
OPTIONAL_ASSETS = ["likelihood"]
COLUMNS = ["unit_id", "item_id", "datetime", "flooded_pixels", "kept"] + [
    f"{role}_path" for role in REQUIRED_ASSETS + OPTIONAL_ASSETS
]


def gdal_path(href):
    if href.startswith(("http://", "https://")):
        return f"/vsicurl/{href}"
    if href.startswith("file://"):
        # file:///C:/path on Windows, file:///path elsewhere
        return url2pathname(urlparse(href).path)
    if href.startswith("s3://"):
        return f"/vsis3/{href[len('s3://'):]}"
    return href


def screen_item(job):
    """
    Reads the flood layer of one scene once, at the screening resolution, over the windows of all units it covers,
    counts the flooded working cells in each unit's domain and copies the unit windows of all layers where there
    are more than min_valid_pixels. Runs in a thread: GDAL releases the GIL while it waits for the network.
    :param job: dict with item_id, datetime, units, grids, assets (role -> href) and params
    :return: one screening row per unit; empty if the scene could not be read (it is retried in the next run)
    """
    p = job["params"]
    screen = p["screen"]
    list_grid = [job["grids"][unit_id] for unit_id in job["units"]]
    xmin = min(g["x0"] for g in list_grid)
    ymax = max(g["y1"] for g in list_grid)
    xmax = max(g["x0"] + g["mask"].shape[1] * screen for g in list_grid)
    ymin = min(g["y1"] - g["mask"].shape[0] * screen for g in list_grid)
    safe_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", job["item_id"])
    list_rows = []
    try:
        # Mean of the 0/1 flood layer per screening cell = flooded share (nodata excluded).
        share, _ = gis.Raster.read_window(
            gdal_path(job["assets"]["flood"]),
            (xmin, ymin, xmax, ymax),
            screen,
            "average",
            p["warp_memory_mb"],
            dst_wkt=p["wkt"],
            src_nodata=p["nodata"],
        )
        for unit_id, grid in zip(job["units"], list_grid):
            r0 = int(round((ymax - grid["y1"]) / screen))
            c0 = int(round((grid["x0"] - xmin) / screen))
            part = share[
                r0 : r0 + grid["mask"].shape[0], c0 : c0 + grid["mask"].shape[1]
            ]
            flooded = (
                float(np.nansum(np.where(grid["mask"], part, 0.0)))
                * screen**2
                / p["resolution"] ** 2
            )
            row = {
                "unit_id": unit_id,
                "item_id": job["item_id"],
                "datetime": job["datetime"],
                "flooded_pixels": round(flooded),
                "kept": flooded > p["min_valid_pixels"],
            }
            for role in REQUIRED_ASSETS + OPTIONAL_ASSETS:
                row[f"{role}_path"] = ""
                if not row["kept"] or role not in job["assets"]:
                    continue
                target = Path(p["gfm_dir"]) / unit_id / f"{safe_id}_{role}.tif"
                target.parent.mkdir(parents=True, exist_ok=True)
                window = grid["window"]
                gdal.Translate(
                    str(target),
                    gdal_path(job["assets"][role]),
                    projWin=[window[0], window[3], window[2], window[1]],
                    projWinSRS=p["wkt"],
                    creationOptions=[
                        "TILED=YES",
                        f"COMPRESS={gis.Raster.best_compression()}",
                        "PREDICTOR=2",
                    ],
                )
                row[f"{role}_path"] = str(target)
            list_rows.append(row)
    except RuntimeError as error:
        logging.getLogger("download_gfm").warning(
            "Item %s unreadable (%s); it is retried in the next run.",
            job["item_id"],
            error,
        )
        return []
    time.sleep(p["pause"])
    return list_rows


if __name__ == "__main__":
    run = gis.Run.env_overrides(
        CONFIG_PATH=CONFIG_PATH, COUNTRY=COUNTRY, FORCE=FORCE, UNITS=UNITS
    )
    cfg = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], KEYS)
    general, country, sources = cfg["general"], cfg["country"], cfg["sources"]
    work_dir = Path(general["work_dir"]) / run["COUNTRY"]
    out_dir = Path(general["output_dir"]) / run["COUNTRY"]
    log = gis.Run.logger("download_gfm", out_dir / "logs", general["log_level"])
    gis.Config.warn_placeholders(cfg)
    assets = sources["gfm_assets"]
    list_miss = [k for k in ["gfm_stac_url", "gfm_collection"] if not sources[k]]
    list_miss += [f"gfm_assets.{k}" for k in REQUIRED_ASSETS if not assets.get(k)]
    if list_miss:
        raise SystemExit(
            f"GFM sources not configured: {list_miss}. The STAC endpoint, collection and asset names have not been "
            "verified against the GFM documentation; fill them in [sources]."
        )
    list_need = [
        work_dir / "dtm" / "dtm.vrt",
        out_dir / "units.csv",
        work_dir / "units_domain.gpkg",
    ]
    if not all(path.exists() for path in list_need):
        raise SystemExit(
            "dtm.vrt, units.csv or units_domain.gpkg is missing; run prepare_dtm.py and build_units.py."
        )
    wkt = gdal.Open(str(list_need[0])).GetProjection()
    screen = sources["gfm_screen_resolution_m"]
    pause = sources["request_pause_s"]
    headers = {}
    if sources["gfm_token_env"]:
        token = os.environ.get(sources["gfm_token_env"], "")
        if not token:
            raise SystemExit(
                f"Environment variable {sources['gfm_token_env']} with the GFM token is not set."
            )
        headers["Authorization"] = f"Bearer {token}"
        gdal.SetConfigOption("GDAL_HTTP_HEADERS", f"Authorization: Bearer {token}")
    gdal.SetConfigOption("GDAL_HTTP_MAX_RETRY", str(sources["max_retries"]))
    gdal.SetConfigOption("GDAL_HTTP_RETRY_DELAY", str(pause))
    gdal.SetConfigOption("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
    date_to = sources["gfm_date_to"] or datetime.date.today().isoformat()
    period = f"{sources['gfm_date_from']}T00:00:00Z/{date_to}T23:59:59Z"

    df_all = pd.read_csv(out_dir / "units.csv", dtype={"unit_id": str})
    df_units = df_all[df_all["unit_id"].isin(run["UNITS"])] if run["UNITS"] else df_all
    gdf_dom = gpd.read_file(work_dir / "units_domain.gpkg").set_index("unit_id")
    gfm_dir = work_dir / "gfm"
    gfm_dir.mkdir(parents=True, exist_ok=True)
    # Checkpoint: every screened pair of scene and unit, kept or not; a rerun screens only new pairs.
    scr_path = gfm_dir / "screening.csv"
    df_scr = (
        pd.read_csv(scr_path, dtype={"unit_id": str}, keep_default_na=False)
        if scr_path.exists()
        else pd.DataFrame(columns=COLUMNS)
    )
    if run["FORCE"]:
        df_scr = df_scr[~df_scr["unit_id"].isin(df_units["unit_id"])]
        df_scr.to_csv(scr_path, index=False)
    done = set(zip(df_scr["item_id"], df_scr["unit_id"]))

    # --- Screening grid per unit: its window on a grid of screening cells shared by all units, domain mask ---
    grids, list_core = {}, []
    for unit in df_units.itertuples():
        x0 = float(np.floor(unit.win_xmin / screen) * screen)
        x1 = float(np.ceil(unit.win_xmax / screen) * screen)
        y0 = float(np.floor(unit.win_ymin / screen) * screen)
        y1 = float(np.ceil(unit.win_ymax / screen) * screen)
        grid_shape = (int(round((y1 - y0) / screen)), int(round((x1 - x0) / screen)))
        mask = (
            gis.Raster.rasterize_geometries(
                [gdf_dom.geometry[unit.unit_id]],
                [1],
                grid_shape,
                (x0, screen, 0.0, y1, 0.0, -screen),
                wkt,
                True,
                np.uint8,
            )
            > 0
        )
        grids[unit.unit_id] = {
            "x0": x0,
            "y1": y1,
            "mask": mask,
            "window": (unit.win_xmin, unit.win_ymin, unit.win_xmax, unit.win_ymax),
        }
        list_core.append(
            shapely.segmentize(
                shapely.box(
                    unit.core_xmin, unit.core_ymin, unit.core_xmax, unit.core_ymax
                ),
                500.0,
            )
        )
    cores = gpd.GeoSeries(
        list_core, index=df_units["unit_id"].to_numpy(), crs=wkt
    ).to_crs(4326)

    # --- One STAC item search over all selected units, with paging ---
    items = {}
    response = gis.Http.get_json(
        f"{sources['gfm_stac_url'].rstrip('/')}/search",
        params={
            "collections": [sources["gfm_collection"]],
            "bbox": [round(float(v), 6) for v in cores.total_bounds],
            "datetime": period,
            "limit": sources["gfm_page_limit"],
        },
        retries=sources["max_retries"],
        pause_s=pause,
        headers=headers,
    )
    while True:
        for feature in response.get("features", []):
            items[feature["id"]] = feature
        list_next = [
            link for link in response.get("links", []) if link.get("rel") == "next"
        ]
        if not list_next:
            break
        time.sleep(pause)
        link = list_next[0]
        response = gis.Http.get_json(
            link["href"],
            body=(
                link.get("body")
                if str(link.get("method", "GET")).upper() == "POST"
                else None
            ),
            retries=sources["max_retries"],
            pause_s=pause,
            headers=headers,
        )

    # --- Scenes and the units whose core their footprint touches ---
    params = {
        "screen": screen,
        "resolution": country["working_resolution_m"],
        "wkt": wkt,
        "nodata": (
            sources["gfm_nodata_values"][0] if sources["gfm_nodata_values"] else None
        ),
        "warp_memory_mb": general["warp_memory_limit_mb"],
        "min_valid_pixels": cfg["evidence"]["min_valid_pixels"],
        "gfm_dir": str(gfm_dir),
        "pause": pause,
    }
    list_job = []
    for item_id, item in sorted(
        items.items(), key=lambda kv: kv[1]["properties"].get("datetime") or ""
    ):
        available = item.get("assets", {})
        list_abs = [assets[r] for r in REQUIRED_ASSETS if assets[r] not in available]
        if list_abs:
            log.warning("Item %s lacks assets %s; skipped.", item_id, list_abs)
            continue
        footprint = (
            shape(item["geometry"])
            if item.get("geometry")
            else shapely.box(*item["bbox"])
        )
        list_mem = [
            u
            for u in cores.index[cores.intersects(footprint)]
            if (item_id, u) not in done
        ]
        if not list_mem:
            continue
        list_job.append(
            {
                "item_id": item_id,
                "datetime": item["properties"].get("datetime")
                or item["properties"].get("start_datetime"),
                "units": list_mem,
                "grids": {u: grids[u] for u in list_mem},
                "assets": {
                    role: available[assets[role]]["href"]
                    for role in REQUIRED_ASSETS + OPTIONAL_ASSETS
                    if assets.get(role) and assets[role] in available
                },
                "params": params,
            }
        )
    log.info(
        "%d GFM items found for %d units; %d items to screen (%d scene-unit pairs screened in earlier runs).",
        len(items),
        len(df_units),
        len(list_job),
        len(done),
    )

    # --- Screening and copies, several scenes at a time ---
    start, kept = time.perf_counter(), 0
    with ThreadPoolExecutor(max_workers=sources["gfm_download_threads"]) as pool:
        for k, list_rows in enumerate(pool.map(screen_item, list_job), start=1):
            if list_rows:
                pd.DataFrame(list_rows, columns=COLUMNS).to_csv(
                    scr_path,
                    mode="a",
                    header=not scr_path.exists(),
                    index=False,
                )
                kept += sum(bool(r["kept"]) for r in list_rows)
            if k % 50 == 0 or k == len(list_job):
                log.info(
                    "Screened %d of %d items in %s; %d unit windows kept.",
                    k,
                    len(list_job),
                    gis.Run.format_seconds(time.perf_counter() - start),
                    kept,
                )

    # --- Index of the kept unit windows (all runs, all units) ---
    df_scr = (
        pd.read_csv(scr_path, dtype={"unit_id": str}, keep_default_na=False)
        if scr_path.exists()
        else pd.DataFrame(columns=COLUMNS)
    )
    keep = df_scr["kept"].astype(str).str.lower().isin(["true", "1"])
    df_idx = df_scr[keep & df_scr["unit_id"].isin(df_all["unit_id"])]
    df_idx = df_idx.drop(columns="kept").sort_values(["unit_id", "datetime"])
    df_idx.to_csv(out_dir / "gfm_index.csv", index=False)
    log.info(
        "gfm_index.csv: %d scene windows in %d units (screened pairs: %d).",
        len(df_idx),
        df_idx["unit_id"].nunique() if len(df_idx) else 0,
        len(df_scr),
    )
