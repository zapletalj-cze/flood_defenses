# Author: Jakub Zapletal
# Date: 2026-11-05
# Stage F: downloads GloFAS and EFAS historical discharge from the Copernicus EWDS (cdsapi) for the country bounding
# box, in requests of at most 500 layers, in parallel with resume, and copies the static upstream area files.

import json
import os
import shutil
import threading
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import url2pathname

import geopandas as gpd
import pandas as pd
from osgeo import gdal

import gis

gdal.UseExceptions()

# --- Run variables (no command-line parameters; DEFENSES_<NAME> environment variables override them) ---
CONFIG_PATH = "config/defenses.toml"
COUNTRY = "XX"

KEYS = [
    "general.work_dir",
    "general.output_dir",
    "general.log_level",
    "discharge.sources",
    "discharge.ewds_api_url",
    "discharge.raw_dir",
    "discharge.max_parallel_requests",
    "discharge.bbox_buffer_deg",
    *[
        f"discharge.{s}.{k}"
        for s in ("glofas", "efas")
        for k in (
            "dataset_id",
            "request",
            "data_format",
            "download_format",
            "first_year",
            "last_year",
            "request_unit",
            "time_keys",
            "area_key",
            "step_hours",
            "upstream_area_file",
            "netcdf_variable",
        )
    ],
]
MAX_LAYERS = 500  # verified: one EWDS request returns at most 500 layers
RETRIES = 4  # attempts per request; waits of 30 s, 60 s, 120 s between them
EXTENSIONS = {"grib": "grib", "grib2": "grib", "netcdf": "nc", "netcdf4": "nc"}


if __name__ == "__main__":
    run = gis.Run.env_overrides(CONFIG_PATH=CONFIG_PATH, COUNTRY=COUNTRY)
    cfg = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], KEYS)
    general, dis = cfg["general"], cfg["discharge"]
    work_dir = Path(general["work_dir"]) / run["COUNTRY"]
    out_dir = Path(general["output_dir"]) / run["COUNTRY"]
    log = gis.Run.logger("download_discharge", out_dir / "logs", general["log_level"])
    gis.Config.warn_placeholders(cfg)
    raw_dir = Path(dis["raw_dir"]) if dis["raw_dir"] else work_dir / "discharge" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()

    # --- Area: country bounding box plus bbox_buffer_deg, as north, west, south, east in degrees ---
    fp_path = out_dir / "dtm_footprint.gpkg"
    if not fp_path.exists():
        raise SystemExit("dtm_footprint.gpkg missing; run prepare_dtm.py first.")
    west, south, east, north = gpd.read_file(fp_path).to_crs(4326).total_bounds
    b = dis["bbox_buffer_deg"]
    area = [
        round(north + b, 4),
        round(west - b, 4),
        round(south - b, 4),
        round(east + b, 4),
    ]
    log.info("Area (north, west, south, east): %s.", area)

    # --- Requests per source, split by request_unit so that none exceeds 500 layers ---
    list_job, list_skip_src = [], []
    for src in dis["sources"]:
        c = dis[src]
        list_miss = [
            k
            for k, ok in (
                ("dataset_id", bool(c["dataset_id"])),
                ("request", bool(c["request"])),
                ("data_format", c["data_format"] in EXTENSIONS),
                ("first_year", c["first_year"] > 0),
                ("last_year", c["last_year"] >= c["first_year"] > 0),
                (
                    "time_keys (year, month, day)",
                    {"year", "month", "day"} <= set(c["time_keys"]),
                ),
                ("area_key", bool(c["area_key"])),
            )
            if not ok
        ]
        if list_miss:
            # Values come from the EWDS download form ("Show API request") and are never guessed.
            log.error(
                "%s: configuration incomplete %s; copy them from the EWDS download form of %s. Source skipped.",
                src,
                list_miss,
                c["dataset_id"] or src,
            )
            list_skip_src.append(src)
            continue
        if c["dataset_id"] == "cems-glofas-historical" and c["request"].get(
            "product_type"
        ) not in (
            None,
            ["consolidated"],
            "consolidated",
        ):
            raise SystemExit(
                "GloFAS: product_type must be consolidated only (no mixing with intermediate data)."
            )
        steps = 24 // c["step_hours"]
        tk = c["time_keys"]
        for year in range(c["first_year"], c["last_year"] + 1):
            months = [None] if c["request_unit"] == "year" else list(range(1, 13))
            for month in months:
                month_list = list(range(1, 13)) if month is None else [month]
                days = pd.date_range(
                    f"{year}-{month_list[0]:02d}-01",
                    pd.Timestamp(f"{year}-{month_list[-1]:02d}-01")
                    + pd.offsets.MonthEnd(0),
                    freq="D",
                )
                layers = len(days) * steps
                if layers > MAX_LAYERS:
                    raise SystemExit(
                        f"{src}: {layers} layers per request exceed {MAX_LAYERS}; use request_unit = month."
                    )
                request = {
                    **c["request"],
                    tk["year"]: [str(year)],
                    tk["month"]: [f"{m:02d}" for m in month_list],
                    tk["day"]: [f"{d:02d}" for d in range(1, 32)],
                    c["area_key"]: area,
                    "data_format": c["data_format"],
                }
                if c["download_format"]:
                    request["download_format"] = c["download_format"]
                name = f"{year}" if month is None else f"{year}_{month:02d}"
                list_job.append(
                    {
                        "source": src,
                        "dataset": c["dataset_id"],
                        "request": request,
                        "target": raw_dir
                        / src
                        / f"{name}.{EXTENSIONS[c['data_format']]}",
                        "layers": layers,
                        "variable": c["netcdf_variable"],
                    }
                )
    if not list_job:
        raise SystemExit(f"No source to download (skipped: {list_skip_src}).")

    # --- Client: cdsapi against the EWDS, or a local folder (file://) holding the same file names ---
    local = dis["ewds_api_url"].startswith("file://")
    if local:
        # file:///C:/folder on Windows, file:///folder elsewhere
        local_root = Path(url2pathname(urlparse(dis["ewds_api_url"]).path))
        log.info(
            "ewds_api_url is a local folder: files are copied from %s.", local_root
        )
    else:
        try:
            import cdsapi
        except ImportError:
            raise SystemExit(
                "cdsapi is not installed (pixi.toml); it is needed for the EWDS downloads."
            )
        # The personal key comes from ~/.cdsapirc or the CDSAPI_KEY environment variable, never from this repository.
        client = cdsapi.Client(url=dis["ewds_api_url"], quiet=True, progress=False)

    def layer_count(path, variable):
        # Integrity check: the file (its discharge variable for NetCDF) opens with GDAL and holds the requested
        # number of time steps.
        try:
            name = (
                f'NETCDF:"{path}":{variable}'
                if variable and path.suffix == ".nc"
                else str(path)
            )
            return gdal.Open(name).RasterCount
        except Exception:  # noqa: BLE001 - any unreadable file counts as incomplete
            return -1

    lock = threading.Lock()
    list_log = []

    def fetch(job):
        target = job["target"]
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists() and layer_count(target, job["variable"]) == job["layers"]:
            return {
                **job,
                "status": "complete, skipped",
                "attempts": 0,
                "seconds": 0.0,
                "error": "",
            }
        begin, error = time.perf_counter(), ""
        for attempt in range(1, RETRIES + 1):
            tmp = target.with_suffix(target.suffix + ".part")
            try:
                if local:
                    shutil.copy(local_root / job["source"] / target.name, tmp)
                else:
                    client.retrieve(job["dataset"], job["request"], str(tmp))
                if zipfile.is_zipfile(tmp):
                    # download_format = zip: the archive holds one file
                    with zipfile.ZipFile(tmp) as archive, archive.open(
                        archive.namelist()[0]
                    ) as member:
                        with open(target, "wb") as f:
                            shutil.copyfileobj(member, f)
                    tmp.unlink()
                else:
                    os.replace(tmp, target)
                found = layer_count(target, job["variable"])
                if found != job["layers"]:
                    raise ValueError(f"{found} time steps instead of {job['layers']}")
                return {
                    **job,
                    "status": "downloaded",
                    "attempts": attempt,
                    "seconds": time.perf_counter() - begin,
                    "error": "",
                }
            except Exception as exc:  # noqa: BLE001 - retried, then logged
                error = str(exc)[:300]
                tmp.unlink(missing_ok=True)
                if attempt < RETRIES:
                    time.sleep(30 * 2 ** (attempt - 1) if not local else 0)
        return {
            **job,
            "status": "failed",
            "attempts": RETRIES,
            "seconds": time.perf_counter() - begin,
            "error": error,
        }

    with ThreadPoolExecutor(max_workers=dis["max_parallel_requests"]) as pool:
        for k, res in enumerate(pool.map(fetch, list_job), start=1):
            size = res["target"].stat().st_size / 1e6 if res["target"].exists() else 0.0
            with lock:
                list_log.append(
                    {
                        "source": res["source"],
                        "file": res["target"].name,
                        "request": json.dumps(res["request"]),
                        "status": res["status"],
                        "size_mb": round(size, 3),
                        "seconds": round(res["seconds"], 1),
                        "attempts": res["attempts"],
                        "error": res["error"],
                    }
                )
            log.info(
                "[%d/%d] %s %s: %s (%.1f MB).",
                k,
                len(list_job),
                res["source"],
                res["target"].name,
                res["status"],
                size,
            )
    df_log = pd.DataFrame(list_log)
    log_path = raw_dir / "download_log.csv"
    df_log.to_csv(log_path, mode="a", header=not log_path.exists(), index=False)

    # --- Static upstream area files: copied next to the discharge with their source recorded ---
    sources = {}
    for src in dis["sources"]:
        path = dis[src]["upstream_area_file"]
        if not path:
            log.warning(
                "%s: upstream_area_file not set; discharge_prepare.py needs it.", src
            )
            continue
        target = raw_dir / src / f"upstream_area{Path(path).suffix}"
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(path, target)
        sources[src] = {
            "upstream_area_file": str(path),
            "copied_to": str(target),
            "dataset_id": dis[src]["dataset_id"],
        }
    (raw_dir / "discharge_sources.json").write_text(json.dumps(sources, indent=1))
    failed = df_log[df_log["status"] == "failed"]
    log.info(
        "Discharge download: %d requests, %d downloaded, %d complete and skipped, %d failed; %.1f MB; %.1f s.",
        len(df_log),
        int((df_log["status"] == "downloaded").sum()),
        int((df_log["status"] == "complete, skipped").sum()),
        len(failed),
        df_log["size_mb"].sum(),
        time.perf_counter() - start,
    )
    for src in sorted(set(failed["source"])):
        log.error(
            "%s: %d requests failed (first error: %s). %s",
            src,
            int((failed["source"] == src).sum()),
            failed[failed["source"] == src]["error"].iloc[0],
            "The workflow continues with GloFAS only." if src == "efas" else "",
        )
    if list_skip_src or len(failed):
        print(f"Skipped sources: {list_skip_src}; failed requests: {len(failed)}.")
