# Author: Jakub Zapletal
# Date: 2026-11-04
# Downloads the GHSL built-up volume tiles (total and non-residential, 100 m) intersecting the DTM footprint for the
# configured or latest available epoch and mosaics each product as a VRT.

import datetime
import json
import urllib.error
import urllib.request
from pathlib import Path

import geopandas as gpd
import shapely
from osgeo import gdal
from shapely.geometry import box

import gis

# --- Run variables (no command-line parameters; DEFENSES_<NAME> environment variables override them) ---
CONFIG_PATH = "config/defenses.toml"
COUNTRY = "XX"

KEYS = [
    "general.work_dir",
    "general.output_dir",
    "general.log_level",
    "sources.ghsl_url_templates",
    "sources.ghsl_tile_schema",
    "sources.ghsl_tile_field",
    "sources.ghsl_epoch",
    "sources.request_pause_s",
    "sources.max_retries",
]
PRODUCTS = ["total", "non_residential"]
FIRST_EPOCH = 1975  # first GHSL epoch; epochs follow in 5-year steps

if __name__ == "__main__":
    run = gis.Run.env_overrides(CONFIG_PATH=CONFIG_PATH, COUNTRY=COUNTRY)
    cfg = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], KEYS)
    general, sources = cfg["general"], cfg["sources"]
    work_dir = Path(general["work_dir"]) / run["COUNTRY"] / "ghsl"
    out_dir = Path(general["output_dir"]) / run["COUNTRY"]
    log = gis.Run.logger("download_ghsl", out_dir / "logs", general["log_level"])
    gis.Config.warn_placeholders(cfg)
    templates = sources["ghsl_url_templates"]
    missing = [p for p in PRODUCTS if not templates.get(p)]
    if missing or not sources["ghsl_tile_schema"] or not sources["ghsl_tile_field"]:
        raise SystemExit(
            f"GHSL sources are not configured (empty: {missing or ''} {'ghsl_tile_schema' if not sources['ghsl_tile_schema'] else ''} "
            f"{'ghsl_tile_field' if not sources['ghsl_tile_field'] else ''}). The tile URLs and the tiling schema have not "
            "been verified against the JRC GHSL documentation; fill them in [sources]."
        )
    footprint_path = out_dir / "dtm_footprint.gpkg"
    if not footprint_path.exists():
        raise SystemExit("dtm_footprint.gpkg is missing; run prepare_dtm.py first.")
    footprint = gpd.read_file(footprint_path)

    # --- Tiling schema: tiles intersecting the DTM footprint ---
    schema_source = sources["ghsl_tile_schema"]
    if "://" in schema_source:
        schema_source = gis.Http.download(
            schema_source,
            work_dir / Path(schema_source).name,
            sources["max_retries"],
            sources["request_pause_s"],
        )
    if str(schema_source).endswith(".zip"):
        shapefiles = [
            n
            for n in gdal.ReadDir(f"/vsizip/{schema_source}") or []
            if n.endswith(".shp")
        ]
        if len(shapefiles) != 1:
            raise SystemExit(
                f"Expected one shapefile in {schema_source}, found {shapefiles}."
            )
        schema_source = f"/vsizip/{schema_source}/{shapefiles[0]}"
    gis.Vector.check_single_layer(schema_source)
    schema = gpd.read_file(schema_source)
    if sources["ghsl_tile_field"] not in schema.columns:
        raise SystemExit(
            f"Tile field {sources['ghsl_tile_field']!r} not in the tiling schema; columns: {list(schema.columns)}"
        )
    area = gpd.GeoSeries([box(*footprint.total_bounds)], crs=footprint.crs).to_crs(
        schema.crs
    )
    area = shapely.segmentize(area.iloc[0], max(area.iloc[0].length / 400, 1.0))
    tiles = sorted(
        schema.loc[schema.intersects(area), sources["ghsl_tile_field"]].astype(str)
    )
    if not tiles:
        raise SystemExit("No GHSL tile intersects the DTM footprint.")
    log.info(
        "GHSL tiles intersecting the DTM footprint: %d (%s).",
        len(tiles),
        ", ".join(tiles),
    )

    # --- Epoch: configured, otherwise the latest one the server provides (not after the current year) ---
    epoch = sources["ghsl_epoch"]
    if not epoch:
        year = datetime.date.today().year
        for candidate in range(year - year % 5, FIRST_EPOCH - 1, -5):
            url = templates["total"].format(epoch=candidate, tile=tiles[0])
            try:
                with urllib.request.urlopen(
                    urllib.request.Request(url, method="HEAD"), timeout=60
                ):
                    epoch = candidate
                    break
            except (urllib.error.URLError, OSError):
                continue
        if not epoch:
            raise SystemExit(
                "No GHSL epoch found on the server for the configured templates."
            )
        log.info("Latest available GHSL epoch: %d.", epoch)
    if epoch > datetime.date.today().year - 5:
        log.warning(
            "GHSL epoch %d is modelled, not observed (see the GHSL documentation).",
            epoch,
        )

    # --- Tiles and VRT mosaics ---
    record = {"epoch": epoch, "tiles": tiles, "products": {}}
    for product in PRODUCTS:
        paths = []
        for tile in tiles:
            url = templates[product].format(epoch=epoch, tile=tile)
            target = work_dir / product / url.rsplit("/", 1)[-1]
            gis.Http.download(
                url, target, sources["max_retries"], sources["request_pause_s"]
            )
            if target.suffix == ".zip":
                rasters = [
                    n
                    for n in gdal.ReadDir(f"/vsizip/{target}") or []
                    if n.lower().endswith(".tif")
                ]
                if len(rasters) != 1:
                    raise SystemExit(
                        f"Expected one GeoTIFF in {target}, found {rasters}."
                    )
                paths.append(f"/vsizip/{target}/{rasters[0]}")
            else:
                paths.append(str(target))
        vrt = gis.Raster.build_vrt(paths, out_dir / f"ghsl_{product}.vrt")
        record["products"][product] = {
            "vrt": vrt,
            "sources": paths,
            "url_template": templates[product],
        }
        info = gis.Raster.inspect(vrt)
        log.info(
            "GHSL %s: %d tiles, %d x %d cells of %.0f m -> %s",
            product,
            len(paths),
            info["width"],
            info["height"],
            info["pixel_x"],
            vrt,
        )
    (out_dir / "ghsl_sources.json").write_text(json.dumps(record, indent=2))
