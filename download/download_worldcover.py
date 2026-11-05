# Author: Jakub Zapletal
# Date: 2026-11-04
# Downloads the ESA WorldCover 10 m map tiles that intersect the river corridor (official 3x3 degree tile grid) and
# mosaics them as worldcover.vrt for HAND_YZ.py; records version, tiles and licence in worldcover_sources.json.

import json
from pathlib import Path

import geopandas as gpd
import shapely
from osgeo import gdal

import gis

# --- Run variables (no command-line parameters; DEFENSES_<NAME> environment variables override them) ---
CONFIG_PATH = "config/defenses.toml"
COUNTRY = "NOR"

KEYS = [
    "general.work_dir",
    "general.output_dir",
    "general.log_level",
    "sources.worldcover_url_template",
    "sources.worldcover_grid_url",
    "sources.worldcover_tile_field",
    "sources.worldcover_version",
    "sources.worldcover_year",
    "sources.request_pause_s",
    "sources.max_retries",
]

if __name__ == "__main__":
    run = gis.Run.env_overrides(CONFIG_PATH=CONFIG_PATH, COUNTRY=COUNTRY)
    cfg = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], KEYS)
    general, sources = cfg["general"], cfg["sources"]
    work_dir = Path(general["work_dir"]) / run["COUNTRY"] / "worldcover"
    out_dir = Path(general["output_dir"]) / run["COUNTRY"]
    log = gis.Run.logger("download_worldcover", out_dir / "logs", general["log_level"])
    gis.Config.warn_placeholders(cfg)
    if not sources["worldcover_url_template"] or not sources["worldcover_grid_url"]:
        raise SystemExit(
            "sources.worldcover_url_template or worldcover_grid_url is empty; set them or list local tiles "
            "in countries.<code>.lulc.worldcover."
        )
    # The corridor holds every compartment; the DTM footprint is the fallback before build_units.py ran.
    area_path = next(
        (
            p
            for p in [out_dir / "corridor.gpkg", out_dir / "dtm_footprint.gpkg"]
            if p.exists()
        ),
        None,
    )
    if area_path is None:
        raise SystemExit(
            "corridor.gpkg and dtm_footprint.gpkg are missing; run prepare_dtm.py and build_units.py first."
        )
    gis.Vector.check_single_layer(area_path)
    gdf_area = gpd.read_file(area_path)
    area = gpd.GeoSeries(
        [shapely.segmentize(gdf_area.union_all(), 1000.0)], crs=gdf_area.crs
    ).to_crs(4326)

    # --- Tiles of the official grid that intersect the area ---
    pause, retries = sources["request_pause_s"], sources["max_retries"]
    grid_path = gis.Http.download(
        sources["worldcover_grid_url"],
        work_dir / Path(sources["worldcover_grid_url"]).name,
        retries,
        pause,
    )
    gis.Vector.check_single_layer(grid_path)
    gdf_grid = gpd.read_file(grid_path).to_crs(4326)
    field = sources["worldcover_tile_field"]
    if field not in gdf_grid.columns:
        raise SystemExit(
            f"Field {field!r} not found in the WorldCover grid; columns: {list(gdf_grid.columns)}"
        )
    list_tile = sorted(
        gdf_grid.loc[gdf_grid.intersects(area.iloc[0]), field].astype(str)
    )
    if not list_tile:
        raise SystemExit("No WorldCover tile intersects the corridor.")
    log.info(
        "WorldCover %s (%d): %d tiles intersect %s: %s.",
        sources["worldcover_version"],
        sources["worldcover_year"],
        len(list_tile),
        area_path.name,
        ", ".join(list_tile),
    )

    # --- Download (skipped when the local file has the size the server reports), mosaic ---
    list_path, list_url = [], []
    for tile in list_tile:
        url = sources["worldcover_url_template"].format(
            version=sources["worldcover_version"],
            year=sources["worldcover_year"],
            tile=tile,
        )
        list_path.append(
            gis.Http.download(url, work_dir / Path(url).name, retries, pause)
        )
        list_url.append(url)
    vrt = gis.Raster.build_vrt(list_path, out_dir / "worldcover.vrt")
    metadata = gdal.Open(list_path[0]).GetMetadata()
    record = {
        "version": sources["worldcover_version"],
        "year": sources["worldcover_year"],
        "tiles": list_tile,
        "urls": list_url,
        "files": list_path,
        "licence": metadata.get("license", ""),
        "attribution": metadata.get("copyright", ""),
        "legend": metadata.get("legend", ""),
    }
    (out_dir / "worldcover_sources.json").write_text(json.dumps(record, indent=2))
    log.info(
        "worldcover.vrt over %d tiles; licence: %s.",
        len(list_path),
        record["licence"] or "not stated in the file",
    )
