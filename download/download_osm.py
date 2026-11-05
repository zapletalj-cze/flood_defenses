# Author: Jakub Zapletal
# Date: 2026-11-04
# Downloads the country OpenStreetMap extract (Geofabrik .osm.pbf) and extracts hydraulic structures (as points),
# drains, buildings, water polygons, road and railway crossings, flood walls and industrial areas inside the DTM
# footprint into single-layer GeoPackages with GDAL; logs counts of candidate wall tags.

from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

import geopandas as gpd
from osgeo import gdal
from shapely.geometry import box

import gis

# --- Run variables (no command-line parameters; DEFENSES_<NAME> environment variables override them) ---
CONFIG_PATH = "config/defenses.toml"
COUNTRY = "XX"
FORCE = False  # extract again although the outputs exist

KEYS = [
    "general.work_dir",
    "general.output_dir",
    "general.gdal_threads_per_worker",
    "general.log_level",
    "country.geofabrik_region",
    "country.osm_pbf",
    "sources.geofabrik_url_template",
    "sources.osm_tags_hydraulic",
    "sources.osm_tags_drains",
    "sources.osm_tags_water",
    "country.building_footprints",
    "country.water_polygons",
    "country.context.osm_industrial",
    "barriers.osm_tags_walls",
    "barriers.osm_tags_crossings",
    "barriers.osm_wall_candidates",
    "sources.request_pause_s",
    "sources.max_retries",
]
OSM_CONFIG = Path(__file__).resolve().parents[1] / "config" / "osmconf.ini"
# Columns the custom osmconf.ini exposes per OSM layer (checked against the tag lists of the config).
EXPOSED = {
    "points": [
        "name",
        "man_made",
        "waterway",
        "lock",
        "building",
        "barrier",
        "highway",
    ],
    "lines": [
        "name",
        "man_made",
        "waterway",
        "building",
        "intermittent",
        "highway",
        "railway",
        "barrier",
        "wall",
    ],
    "multipolygons": [
        "name",
        "man_made",
        "waterway",
        "building",
        "landuse",
        "natural",
        "water",
        "barrier",
        "wall",
    ],
}
# Layer each tag table is matched in (the hydraulic table is matched in all three layers).
LAYER_OF = {
    "drains": "lines",
    "water": "multipolygons",
    "crossings": "lines",
    "walls": "lines",
    "wall_candidates": "lines",
    "industrial": "multipolygons",
}

if __name__ == "__main__":
    run = gis.Run.env_overrides(CONFIG_PATH=CONFIG_PATH, COUNTRY=COUNTRY, FORCE=FORCE)
    cfg = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], KEYS)
    general, country, sources = cfg["general"], cfg["country"], cfg["sources"]
    work_dir = Path(general["work_dir"]) / run["COUNTRY"]
    out_dir = Path(general["output_dir"]) / run["COUNTRY"]
    log = gis.Run.logger("download_osm", out_dir / "logs", general["log_level"])
    gis.Config.warn_placeholders(cfg)

    tags = {
        "hydraulic": sources["osm_tags_hydraulic"],
        "drains": sources["osm_tags_drains"],
        "water": sources["osm_tags_water"],
        "crossings": cfg["barriers"]["osm_tags_crossings"],
        "walls": cfg["barriers"]["osm_tags_walls"],
        "wall_candidates": cfg["barriers"]["osm_wall_candidates"],
        "industrial": country["context"]["osm_industrial"],
    }
    for name, table in tags.items():
        for key, list_val in table.items():
            if not isinstance(list_val, list) or not all(
                isinstance(v, str) for v in list_val
            ):
                raise SystemExit(
                    f"sources.osm_tags_{name}.{key} must be a list of tag values."
                )
    list_miss = sorted(
        {
            k
            for k in tags["hydraulic"]
            if any(k not in cols for cols in EXPOSED.values())
        }
        | {
            k
            for name, layer in LAYER_OF.items()
            for k in tags[name]
            if k not in EXPOSED[layer]
        }
    )
    if list_miss:
        raise SystemExit(
            f"OSM keys {list_miss} are not exposed as columns by {OSM_CONFIG}; add them to the "
            "attributes of the [points], [lines] and [multipolygons] sections and to EXPOSED."
        )

    # --- Country extract ---
    if country["osm_pbf"]:
        pbf = Path(country["osm_pbf"])
        log.info("Using the local OSM extract %s.", pbf)
    else:
        template = sources["geofabrik_url_template"]
        if not template or not country["geofabrik_region"]:
            raise SystemExit(
                "sources.geofabrik_url_template or countries.<code>.geofabrik_region is empty. The URL "
                "pattern has not been verified against Geofabrik's documentation; fill it in (e.g. with "
                "a {region} placeholder) or set countries.<code>.osm_pbf to a downloaded extract."
            )
        url = template.format(region=country["geofabrik_region"])
        pbf = work_dir / "osm" / url.rsplit("/", 1)[-1]
        gis.Http.download(url, pbf, sources["max_retries"], sources["request_pause_s"])

    # --- Extraction inside the DTM footprint (one VectorTranslate call per output) ---
    fp_path = out_dir / "dtm_footprint.gpkg"
    vrt_path = work_dir / "dtm" / "dtm.vrt"
    if not fp_path.exists() or not vrt_path.exists():
        raise SystemExit(
            "dtm_footprint.gpkg or dtm.vrt is missing; run prepare_dtm.py first."
        )
    wkt = gdal.Open(str(vrt_path)).GetProjection()
    # OSM data are in WGS 84; GDAL accepts no spatial-filter CRS together with SQL, so the footprint's
    # bounding box (densified before reprojection) is given in WGS 84.
    extent = gpd.GeoSeries(
        [box(*gpd.read_file(fp_path).to_crs(wkt).total_bounds)], crs=wkt
    )
    list_bnd = list(
        extent.segmentize(extent.iloc[0].length / 400).to_crs(4326).total_bounds
    )
    gdal.SetConfigOption("OSM_CONFIG_FILE", str(OSM_CONFIG))
    # The OSM driver keeps a temporary index; place it on the work disk.
    gdal.SetConfigOption("CPL_TMPDIR", str(work_dir / "osm"))
    (work_dir / "osm").mkdir(parents=True, exist_ok=True)

    def condition(table, layer):
        list_part = [
            f"{key} IN ({', '.join(repr(v) for v in values)})"
            for key, values in table.items()
            if key in EXPOSED[layer] and values
        ]
        return " OR ".join(list_part) if list_part else "0"

    list_col = sorted(set(tags["hydraulic"]) | {"name"})
    select = ", ".join(list_col)
    outputs = {
        "osm_hydraulic_points": (
            f"SELECT osm_id AS osm_id, 'node' AS osm_element, NULL AS osm_way_id, {select} FROM points "
            f"WHERE {condition(tags['hydraulic'], 'points')} "
            f"UNION ALL SELECT osm_id AS osm_id, 'way' AS osm_element, NULL AS osm_way_id, {select} FROM lines "
            f"WHERE {condition(tags['hydraulic'], 'lines')} "
            f"UNION ALL SELECT osm_id AS osm_id, 'area' AS osm_element, "
            f"osm_way_id AS osm_way_id, {select} FROM multipolygons "
            f"WHERE {condition(tags['hydraulic'], 'multipolygons')}",
            "GEOMETRY",
        ),
        "osm_drains": (
            f"SELECT osm_id, {', '.join(sorted(set(tags['drains']) | {'name'}))} FROM lines "
            f"WHERE {condition(tags['drains'], 'lines')}",
            "LINESTRING",
        ),
        # Stage D: crossings at gaps of protection lines, industrial areas for the context zones.
        "osm_crossings": (
            f"SELECT osm_id, {', '.join(sorted(set(tags['crossings']) | {'name'}))} FROM lines "
            f"WHERE {condition(tags['crossings'], 'lines')}",
            "LINESTRING",
        ),
        "osm_industrial": (
            f"SELECT osm_id, osm_way_id, {', '.join(sorted(set(tags['industrial']) | {'name'}))} "
            f"FROM multipolygons "
            f"WHERE {condition(tags['industrial'], 'multipolygons')}",
            "MULTIPOLYGON",
        ),
        # Counts only: lines with any candidate wall tag (see barriers.osm_wall_candidates).
        "osm_wall_candidates": (
            f"SELECT osm_id, {', '.join(sorted(set(tags['wall_candidates']) | {'name'}))} "
            f"FROM lines WHERE {condition(tags['wall_candidates'], 'lines')}",
            "LINESTRING",
        ),
    }
    # Always extracted: Stage E combines OSM buildings with national footprints; man_made serves its exclusions.
    outputs["osm_buildings"] = (
        "SELECT osm_id, osm_way_id, building, man_made FROM multipolygons "
        "WHERE building IS NOT NULL",
        "MULTIPOLYGON",
    )
    if not country["water_polygons"]:
        outputs["osm_water"] = (
            f"SELECT osm_id, osm_way_id, {', '.join(sorted(set(tags['water']) | {'name'}))} "
            f"FROM multipolygons "
            f"WHERE {condition(tags['water'], 'multipolygons')}",
            "MULTIPOLYGON",
        )
    if tags["walls"]:
        outputs["osm_walls"] = (
            f"SELECT osm_id, {', '.join(sorted(set(tags['walls']) | {'name'}))} FROM lines "
            f"WHERE {condition(tags['walls'], 'lines')}",
            "LINESTRING",
        )
    else:
        log.info(
            "barriers.osm_tags_walls is empty: no osm_walls.gpkg; choose the tags from the candidate counts below."
        )
    list_extract = []
    for name, (sql, geometry_type) in outputs.items():
        path = (
            work_dir / "osm" / f"{name}.gpkg"
            if name == "osm_wall_candidates"
            else out_dir / f"{name}.gpkg"
        )
        if path.exists() and not run["FORCE"]:
            try:
                gis.Vector.check_single_layer(path)
            except (RuntimeError, SystemExit):
                log.warning("%s is incomplete; re-extracting.", path.name)
                path.unlink()
        if path.exists() and not run["FORCE"]:
            log.info("%s exists; extraction skipped (set FORCE to redo).", path.name)
        else:
            list_extract.append((path, name, sql, geometry_type))

    def translate_output(task):
        path, name, sql, geometry_type = task
        path.unlink(missing_ok=True)
        source = None
        try:
            source = gdal.OpenEx(
                str(pbf),
                gdal.OF_VECTOR | gdal.OF_READONLY,
                open_options=["INTERLEAVED_READING=YES"],
            )
            if source is None:
                raise RuntimeError(f"GDAL could not open {pbf} as a vector dataset.")
            translated = gdal.VectorTranslate(
                str(path),
                source,
                format="GPKG",
                SQLStatement=sql,
                SQLDialect="OGRSQL",
                layerName=name,
                spatFilter=list_bnd,
                dstSRS=wkt,
                geometryType=geometry_type,
                datasetCreationOptions=["ADD_GPKG_OGR_CONTENTS=NO"],
            )
            if translated is None:
                raise RuntimeError(f"GDAL did not create {path.name}.")
            translated = None
            source = None
            gis.Vector.check_single_layer(path)
        except Exception:
            source = None
            path.unlink(missing_ok=True)
            raise

    if list_extract:
        n_threads = max(
            1, min(general["gdal_threads_per_worker"], len(list_extract))
        )
        log.info(
            "Extracting %d OSM layers with %d parallel workers.",
            len(list_extract),
            n_threads,
        )
        with ThreadPoolExecutor(max_workers=n_threads) as executor:
            list(executor.map(translate_output, list_extract))

    for name in outputs:
        path = (
            work_dir / "osm" / f"{name}.gpkg"
            if name == "osm_wall_candidates"
            else out_dir / f"{name}.gpkg"
        )
        gis.Vector.check_single_layer(path)
        gdf_ext = gpd.read_file(path)
        if name == "osm_hydraulic_points":
            gdf_ext["osm_id"] = gdf_ext["osm_id"].fillna(gdf_ext["osm_way_id"])
            gdf_ext = gdf_ext.drop(columns="osm_way_id")
            gdf_ext = gdf_ext.set_geometry(gdf_ext.geometry.representative_point())
            gdf_ext.to_file(path, driver="GPKG")
        log.info("%s: %d features.", path.name, len(gdf_ext))
        # Counts per tag value, to adjust the tag lists in the configuration.
        table = {
            "osm_hydraulic_points": tags["hydraulic"],
            "osm_drains": tags["drains"],
            "osm_crossings": tags["crossings"],
            "osm_industrial": tags["industrial"],
            "osm_wall_candidates": tags["wall_candidates"],
            "osm_water": tags["water"],
            "osm_walls": tags["walls"],
        }.get(name, {})
        for key in table:
            if key in gdf_ext.columns:
                counts = gdf_ext[key].dropna().value_counts()
                for value, count in counts.items():
                    log.info("  %s=%s: %d", key, value, count)
