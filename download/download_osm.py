# Author: Jakub Zapletal
# Date: 2026-11-04
# Downloads the country OpenStreetMap extract (Geofabrik .osm.pbf) and extracts hydraulic structures (as points)
# and drain lines inside the DTM footprint into single-layer GeoPackages with GDAL.

from pathlib import Path

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
    "general.log_level",
    "country.geofabrik_region",
    "country.osm_pbf",
    "sources.geofabrik_url_template",
    "sources.osm_tags_hydraulic",
    "sources.osm_tags_drains",
    "sources.request_pause_s",
    "sources.max_retries",
]
OSM_CONFIG = Path(__file__).resolve().parents[1] / "config" / "osmconf.ini"
# Columns the custom osmconf.ini exposes per OSM layer (checked against the tag lists of the config).
EXPOSED = {
    "points": ["name", "man_made", "waterway", "lock", "building"],
    "lines": ["name", "man_made", "waterway", "building", "intermittent"],
    "multipolygons": ["name", "man_made", "waterway", "building"],
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
    }
    for name, table in tags.items():
        for key, values in table.items():
            if not isinstance(values, list) or not all(
                isinstance(v, str) for v in values
            ):
                raise SystemExit(
                    f"sources.osm_tags_{name}.{key} must be a list of tag values."
                )
    missing = sorted(
        {
            k
            for k in tags["hydraulic"]
            if any(k not in cols for cols in EXPOSED.values())
        }
        | {k for k in tags["drains"] if k not in EXPOSED["lines"]}
    )
    if missing:
        raise SystemExit(
            f"OSM keys {missing} are not exposed as columns by {OSM_CONFIG}; add them to the "
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
    footprint_path = out_dir / "dtm_footprint.gpkg"
    vrt_path = work_dir / "dtm" / "dtm.vrt"
    if not footprint_path.exists() or not vrt_path.exists():
        raise SystemExit(
            "dtm_footprint.gpkg or dtm.vrt is missing; run prepare_dtm.py first."
        )
    wkt = gdal.Open(str(vrt_path)).GetProjection()
    # OSM data are in WGS 84; GDAL accepts no spatial-filter CRS together with SQL, so the footprint's
    # bounding box (densified before reprojection) is given in WGS 84.
    extent = gpd.GeoSeries(
        [box(*gpd.read_file(footprint_path).to_crs(wkt).total_bounds)], crs=wkt
    )
    bounds = list(
        extent.segmentize(extent.iloc[0].length / 400).to_crs(4326).total_bounds
    )
    gdal.SetConfigOption("OSM_CONFIG_FILE", str(OSM_CONFIG))
    # The OSM driver keeps a temporary index; place it on the work disk.
    gdal.SetConfigOption("CPL_TMPDIR", str(work_dir / "osm"))
    (work_dir / "osm").mkdir(parents=True, exist_ok=True)

    def condition(table, layer):
        parts = [
            f"{key} IN ({', '.join(repr(v) for v in values)})"
            for key, values in table.items()
            if key in EXPOSED[layer] and values
        ]
        return " OR ".join(parts) if parts else "0"

    columns = sorted(set(tags["hydraulic"]) | {"name"})
    select = ", ".join(columns)
    outputs = {
        "osm_hydraulic_points": (
            f"SELECT osm_id, 'node' AS osm_element, {select}, geometry FROM points "
            f"WHERE {condition(tags['hydraulic'], 'points')} "
            f"UNION ALL SELECT osm_id, 'way', {select}, ST_PointOnSurface(geometry) FROM lines "
            f"WHERE {condition(tags['hydraulic'], 'lines')} "
            f"UNION ALL SELECT COALESCE(osm_id, osm_way_id), 'area', {select}, "
            f"ST_PointOnSurface(geometry) FROM multipolygons "
            f"WHERE {condition(tags['hydraulic'], 'multipolygons')}",
            "POINT",
        ),
        "osm_drains": (
            f"SELECT osm_id, {', '.join(sorted(set(tags['drains']) | {'name'}))}, geometry FROM lines "
            f"WHERE {condition(tags['drains'], 'lines')}",
            "LINESTRING",
        ),
    }
    for name, (sql, geometry_type) in outputs.items():
        path = out_dir / f"{name}.gpkg"
        if path.exists() and not run["FORCE"]:
            log.info("%s exists; extraction skipped (set FORCE to redo).", path.name)
        else:
            path.unlink(missing_ok=True)
            gdal.VectorTranslate(
                str(path),
                str(pbf),
                format="GPKG",
                SQLStatement=sql,
                SQLDialect="SQLITE",
                layerName=name,
                spatFilter=bounds,
                dstSRS=wkt,
                geometryType=geometry_type,
            )
        gis.Vector.check_single_layer(path)
        extracted = gpd.read_file(path)
        log.info("%s: %d features.", path.name, len(extracted))
        # Counts per tag value, to adjust the tag lists in the configuration.
        table = tags["hydraulic"] if name == "osm_hydraulic_points" else tags["drains"]
        for key in table:
            if key in extracted.columns:
                counts = extracted[key].dropna().value_counts()
                for value, count in counts.items():
                    log.info("  %s=%s: %d", key, value, count)
