# Author: Jakub Zapletal
# Date: 2026-11-04
# Decides which river features are water sources (closed by a gate at a levee crossing or by a pumping station at
# their downstream end means not) and selects the pumping stations next to levees or rivers, from OSM structures.

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from osgeo import gdal

import gis

# --- Run variables (no command-line parameters; DEFENSES_<NAME> environment variables override them) ---
CONFIG_PATH = "config/defenses.toml"
COUNTRY = "XX"

KEYS = [
    "general.work_dir",
    "general.output_dir",
    "general.log_level",
    "country.levees",
    "country.rivers",
    "country.endpoint_elevation_radius_m",
    "country.compartments.gap_tolerance_m",
    "country.rp.levee_adjacent_distance_m",
    "sources.osm_tags_hydraulic",
    "sources.osm_tags_closing",
    "sources.osm_tags_pumping",
]

if __name__ == "__main__":
    run = gis.Run.env_overrides(CONFIG_PATH=CONFIG_PATH, COUNTRY=COUNTRY)
    cfg = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], KEYS)
    general, country, sources = cfg["general"], cfg["country"], cfg["sources"]
    work_dir = Path(general["work_dir"]) / run["COUNTRY"]
    out_dir = Path(general["output_dir"]) / run["COUNTRY"]
    log = gis.Run.logger("osm_hydraulics", out_dir / "logs", general["log_level"])
    gis.Config.warn_placeholders(cfg)
    vrt_path = work_dir / "dtm" / "dtm.vrt"
    points_path = out_dir / "osm_hydraulic_points.gpkg"
    if not vrt_path.exists() or not points_path.exists():
        raise SystemExit(
            "dtm.vrt or osm_hydraulic_points.gpkg is missing; run prepare_dtm.py and download_osm.py first."
        )
    dtm = gdal.Open(str(vrt_path))
    wkt = dtm.GetProjection()
    gap = country["compartments"]["gap_tolerance_m"]
    near = country["rp"]["levee_adjacent_distance_m"]

    rivers = gis.Vector.read_merged(country["rivers"], wkt, keep_columns=True)
    rivers.insert(0, "river_id", np.arange(len(rivers)))
    levees = gis.Vector.read_merged(country["levees"], wkt)
    gis.Vector.check_single_layer(points_path)
    structures = gpd.read_file(points_path).to_crs(wkt)

    def matching(frame, table):
        hit = np.zeros(len(frame), bool)
        for key, values in table.items():
            if key in frame.columns:
                hit |= frame[key].isin(values).to_numpy()
        return hit

    def describe(frame):
        keys = [k for k in sources["osm_tags_hydraulic"] if k in frame.columns]
        return frame[keys].apply(
            lambda r: ", ".join(f"{k}={v}" for k, v in r.items() if pd.notna(v)), axis=1
        )

    closing = structures[matching(structures, sources["osm_tags_closing"])].reset_index(
        drop=True
    )
    pumping = structures[matching(structures, sources["osm_tags_pumping"])].reset_index(
        drop=True
    )
    closing["tags"] = describe(closing) if len(closing) else []
    log.info(
        "OSM structures: %d in total, %d closing structures, %d pumping stations.",
        len(structures),
        len(closing),
        len(pumping),
    )

    reasons = pd.Series("", index=rivers.index, dtype=object)

    # --- Closing structures at river-levee crossings ---
    if len(levees) and len(closing):
        river_idx, levee_idx = levees.sindex.query(
            rivers.geometry.values, predicate="intersects"
        )
        crossings = shapely.intersection(
            rivers.geometry.values[river_idx], levees.geometry.values[levee_idx]
        )
        parts = shapely.get_parts(crossings, return_index=True)
        points = shapely.point_on_surface(parts[0])
        owner = river_idx[parts[1]]
        cross_idx, structure_idx = closing.sindex.query(
            points, predicate="dwithin", distance=gap
        )
        if len(cross_idx):
            hits = (
                pd.DataFrame(
                    {
                        "river": owner[cross_idx],
                        "tags": closing["tags"].to_numpy()[structure_idx],
                        "distance": shapely.distance(
                            points[cross_idx], closing.geometry.values[structure_idx]
                        ),
                    }
                )
                .sort_values("distance")
                .drop_duplicates("river")
            )
            for river, tags, distance in hits.itertuples(index=False):
                reasons[river] = (
                    f"closing structure ({tags}) {distance:.1f} m from a levee crossing"
                )
        log.info(
            "River-levee crossings: %d; closed by a nearby closing structure: %d river features.",
            len(points),
            int((reasons != "").sum()),
        )

    # --- Downstream ends from DTM elevations at both line ends ---
    coords = [shapely.get_coordinates(g) for g in rivers.geometry]
    first = np.array([c[0] for c in coords])
    last = np.array([c[-1] for c in coords])
    radius = country["endpoint_elevation_radius_m"]
    elevation = np.full((len(rivers), 2), np.nan)
    for k, (start, end) in enumerate(zip(first, last)):
        for j, (x, y) in enumerate((start, end)):
            window, _ = gis.Raster.read_native(
                dtm, (x - radius, y - radius, x + radius, y + radius)
            )
            if window.size and np.isfinite(window).any():
                elevation[k, j] = np.nanmedian(window)
        if (k + 1) % 5000 == 0:
            log.info(
                "  endpoint elevations: %d of %d river features", k + 1, len(rivers)
            )
    by_dtm = np.isfinite(elevation).all(axis=1) & (
        np.abs(elevation[:, 0] - elevation[:, 1]) > 1e-3
    )
    start_lower = by_dtm & (elevation[:, 0] < elevation[:, 1])
    downstream = np.where(start_lower[:, None], first, last)
    rivers["downstream_x"], rivers["downstream_y"] = downstream[:, 0], downstream[:, 1]
    rivers["downstream_method"] = np.where(by_dtm, "dtm elevation", "line direction")
    log.info(
        "Downstream ends: %d from DTM elevations, %d from the line direction (equal or missing elevations).",
        int(by_dtm.sum()),
        int((~by_dtm).sum()),
    )

    # --- Pumping stations at the downstream end ---
    if len(pumping):
        end_idx, pump_idx = pumping.sindex.query(
            shapely.points(downstream), predicate="dwithin", distance=gap
        )
        for river, pump in zip(end_idx, pump_idx):
            if not reasons[river]:
                distance = shapely.distance(
                    shapely.points(downstream[river]), pumping.geometry.values[pump]
                )
                reasons[river] = (
                    f"pumping station {distance:.1f} m from the downstream end"
                )

    rivers["is_source"] = reasons == ""
    rivers["reason"] = reasons
    gis.Vector.write_gpkg(rivers, out_dir / "rivers_sources.gpkg")
    log.info(
        "rivers_sources.gpkg: %d river features, %d closed during floods (not water sources).",
        len(rivers),
        int((~rivers.is_source).sum()),
    )

    # --- Pumping stations next to levees or rivers ---
    near_levee = np.zeros(len(pumping), bool)
    near_river = np.zeros(len(pumping), bool)
    if len(pumping) and len(levees):
        near_levee[
            np.unique(
                levees.sindex.query(
                    pumping.geometry.values, predicate="dwithin", distance=near
                )[0]
            )
        ] = True
    if len(pumping) and len(rivers):
        near_river[
            np.unique(
                rivers.sindex.query(
                    pumping.geometry.values, predicate="dwithin", distance=near
                )[0]
            )
        ] = True
    selected = pumping.assign(near_levee=near_levee, near_river=near_river)[
        near_levee | near_river
    ]
    gis.Vector.write_gpkg(selected.reset_index(drop=True), out_dir / "osm_pumping.gpkg")
    log.info(
        "osm_pumping.gpkg: %d of %d pumping stations within %.0f m of a levee or river.",
        len(selected),
        len(pumping),
        near,
    )
