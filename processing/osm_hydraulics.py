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

    gdf_rn = gis.Vector.read_merged(country["rivers"], wkt, keep_columns=True)
    gdf_rn.insert(0, "river_id", np.arange(len(gdf_rn)))
    gdf_lv = gis.Vector.read_merged(country["levees"], wkt)
    gis.Vector.check_single_layer(points_path)
    gdf_st = gpd.read_file(points_path).to_crs(wkt)

    def matching(frame, table):
        hit = np.zeros(len(frame), bool)
        for key, list_val in table.items():
            if key in frame.columns:
                hit |= frame[key].isin(list_val).to_numpy()
        return hit

    def describe(frame):
        list_key = [k for k in sources["osm_tags_hydraulic"] if k in frame.columns]
        return frame[list_key].apply(
            lambda r: ", ".join(f"{k}={v}" for k, v in r.items() if pd.notna(v)), axis=1
        )

    gdf_close = gdf_st[matching(gdf_st, sources["osm_tags_closing"])].reset_index(
        drop=True
    )
    gdf_pump = gdf_st[matching(gdf_st, sources["osm_tags_pumping"])].reset_index(
        drop=True
    )
    gdf_close["tags"] = describe(gdf_close) if len(gdf_close) else []
    log.info(
        "OSM structures: %d in total, %d closing structures, %d pumping stations.",
        len(gdf_st),
        len(gdf_close),
        len(gdf_pump),
    )

    reasons = pd.Series("", index=gdf_rn.index, dtype=object)

    # --- Closing structures at river-levee crossings ---
    if len(gdf_lv) and len(gdf_close):
        river_idx, levee_idx = gdf_lv.sindex.query(
            gdf_rn.geometry.values, predicate="intersects"
        )
        crossings = shapely.intersection(
            gdf_rn.geometry.values[river_idx], gdf_lv.geometry.values[levee_idx]
        )
        parts = shapely.get_parts(crossings, return_index=True)
        points = shapely.point_on_surface(parts[0])
        owner = river_idx[parts[1]]
        cross_idx, structure_idx = gdf_close.sindex.query(
            points, predicate="dwithin", distance=gap
        )
        if len(cross_idx):
            df_hits = (
                pd.DataFrame(
                    {
                        "river": owner[cross_idx],
                        "tags": gdf_close["tags"].to_numpy()[structure_idx],
                        "distance": shapely.distance(
                            points[cross_idx], gdf_close.geometry.values[structure_idx]
                        ),
                    }
                )
                .sort_values("distance")
                .drop_duplicates("river")
            )
            for river, tags, distance in df_hits.itertuples(index=False):
                reasons[river] = (
                    f"closing structure ({tags}) {distance:.1f} m from a levee crossing"
                )
        log.info(
            "River-levee crossings: %d; closed by a nearby closing structure: %d river features.",
            len(points),
            int((reasons != "").sum()),
        )

    # --- Downstream ends from DTM elevations at both line ends ---
    list_xy = [shapely.get_coordinates(g) for g in gdf_rn.geometry]
    first = np.array([c[0] for c in list_xy])
    last = np.array([c[-1] for c in list_xy])
    radius = country["endpoint_elevation_radius_m"]
    elevation = np.full((len(gdf_rn), 2), np.nan)
    for k, (start, end) in enumerate(zip(first, last)):
        for j, (x, y) in enumerate((start, end)):
            window, _ = gis.Raster.read_native(
                dtm, (x - radius, y - radius, x + radius, y + radius)
            )
            if window.size and np.isfinite(window).any():
                elevation[k, j] = np.nanmedian(window)
        if (k + 1) % 5000 == 0:
            log.info(
                "  endpoint elevations: %d of %d river features", k + 1, len(gdf_rn)
            )
    by_dtm = np.isfinite(elevation).all(axis=1) & (
        np.abs(elevation[:, 0] - elevation[:, 1]) > 1e-3
    )
    start_lower = by_dtm & (elevation[:, 0] < elevation[:, 1])
    downstream = np.where(start_lower[:, None], first, last)
    gdf_rn["downstream_x"], gdf_rn["downstream_y"] = downstream[:, 0], downstream[:, 1]
    gdf_rn["downstream_method"] = np.where(by_dtm, "dtm elevation", "line direction")
    log.info(
        "Downstream ends: %d from DTM elevations, %d from the line direction (equal or missing elevations).",
        int(by_dtm.sum()),
        int((~by_dtm).sum()),
    )

    # --- Pumping stations at the downstream end ---
    if len(gdf_pump):
        end_idx, pump_idx = gdf_pump.sindex.query(
            shapely.points(downstream), predicate="dwithin", distance=gap
        )
        for river, pump in zip(end_idx, pump_idx):
            if not reasons[river]:
                distance = shapely.distance(
                    shapely.points(downstream[river]), gdf_pump.geometry.values[pump]
                )
                reasons[river] = (
                    f"pumping station {distance:.1f} m from the downstream end"
                )

    gdf_rn["is_source"] = reasons == ""
    gdf_rn["reason"] = reasons
    gis.Vector.write_gpkg(gdf_rn, out_dir / "rivers_sources.gpkg")
    log.info(
        "rivers_sources.gpkg: %d river features, %d closed during floods (not water sources).",
        len(gdf_rn),
        int((~gdf_rn.is_source).sum()),
    )

    # --- Pumping stations next to levees or rivers ---
    near_levee = np.zeros(len(gdf_pump), bool)
    near_river = np.zeros(len(gdf_pump), bool)
    if len(gdf_pump) and len(gdf_lv):
        near_levee[
            np.unique(
                gdf_lv.sindex.query(
                    gdf_pump.geometry.values, predicate="dwithin", distance=near
                )[0]
            )
        ] = True
    if len(gdf_pump) and len(gdf_rn):
        near_river[
            np.unique(
                gdf_rn.sindex.query(
                    gdf_pump.geometry.values, predicate="dwithin", distance=near
                )[0]
            )
        ] = True
    gdf_sel = gdf_pump.assign(near_levee=near_levee, near_river=near_river)[
        near_levee | near_river
    ]
    gis.Vector.write_gpkg(gdf_sel.reset_index(drop=True), out_dir / "osm_pumping.gpkg")
    log.info(
        "osm_pumping.gpkg: %d of %d pumping stations within %.0f m of a levee or river.",
        len(gdf_sel),
        len(gdf_pump),
        near,
    )
