# Author: Jakub Zapletal
# Date: 2026-11-04
# Builds the context zones of the river corridor (industrial, dense urban, urban, rural) from CORINE, Urban Atlas, the
# GHSL settlement model; earlier contexts of context.order win where sources overlap.

import time
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
    "country.context.order",
    "country.context.corine_industrial_codes",
    "country.context.corine_dense_urban_codes",
    "country.context.urban_atlas_industrial_codes",
    "country.context.urban_atlas_dense_urban_codes",
    "country.context.ghsl_smod_dense_urban",
    "country.context.ghsl_smod_urban",
    "country.context.min_zone_area_m2",
]

if __name__ == "__main__":
    run = gis.Run.env_overrides(CONFIG_PATH=CONFIG_PATH, COUNTRY=COUNTRY)
    cfg = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], KEYS)
    general, context = cfg["general"], cfg["country"]["context"]
    work_dir = Path(general["work_dir"]) / run["COUNTRY"]
    out_dir = Path(general["output_dir"]) / run["COUNTRY"]
    log = gis.Run.logger("context_zones", out_dir / "logs", general["log_level"])
    gis.Config.warn_placeholders(cfg)
    corridor_path, vrt_path = out_dir / "corridor.gpkg", work_dir / "dtm" / "dtm.vrt"
    if not corridor_path.exists() or not vrt_path.exists():
        raise SystemExit(
            "corridor.gpkg or dtm.vrt is missing; run prepare_dtm.py and build_units.py first."
        )
    start = time.perf_counter()
    wkt = gdal.Open(str(vrt_path)).GetProjection()
    corridor = shapely.union_all(
        gpd.read_file(corridor_path).to_crs(wkt).geometry.values
    )

    # --- Source polygons per context, in the DTM CRS ---
    pieces = {name: [] for name in context["order"]}
    list_rsrc = [
        (
            "CORINE",
            out_dir / "corine_codes.tif",
            "industrial",
            context["corine_industrial_codes"],
        ),
        (
            "CORINE",
            out_dir / "corine_codes.tif",
            "dense_urban",
            context["corine_dense_urban_codes"],
        ),
        (
            "GHSL SMOD",
            out_dir / "ghsl_smod.vrt",
            "dense_urban",
            context["ghsl_smod_dense_urban"],
        ),
        ("GHSL SMOD", out_dir / "ghsl_smod.vrt", "urban", context["ghsl_smod_urban"]),
    ]
    for label, path, name, list_code in list_rsrc:
        if not list_code:
            continue
        if not path.exists():
            log.warning(
                "%s (%s) not found: no %s zones from it.", label, path.name, name
            )
            continue
        ds = gdal.Open(str(path))
        raster_wkt = ds.GetProjection()
        area = gpd.GeoSeries([shapely.segmentize(corridor, 1000.0)], crs=wkt).to_crs(
            raster_wkt
        )
        values, window_gt = gis.Raster.read_native(ds, tuple(area.total_bounds))
        mask = np.isin(values, np.asarray(list_code, np.float32))
        if not mask.any():
            log.info("%s: no cell with codes %s in the corridor.", label, list_code)
            continue
        gdf_poly = gis.Raster.polygonize_labels(
            mask.astype(np.int32), window_gt, raster_wkt
        )
        gdf_poly = gdf_poly[gdf_poly["label"] == 1].to_crs(wkt)
        pieces[name].append((label, shapely.union_all(gdf_poly.geometry.values)))
        log.info(
            "%s codes %s -> %s: %.2f km2 in the raster window.",
            label,
            list_code,
            name,
            gdf_poly.area.sum() / 1e6,
        )
    list_vsrc = [
        (
            "Urban Atlas",
            out_dir / "urban_atlas.gpkg",
            "industrial",
            context["urban_atlas_industrial_codes"],
        ),
        (
            "Urban Atlas",
            out_dir / "urban_atlas.gpkg",
            "dense_urban",
            context["urban_atlas_dense_urban_codes"],
        ),
    ]
    for label, path, name, list_code in list_vsrc:
        if list_code is not None and not list_code:
            continue
        if not path.exists():
            log.warning(
                "%s (%s) not found: no %s zones from it.", label, path.name, name
            )
            continue
        gis.Vector.check_single_layer(path)
        gdf_vec = gpd.read_file(path).to_crs(wkt)
        if list_code is not None:
            gdf_vec = gdf_vec[gdf_vec["ua_code"].isin(list_code)]
        gdf_vec = gdf_vec[gdf_vec.geometry.notna()]
        # Urban Atlas polygons can be invalid (self-touching rings); repair instead of dropping them.
        gdf_vec = gdf_vec.set_geometry(shapely.make_valid(gdf_vec.geometry.values))
        if gdf_vec.empty:
            log.info("%s: no polygon for %s.", label, name)
            continue
        pieces[name].append((label, shapely.union_all(gdf_vec.geometry.values)))
        log.info(
            "%s -> %s: %d polygons, %.2f km2.",
            label,
            name,
            len(gdf_vec),
            gdf_vec.area.sum() / 1e6,
        )

    # --- Zones in context order: earlier contexts win; small parts fall through to later contexts ---
    list_rows, taken = [], shapely.Polygon()
    for name in context["order"]:
        if name == "rural":
            geometry = shapely.difference(corridor, taken)
            labels = "remainder of the corridor"
            list_part = [g for g in shapely.get_parts(geometry) if not g.is_empty]
        else:
            if not pieces[name]:
                continue
            geometry = shapely.difference(
                shapely.intersection(
                    shapely.union_all([g for _, g in pieces[name]]), corridor
                ),
                taken,
            )
            labels = ", ".join(sorted({label for label, _ in pieces[name]}))
            list_part = [
                g
                for g in shapely.get_parts(geometry)
                if not g.is_empty and g.area >= context["min_zone_area_m2"]
            ]
            dropped = len(shapely.get_parts(geometry)) - len(list_part)
            if dropped:
                log.info(
                    "%s: %d parts smaller than %.0f m2 dropped; later contexts take their area.",
                    name,
                    dropped,
                    context["min_zone_area_m2"],
                )
        kept = shapely.union_all(list_part) if list_part else shapely.Polygon()
        taken = shapely.union_all([taken, kept])
        list_rows += [
            {"context": name, "sources": labels, "geometry": g} for g in list_part
        ]
    gdf_zone = (
        gpd.GeoDataFrame(list_rows, geometry="geometry", crs=wkt)
        if list_rows
        else gpd.GeoDataFrame(
            {"context": pd.Series(dtype=str), "sources": pd.Series(dtype=str)},
            geometry=gpd.GeoSeries([], crs=wkt),
        )
    )
    gdf_zone["area_m2"] = gdf_zone.area
    gdf_zone.insert(0, "zone_id", np.arange(len(gdf_zone)))
    gis.Vector.write_gpkg(gdf_zone, out_dir / "zones.gpkg")
    df_sum = gdf_zone.groupby("context")["area_m2"].agg(["count", "sum"])
    log.info(
        "zones.gpkg in %.1f s: %s (corridor %.1f km2).",
        time.perf_counter() - start,
        "; ".join(
            f"{name} {int(df_sum.loc[name, 'count'])} zones {df_sum.loc[name, 'sum'] / 1e6:.2f} km2"
            for name in context["order"]
            if name in df_sum.index
        ),
        corridor.area / 1e6,
    )
