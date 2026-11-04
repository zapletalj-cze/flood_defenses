# Author: Jakub Zapletal
# Date: 2026-11-04
# HAND_YZ: land use / land cover exposure of each compartment (ESA WorldCover; Urban Atlas where available), valued
# with category damages at the HAND depth of the protection level. Supplements the GHSL value of rp_index.py.

import hashlib
import json
import logging
import os
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from osgeo import gdal

import gis

# --- Run variables (no command-line parameters; DEFENSES_<NAME> environment variables override them) ---
CONFIG_PATH = "config/defenses.toml"
COUNTRY = "XX"
FORCE = False  # recompute units whose land cover exposure is up to date
UNITS = []  # unit ids to process; empty = all units with compartments

KEYS = [
    "general.work_dir",
    "general.output_dir",
    "general.n_workers",
    "general.gdal_cachemax_mb",
    "general.warp_memory_limit_mb",
    "general.gdal_threads_per_worker",
    "general.log_level",
    "country.lulc.worldcover",
    "country.lulc.urban_atlas",
    "country.lulc.urban_atlas_layer",
    "country.lulc.urban_atlas_class_field",
    "country.lulc.worldcover_categories",
    "country.lulc.urban_atlas_categories",
    "country.lulc.max_damage_csv",
    "country.lulc.depth_damage_csv",
]


def exposure_unit(task):
    start = time.perf_counter()
    unit_dir = Path(task["unit_dir"])
    labels_ds = gdal.Open(str(unit_dir / "labels.tif"))
    gt, wkt = labels_ds.GetGeoTransform(), labels_ds.GetProjection()
    labels = labels_ds.ReadAsArray()
    # Only the bounding box of the unit's kept compartments is valued.
    kept = np.isin(labels, np.array(list(task["levels"]), dtype=labels.dtype))
    rows_any, cols_any = np.flatnonzero(kept.any(axis=1)), np.flatnonzero(
        kept.any(axis=0)
    )
    r0, r1, c0, c1 = rows_any[0], rows_any[-1] + 1, cols_any[0], cols_any[-1] + 1
    labels = np.where(kept, labels, 0)[r0:r1, c0:c1]
    del kept
    hand = (
        gdal.Open(str(unit_dir / "hand_undef.tif"))
        .ReadAsArray(int(c0), int(r0), int(c1 - c0), int(r1 - r0))
        .astype(np.float32)
    )
    res = gt[1]
    gt = (gt[0] + c0 * res, res, 0.0, gt[3] - r0 * res, 0.0, -res)
    bounds = (
        gt[0],
        gt[3] - labels.shape[0] * res,
        gt[0] + labels.shape[1] * res,
        gt[3],
    )
    # Water depth when the water stands at the protection level (the design event), as in rp_index.py.
    level = np.full(labels.max() + 1, np.nan, np.float32)
    level[np.array(list(task["levels"]), dtype=np.int64)] = list(
        task["levels"].values()
    )
    depth = np.clip(level[labels] - hand, 0.0, None)
    inside = np.isfinite(depth)

    # Category per cell: WorldCover classes through the lookup table; Urban Atlas replaces them where it exists,
    # also with classes that carry no value.
    classes, _ = gis.Raster.read_window(
        task["worldcover"], bounds, res, "near", task["warp_memory_mb"], dst_wkt=wkt
    )
    lookup = np.asarray(task["worldcover_lookup"], np.int16)
    category = np.zeros(labels.shape, np.int16)
    known = np.isfinite(classes) & (classes >= 0) & (classes < len(lookup))
    category[known] = lookup[classes[known].astype(np.int64)]
    del classes, known
    ua_cells = 0
    if task["urban_atlas"]:
        gdf_poly = gpd.read_file(task["urban_atlas"], bbox=bounds)
        if len(gdf_poly):
            covered = (
                gis.Raster.rasterize_geometries(
                    gdf_poly.geometry,
                    np.ones(len(gdf_poly)),
                    labels.shape,
                    gt,
                    wkt,
                    False,
                    np.uint8,
                )
                > 0
            )
            burnt = gis.Raster.rasterize_geometries(
                gdf_poly.geometry,
                gdf_poly["category_id"].to_numpy(),
                labels.shape,
                gt,
                wkt,
                False,
                np.int16,
            )
            category[covered] = burnt[covered]
            ua_cells = int((covered & inside).sum())

    # Damage per cell = value of the land (per m2 times cell area) times the category's damage fraction at depth.
    cell_area = res * res
    n_labels = labels.max() + 1
    per_category = {}
    exposed = np.zeros(n_labels)
    for k, (name, value_m2, list_frac) in enumerate(task["categories"], start=1):
        cells = inside & (category == k)
        damage = (
            value_m2 * cell_area * np.interp(depth[cells], task["depth_m"], list_frac)
        )
        per_category[name] = np.bincount(
            labels[cells], weights=damage, minlength=n_labels
        )
        exposed += np.bincount(
            labels[cells], weights=(damage > 0) * cell_area, minlength=n_labels
        )
    list_id = list(task["levels"])
    df_res = pd.DataFrame(
        {
            "unit_label": list_id,
            "value_lulc": np.sum([v[list_id] for v in per_category.values()], axis=0),
            "lulc_exposed_area_m2": exposed[list_id],
            **{
                f"lulc_{name}": values[list_id] for name, values in per_category.items()
            },
        }
    )
    df_res.to_csv(unit_dir / "hand_yz.csv", index=False)
    seconds = time.perf_counter() - start
    done = {
        "unit_id": task["unit_id"],
        "parameters": task["parameters_hash"],
        "seconds": seconds,
        "peak_rss_mb": gis.Run.peak_rss_mb(),
        "cells": int(labels.size),
        "urban_atlas_cells": ua_cells,
        "summary": (
            f"unit {task['unit_id']}: {len(df_res)} compartments, land cover value "
            f"{df_res['value_lulc'].sum():.3g}, Urban Atlas cells {ua_cells}, {seconds:.1f} s"
        ),
    }
    marker = unit_dir / "hand_yz.json.tmp"
    marker.write_text(json.dumps(done, indent=2))
    os.replace(marker, unit_dir / "hand_yz.json")
    logging.getLogger("HAND_YZ").info(done["summary"])
    return done


if __name__ == "__main__":
    run = gis.Run.env_overrides(
        CONFIG_PATH=CONFIG_PATH, COUNTRY=COUNTRY, FORCE=FORCE, UNITS=UNITS
    )
    cfg = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], KEYS)
    general, lulc = cfg["general"], cfg["country"]["lulc"]
    work_dir = Path(general["work_dir"]) / run["COUNTRY"]
    out_dir = Path(general["output_dir"]) / run["COUNTRY"]
    log = gis.Run.logger("HAND_YZ", out_dir / "logs", general["log_level"])
    gis.Config.warn_placeholders(cfg)
    gis.Run.configure_gdal(
        general["gdal_cachemax_mb"], general["gdal_threads_per_worker"]
    )
    if not (out_dir / "compartments.gpkg").exists():
        raise SystemExit("compartments.gpkg is missing; run compartments.py first.")
    if not lulc["max_damage_csv"] or not lulc["depth_damage_csv"]:
        raise SystemExit(
            "Set countries.<code>.lulc.max_damage_csv and depth_damage_csv (values per category)."
        )
    if lulc["worldcover"]:
        worldcover = gis.Raster.build_vrt(
            lulc["worldcover"], work_dir / "lulc" / "worldcover.vrt"
        )
    elif (out_dir / "worldcover.vrt").exists():
        worldcover = str(out_dir / "worldcover.vrt")
    else:
        raise SystemExit(
            "No WorldCover data: run download_worldcover.py or list tiles in countries.<code>.lulc.worldcover."
        )

    # --- Categories, values per m2 of land and damage curves ---
    list_cat = sorted(
        set(lulc["worldcover_categories"].values())
        | set(lulc["urban_atlas_categories"].values())
    )
    if not list_cat:
        raise SystemExit(
            "No land cover categories: fill lulc.worldcover_categories or lulc.urban_atlas_categories."
        )
    index = {name: k for k, name in enumerate(list_cat, start=1)}
    df_maxd = pd.read_csv(lulc["max_damage_csv"])
    df_curve = pd.read_csv(lulc["depth_damage_csv"]).sort_values("depth_m")
    list_miss = [c for c in ["country"] + list_cat if c not in df_maxd.columns]
    list_miss += [
        f"{c} (curve)" for c in ["depth_m"] + list_cat if c not in df_curve.columns
    ]
    if list_miss:
        raise SystemExit(
            f"Columns missing in the land cover damage tables: {list_miss}."
        )
    df_row = df_maxd[
        df_maxd["country"].astype(str).str.upper() == run["COUNTRY"].upper()
    ]
    if len(df_row) != 1:
        raise SystemExit(
            f"{lulc['max_damage_csv']} needs exactly one row with country = {run['COUNTRY']}."
        )
    if (
        not df_curve["depth_m"].is_unique
        or ((df_curve[list_cat] < 0) | (df_curve[list_cat] > 1)).any().any()
    ):
        raise SystemExit(
            f"{lulc['depth_damage_csv']}: depths must be unique and fractions within 0..1."
        )
    values = {name: float(df_row[name].iloc[0]) for name in list_cat}
    lookup = np.zeros(256, np.int16)
    for code, name in lulc["worldcover_categories"].items():
        lookup[int(code)] = index[name]
    log.info(
        "Land cover categories (value per m2): %s.",
        ", ".join(f"{k} {v:g}" for k, v in values.items()),
    )

    # --- Urban Atlas: single-layer copies, class codes to category ids, only polygons touching compartments ---
    gdf_comp = gpd.read_file(out_dir / "compartments.gpkg")
    wkt = gdf_comp.crs.to_wkt()
    urban_atlas, list_ua = "", []
    if lulc["urban_atlas"]:
        list_path = []
        for path in gis.Config.as_list(lulc["urban_atlas"]):
            list_ua.append((str(path), os.path.getmtime(path)))
            list_name = gis.Vector.layer_names(path)
            if len(list_name) == 1:
                list_path.append(path)
            elif lulc["urban_atlas_layer"]:
                list_path.append(
                    gis.Vector.single_layer_copy(
                        path,
                        lulc["urban_atlas_layer"],
                        work_dir / "lulc" / f"{Path(path).stem}_layer.gpkg",
                    )
                )
                log.info(
                    "%s: layer %s copied (layers %s).",
                    Path(path).name,
                    lulc["urban_atlas_layer"],
                    list_name,
                )
            else:
                raise SystemExit(
                    f"{path} has layers {list_name}; set countries.<code>.lulc.urban_atlas_layer."
                )
        field = lulc["urban_atlas_class_field"]
        gdf_poly = gis.Vector.read_merged(
            list_path, wkt, bbox=tuple(gdf_comp.total_bounds), keep_columns=[field]
        )
        if field not in gdf_poly.columns:
            raise SystemExit(f"Field {field!r} not found in the Urban Atlas files.")
        touching, _ = gdf_comp.sindex.query(
            gdf_poly.geometry.values, predicate="intersects"
        )
        gdf_poly = gdf_poly.iloc[np.unique(touching)].reset_index(drop=True)
        codes = (
            gdf_poly[field].astype(str).str.strip().str.replace(r"\.0$", "", regex=True)
        )
        mapping = {
            str(code): index[name]
            for code, name in lulc["urban_atlas_categories"].items()
        }
        gdf_poly["category_id"] = codes.map(mapping).fillna(0).astype(int)
        list_unmap = sorted(set(codes[gdf_poly["category_id"] == 0]))
        log.info(
            "Urban Atlas: %d polygons touch compartments; classes without a category (no value): %s.",
            len(gdf_poly),
            ", ".join(list_unmap) or "none",
        )
        urban_atlas = str(work_dir / "lulc" / "urban_atlas_categories.gpkg")
        gis.Vector.write_gpkg(gdf_poly[["category_id", "geometry"]], urban_atlas)

    # --- Exposure per unit (parallel, checkpointed with a hash of the inputs) ---
    parameters = {
        "worldcover": worldcover,
        "worldcover_lookup": lookup.tolist(),
        "urban_atlas": urban_atlas,
        "categories": [(n, values[n], df_curve[n].tolist()) for n in list_cat],
        "depth_m": df_curve["depth_m"].tolist(),
    }
    # The hash covers the parameters and the source files (VRT members, Urban Atlas files, class mapping).
    stamp = {
        "worldcover_files": [
            (f, os.path.getmtime(f)) for f in gdal.Open(worldcover).GetFileList()
        ],
        "urban_atlas_files": list_ua,
        "urban_atlas_categories": lulc["urban_atlas_categories"],
        "urban_atlas_field": lulc["urban_atlas_class_field"],
    }
    phash = hashlib.sha256(
        json.dumps([parameters, stamp], sort_keys=True, default=str).encode()
    ).hexdigest()
    list_unit = sorted(gdf_comp["unit_id"].unique())
    if run["UNITS"]:
        list_unit = [u for u in list_unit if u in run["UNITS"]]
    list_task, list_skip = [], []
    for unit_id in list_unit:
        unit_dir = work_dir / "units" / unit_id
        marker = unit_dir / "hand_yz.json"
        if (
            marker.exists()
            and not run["FORCE"]
            and json.loads(marker.read_text()).get("parameters") == phash
            and marker.stat().st_mtime >= (unit_dir / "done.json").stat().st_mtime
        ):
            list_skip.append(unit_id)
            continue
        gdf_mem = gdf_comp[gdf_comp["unit_id"] == unit_id]
        list_task.append(
            {
                "unit_id": unit_id,
                "unit_dir": str(unit_dir),
                "levels": dict(
                    zip(
                        gdf_mem["unit_label"].astype(int),
                        gdf_mem["protection_level_m"].astype(float),
                    )
                ),
                "warp_memory_mb": general["warp_memory_limit_mb"],
                "parameters_hash": phash,
                **parameters,
            }
        )
    log.info(
        "Units with compartments: %d to process, %d up to date and skipped.",
        len(list_task),
        len(list_skip),
    )
    list_res, list_fail = gis.Run.run_units(
        exposure_unit,
        list_task,
        general["n_workers"],
        general["gdal_cachemax_mb"],
        general["gdal_threads_per_worker"],
        [1] * len(list_task),
        log,
    )
    if list_fail:
        raise SystemExit(f"Units failed: {list_fail}; rerun to retry them.")

    list_exp = [
        pd.read_csv(work_dir / "units" / u / "hand_yz.csv").assign(unit_id=u)
        for u in list_unit
        if (work_dir / "units" / u / "hand_yz.csv").exists()
    ]
    df_exp = gdf_comp[["compartment_id", "unit_id", "unit_label"]].merge(
        pd.concat(list_exp, ignore_index=True),
        on=["unit_id", "unit_label"],
        how="inner",
    )
    df_exp.to_csv(out_dir / "exposure_lulc.csv", index=False)
    log.info(
        "exposure_lulc.csv: %d compartments, land cover value %.4g (%s).",
        len(df_exp),
        df_exp["value_lulc"].sum(),
        ", ".join(f"{n} {df_exp[f'lulc_{n}'].sum():.3g}" for n in list_cat),
    )
