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
        polygons = gpd.read_file(task["urban_atlas"], bbox=bounds)
        if len(polygons):
            covered = (
                gis.Raster.rasterize_geometries(
                    polygons.geometry,
                    np.ones(len(polygons)),
                    labels.shape,
                    gt,
                    wkt,
                    False,
                    np.uint8,
                )
                > 0
            )
            burnt = gis.Raster.rasterize_geometries(
                polygons.geometry,
                polygons["category_id"].to_numpy(),
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
    for k, (name, value_m2, fractions) in enumerate(task["categories"], start=1):
        cells = inside & (category == k)
        damage = (
            value_m2 * cell_area * np.interp(depth[cells], task["depth_m"], fractions)
        )
        per_category[name] = np.bincount(
            labels[cells], weights=damage, minlength=n_labels
        )
        exposed += np.bincount(
            labels[cells], weights=(damage > 0) * cell_area, minlength=n_labels
        )
    ids = list(task["levels"])
    result = pd.DataFrame(
        {
            "unit_label": ids,
            "value_lulc": np.sum([v[ids] for v in per_category.values()], axis=0),
            "lulc_exposed_area_m2": exposed[ids],
            **{f"lulc_{name}": values[ids] for name, values in per_category.items()},
        }
    )
    result.to_csv(unit_dir / "hand_yz.csv", index=False)
    seconds = time.perf_counter() - start
    done = {
        "unit_id": task["unit_id"],
        "parameters": task["parameters_hash"],
        "seconds": seconds,
        "peak_rss_mb": gis.Run.peak_rss_mb(),
        "cells": int(labels.size),
        "urban_atlas_cells": ua_cells,
        "summary": (
            f"unit {task['unit_id']}: {len(result)} compartments, land cover value "
            f"{result['value_lulc'].sum():.3g}, Urban Atlas cells {ua_cells}, {seconds:.1f} s"
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
    categories = sorted(
        set(lulc["worldcover_categories"].values())
        | set(lulc["urban_atlas_categories"].values())
    )
    if not categories:
        raise SystemExit(
            "No land cover categories: fill lulc.worldcover_categories or lulc.urban_atlas_categories."
        )
    index = {name: k for k, name in enumerate(categories, start=1)}
    max_damage = pd.read_csv(lulc["max_damage_csv"])
    curves = pd.read_csv(lulc["depth_damage_csv"]).sort_values("depth_m")
    missing = [c for c in ["country"] + categories if c not in max_damage.columns]
    missing += [
        f"{c} (curve)" for c in ["depth_m"] + categories if c not in curves.columns
    ]
    if missing:
        raise SystemExit(f"Columns missing in the land cover damage tables: {missing}.")
    row = max_damage[
        max_damage["country"].astype(str).str.upper() == run["COUNTRY"].upper()
    ]
    if len(row) != 1:
        raise SystemExit(
            f"{lulc['max_damage_csv']} needs exactly one row with country = {run['COUNTRY']}."
        )
    if (
        not curves["depth_m"].is_unique
        or ((curves[categories] < 0) | (curves[categories] > 1)).any().any()
    ):
        raise SystemExit(
            f"{lulc['depth_damage_csv']}: depths must be unique and fractions within 0..1."
        )
    values = {name: float(row[name].iloc[0]) for name in categories}
    lookup = np.zeros(256, np.int16)
    for code, name in lulc["worldcover_categories"].items():
        lookup[int(code)] = index[name]
    log.info(
        "Land cover categories (value per m2): %s.",
        ", ".join(f"{k} {v:g}" for k, v in values.items()),
    )

    # --- Urban Atlas: single-layer copies, class codes to category ids, only polygons touching compartments ---
    compartments = gpd.read_file(out_dir / "compartments.gpkg")
    wkt = compartments.crs.to_wkt()
    urban_atlas, ua_sources = "", []
    if lulc["urban_atlas"]:
        paths = []
        for path in gis.Config.as_list(lulc["urban_atlas"]):
            ua_sources.append((str(path), os.path.getmtime(path)))
            names = gis.Vector.layer_names(path)
            if len(names) == 1:
                paths.append(path)
            elif lulc["urban_atlas_layer"]:
                paths.append(
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
                    names,
                )
            else:
                raise SystemExit(
                    f"{path} has layers {names}; set countries.<code>.lulc.urban_atlas_layer."
                )
        field = lulc["urban_atlas_class_field"]
        polygons = gis.Vector.read_merged(
            paths, wkt, bbox=tuple(compartments.total_bounds), keep_columns=[field]
        )
        if field not in polygons.columns:
            raise SystemExit(f"Field {field!r} not found in the Urban Atlas files.")
        touching, _ = compartments.sindex.query(
            polygons.geometry.values, predicate="intersects"
        )
        polygons = polygons.iloc[np.unique(touching)].reset_index(drop=True)
        codes = (
            polygons[field].astype(str).str.strip().str.replace(r"\.0$", "", regex=True)
        )
        mapping = {
            str(code): index[name]
            for code, name in lulc["urban_atlas_categories"].items()
        }
        polygons["category_id"] = codes.map(mapping).fillna(0).astype(int)
        unmapped = sorted(set(codes[polygons["category_id"] == 0]))
        log.info(
            "Urban Atlas: %d polygons touch compartments; classes without a category (no value): %s.",
            len(polygons),
            ", ".join(unmapped) or "none",
        )
        urban_atlas = str(work_dir / "lulc" / "urban_atlas_categories.gpkg")
        gis.Vector.write_gpkg(polygons[["category_id", "geometry"]], urban_atlas)

    # --- Exposure per unit (parallel, checkpointed with a hash of the inputs) ---
    parameters = {
        "worldcover": worldcover,
        "worldcover_lookup": lookup.tolist(),
        "urban_atlas": urban_atlas,
        "categories": [(n, values[n], curves[n].tolist()) for n in categories],
        "depth_m": curves["depth_m"].tolist(),
    }
    # The hash covers the parameters and the source files (VRT members, Urban Atlas files, class mapping).
    stamp = {
        "worldcover_files": [
            (f, os.path.getmtime(f)) for f in gdal.Open(worldcover).GetFileList()
        ],
        "urban_atlas_files": ua_sources,
        "urban_atlas_categories": lulc["urban_atlas_categories"],
        "urban_atlas_field": lulc["urban_atlas_class_field"],
    }
    parameters_hash = hashlib.sha256(
        json.dumps([parameters, stamp], sort_keys=True, default=str).encode()
    ).hexdigest()
    units = sorted(compartments["unit_id"].unique())
    if run["UNITS"]:
        units = [u for u in units if u in run["UNITS"]]
    tasks, skipped = [], []
    for unit_id in units:
        unit_dir = work_dir / "units" / unit_id
        marker = unit_dir / "hand_yz.json"
        if (
            marker.exists()
            and not run["FORCE"]
            and json.loads(marker.read_text()).get("parameters") == parameters_hash
            and marker.stat().st_mtime >= (unit_dir / "done.json").stat().st_mtime
        ):
            skipped.append(unit_id)
            continue
        members = compartments[compartments["unit_id"] == unit_id]
        tasks.append(
            {
                "unit_id": unit_id,
                "unit_dir": str(unit_dir),
                "levels": dict(
                    zip(
                        members["unit_label"].astype(int),
                        members["protection_level_m"].astype(float),
                    )
                ),
                "warp_memory_mb": general["warp_memory_limit_mb"],
                "parameters_hash": parameters_hash,
                **parameters,
            }
        )
    log.info(
        "Units with compartments: %d to process, %d up to date and skipped.",
        len(tasks),
        len(skipped),
    )
    results, failed = gis.Run.run_units(
        exposure_unit,
        tasks,
        general["n_workers"],
        general["gdal_cachemax_mb"],
        general["gdal_threads_per_worker"],
        [1] * len(tasks),
        log,
    )
    if failed:
        raise SystemExit(f"Units failed: {failed}; rerun to retry them.")

    tables = [
        pd.read_csv(work_dir / "units" / u / "hand_yz.csv").assign(unit_id=u)
        for u in units
        if (work_dir / "units" / u / "hand_yz.csv").exists()
    ]
    table = compartments[["compartment_id", "unit_id", "unit_label"]].merge(
        pd.concat(tables, ignore_index=True), on=["unit_id", "unit_label"], how="inner"
    )
    table.to_csv(out_dir / "exposure_lulc.csv", index=False)
    log.info(
        "exposure_lulc.csv: %d compartments, land cover value %.4g (%s).",
        len(table),
        table["value_lulc"].sum(),
        ", ".join(f"{n} {table[f'lulc_{n}'].sum():.3g}" for n in categories),
    )
