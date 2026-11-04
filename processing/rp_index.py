# Author: Jakub Zapletal
# Date: 2026-11-04
# Estimates a return period of protection per compartment from exposure (GHSL building volume valued with JRC-based
# maximum damages and depth-damage curves, plus the land cover value of HAND_YZ.py when enabled) per metre of levee
# and from the river size (upstream area).

import json
import logging
import os
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
FORCE = False  # recompute the exposure of units that are already finished
UNITS = []  # unit ids to process; empty = all units with compartments

KEYS = [
    "general.work_dir",
    "general.output_dir",
    "general.n_workers",
    "general.gdal_cachemax_mb",
    "general.warp_memory_limit_mb",
    "general.gdal_threads_per_worker",
    "general.log_level",
    "country.rivers",
    "country.river_upstream_area_field",
    "country.upstream_area_raster",
    "country.overrides",
    "country.rp.classes",
    "country.rp.distribution_classes",
    "country.rp.distribution_shares",
    "country.rp.max_rp",
    "country.rp.cap_levels_m",
    "country.rp.cap_rp",
    "country.rp.score_weight_value",
    "country.rp.score_weight_upstream",
    "country.rp.levee_adjacent_distance_m",
    "country.rp.storey_height_m",
    "country.rp.max_damage_csv",
    "country.rp.depth_damage_csv",
    "country.rp.value_total_tolerance",
    "country.lulc.enabled",
    "sources.ghsl_nodata_as_zero",
]
OBSERVATION_COLUMNS = {  # filled by rp_bounds.py (Stage C)
    "obs_rp_low": np.nan,
    "obs_rp_high": np.nan,
    "n_held": 0,
    "n_flooded": 0,
    "obs_conflict": False,
    "model_obs_conflict": False,
}


def index_unit(task):
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
    hand = (
        gdal.Open(str(unit_dir / "hand_undef.tif"))
        .ReadAsArray(int(c0), int(r0), int(c1 - c0), int(r1 - r0))
        .astype(np.float32)
    )
    del kept
    res = gt[1]
    gt = (gt[0] + c0 * res, res, 0.0, gt[3] - r0 * res, 0.0, -res)
    rows, cols = labels.shape
    bounds = (gt[0], gt[3] - rows * res, gt[0] + cols * res, gt[3])

    # Protection level of each kept compartment as a per-cell raster; depth of water below it.
    level = np.full(labels.max() + 1, np.nan, np.float32)
    level[np.array(list(task["levels"]), dtype=np.int64)] = list(
        task["levels"].values()
    )
    depth = np.clip(level[labels] - hand, 0.0, None)
    inside = np.isfinite(depth)

    # Building volume warped to the unit grid preserving totals: the average of the source cells (volume per
    # source cell) times the ratio of the cell areas is the volume per working cell.
    volume, checks = {}, {}
    for product, vrt in task["ghsl"].items():
        source = gdal.Open(vrt)
        source_gt, source_wkt = source.GetGeoTransform(), source.GetProjection()
        warped, _ = gis.Raster.read_window(
            vrt, bounds, res, "average", task["warp_memory_mb"], dst_wkt=wkt
        )
        if task["nodata_as_zero"]:
            warped = np.nan_to_num(warped, nan=0.0)
        volume[product] = warped * (res * res) / abs(source_gt[1] * source_gt[5])
        # Independent total: source cells weighted by the share of their area inside the window (window
        # reprojected to the source CRS).
        window = gpd.GeoSeries(
            [shapely.segmentize(shapely.box(*bounds), res * 10)], crs=wkt
        )
        polygon = window.to_crs(source_wkt).iloc[0]
        cells, cells_gt = gis.Raster.read_native(source, polygon.bounds)
        if cells.size:
            # Exact share of each source cell inside the window polygon.
            x_edges = cells_gt[0] + np.arange(cells.shape[1] + 1) * cells_gt[1]
            y_edges = cells_gt[3] + np.arange(cells.shape[0] + 1) * cells_gt[5]
            cell_boxes = shapely.box(
                x_edges[None, :-1],
                y_edges[1:, None],
                x_edges[None, 1:],
                y_edges[:-1, None],
            )
            share = shapely.area(shapely.intersection(cell_boxes, polygon)) / abs(
                cells_gt[1] * cells_gt[5]
            )
            expected = float(np.nansum(np.nan_to_num(cells, nan=0.0) * share))
            got = float(volume[product].sum())
            checks[product] = got / expected if expected > 0 else float("nan")

    residential = np.clip(volume["total"] - volume["non_residential"], 0.0, None)
    floors = task["storey_height_m"]
    value_res = residential / floors * task["max_damage"]["residential"]
    value_nonres = (
        volume["non_residential"] / floors * task["max_damage"]["non_residential"]
    )
    curve = task["depth_damage"]
    damage = np.zeros(labels.shape, np.float64)
    damage[inside] = value_res[inside] * np.interp(
        depth[inside], curve["depth_m"], curve["res"]
    ) + value_nonres[inside] * np.interp(
        depth[inside], curve["depth_m"], curve["nonres"]
    )
    exposed = np.bincount(
        labels.ravel(), weights=damage.ravel(), minlength=labels.max() + 1
    )
    result = pd.DataFrame(
        {
            "unit_label": list(task["levels"]),
            "value_exposed": exposed[list(task["levels"])],
        }
    )
    result.to_csv(unit_dir / "rp_index.csv", index=False)
    seconds = time.perf_counter() - start
    done = {
        "unit_id": task["unit_id"],
        "seconds": seconds,
        "total_ratio": checks,
        "peak_rss_mb": gis.Run.peak_rss_mb(),
        "summary": (
            f"unit {task['unit_id']}: {len(result)} compartments valued in {seconds:.1f} s; warped/source "
            "volume " + ", ".join(f"{k} {v:.4f}" for k, v in checks.items())
        ),
    }
    marker = unit_dir / "rp_index.json.tmp"
    marker.write_text(json.dumps(done, indent=2))
    os.replace(marker, unit_dir / "rp_index.json")
    logging.getLogger("rp_index").info(done["summary"])
    return done


if __name__ == "__main__":
    run = gis.Run.env_overrides(
        CONFIG_PATH=CONFIG_PATH, COUNTRY=COUNTRY, FORCE=FORCE, UNITS=UNITS
    )
    cfg = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], KEYS)
    general, country, rp = cfg["general"], cfg["country"], cfg["country"]["rp"]
    work_dir = Path(general["work_dir"]) / run["COUNTRY"]
    out_dir = Path(general["output_dir"]) / run["COUNTRY"]
    log = gis.Run.logger("rp_index", out_dir / "logs", general["log_level"])
    gis.Config.warn_placeholders(cfg)
    gis.Run.configure_gdal(
        general["gdal_cachemax_mb"], general["gdal_threads_per_worker"]
    )
    needed = [
        out_dir / "compartments.gpkg",
        out_dir / "barriers.gpkg",
        out_dir / "ghsl_total.vrt",
        out_dir / "ghsl_non_residential.vrt",
    ]
    if not all(path.exists() for path in needed):
        raise SystemExit(
            f"Missing inputs {[str(p) for p in needed if not p.exists()]}; run compartments.py "
            "and download_ghsl.py first."
        )
    if not country["river_upstream_area_field"] and not country["upstream_area_raster"]:
        raise SystemExit(
            "Set countries.<code>.river_upstream_area_field or upstream_area_raster: the score "
            "needs the upstream area of the source river."
        )

    # --- Damage tables ---
    max_damage = pd.read_csv(rp["max_damage_csv"])
    expected = {"country", "residential", "non_residential"}
    if not expected <= set(max_damage.columns):
        raise SystemExit(
            f"{rp['max_damage_csv']} needs the columns {sorted(expected)}; found "
            f"{list(max_damage.columns)}."
        )
    row = max_damage[
        max_damage["country"].astype(str).str.upper() == run["COUNTRY"].upper()
    ]
    if len(row) != 1:
        raise SystemExit(
            f"{rp['max_damage_csv']} needs exactly one row with country = {run['COUNTRY']}."
        )
    damage_values = {
        k: float(row[k].iloc[0]) for k in ["residential", "non_residential"]
    }
    curve = pd.read_csv(rp["depth_damage_csv"])
    if not {"depth_m", "res", "nonres"} <= set(curve.columns):
        raise SystemExit(
            f"{rp['depth_damage_csv']} needs the columns depth_m, res, nonres."
        )
    curve = curve.sort_values("depth_m")
    if (
        not curve["depth_m"].is_unique
        or ((curve[["res", "nonres"]] < 0) | (curve[["res", "nonres"]] > 1)).any().any()
    ):
        raise SystemExit(
            f"{rp['depth_damage_csv']}: depths must be unique and fractions within 0..1."
        )
    log.info(
        "Maximum damage for %s: residential %.1f, non-residential %.1f per m2 of floor area.",
        run["COUNTRY"],
        damage_values["residential"],
        damage_values["non_residential"],
    )

    # --- Exposure per unit (parallel, checkpointed) ---
    compartments = gpd.read_file(out_dir / "compartments.gpkg")
    wkt = compartments.crs.to_wkt()
    units = sorted(compartments["unit_id"].unique())
    if run["UNITS"]:
        units = [u for u in units if u in run["UNITS"]]
    tasks, skipped = [], []
    for unit_id in units:
        unit_dir = work_dir / "units" / unit_id
        marker = unit_dir / "rp_index.json"
        if (
            marker.exists()
            and not run["FORCE"]
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
                "ghsl": {
                    p: str(out_dir / f"ghsl_{p}.vrt")
                    for p in ["total", "non_residential"]
                },
                "warp_memory_mb": general["warp_memory_limit_mb"],
                "nodata_as_zero": cfg["sources"]["ghsl_nodata_as_zero"],
                "storey_height_m": rp["storey_height_m"],
                "max_damage": damage_values,
                "depth_damage": {
                    c: curve[c].to_list() for c in ["depth_m", "res", "nonres"]
                },
            }
        )
    log.info(
        "Units with compartments: %d to process, %d finished and skipped.",
        len(tasks),
        len(skipped),
    )
    results, failed = gis.Run.run_units(
        index_unit,
        tasks,
        general["n_workers"],
        general["gdal_cachemax_mb"],
        general["gdal_threads_per_worker"],
        [1] * len(tasks),
        log,
    )
    if failed:
        raise SystemExit(f"Units failed: {failed}; rerun to retry them.")
    ratios = [
        r
        for result in results
        for r in result["total_ratio"].values()
        if np.isfinite(r)
    ]
    if ratios:
        worst = max(abs(r - 1) for r in ratios)
        (log.warning if worst > rp["value_total_tolerance"] else log.info)(
            "Value warp check: warped totals differ from the source totals by at most %.3f %% (tolerance %.1f %%).",
            100 * worst,
            100 * rp["value_total_tolerance"],
        )

    values = pd.concat(
        [
            pd.read_csv(work_dir / "units" / u / "rp_index.csv").assign(unit_id=u)
            for u in units
            if (work_dir / "units" / u / "rp_index.csv").exists()
        ],
        ignore_index=True,
    )
    result = compartments.merge(values, on=["unit_id", "unit_label"], how="left")
    if result["value_exposed"].isna().any():
        log.warning(
            "%d compartments have no exposure value (units not processed).",
            int(result["value_exposed"].isna().sum()),
        )
    # Land cover exposure (HAND_YZ.py) supplements the building value.
    result["value_ghsl"] = result["value_exposed"]
    result["value_lulc"] = 0.0
    if country["lulc"]["enabled"]:
        lulc_path = out_dir / "exposure_lulc.csv"
        if not lulc_path.exists():
            raise SystemExit(
                "lulc.enabled is true but exposure_lulc.csv is missing; run HAND_YZ.py first."
            )
        lulc = pd.read_csv(lulc_path).drop(columns=["unit_id", "unit_label"])
        result = result.drop(columns="value_lulc").merge(
            lulc, on="compartment_id", how="left"
        )
        if result["value_lulc"].isna().any():
            log.warning(
                "%d compartments have no land cover value (HAND_YZ.py not run for their units); 0 is used.",
                int(result["value_lulc"].isna().sum()),
            )
        lulc_columns = [c for c in lulc.columns if c != "compartment_id"]
        result[lulc_columns] = result[lulc_columns].fillna(0.0)
        result["value_exposed"] = result["value_ghsl"] + result["value_lulc"]
        log.info(
            "Land cover value added: %.4g on top of %.4g from buildings (%.1f %%).",
            result["value_lulc"].sum(),
            result["value_ghsl"].sum(),
            100 * result["value_lulc"].sum() / max(result["value_ghsl"].sum(), 1e-9),
        )

    # --- Levee length, upstream area ---
    barriers = gpd.read_file(out_dir / "barriers.gpkg").to_crs(wkt)
    result["levee_length_m"] = gis.Vector.length_near_boundary(
        result, barriers, rp["levee_adjacent_distance_m"]
    )
    result["no_levee"] = result["levee_length_m"] <= 0
    sources_points = gpd.GeoDataFrame(
        geometry=gpd.points_from_xy(result["source_x"], result["source_y"]), crs=wkt
    )
    if country["river_upstream_area_field"]:
        field = country["river_upstream_area_field"]
        rivers = gis.Vector.read_merged(country["rivers"], wkt, keep_columns=[field])
        if field not in rivers.columns:
            raise SystemExit(f"Column {field!r} not found in the river files.")
        nearest = gpd.sjoin_nearest(
            sources_points, rivers[[field, "geometry"]], how="left"
        )
        nearest = nearest[~nearest.index.duplicated()]
        result["upstream_area_km2"] = pd.to_numeric(
            nearest[field], errors="coerce"
        ).to_numpy()
    else:
        raster = gdal.Open(country["upstream_area_raster"])
        points = sources_points.to_crs(raster.GetProjection())
        result["upstream_area_km2"] = [
            (
                float(
                    np.nanmean(gis.Raster.read_native(raster, (p.x, p.y, p.x, p.y))[0])
                )
                if not p.is_empty
                else np.nan
            )
            for p in points.geometry
        ]
    missing_area = ~(result["upstream_area_km2"] > 0)
    if missing_area.any():
        log.warning(
            "%d compartments have no upstream area; their score uses the value term only.",
            int(missing_area.sum()),
        )

    # --- Score and quantile mapping (compartments with levees only) ---
    scored = ~result["no_levee"]
    result["v_rel"] = np.where(
        scored, result["value_exposed"] / result["levee_length_m"].where(scored), np.nan
    )
    with np.errstate(divide="ignore"):
        rank_value = gis.Stats.percentile_rank(
            np.where(scored, np.log(result["v_rel"]), np.nan)
        )
        rank_area = gis.Stats.percentile_rank(
            np.where(
                scored & ~missing_area, np.log(result["upstream_area_km2"]), np.nan
            )
        )
    weight_v, weight_a = rp["score_weight_value"], rp["score_weight_upstream"]
    result["score"] = np.where(
        np.isnan(rank_area),
        rank_value,
        (weight_v * rank_value + weight_a * rank_area) / (weight_v + weight_a),
    )
    result.loc[~scored, "score"] = np.nan
    result["rp_model"] = gis.Stats.quantile_classes(
        result["score"], rp["distribution_classes"], rp["distribution_shares"]
    )

    # --- Cap, final RP and interval ---
    cap_index = np.searchsorted(
        rp["cap_levels_m"], result["protection_level_m"].to_numpy(), side="right"
    )
    result["rp_cap"] = np.where(
        cap_index < len(rp["cap_rp"]),
        np.asarray(rp["cap_rp"], float)[np.minimum(cap_index, len(rp["cap_rp"]) - 1)],
        np.nan,
    )
    candidates = np.column_stack(
        [
            result["rp_model"],
            result["rp_cap"],
            np.full(len(result), float(rp["max_rp"])),
        ]
    )
    result["rp"] = np.nanmin(candidates, axis=1)
    lowered = result["rp_model"].isna() | (result["rp"] < result["rp_model"])
    result["rp_source"] = np.where(lowered, "cap", "model")
    classes = np.asarray(rp["classes"], float)
    position = np.searchsorted(classes, result["rp"].to_numpy())
    result["rp_low"] = classes[np.maximum(position - 1, 0)]
    upper = classes[np.minimum(position + 1, len(classes) - 1)]
    bound = np.nanmin(
        np.column_stack(
            [upper, result["rp_cap"], np.full(len(result), float(rp["max_rp"]))]
        ),
        axis=1,
    )
    result["rp_high"] = np.maximum(result["rp"], bound)

    # --- Overrides ---
    if country["overrides"]:
        overrides = gis.Vector.read_merged(
            country["overrides"], wkt, keep_columns=["rp", "note"]
        )
        joined = gpd.sjoin(
            overrides, result[["geometry"]], predicate="within", how="inner"
        )
        conflicts = joined.groupby("index_right")["rp"].nunique()
        if (conflicts > 1).any():
            log.warning(
                "%d compartments contain override points with different rp; the first is used.",
                int((conflicts > 1).sum()),
            )
        first = joined.drop_duplicates("index_right").set_index("index_right")
        result.loc[first.index, ["rp", "rp_low", "rp_high"]] = np.repeat(
            first["rp"].astype(float).to_numpy()[:, None], 3, axis=1
        )
        result.loc[first.index, "rp_source"] = "override"
        log.info("Overrides applied to %d compartments.", len(first))

    for column, value in OBSERVATION_COLUMNS.items():
        result[column] = value
    for column in ["rp_model", "rp_cap", "rp", "rp_low", "rp_high"]:
        result[column] = result[column].astype(float)
    gis.Vector.write_gpkg(result, out_dir / "compartments_rp.gpkg")
    shares = result.loc[scored, "rp_model"].value_counts(normalize=True).sort_index()
    log.info(
        "Model classes (share of compartments with levees): %s; targets %s.",
        ", ".join(f"{int(k)}: {v:.2f}" for k, v in shares.items()),
        ", ".join(
            f"{c}: {s:.2f}"
            for c, s in zip(rp["distribution_classes"], rp["distribution_shares"])
        ),
    )
    log.info(
        "RP sources: %s. compartments_rp.gpkg written (%d compartments).",
        ", ".join(f"{k} {v}" for k, v in result["rp_source"].value_counts().items()),
        len(result),
    )
