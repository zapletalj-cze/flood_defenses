# Author: Jakub Zapletal
# Date: 2026-11-04
# Estimates a return period of protection per protected area from exposure (GHSL building volume valued with JRC-based
# maximum damages and depth-damage curves, plus the land cover value of HAND_YZ.py when enabled) per metre of levee
# and from the river size (upstream area). LAYER "A": protected_areas.gpkg (all water sources together), "B":
# protected_areas_by_group.gpkg (per river group, with the labels and HAND of each group).

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
UNITS = []  # unit ids to process; empty = all units with protected areas
LAYER = "A"  # "A": protected areas of all water sources together, "B": per river group

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
    labels_ds = gdal.Open(task["labels_path"])
    gt, wkt = labels_ds.GetGeoTransform(), labels_ds.GetProjection()
    labels = labels_ds.ReadAsArray()
    # Only the bounding box of the unit's kept protected areas is valued.
    kept = np.isin(labels, np.array(list(task["levels"]), dtype=labels.dtype))
    rows_any, cols_any = np.flatnonzero(kept.any(axis=1)), np.flatnonzero(
        kept.any(axis=0)
    )
    r0, r1, c0, c1 = rows_any[0], rows_any[-1] + 1, cols_any[0], cols_any[-1] + 1
    labels = np.where(kept, labels, 0)[r0:r1, c0:c1]
    hand = (
        gdal.Open(task["hand_path"])
        .ReadAsArray(int(c0), int(r0), int(c1 - c0), int(r1 - r0))
        .astype(np.float32)
    )
    del kept
    res = gt[1]
    gt = (gt[0] + c0 * res, res, 0.0, gt[3] - r0 * res, 0.0, -res)
    rows, cols = labels.shape
    bounds = (gt[0], gt[3] - rows * res, gt[0] + cols * res, gt[3])

    # Protection level of each kept protected area as a per-cell raster; depth of water below it.
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
    df_res = pd.DataFrame(
        {
            "unit_label": list(task["levels"]),
            "value_exposed": exposed[list(task["levels"])],
        }
    )
    df_res.to_csv(task["csv_path"], index=False)
    seconds = time.perf_counter() - start
    done = {
        "unit_id": task["unit_id"],
        "seconds": seconds,
        "total_ratio": checks,
        "peak_rss_mb": gis.Run.peak_rss_mb(),
        "summary": (
            f"unit {task['unit_id']}{task['group_text']}: {len(df_res)} protected areas valued in {seconds:.1f} s; warped/source "
            "volume " + ", ".join(f"{k} {v:.4f}" for k, v in checks.items())
        ),
    }
    marker = task["marker_path"] + ".tmp"
    Path(marker).write_text(json.dumps(done, indent=2))
    os.replace(marker, task["marker_path"])
    logging.getLogger("rp_index").info(done["summary"])
    return done


if __name__ == "__main__":
    run = gis.Run.env_overrides(
        CONFIG_PATH=CONFIG_PATH, COUNTRY=COUNTRY, FORCE=FORCE, UNITS=UNITS, LAYER=LAYER
    )
    cfg = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], KEYS)
    general, country, rp = cfg["general"], cfg["country"], cfg["country"]["rp"]
    files = gis.Run.layer_files(run["LAYER"])
    work_dir = Path(general["work_dir"]) / run["COUNTRY"]
    out_dir = Path(general["output_dir"]) / run["COUNTRY"]
    log = gis.Run.logger("rp_index", out_dir / "logs", general["log_level"])
    gis.Config.warn_placeholders(cfg)
    gis.Run.configure_gdal(
        general["gdal_cachemax_mb"], general["gdal_threads_per_worker"]
    )
    list_need = [
        out_dir / files["areas"],
        out_dir / "barriers.gpkg",
        out_dir / "ghsl_total.vrt",
        out_dir / "ghsl_non_residential.vrt",
    ]
    if not all(path.exists() for path in list_need):
        raise SystemExit(
            f"Missing inputs {[str(p) for p in list_need if not p.exists()]}; run protected_areas.py "
            "and download_ghsl.py first."
        )
    if not country["river_upstream_area_field"] and not country["upstream_area_raster"]:
        raise SystemExit(
            "Set countries.<code>.river_upstream_area_field or upstream_area_raster: the score "
            "needs the upstream area of the source river."
        )

    # --- Damage tables ---
    df_maxd = pd.read_csv(rp["max_damage_csv"])
    expected = {"country", "residential", "non_residential"}
    if not expected <= set(df_maxd.columns):
        raise SystemExit(
            f"{rp['max_damage_csv']} needs the columns {sorted(expected)}; found "
            f"{list(df_maxd.columns)}."
        )
    df_row = df_maxd[
        df_maxd["country"].astype(str).str.upper() == run["COUNTRY"].upper()
    ]
    if len(df_row) != 1:
        raise SystemExit(
            f"{rp['max_damage_csv']} needs exactly one row with country = {run['COUNTRY']}."
        )
    damage_values = {
        k: float(df_row[k].iloc[0]) for k in ["residential", "non_residential"]
    }
    df_curve = pd.read_csv(rp["depth_damage_csv"])
    if not {"depth_m", "res", "nonres"} <= set(df_curve.columns):
        raise SystemExit(
            f"{rp['depth_damage_csv']} needs the columns depth_m, res, nonres."
        )
    df_curve = df_curve.sort_values("depth_m")
    if (
        not df_curve["depth_m"].is_unique
        or ((df_curve[["res", "nonres"]] < 0) | (df_curve[["res", "nonres"]] > 1))
        .any()
        .any()
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
    gdf_pa = gpd.read_file(out_dir / files["areas"])
    wkt = gdf_pa.crs.to_wkt()
    # one task per unit (file A) or per unit and river group (file B)
    keys = ["unit_id", "group_id"] if files["layer"] == "B" else ["unit_id"]
    list_key = sorted(gdf_pa[keys].drop_duplicates().itertuples(index=False, name=None))
    if run["UNITS"]:
        list_key = [k for k in list_key if k[0] in run["UNITS"]]
    list_task, list_skip = [], []
    for key in list_key:
        unit_id, group_id = key[0], (key[1] if len(key) > 1 else None)
        unit_dir = work_dir / "units" / unit_id
        marker = gis.Run.unit_layer_path(unit_dir, "rp_index.json", group_id)
        if (
            Path(marker).exists()
            and not run["FORCE"]
            and Path(marker).stat().st_mtime >= (unit_dir / "done.json").stat().st_mtime
        ):
            list_skip.append(key)
            continue
        gdf_mem = gdf_pa[
            (gdf_pa["unit_id"] == unit_id)
            & (gdf_pa["group_id"] == group_id if group_id is not None else True)
        ]
        list_task.append(
            {
                "unit_id": unit_id,
                "group_text": f" group {group_id}" if group_id is not None else "",
                "labels_path": gis.Run.unit_layer_path(
                    unit_dir, "labels.tif", group_id
                ),
                "hand_path": gis.Run.unit_layer_path(
                    unit_dir, "hand_undef.tif", group_id
                ),
                "csv_path": gis.Run.unit_layer_path(unit_dir, "rp_index.csv", group_id),
                "marker_path": marker,
                "levels": dict(
                    zip(
                        gdf_mem["unit_label"].astype(int),
                        gdf_mem["protection_level_m"].astype(float),
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
                    c: df_curve[c].to_list() for c in ["depth_m", "res", "nonres"]
                },
            }
        )
    log.info(
        "File %s: units%s with protected areas: %d to process, %d finished and skipped.",
        files["layer"],
        " and river groups" if files["layer"] == "B" else "",
        len(list_task),
        len(list_skip),
    )
    list_res, list_fail = gis.Run.run_units(
        index_unit,
        list_task,
        general["n_workers"],
        general["gdal_cachemax_mb"],
        general["gdal_threads_per_worker"],
        [1] * len(list_task),
        log,
    )
    if list_fail:
        raise SystemExit(f"Units failed: {list_fail}; rerun to retry them.")
    list_ratio = [
        r
        for result in list_res
        for r in result["total_ratio"].values()
        if np.isfinite(r)
    ]
    if list_ratio:
        worst = max(abs(r - 1) for r in list_ratio)
        (log.warning if worst > rp["value_total_tolerance"] else log.info)(
            "Value warp check: warped totals differ from the source totals by at most %.3f %% (tolerance %.1f %%).",
            100 * worst,
            100 * rp["value_total_tolerance"],
        )

    list_val = []
    for key in list_key:
        group_id = key[1] if len(key) > 1 else None
        path = gis.Run.unit_layer_path(
            work_dir / "units" / key[0], "rp_index.csv", group_id
        )
        if Path(path).exists():
            df_one = pd.read_csv(path).assign(unit_id=key[0])
            if group_id is not None:
                df_one["group_id"] = group_id
            list_val.append(df_one)
    df_val = (
        pd.concat(list_val, ignore_index=True)
        if list_val
        else pd.DataFrame(columns=keys + ["unit_label", "value_exposed"])
    )
    gdf_rp = gdf_pa.merge(df_val, on=keys + ["unit_label"], how="left")
    if gdf_rp["value_exposed"].isna().any():
        log.warning(
            "%d protected areas have no exposure value (units not processed).",
            int(gdf_rp["value_exposed"].isna().sum()),
        )
    # Land cover exposure (HAND_YZ.py) supplements the building value.
    gdf_rp["value_ghsl"] = gdf_rp["value_exposed"]
    gdf_rp["value_lulc"] = 0.0
    if country["lulc"]["enabled"]:
        lulc_path = out_dir / f"exposure_lulc{files['suffix']}.csv"
        if not lulc_path.exists():
            raise SystemExit(
                f"lulc.enabled is true but {lulc_path.name} is missing; run HAND_YZ.py with LAYER = "
                f"{files['layer']} first."
            )
        df_lulc = pd.read_csv(lulc_path).drop(
            columns=["unit_id", "unit_label", "group_id"], errors="ignore"
        )
        gdf_rp = gdf_rp.drop(columns="value_lulc").merge(
            df_lulc, on="protected_area_id", how="left"
        )
        if gdf_rp["value_lulc"].isna().any():
            log.warning(
                "%d protected areas have no land cover value (HAND_YZ.py not run for their units); 0 is used.",
                int(gdf_rp["value_lulc"].isna().sum()),
            )
        list_lulc_col = [c for c in df_lulc.columns if c != "protected_area_id"]
        gdf_rp[list_lulc_col] = gdf_rp[list_lulc_col].fillna(0.0)
        gdf_rp["value_exposed"] = gdf_rp["value_ghsl"] + gdf_rp["value_lulc"]
        log.info(
            "Land cover value added: %.4g on top of %.4g from buildings (%.1f %%).",
            gdf_rp["value_lulc"].sum(),
            gdf_rp["value_ghsl"].sum(),
            100 * gdf_rp["value_lulc"].sum() / max(gdf_rp["value_ghsl"].sum(), 1e-9),
        )

    # --- Levee length, upstream area ---
    gdf_bar = gpd.read_file(out_dir / "barriers.gpkg").to_crs(wkt)
    gdf_rp["levee_length_m"] = gis.Vector.length_near_boundary(
        gdf_rp, gdf_bar, rp["levee_adjacent_distance_m"]
    )
    gdf_rp["no_levee"] = gdf_rp["levee_length_m"] <= 0
    gdf_src = gpd.GeoDataFrame(
        geometry=gpd.points_from_xy(gdf_rp["source_x"], gdf_rp["source_y"]), crs=wkt
    )
    if country["river_upstream_area_field"]:
        field = country["river_upstream_area_field"]
        gdf_rn = gis.Vector.read_merged(country["rivers"], wkt, keep_columns=[field])
        if field not in gdf_rn.columns:
            raise SystemExit(f"Column {field!r} not found in the river files.")
        gdf_sjoin_src_rn = gpd.sjoin_nearest(
            gdf_src, gdf_rn[[field, "geometry"]], how="left"
        )
        gdf_sjoin_src_rn = gdf_sjoin_src_rn[~gdf_sjoin_src_rn.index.duplicated()]
        gdf_rp["upstream_area_km2"] = pd.to_numeric(
            gdf_sjoin_src_rn[field], errors="coerce"
        ).to_numpy()
    else:
        raster = gdal.Open(country["upstream_area_raster"])
        gdf_pts = gdf_src.to_crs(raster.GetProjection())
        gdf_rp["upstream_area_km2"] = [
            (
                float(
                    np.nanmean(gis.Raster.read_native(raster, (p.x, p.y, p.x, p.y))[0])
                )
                if not p.is_empty
                else np.nan
            )
            for p in gdf_pts.geometry
        ]
    missing_area = ~(gdf_rp["upstream_area_km2"] > 0)
    if missing_area.any():
        log.warning(
            "%d protected areas have no upstream area; their score uses the value term only.",
            int(missing_area.sum()),
        )

    # --- Score and quantile mapping (protected areas with levees only) ---
    scored = ~gdf_rp["no_levee"]
    gdf_rp["v_rel"] = np.where(
        scored, gdf_rp["value_exposed"] / gdf_rp["levee_length_m"].where(scored), np.nan
    )
    with np.errstate(divide="ignore"):
        rank_value = gis.Stats.percentile_rank(
            np.where(scored, np.log(gdf_rp["v_rel"]), np.nan)
        )
        rank_area = gis.Stats.percentile_rank(
            np.where(
                scored & ~missing_area, np.log(gdf_rp["upstream_area_km2"]), np.nan
            )
        )
    weight_v, weight_a = rp["score_weight_value"], rp["score_weight_upstream"]
    gdf_rp["score"] = np.where(
        np.isnan(rank_area),
        rank_value,
        (weight_v * rank_value + weight_a * rank_area) / (weight_v + weight_a),
    )
    gdf_rp.loc[~scored, "score"] = np.nan
    gdf_rp["rp_model"] = gis.Stats.quantile_classes(
        gdf_rp["score"], rp["distribution_classes"], rp["distribution_shares"]
    )

    # --- Cap, final RP and interval ---
    cap_index = np.searchsorted(
        rp["cap_levels_m"], gdf_rp["protection_level_m"].to_numpy(), side="right"
    )
    gdf_rp["rp_cap"] = np.where(
        cap_index < len(rp["cap_rp"]),
        np.asarray(rp["cap_rp"], float)[np.minimum(cap_index, len(rp["cap_rp"]) - 1)],
        np.nan,
    )
    candidates = np.column_stack(
        [
            gdf_rp["rp_model"],
            gdf_rp["rp_cap"],
            np.full(len(gdf_rp), float(rp["max_rp"])),
        ]
    )
    gdf_rp["rp"] = np.nanmin(candidates, axis=1)
    lowered = gdf_rp["rp_model"].isna() | (gdf_rp["rp"] < gdf_rp["rp_model"])
    gdf_rp["rp_source"] = np.where(lowered, "cap", "model")
    classes = np.asarray(rp["classes"], float)
    position = np.searchsorted(classes, gdf_rp["rp"].to_numpy())
    gdf_rp["rp_low"] = classes[np.maximum(position - 1, 0)]
    upper = classes[np.minimum(position + 1, len(classes) - 1)]
    bound = np.nanmin(
        np.column_stack(
            [upper, gdf_rp["rp_cap"], np.full(len(gdf_rp), float(rp["max_rp"]))]
        ),
        axis=1,
    )
    gdf_rp["rp_high"] = np.maximum(gdf_rp["rp"], bound)

    # --- Overrides ---
    if country["overrides"]:
        gdf_ovr = gis.Vector.read_merged(
            country["overrides"], wkt, keep_columns=["rp", "note"]
        )
        gdf_sjoin_ovr_rp = gpd.sjoin(
            gdf_ovr, gdf_rp[["geometry"]], predicate="within", how="inner"
        )
        conflicts = gdf_sjoin_ovr_rp.groupby("index_right")["rp"].nunique()
        if (conflicts > 1).any():
            log.warning(
                "%d protected areas contain override points with different rp; the first is used.",
                int((conflicts > 1).sum()),
            )
        gdf_first = gdf_sjoin_ovr_rp.drop_duplicates("index_right").set_index(
            "index_right"
        )
        gdf_rp.loc[gdf_first.index, ["rp", "rp_low", "rp_high"]] = np.repeat(
            gdf_first["rp"].astype(float).to_numpy()[:, None], 3, axis=1
        )
        gdf_rp.loc[gdf_first.index, "rp_source"] = "override"
        log.info("Overrides applied to %d protected areas.", len(gdf_first))

    for column, value in OBSERVATION_COLUMNS.items():
        gdf_rp[column] = value
    for column in ["rp_model", "rp_cap", "rp", "rp_low", "rp_high"]:
        gdf_rp[column] = gdf_rp[column].astype(float)
    gis.Vector.write_gpkg(gdf_rp, out_dir / files["rp"])
    shares = gdf_rp.loc[scored, "rp_model"].value_counts(normalize=True).sort_index()
    log.info(
        "Model classes (share of protected areas with levees): %s; targets %s.",
        ", ".join(f"{int(k)}: {v:.2f}" for k, v in shares.items()),
        ", ".join(
            f"{c}: {s:.2f}"
            for c, s in zip(rp["distribution_classes"], rp["distribution_shares"])
        ),
    )
    log.info(
        "RP sources: %s. %s written (%d protected areas).",
        ", ".join(f"{k} {v}" for k, v in gdf_rp["rp_source"].value_counts().items()),
        files["rp"],
        len(gdf_rp),
    )
