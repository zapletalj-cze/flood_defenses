# Author: Jakub Zapletal
# Date: 2026-11-05
# Stage E, step E5: sensitivity of the implied protection to level method, built_share_min and subsegment length;
# recall of built area in reference protected areas and wet built area inside implied protected areas per run.

import json
import os
import subprocess
import sys
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from osgeo import gdal

import gis
from plots import INK, Plots

# --- Run variables (no command-line parameters; DEFENSES_<NAME> environment variables override them) ---
CONFIG_PATH = "config/defenses.toml"
COUNTRY = "XX"
BASE_RUN = "base"  # run of implied_protection.py whose curves (base length) every sensitivity run reuses
EVAL_ID = "sensitivity"  # plots in <output_dir>/<country>/plots/implied/<EVAL_ID>
SHOW_PLOTS = False  # show up to implied.plots.show_max figures per type while running (needs a display)
FORCE = False  # recompute runs whose implied_protected_areas.gpkg exists

KEYS = [
    "general.work_dir",
    "general.output_dir",
    "general.log_level",
    "general.n_workers",
    "general.gdal_cachemax_mb",
    "general.warp_memory_limit_mb",
    "general.gdal_threads_per_worker",
    "country.reference_protected_areas",
    "implied.units",
    "implied.recompute_hand",
    "implied.quantiles",
    "implied.subsegment_length_m",
    "implied.subsegment_length_grid",
    "implied.built_share_min_grid",
    *[
        f"implied.plots.{k}"
        for k in (
            "enabled",
            "dpi",
            "max_per_type",
            "sample_seed",
            "show_max",
            "plot_all_in_units",
        )
    ],
]
GFM_KEYS = [
    "sources.gfm_flood_values",
    "sources.gfm_exclusion_values",
    "sources.gfm_reference_water_values",
    "sources.gfm_nodata_values",
]


def evaluate_unit(task):
    # Worker: per unit (core cells), built area of the reference areas covered by explicit and by explicit plus implied
    # protected areas of every run; wet built area per (run, implied protected area, scene).
    start = time.perf_counter()
    unit_dir = Path(task["unit_dir"])
    arrays = {}
    for name in ("implied_built.tif", task["hand_name"]):
        ds = gdal.Open(str(unit_dir / name))
        arrays[name] = ds.GetRasterBand(1).ReadAsArray().astype(np.float32, copy=False)
        gt, wkt = ds.GetGeoTransform(), ds.GetProjection()
        ds = None
    built, hand = arrays["implied_built.tif"], arrays[task["hand_name"]]
    shape, res = built.shape, gt[1]
    bounds = (gt[0], gt[3] - shape[0] * res, gt[0] + shape[1] * res, gt[3])
    core = np.zeros(shape, bool)
    r0 = int(round((gt[3] - task["core"][3]) / res))
    r1 = int(round((gt[3] - task["core"][1]) / res))
    c0 = int(round((task["core"][0] - gt[0]) / res))
    c1 = int(round((task["core"][2] - gt[0]) / res))
    core[max(r0, 0) : r1, max(c0, 0) : c1] = True
    area = np.where(core, built, 0.0) * res * res

    def burn(wkbs, values, dtype):
        if not wkbs:
            return np.zeros(shape, dtype)
        return gis.Raster.rasterize_geometries(
            shapely.from_wkb(wkbs), values, shape, gt, wkt, False, dtype
        )

    ref = burn(task["reference"], np.ones(len(task["reference"])), np.uint8) > 0
    expl = burn(task["explicit"], np.ones(len(task["explicit"])), np.uint8) > 0
    out = {
        "unit_id": task["unit_id"],
        "b_ref": float(area[ref].sum()),
        "b_expl": float(area[ref & expl].sum()),
        "b_cov": {},
        "wet": [],
    }
    labels = {}
    for run_id, polys in task["runs"].items():
        lab = burn(polys, np.arange(1, len(polys) + 1), np.int32)
        labels[run_id] = lab
        out["b_cov"][run_id] = float(area[ref & (expl | (lab > 0))].sum())
    for scene in task["scenes"]:
        _, wet, _ = gis.Evidence.scene_masks(
            scene, bounds, res, wkt, hand, task["gfm"], task["warp_mb"]
        )
        wet_area = np.where(wet, area, 0.0)
        for run_id, lab in labels.items():
            sums = np.bincount(
                lab.ravel(),
                weights=wet_area.ravel(),
                minlength=len(task["runs"][run_id]) + 1,
            )
            for j in np.flatnonzero(sums[1:] > 0):
                out["wet"].append(
                    {
                        "run_id": run_id,
                        "index": int(j),
                        "item_id": scene["item_id"],
                        "wet_built_m2": float(sums[j + 1]),
                    }
                )
    out["summary"] = (
        f"unit {task['unit_id']}: {len(task['runs'])} runs, {len(task['scenes'])} scenes, "
        f"{time.perf_counter() - start:.1f} s"
    )
    out["peak_rss_mb"] = gis.Run.peak_rss_mb()
    return out


if __name__ == "__main__":
    run = gis.Run.env_overrides(
        CONFIG_PATH=CONFIG_PATH,
        COUNTRY=COUNTRY,
        BASE_RUN=BASE_RUN,
        EVAL_ID=EVAL_ID,
        SHOW_PLOTS=SHOW_PLOTS,
        FORCE=FORCE,
    )
    cfg = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], KEYS)
    general, country, imp = cfg["general"], cfg["country"], cfg["implied"]
    work_dir = Path(general["work_dir"]) / run["COUNTRY"]
    out_dir = Path(general["output_dir"]) / run["COUNTRY"]
    log = gis.Run.logger("implied_evaluate", out_dir / "logs", general["log_level"])
    gis.Config.warn_placeholders(cfg)
    start = time.perf_counter()
    base_dir = out_dir / "implied" / run["BASE_RUN"]
    if not (base_dir / "curves.npz").exists():
        raise SystemExit(
            f"implied/{run['BASE_RUN']}/curves.npz missing; run implied_protection.py first."
        )

    # --- Sensitivity grid: each combination is a run of implied_protection.py reusing the base curves ---
    methods = ["onset"] + [f"q{q:.2f}" for q in imp["quantiles"]]
    list_run = [
        {
            "run_id": f"{m}_bs{b:g}_L{length:g}",
            "method": m,
            "built_share_min": b,
            "length_m": length,
        }
        for length in imp["subsegment_length_grid"]
        for b in imp["built_share_min_grid"]
        for m in methods
    ]
    log.info(
        "Sensitivity grid: %d runs (%d methods x %d shares x %d lengths) on %s.",
        len(list_run),
        len(methods),
        len(imp["built_share_min_grid"]),
        len(imp["subsegment_length_grid"]),
        f"pilot units {imp['units']}" if imp["units"] else "all units",
    )
    t = time.perf_counter()
    for r in list_run:
        target = out_dir / "implied" / r["run_id"] / "implied_protected_areas.gpkg"
        if target.exists() and not run["FORCE"]:
            continue
        env = dict(
            os.environ,
            DEFENSES_CONFIG_PATH=str(run["CONFIG_PATH"]),
            DEFENSES_COUNTRY=run["COUNTRY"],
            DEFENSES_RUN_ID=r["run_id"],
            DEFENSES_CURVES_FROM=run["BASE_RUN"],
            DEFENSES_SUBSEGMENT_LENGTH_M=str(r["length_m"]),
            DEFENSES_LEVEL_METHOD=r["method"],
            DEFENSES_BUILT_SHARE_MIN=str(r["built_share_min"]),
            DEFENSES_PLOTS="0",
        )
        proc = subprocess.run(
            [sys.executable, "-m", "processing.implied_protection"],
            env=env,
            capture_output=True,
            text=True,
        )
        (out_dir / "implied" / r["run_id"]).mkdir(parents=True, exist_ok=True)
        (out_dir / "implied" / r["run_id"] / "implied_protection.log").write_text(
            proc.stdout + proc.stderr
        )
        if proc.returncode != 0:
            raise SystemExit(f"Run {r['run_id']} failed:\n{proc.stderr[-2000:]}")
        log.info("Run %s done.", r["run_id"])
    runs_s = time.perf_counter() - t

    # --- Per-unit evaluation of all runs ---
    wkt = gdal.Open(str(work_dir / "dtm" / "dtm.vrt")).GetProjection()
    gdf_ref = (
        gis.Vector.read_merged(country["reference_protected_areas"], wkt)
        if country["reference_protected_areas"]
        else None
    )
    if gdf_ref is None:
        log.warning(
            "No reference_protected_areas: the positive metric is not computed."
        )
    explicit_path = out_dir / "protected_areas_rp.gpkg"
    if not explicit_path.exists():
        explicit_path = out_dir / "protected_areas.gpkg"
    gdf_expl = gpd.read_file(explicit_path).to_crs(wkt)
    gdf_runs = {
        r["run_id"]: gpd.read_file(
            out_dir / "implied" / r["run_id"] / "implied_protected_areas.gpkg"
        ).to_crs(wkt)
        for r in list_run
    }
    scenes_by_unit, gfm = {}, None
    if (out_dir / "gfm_index.csv").exists():
        cfg_gfm = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], GFM_KEYS)
        src = cfg_gfm["sources"]
        gfm = {
            "flood_values": src["gfm_flood_values"],
            "exclusion_values": src["gfm_exclusion_values"],
            "reference_water_values": src["gfm_reference_water_values"],
            "nodata_values": src["gfm_nodata_values"],
        }
        df_idx = pd.read_csv(out_dir / "gfm_index.csv", dtype={"unit_id": str}).fillna(
            ""
        )
        for unit_id, df_scn in df_idx.groupby("unit_id"):
            scenes_by_unit[unit_id] = df_scn.to_dict("records")
    else:
        log.warning(
            "gfm_index.csv missing: the negative metric is not computed (no scenes)."
        )
    log.info(
        "EMS flood extents are not read yet (Stage D5); the negative metric uses GFM scenes only."
    )
    df_units = pd.read_csv(out_dir / "units.csv", dtype={"unit_id": str})
    if imp["units"]:
        df_units = df_units[df_units["unit_id"].isin(imp["units"])]
    df_units = df_units[
        df_units["unit_id"].map(
            lambda u: (work_dir / "units" / u / "implied_built.tif").exists()
        )
    ]
    hand_name = "hand_undef_snapped.tif" if imp["recompute_hand"] else "hand_undef.tif"
    ref_tree = shapely.STRtree(gdf_ref.geometry.values) if gdf_ref is not None else None
    expl_tree = shapely.STRtree(gdf_expl.geometry.values)
    run_trees = {k: shapely.STRtree(g.geometry.values) for k, g in gdf_runs.items()}
    list_task, list_w = [], []
    for u in df_units.itertuples():
        window = shapely.box(u.win_xmin, u.win_ymin, u.win_xmax, u.win_ymax)
        list_task.append(
            {
                "unit_id": u.unit_id,
                "unit_dir": str(work_dir / "units" / u.unit_id),
                "hand_name": hand_name,
                "core": (u.core_xmin, u.core_ymin, u.core_xmax, u.core_ymax),
                "reference": (
                    list(
                        shapely.to_wkb(
                            gdf_ref.geometry.values[
                                ref_tree.query(window, predicate="intersects")
                            ]
                        )
                    )
                    if ref_tree is not None
                    else []
                ),
                "explicit": list(
                    shapely.to_wkb(
                        gdf_expl.geometry.values[
                            expl_tree.query(window, predicate="intersects")
                        ]
                    )
                ),
                # all protected areas of every run, in their order, so the burnt value maps back to the row
                "runs": {
                    k: list(shapely.to_wkb(g.geometry.values))
                    for k, g in gdf_runs.items()
                },
                "scenes": scenes_by_unit.get(u.unit_id, []) if gfm else [],
                "gfm": gfm,
                "warp_mb": general["warp_memory_limit_mb"],
            }
        )
        list_w.append(float(u.est_window_cells))
    t = time.perf_counter()
    list_res, list_fail = gis.Run.run_units(
        evaluate_unit,
        list_task,
        general["n_workers"],
        general["gdal_cachemax_mb"],
        general["gdal_threads_per_worker"],
        list_w,
        log,
    )
    eval_s = time.perf_counter() - t
    if list_fail:
        raise SystemExit(f"Failed units: {list_fail}")

    # --- Metrics per run ---
    # one event RP per test, combined over the sources by source_mode (Stage F)
    event = (
        pd.read_csv(out_dir / "event_rp_combined.csv")
        if (out_dir / "event_rp_combined.csv").exists()
        else None
    )
    if event is None:
        log.warning(
            "event_rp_combined.csv (Stage F) missing: every wet built cell inside implied protected areas counts as potential."
        )
    df_wet = pd.DataFrame(
        [w for r in list_res for w in r["wet"]],
        columns=["run_id", "index", "item_id", "wet_built_m2"],
    )
    df_wet = df_wet.groupby(["run_id", "index", "item_id"], as_index=False)[
        "wet_built_m2"
    ].sum()
    list_test, list_row = [], []
    b_ref = sum(r["b_ref"] for r in list_res)
    b_expl = sum(r["b_expl"] for r in list_res)
    for r in list_run:
        gdf_r = gdf_runs[r["run_id"]]
        b_cov = sum(res_u["b_cov"][r["run_id"]] for res_u in list_res)
        wet_r = df_wet[df_wet["run_id"] == r["run_id"]]
        tot = {"error": 0.0, "ok": 0.0, "potential": 0.0}
        for w in wet_r.itertuples():
            pa = gdf_r.iloc[w.index]
            rp_event = np.nan
            if event is not None:
                hit = event[
                    (event["protected_area_id"] == pa["implied_id"])
                    & (event["item_id"] == w.item_id)
                ]
                if hit.empty:
                    # an event RP is a property of the flood at the river: same river and scene
                    river = event["protected_area_id"].astype(str).str.split("_").str[1]
                    hit = event[
                        event["protected_area_id"]
                        .astype(str)
                        .str.startswith("implied_")
                        & (river == str(pa["river_id"]))
                        & (event["item_id"] == w.item_id)
                    ]
                if hit["rp_combined"].notna().any():
                    rp_event = float(hit["rp_combined"].min())
            status = (
                "potential"
                if np.isnan(rp_event)
                else ("error" if rp_event <= pa["rp"] else "ok")
            )
            tot[status] += w.wet_built_m2
            list_test.append(
                {
                    "run_id": r["run_id"],
                    "implied_id": pa["implied_id"],
                    "item_id": w.item_id,
                    "wet_built_m2": w.wet_built_m2,
                    "rp": pa["rp"],
                    "rp_event": rp_event,
                    "status": status,
                }
            )
        wet_all = sum(tot.values())
        area = gdf_r.area.sum() if len(gdf_r) else 0.0
        list_row.append(
            {
                **r,
                "implied_protected_areas": len(gdf_r),
                "implied_area_km2": area / 1e6,
                "recall_built": b_cov / b_ref if b_ref > 0 else np.nan,
                "recall_built_explicit": b_expl / b_ref if b_ref > 0 else np.nan,
                "implied_gain": (b_cov - b_expl) / b_ref if b_ref > 0 else np.nan,
                "wet_built_m2": wet_all,
                "wrongly_protected_share": (
                    tot["error"] / wet_all if wet_all > 0 else np.nan
                ),
                "potential_share": (
                    tot["potential"] / wet_all if wet_all > 0 else np.nan
                ),
                "historic_low_core_area_share": (
                    gdf_r.area[gdf_r["historic_low_core"].astype(bool)].sum() / area
                    if area > 0
                    else np.nan
                ),
                **{
                    f"confidence_{c}": (
                        int((gdf_r["confidence"] == c).sum()) if len(gdf_r) else 0
                    )
                    for c in ("high", "medium", "low")
                },
                "by_context": (
                    json.dumps(gdf_r["context"].value_counts().to_dict())
                    if len(gdf_r)
                    else "{}"
                ),
            }
        )
    df_sens = pd.DataFrame(list_row)
    df_sens.to_csv(out_dir / "sensitivity.csv", index=False)
    pd.DataFrame(
        list_test,
        columns=[
            "run_id",
            "implied_id",
            "item_id",
            "wet_built_m2",
            "rp",
            "rp_event",
            "status",
        ],
    ).to_csv(out_dir / "sensitivity_tests.csv", index=False)
    log.info(
        "sensitivity.csv: %d runs; recall_built %s; wrongly_protected_share %s.",
        len(df_sens),
        (
            f"{df_sens['recall_built'].min():.3f}-{df_sens['recall_built'].max():.3f}"
            if b_ref > 0
            else "n/a"
        ),
        (
            f"{df_sens['wrongly_protected_share'].min():.3f}-{df_sens['wrongly_protected_share'].max():.3f}"
            if df_sens["wrongly_protected_share"].notna().any()
            else "n/a"
        ),
    )

    # --- Plots: recall_built against wrongly_protected_share, small multiples by subsegment length ---
    t = time.perf_counter()
    if imp["plots"]["enabled"]:
        plots = Plots(
            out_dir / "plots" / "implied" / run["EVAL_ID"],
            imp["plots"],
            run["SHOW_PLOTS"],
            log,
        )
        plt = plots.plt
        styles = Plots.method_styles(methods)
        plots.select("sensitivity", ["sensitivity"])
        lengths = list(imp["subsegment_length_grid"])
        fig, axes = plt.subplots(
            1,
            len(lengths),
            figsize=(4.2 * len(lengths) + 1, 4.6),
            sharey=True,
            squeeze=False,
            layout="constrained",
        )
        y_metric = "recall_built"
        x_metric = (
            "wrongly_protected_share"
            if df_sens["wrongly_protected_share"].notna().any()
            else "potential_share"
        )
        for ax, length in zip(axes[0], lengths):
            sub = df_sens[df_sens["length_m"] == length]
            for m in methods:
                sm = sub[sub["method"] == m]
                ax.plot(
                    sm[x_metric],
                    sm[y_metric],
                    "o",
                    color=styles[m]["color"],
                    ms=7,
                    label=m,
                )
                for row in sm.itertuples():
                    ax.annotate(
                        f"{row.built_share_min:g}",
                        (getattr(row, x_metric), getattr(row, y_metric)),
                        fontsize=7,
                        color=INK["secondary"],
                        xytext=(4, 3),
                        textcoords="offset points",
                    )
            ax.set_title(f"Subsegments of {length:g} m")
            ax.set_xlabel(x_metric.replace("_", " "))
        axes[0, 0].set_ylabel("recall built (explicit + implied)")
        handles, labels_ = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels_, loc="outside lower center", ncols=len(methods))
        fig.suptitle(
            f"sensitivity · {len(df_sens)} runs · points labelled with built_share_min",
            x=0.01,
            ha="left",
        )
        plots.save(fig, "sensitivity", "sensitivity", {"runs": len(df_sens)})
        log.info("Plots: %s", plots.finish())
    perf = {
        "runs": len(list_run),
        "runs_s": round(runs_s, 2),
        "evaluation_s": round(eval_s, 2),
        "plots_s": round(time.perf_counter() - t, 2),
        "total_s": round(time.perf_counter() - start, 2),
        "peak_rss_mb_worker": (
            round(max(r["peak_rss_mb"] for r in list_res)) if list_res else 0
        ),
    }
    (out_dir / "implied_evaluate_performance.json").write_text(
        json.dumps(perf, indent=1)
    )
    log.info("Timings: %s", ", ".join(f"{k} {v}" for k, v in perf.items()))
