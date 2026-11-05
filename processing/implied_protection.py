# Author: Jakub Zapletal
# Date: 2026-11-05
# Stage E, steps E3-E4: implied protection in dense built-up reaches. Built-up fraction, HAND curves and levels per
# bank subsegment (A-B); implied cells, merging, precedence, temporal signal, corroboration and RP (C-F).

import json
import shutil
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from osgeo import gdal

import gis
from plots import COLOURS, INK, Plots

# --- Run variables (no command-line parameters; DEFENSES_<NAME> environment variables override them) ---
CONFIG_PATH = "config/defenses.toml"
COUNTRY = "XX"
RUN_ID = "base"  # results in <output_dir>/<country>/implied/<RUN_ID>, plots in .../plots/implied/<RUN_ID>
SHOW_PLOTS = False  # show up to implied.plots.show_max figures per type while running (needs a display)

KEYS = [
    "general.work_dir",
    "general.output_dir",
    "general.log_level",
    "general.n_workers",
    "general.gdal_cachemax_mb",
    "general.warp_memory_limit_mb",
    "general.gdal_threads_per_worker",
    "country.building_footprints",
    "country.working_resolution_m",
    "country.context.order",
    *[
        f"implied.{k}"
        for k in (
            "units",
            "recompute_hand",
            "envelope_hand_m",
            "subsegment_length_m",
            "built_source",
            "exclude_building_tags",
            "built_share_min",
            "built_area_min_m2",
            "hand_bin_m",
            "min_bin_area_m2",
            "smooth_bins",
            "onset_fraction",
            "onset_persist_bins",
            "quantiles",
            "primary_level_method",
            "sharp_window_m",
            "sharp_high",
            "ghsl_epochs",
        )
    ],
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
TILE = 1024  # working cells per side of a footprint tile (4096 fine cells at a quarter of the cell size)


def curves_unit(task):
    # Worker, pass 1: built-up fraction of the unit grid and envelope / built area per (subsegment, HAND bin).
    start = time.perf_counter()
    unit_dir = Path(task["unit_dir"])
    secs = {}
    t = time.perf_counter()
    ds = gdal.Open(str(unit_dir / "subsegments.tif"))
    gt, wkt = ds.GetGeoTransform(), ds.GetProjection()
    labels = ds.GetRasterBand(1).ReadAsArray()
    ds = None
    shape, res = labels.shape, gt[1]
    ds = gdal.Open(str(unit_dir / task["hand_name"]))
    hand = ds.GetRasterBand(1).ReadAsArray().astype(np.float32)
    ds = None
    df_lut = pd.read_csv(unit_dir / "subsegments.csv")
    n_lab = len(df_lut) + 1
    secs["read"] = time.perf_counter() - t

    # --- A. Built-up fraction ---
    t = time.perf_counter()
    rows_any = np.flatnonzero((labels > 0).any(axis=1))
    cols_any = np.flatnonzero((labels > 0).any(axis=0))
    built_fp = np.zeros(shape, np.float32)
    if task["buildings"] and len(rows_any):
        geoms = shapely.from_wkb(task["buildings"])
        tree = shapely.STRtree(geoms)
        fine = res / 4
        for r0 in range(rows_any[0], rows_any[-1] + 1, TILE):
            for c0 in range(cols_any[0], cols_any[-1] + 1, TILE):
                r1, c1 = min(r0 + TILE, shape[0]), min(c0 + TILE, shape[1])
                if not (labels[r0:r1, c0:c1] > 0).any():
                    continue
                x0, y1 = gt[0] + c0 * res, gt[3] - r0 * res
                tile_box = shapely.box(x0, gt[3] - r1 * res, gt[0] + c1 * res, y1)
                hits = tree.query(tile_box, predicate="intersects")
                if not len(hits):
                    continue
                burnt = gis.Raster.rasterize_geometries(
                    geoms[hits],
                    np.ones(len(hits)),
                    ((r1 - r0) * 4, (c1 - c0) * 4),
                    (x0, fine, 0.0, y1, 0.0, -fine),
                    wkt,
                    False,
                    np.uint8,
                )
                built_fp[r0:r1, c0:c1] = burnt.reshape(r1 - r0, 4, c1 - c0, 4).mean(
                    axis=(1, 3)
                )
    bounds = (gt[0], gt[3] - shape[0] * res, gt[0] + shape[1] * res, gt[3])
    list_ghsl = []
    for path in [task["ghsl_fraction"]] + task["ghsl_epochs"]:
        if not path:
            list_ghsl.append(None)
            continue
        source = gdal.Open(path)
        sgt = source.GetGeoTransform()
        warped, _ = gis.Raster.read_window(
            source, bounds, res, "average", task["warp_mb"], dst_wkt=wkt
        )
        # built-up m2 per source cell, averaged onto the grid, as a share of the source cell
        list_ghsl.append(
            np.clip(np.nan_to_num(warped, nan=0.0) / abs(sgt[1] * sgt[5]), 0.0, 1.0)
        )
    if task["built_source"] == "footprints":
        built = built_fp
    elif task["built_source"] == "ghsl":
        built = list_ghsl[0]
    else:
        built = built_fp if list_ghsl[0] is None else np.maximum(built_fp, list_ghsl[0])
    built = built.astype(np.float32)
    gis.Raster.write_array(unit_dir / "implied_built.tif", built, gt, wkt)
    secs["built"] = time.perf_counter() - t

    # --- B. Envelope and built area per (subsegment, HAND bin), per GHSL epoch, and context ---
    t = time.perf_counter()
    n_bins = task["n_bins"]
    inside = labels > 0
    lab = labels[inside].astype(np.int64)
    b = np.clip(
        np.floor(np.nan_to_num(hand[inside], nan=0.0) / task["hand_bin_m"]),
        0,
        n_bins - 1,
    ).astype(np.int64)
    key = lab * n_bins + b
    size = n_lab * n_bins
    cell = res * res
    env = np.bincount(key, minlength=size).reshape(n_lab, n_bins) * cell
    blt = (np.bincount(key, weights=built[inside], minlength=size) * cell).reshape(
        n_lab, n_bins
    )
    epochs = (
        np.stack(
            [
                (
                    np.bincount(key, weights=g[inside], minlength=size) * cell
                    if g is not None
                    else np.zeros(size)
                ).reshape(n_lab, n_bins)
                for g in list_ghsl[1:]
            ]
        )
        if task["ghsl_epochs"]
        else np.zeros((0, n_lab, n_bins))
    )
    n_ctx = len(task["contexts"]) + 1
    if task["zones"]:
        gdf_zone = gpd.GeoDataFrame(
            {"context": task["zone_context"]},
            geometry=shapely.from_wkb(task["zones"]),
            crs=wkt,
        )
        code = gis.Vector.zones_to_grid(gdf_zone, task["contexts"], shape, gt, wkt)
        ctx = (
            np.bincount(lab * n_ctx + code[inside], minlength=n_lab * n_ctx).reshape(
                n_lab, n_ctx
            )
            * cell
        )
    else:
        ctx = np.zeros((n_lab, n_ctx))
    secs["curves"] = time.perf_counter() - t
    ids = ["none"] + [
        f"{r}_{s}_{k}"
        for r, s, k in zip(df_lut["river_id"], df_lut["side"], df_lut["k"])
    ]
    return {
        "unit_id": task["unit_id"],
        "summary": (
            f"unit {task['unit_id']}: {n_lab - 1} subsegments, built "
            f"{blt[1:].sum() / 1e6:.3f} km2 of {env[1:].sum() / 1e6:.2f} km2 envelope, "
            f"{time.perf_counter() - start:.1f} s"
        ),
        "ids": ids[1:],
        "env": env[1:],
        "built": blt[1:],
        "epochs": epochs[:, 1:],
        "ctx": ctx[1:],
        "secs": secs,
        "peak_rss_mb": gis.Run.peak_rss_mb(),
    }


if __name__ == "__main__":
    run = gis.Run.env_overrides(
        CONFIG_PATH=CONFIG_PATH, COUNTRY=COUNTRY, RUN_ID=RUN_ID, SHOW_PLOTS=SHOW_PLOTS
    )
    cfg = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], KEYS)
    general, country, imp = cfg["general"], cfg["country"], cfg["implied"]
    work_dir = Path(general["work_dir"]) / run["COUNTRY"]
    out_dir = Path(general["output_dir"]) / run["COUNTRY"]
    run_dir = out_dir / "implied" / run["RUN_ID"]
    run_dir.mkdir(parents=True, exist_ok=True)
    log = gis.Run.logger("implied_protection", out_dir / "logs", general["log_level"])
    gis.Config.warn_placeholders(cfg)
    gis.Run.configure_gdal(
        general["gdal_cachemax_mb"], general["gdal_threads_per_worker"]
    )
    shutil.copy(run["CONFIG_PATH"], run_dir / "run_config.toml")
    start = time.perf_counter()

    # --- Inputs: required stops, optional ones switch features off ---
    seg_path = out_dir / "subsegments.gpkg"
    if not seg_path.exists() or not (out_dir / "units.csv").exists():
        raise SystemExit(
            "subsegments.gpkg or units.csv missing; run implied_segments.py first."
        )
    df_units = pd.read_csv(out_dir / "units.csv", dtype={"unit_id": str})
    if imp["units"]:
        df_units = df_units[df_units["unit_id"].isin(imp["units"])]
    df_units = df_units[
        df_units["unit_id"].map(
            lambda u: (work_dir / "units" / u / "subsegments.csv").exists()
        )
    ]
    if df_units.empty:
        raise SystemExit("No unit with subsegments; run implied_segments.py first.")
    hand_name = "hand_undef_snapped.tif" if imp["recompute_hand"] else "hand_undef.tif"
    wkt = gdal.Open(str(work_dir / "dtm" / "dtm.vrt")).GetProjection()
    list_bld = []
    if country["building_footprints"]:
        gdf_bld = gis.Vector.read_merged(
            country["building_footprints"], wkt, keep_columns=True
        )
        list_bld.append(("footprints", gdf_bld))
    if (out_dir / "osm_buildings.gpkg").exists():
        gis.Vector.check_single_layer(out_dir / "osm_buildings.gpkg")
        list_bld.append(
            ("OSM", gpd.read_file(out_dir / "osm_buildings.gpkg").to_crs(wkt))
        )
    list_geom = []
    for label, gdf in list_bld:
        drop = np.zeros(len(gdf), bool)
        for tag, values in imp["exclude_building_tags"].items():
            if tag not in gdf.columns:
                log.info(
                    "%s buildings have no %s attribute; its exclusions do not apply.",
                    label,
                    tag,
                )
                continue
            hit = gdf[tag].astype(str).isin([str(v) for v in values]).to_numpy()
            log.info(
                "%s buildings: %d with %s in %s excluded.",
                label,
                int(hit.sum()),
                tag,
                values,
            )
            drop |= hit
        list_geom.append(gdf.geometry.values[~drop])
        log.info("%s buildings: %d used.", label, int((~drop).sum()))
    buildings = (
        shapely.make_valid(np.concatenate(list_geom))
        if list_geom
        else np.empty(0, object)
    )
    ghsl_fraction = ""
    for name in ["ghsl_built_s_10m.vrt"] + [
        f"ghsl_built_s_{e}.vrt" for e in sorted(imp["ghsl_epochs"], reverse=True)
    ]:
        if (out_dir / name).exists():
            ghsl_fraction = str(out_dir / name)
            break
    ghsl_epochs = [
        str(out_dir / f"ghsl_built_s_{e}.vrt") for e in sorted(imp["ghsl_epochs"])
    ]
    missing_epochs = [p for p in ghsl_epochs if not Path(p).exists()]
    if missing_epochs or not ghsl_epochs:
        log.warning(
            "GHSL built-up surface per epoch missing (%s): the temporal signal is switched off.",
            [Path(p).name for p in missing_epochs] or "implied.ghsl_epochs is empty",
        )
        ghsl_epochs = []
    if imp["built_source"] in ("footprints", "max") and not len(buildings):
        if imp["built_source"] == "footprints":
            raise SystemExit(
                "built_source = footprints, but no building footprints or OSM buildings exist."
            )
        log.warning(
            "No buildings: built_source = max uses the GHSL built-up surface only."
        )
    if imp["built_source"] in ("ghsl", "max") and not ghsl_fraction:
        if imp["built_source"] == "ghsl":
            raise SystemExit(
                "built_source = ghsl, but no GHSL built-up surface exists (download_ghsl.py)."
            )
        log.warning(
            "No GHSL built-up surface: built_source = max uses the footprints only."
        )
    zones_path = out_dir / "zones.gpkg"
    gdf_zone = gpd.read_file(zones_path).to_crs(wkt) if zones_path.exists() else None
    if gdf_zone is None:
        log.warning(
            "zones.gpkg missing: no context labels; the rural score weights apply."
        )

    # --- Pass 1 per unit: built-up fraction and curves ---
    n_bins = int(np.ceil(imp["envelope_hand_m"] / imp["hand_bin_m"] - 1e-9))
    edges = np.arange(n_bins + 1) * imp["hand_bin_m"]
    bld_tree = shapely.STRtree(buildings)
    zone_tree = (
        shapely.STRtree(gdf_zone.geometry.values) if gdf_zone is not None else None
    )
    list_task, list_w = [], []
    for u in df_units.itertuples():
        core = shapely.box(u.core_xmin, u.core_ymin, u.core_xmax, u.core_ymax)
        sel_z = (
            zone_tree.query(core, predicate="intersects")
            if zone_tree is not None
            else []
        )
        list_task.append(
            {
                "unit_id": u.unit_id,
                "unit_dir": str(work_dir / "units" / u.unit_id),
                "hand_name": hand_name,
                "buildings": list(
                    shapely.to_wkb(
                        buildings[bld_tree.query(core, predicate="intersects")]
                    )
                ),
                "ghsl_fraction": (
                    ghsl_fraction if imp["built_source"] != "footprints" else ""
                ),
                "ghsl_epochs": ghsl_epochs,
                "built_source": imp["built_source"],
                "zones": (
                    list(shapely.to_wkb(gdf_zone.geometry.values[sel_z]))
                    if len(sel_z)
                    else []
                ),
                "zone_context": (
                    list(gdf_zone["context"].to_numpy()[sel_z]) if len(sel_z) else []
                ),
                "contexts": country["context"]["order"],
                "n_bins": n_bins,
                "hand_bin_m": imp["hand_bin_m"],
                "warp_mb": general["warp_memory_limit_mb"],
            }
        )
        list_w.append(float(u.est_window_cells))
    t = time.perf_counter()
    list_res, list_fail = gis.Run.run_units(
        curves_unit,
        list_task,
        general["n_workers"],
        general["gdal_cachemax_mb"],
        general["gdal_threads_per_worker"],
        list_w,
        log,
    )
    pass1_s = time.perf_counter() - t
    if list_fail:
        raise SystemExit(f"Failed units in pass 1: {list_fail}")

    # --- B. Curves summed over units per subsegment, levels and eligibility ---
    t = time.perf_counter()
    gdf_seg = gpd.read_file(seg_path)
    pos = {s: i for i, s in enumerate(gdf_seg["subsegment_id"])}
    n_seg, n_ep = len(gdf_seg), len(ghsl_epochs)
    env = np.zeros((n_seg, n_bins))
    blt = np.zeros((n_seg, n_bins))
    ep = np.zeros((n_ep, n_seg, n_bins))
    ctx = np.zeros((n_seg, len(country["context"]["order"]) + 1))
    for res_u in list_res:
        idx = np.array([pos[i] for i in res_u["ids"]], dtype=np.int64)
        np.add.at(env, idx, res_u["env"])
        np.add.at(blt, idx, res_u["built"])
        np.add.at(ctx, idx, res_u["ctx"])
        if n_ep:
            for e in range(n_ep):
                np.add.at(ep[e], idx, res_u["epochs"][e])
    share = np.where(
        env >= imp["min_bin_area_m2"],
        np.divide(blt, env, out=np.zeros_like(blt), where=env > 0),
        np.nan,
    )
    onset_bin, smoothed = gis.Stats.onset(
        share, imp["smooth_bins"], imp["onset_fraction"], imp["onset_persist_bins"]
    )
    window = max(1, int(round(imp["sharp_window_m"] / imp["hand_bin_m"])))
    sharp = gis.Stats.sharpness(smoothed, onset_bin, window)
    quant = gis.Stats.binned_quantiles(blt, edges, imp["quantiles"])
    methods = ["onset"] + [f"q{q:.2f}" for q in imp["quantiles"]]
    styles = Plots.method_styles(methods)
    levels = np.c_[
        np.where(onset_bin >= 0, edges[np.maximum(onset_bin, 0)], np.nan), quant
    ]
    built_total, env_total = blt.sum(axis=1), env.sum(axis=1)
    built_share = np.divide(
        built_total, env_total, out=np.zeros(n_seg), where=env_total > 0
    )
    primary = levels[:, methods.index(imp["primary_level_method"])]
    reason = np.where(
        built_share < imp["built_share_min"],
        "built share below built_share_min",
        np.where(
            built_total < imp["built_area_min_m2"],
            "built area below built_area_min_m2",
            np.where(
                onset_bin < 0, "no onset", np.where(np.isnan(primary), "no level", "")
            ),
        ),
    )
    eligible = reason == ""
    contexts = ["none"] + list(country["context"]["order"])
    context = np.where(
        ctx.sum(axis=1) > 0, np.array(contexts)[np.argmax(ctx, axis=1)], "none"
    )
    gdf_seg["built_area_m2"] = built_total
    gdf_seg["envelope_area_m2"] = env_total
    gdf_seg["built_share_total"] = built_share
    gdf_seg["eligible"] = eligible
    gdf_seg["reason"] = reason
    gdf_seg["sharpness"] = sharp
    for j, m in enumerate(methods):
        gdf_seg[f"level_{m.replace('.', '_')}_m"] = levels[:, j]
    gdf_seg["level_m"] = primary
    gdf_seg["context"] = context
    gis.Vector.write_gpkg(gdf_seg, run_dir / "implied_subsegments.gpkg")
    np.savez_compressed(
        run_dir / "curves.npz",
        subsegment_id=gdf_seg["subsegment_id"].to_numpy().astype(str),
        edges=edges,
        envelope_area=env.astype(np.float32),
        built_area=blt.astype(np.float32),
        epochs=np.array(sorted(imp["ghsl_epochs"]) if n_ep else [], dtype=np.int64),
        built_area_epochs=ep.astype(np.float32),
    )
    curves_s = time.perf_counter() - t
    log.info(
        "implied_subsegments.gpkg: %d subsegments, %d eligible (%s); primary level %s: median %.2f m.",
        n_seg,
        int(eligible.sum()),
        ", ".join(
            f"{r or 'eligible'} {int((reason == r).sum())}" for r in np.unique(reason)
        ),
        imp["primary_level_method"],
        float(np.nanmedian(primary[eligible])) if eligible.any() else float("nan"),
    )

    # --- Plots: curve per subsegment, summary ---
    t = time.perf_counter()
    if imp["plots"]["enabled"]:
        plots = Plots(
            out_dir / "plots" / "implied" / run["RUN_ID"],
            imp["plots"],
            run["SHOW_PLOTS"],
            log,
        )
        plt = plots.plt
        sid = gdf_seg["subsegment_id"].to_numpy()
        gap = np.abs(
            levels[:, 0] - levels[:, methods.index(f"q{min(imp['quantiles']):.2f}")]
        )
        flagged = sid[
            (eligible & (sharp < imp["sharp_high"]))
            | (eligible & (gap > imp["sharp_window_m"]))
            | (~eligible & (onset_bin >= 0))
        ]
        units_of = {}
        if imp["plots"]["plot_all_in_units"]:
            for res_u in list_res:
                for i in res_u["ids"]:
                    units_of.setdefault(i, []).append(res_u["unit_id"])
        centres = (edges[:-1] + edges[1:]) / 2
        for name in plots.select("curve", list(sid), flagged, units_of):
            i = pos[name]
            fig, (ax1, ax2) = plt.subplots(
                2,
                1,
                figsize=(9, 7.5),
                sharex=True,
                height_ratios=[1.6, 1],
                layout="constrained",
            )
            ax1.bar(
                edges[:-1],
                env[i] / 1e4,
                width=imp["hand_bin_m"],
                align="edge",
                color=INK["grid"],
                label="envelope area",
            )
            ax1.set_ylabel("Envelope area per bin (ha)")
            ax1b = ax1.twinx()
            ax1b.plot(
                centres,
                share[i],
                color=INK["muted"],
                lw=0.9,
                marker="o",
                ms=2.5,
                label="built share",
            )
            ax1b.plot(
                centres, smoothed[i], color=INK["primary"], label="smoothed built share"
            )
            ax1b.set_ylabel("Built share")
            ax1b.set_ylim(0, max(1e-3, float(np.nanmax(np.r_[share[i], 0.0])) * 1.15))
            ax1b.grid(False)
            ax1b.spines["right"].set_visible(True)
            cum = np.cumsum(blt[i]) / max(built_total[i], 1e-9)
            ax2.plot(
                edges[1:], cum, color=INK["secondary"], label="cumulative built area"
            )
            ax2.set_ylim(0, 1.02)
            ax2.set_xlabel("HAND (m)")
            ax2.set_ylabel("Cumulative built fraction")
            for j, m in enumerate(methods):
                if np.isfinite(levels[i, j]):
                    for ax in (ax1, ax2):
                        ax.axvline(
                            levels[i, j],
                            **styles[m],
                            label=f"{m} {levels[i, j]:.2f} m" if ax is ax1 else None,
                        )
            handles = sum((a.get_legend_handles_labels()[0] for a in (ax1, ax1b)), [])
            labels_ = sum((a.get_legend_handles_labels()[1] for a in (ax1, ax1b)), [])
            ax2.legend(handles, labels_, loc="lower right", ncols=2, fontsize=8)
            ax1.set_title(
                f"curve · subsegment {name} · run {run['RUN_ID']}\n"
                f"built share {built_share[i]:.2f}, sharpness {sharp[i]:.2f}, context {context[i]}, "
                f"{'eligible' if eligible[i] else 'not eligible: ' + reason[i]}"
            )
            plots.save(
                fig,
                "curve",
                name,
                {
                    "eligible": bool(eligible[i]),
                    "reason": reason[i],
                    "built_share_total": float(built_share[i]),
                    "sharpness": float(sharp[i]),
                    "context": context[i],
                    **{
                        f"level_{m}": float(levels[i, j]) for j, m in enumerate(methods)
                    },
                },
            )
        plots.select("summary", ["summary"])
        fig, axes = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")
        order = [c for c in contexts if (context == c).any()]
        ctx_colours = [COLOURS["context"][c] for c in order]
        for ax, values, title, unit in (
            (
                axes[0, 0],
                primary,
                f"Primary level ({imp['primary_level_method']}), eligible",
                "m",
            ),
            (axes[0, 1], sharp, "Sharpness, eligible", ""),
        ):
            data = [values[eligible & (context == c)] for c in order]
            data = [d[np.isfinite(d)] for d in data]
            if sum(len(d) for d in data):
                ax.hist(
                    data,
                    bins=20,
                    stacked=True,
                    color=ctx_colours,
                    label=order,
                    edgecolor=INK["surface"],
                    linewidth=0.5,
                )
            ax.set_title(title)
            ax.set_xlabel(f"{title.split(',')[0]}{' (' + unit + ')' if unit else ''}")
            ax.set_ylabel("Subsegments")
        data = [built_share[context == c] for c in order]
        axes[1, 0].hist(
            data,
            bins=20,
            stacked=True,
            color=ctx_colours,
            label=order,
            edgecolor=INK["surface"],
            linewidth=0.5,
        )
        axes[1, 0].axvline(
            imp["built_share_min"],
            color=INK["primary"],
            lw=0.8,
            ls="--",
            label="built_share_min",
        )
        axes[1, 0].set_title("Built share of the envelope, all subsegments")
        axes[1, 0].set_xlabel("Built share")
        axes[1, 0].set_ylabel("Subsegments")
        counts = pd.crosstab(
            pd.Series(context, name="context"),
            pd.Series(np.where(eligible, "eligible", "not eligible"), name="status"),
        )
        counts = counts.reindex(order)
        x = np.arange(len(order))
        bottom = np.zeros(len(order))
        for status, colour in (
            ("eligible", INK["primary"]),
            ("not eligible", INK["grid"]),
        ):
            vals = (
                counts[status].to_numpy() if status in counts else np.zeros(len(order))
            )
            axes[1, 1].bar(
                x, vals, bottom=bottom, color=colour, width=0.6, label=status
            )
            bottom += vals
        axes[1, 1].set_xticks(x, order)
        axes[1, 1].set_title("Subsegments per context")
        axes[1, 1].set_ylabel("Subsegments")
        axes[0, 0].legend(fontsize=8)
        axes[1, 1].legend(fontsize=8)
        fig.suptitle(
            f"summary · run {run['RUN_ID']} · {n_seg} subsegments, {int(eligible.sum())} eligible",
            x=0.01,
            ha="left",
        )
        plots.save(
            fig,
            "summary",
            "summary",
            {
                "subsegments": n_seg,
                "eligible": int(eligible.sum()),
                "median_level_m": (
                    float(np.nanmedian(primary[eligible]))
                    if eligible.any()
                    else float("nan")
                ),
            },
        )
        log.info("Plots: %s", plots.finish())
    plots_s = time.perf_counter() - t
    perf = {
        "units": len(list_task),
        "subsegments": n_seg,
        **{
            f"pass1_{k}_s": round(sum(r["secs"][k] for r in list_res), 2)
            for k in ("read", "built", "curves")
        },
        "pass1_wall_s": round(pass1_s, 2),
        "levels_s": round(curves_s, 2),
        "plots_s": round(plots_s, 2),
        "total_s": round(time.perf_counter() - start, 2),
        "peak_rss_mb_worker": round(max(r["peak_rss_mb"] for r in list_res)),
        "peak_rss_mb_main": round(gis.Run.peak_rss_mb()),
    }
    (run_dir / "implied_protection_performance.json").write_text(
        json.dumps(perf, indent=1)
    )
    log.info("Timings: %s", ", ".join(f"{k} {v}" for k, v in perf.items()))
