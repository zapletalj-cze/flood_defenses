# Author: Jakub Zapletal
# Date: 2026-11-13
# Stage H step 3: composites of the GFM scenes of every flood event. For an event (events.gpkg) and every unit its
# river group crosses, the scene windows of gfm_index.csv inside the event window are brought to the unit grid
# (observations only: excluded cells and reference water are not valid) and combined per cell: maximum extent,
# number of valid and of wet observations, first and last wet day. The per-scene masks stay the clipped scene
# layers themselves. Outputs: rasters in work/<country>/events/<event_id>/, event_composites.csv, the comparison
# with official maximum flood extent files (evidence.max_extent_files) and maps in plots/events/<RUN_ID>.
# Run: pixi run event-composites

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
from plots import CATEGORICAL, INK, Plots

# --- Run variables (no command-line parameters; DEFENSES_<NAME> environment variables override them) ---
CONFIG_PATH = "config/defenses.toml"
COUNTRY = "XX"
FORCE = False  # recompute composites that already exist
UNITS = []  # unit ids; empty = all units with scenes
EVENTS = []  # event ids of events.gpkg; empty = all events with scenes
RUN_ID = "base"  # plots in <output_dir>/<country>/plots/events/<RUN_ID>
SHOW_PLOTS = False  # show up to implied.plots.show_max figures per type while running (needs a display)

KEYS = [
    "general.work_dir",
    "general.output_dir",
    "general.n_workers",
    "general.gdal_cachemax_mb",
    "general.warp_memory_limit_mb",
    "general.gdal_threads_per_worker",
    "general.log_level",
    "country.corridor_buffer_m",
    "evidence.use_advisory_flags",
    "evidence.max_extent_files",
    "sources.gfm_flood_values",
    "sources.gfm_exclusion_values",
    "sources.gfm_reference_water_values",
    "sources.gfm_nodata_values",
    "sources.gfm_advisory_values",
    "implied.plots.enabled",
    "implied.plots.dpi",
    "implied.plots.max_per_type",
    "implied.plots.sample_seed",
    "implied.plots.show_max",
    "implied.plots.plot_all_in_units",
]
LAYERS = ["max_extent", "n_valid", "n_wet", "first_wet", "last_wet"]


def composite_unit(task):
    start = time.perf_counter()
    unit_dir, ev_dir = Path(task["unit_dir"]), Path(task["event_dir"])
    ev_dir.mkdir(parents=True, exist_ok=True)
    hand_ds = gdal.Open(str(unit_dir / "hand_undef.tif"))
    gt, wkt = hand_ds.GetGeoTransform(), hand_ds.GetProjection()
    hand = hand_ds.ReadAsArray().astype(np.float32)
    nodata = hand_ds.GetRasterBand(1).GetNoDataValue()
    if nodata is not None:
        hand[hand == nodata] = np.nan
    shape, res = hand.shape, gt[1]
    bounds = (gt[0], gt[3] - shape[0] * res, gt[0] + shape[1] * res, gt[3])
    n_valid = np.zeros(shape, np.uint16)
    n_wet = np.zeros(shape, np.uint16)
    first = np.full(shape, -1, np.int16)
    last = np.full(shape, -1, np.int16)
    origin = np.datetime64(task["window_start"])
    for scene in task["scenes"]:
        valid, wet, _ = gis.Evidence.scene_masks(
            scene, bounds, res, wkt, hand, task["params"], task["warp_memory_mb"]
        )
        day = int((np.datetime64(scene["datetime"][:10]) - origin).astype(int))
        n_valid += valid
        n_wet += wet
        first[wet & (first < 0)] = day
        last[wet] = np.maximum(last[wet], day)
    arrays = {
        "max_extent": np.where(n_valid > 0, (n_wet > 0).astype(np.uint8), 255).astype(
            np.uint8
        ),
        "n_valid": n_valid,
        "n_wet": n_wet,
        "first_wet": first,
        "last_wet": last,
    }
    for name, array in arrays.items():
        gis.Raster.write_array(
            ev_dir / f"{task['unit_id']}_{name}.tif",
            array,
            gt,
            wkt,
            nodata={"max_extent": 255, "first_wet": -1, "last_wet": -1}.get(name),
        )
    days = sorted({s["datetime"][:10] for s in task["scenes"]})
    peak = np.datetime64(task["peak_date"])
    seconds = time.perf_counter() - start
    done = {
        "event_id": task["event_id"],
        "unit_id": task["unit_id"],
        "scenes": len(task["scenes"]),
        "scene_days": len(days),
        "first_scene_day": int((np.datetime64(days[0]) - peak).astype(int)),
        "last_scene_day": int((np.datetime64(days[-1]) - peak).astype(int)),
        "scene_days_after_peak": int(sum(np.datetime64(d) > peak for d in days)),
        "wet_km2": float((n_wet > 0).sum() * res * res / 1e6),
        "observed_km2": float((n_valid > 0).sum() * res * res / 1e6),
        "unit_km2": float(np.isfinite(hand).sum() * res * res / 1e6),
        "valid_observations_mean": (
            float(n_valid[np.isfinite(hand)].mean()) if np.isfinite(hand).any() else 0.0
        ),
        "seconds": seconds,
        "peak_rss_mb": gis.Run.peak_rss_mb(),
        "summary": (
            f"event {task['event_id']} unit {task['unit_id']}: {len(task['scenes'])} scenes on {len(days)} days, wet "
            f"{(n_wet > 0).sum() * res * res / 1e6:.2f} km2, observed {(n_valid > 0).sum() * res * res / 1e6:.1f} km2 "
            f"in {seconds:.1f} s"
        ),
    }
    marker = ev_dir / f"{task['unit_id']}.json.tmp"
    marker.write_text(json.dumps(done, indent=1))
    os.replace(marker, ev_dir / f"{task['unit_id']}.json")
    logging.getLogger("event_composites").info(done["summary"])
    return done


if __name__ == "__main__":
    run = gis.Run.env_overrides(
        CONFIG_PATH=CONFIG_PATH,
        COUNTRY=COUNTRY,
        FORCE=FORCE,
        UNITS=UNITS,
        EVENTS=EVENTS,
        RUN_ID=RUN_ID,
        SHOW_PLOTS=SHOW_PLOTS,
    )
    cfg = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], KEYS)
    general, evidence, sources = cfg["general"], cfg["evidence"], cfg["sources"]
    work_dir = Path(general["work_dir"]) / run["COUNTRY"]
    out_dir = Path(general["output_dir"]) / run["COUNTRY"]
    log = gis.Run.logger("event_composites", out_dir / "logs", general["log_level"])
    gis.Config.warn_placeholders(cfg)
    gis.Run.configure_gdal(
        general["gdal_cachemax_mb"], general["gdal_threads_per_worker"]
    )
    list_need = [
        out_dir / "events.gpkg",
        out_dir / "gfm_index.csv",
        out_dir / "units_core.gpkg",
    ]
    if not all(p.exists() for p in list_need):
        raise SystemExit(
            f"Missing inputs {[str(p) for p in list_need if not p.exists()]}; run flood_events.py and download_gfm.py."
        )
    start = time.perf_counter()
    gdf_ev = gpd.read_file(out_dir / "events.gpkg")
    wkt = gdf_ev.crs.to_wkt()
    if run["EVENTS"]:
        gdf_ev = gdf_ev[gdf_ev["event_id"].isin(run["EVENTS"])]
    df_idx = pd.read_csv(
        out_dir / "gfm_index.csv", dtype={"unit_id": str}, keep_default_na=False
    )
    if run["UNITS"]:
        df_idx = df_idx[df_idx["unit_id"].isin(run["UNITS"])]
    df_idx["day"] = (
        pd.to_datetime(df_idx["datetime"], utc=True, format="ISO8601")
        .dt.tz_localize(None)
        .dt.normalize()
    )
    gdf_core = gpd.read_file(out_dir / "units_core.gpkg").to_crs(wkt)
    gdf_core = gdf_core[gdf_core["unit_id"].isin(set(df_idx["unit_id"]))]
    params = {
        "flood_values": sources["gfm_flood_values"],
        "exclusion_values": sources["gfm_exclusion_values"],
        "reference_water_values": sources["gfm_reference_water_values"],
        "nodata_values": sources["gfm_nodata_values"],
        "advisory_values": (
            sources["gfm_advisory_values"] if evidence["use_advisory_flags"] else []
        ),
    }

    # --- Events with scenes: units the buffered lines of the group cross, scene windows inside the event window ---
    gdf_buf = gdf_ev.set_geometry(
        gdf_ev.geometry.buffer(cfg["country"]["corridor_buffer_m"])
    )
    df_eu = gpd.sjoin(
        gdf_buf, gdf_core[["unit_id", "geometry"]], predicate="intersects"
    )
    by_unit = {u: d.sort_values("datetime") for u, d in df_idx.groupby("unit_id")}
    list_task, list_skip, n_without = [], 0, 0
    for e in df_eu.itertuples():
        d = by_unit[e.unit_id]
        sel = d[
            (d["day"] >= pd.Timestamp(e.window_start))
            & (d["day"] <= pd.Timestamp(e.window_end))
        ]
        if sel.empty:
            n_without += 1
            continue
        ev_dir = work_dir / "events" / e.event_id
        if (ev_dir / f"{e.unit_id}.json").exists() and not run["FORCE"]:
            list_skip += 1
            continue
        list_task.append(
            {
                "unit_id": e.unit_id,
                "event_id": e.event_id,
                "peak_date": e.peak_date,
                "window_start": e.window_start,
                "unit_dir": str(work_dir / "units" / e.unit_id),
                "event_dir": str(ev_dir),
                "scenes": sel.drop(columns="day").to_dict("records"),
                "warp_memory_mb": general["warp_memory_limit_mb"],
                "params": params,
            }
        )
    log.info(
        "Events: %d (of %d in events.gpkg) touch units with GFM scene windows; %d event-unit pairs to compute, %d "
        "finished and skipped, %d without a scene in their window.",
        df_eu["event_id"].nunique(),
        len(gdf_ev),
        len(list_task),
        list_skip,
        n_without,
    )
    list_res, list_fail = gis.Run.run_units(
        composite_unit,
        list_task,
        general["n_workers"],
        general["gdal_cachemax_mb"],
        general["gdal_threads_per_worker"],
        [len(t["scenes"]) for t in list_task],
        log,
    )
    if list_fail:
        log.error(
            "%d event-unit pairs failed in units %s; rerun to retry them.",
            len(list_fail),
            sorted(set(list_fail)),
        )

    # --- Table of all composites ---
    list_done = [
        json.loads(p.read_text())
        for p in sorted((work_dir / "events").glob("*/*.json"))
    ]
    df_comp = pd.DataFrame(list_done).drop(columns=["summary"], errors="ignore")
    df_comp.to_csv(out_dir / "event_composites.csv", index=False)
    if len(df_comp):
        per_event = df_comp.groupby("event_id").agg(
            units=("unit_id", "size"),
            scenes=("scenes", "sum"),
            wet_km2=("wet_km2", "sum"),
            days=("scene_days", "max"),
            after=("scene_days_after_peak", "max"),
        )
        log.info(
            "event_composites.csv: %d composites of %d events; per event: scene days median %.0f, days after the peak "
            "median %.0f, wet area km2 p50 %.2f, max %.1f; events without a scene day after the peak: %d. %.2f s per "
            "composite, peak RSS %.2f GB.",
            len(df_comp),
            len(per_event),
            per_event["days"].median(),
            per_event["after"].median(),
            per_event["wet_km2"].median(),
            per_event["wet_km2"].max(),
            int((per_event["after"] == 0).sum()),
            df_comp["seconds"].mean(),
            df_comp["peak_rss_mb"].max() / 1024,
        )

    # --- Comparison with the official maximum flood extent (files ordered in the GFM portal) ---
    list_cmp = []
    for path in gis.Config.as_list(evidence["max_extent_files"]):
        ds = gdal.Open(str(path))
        if ds is None:
            log.warning("Maximum flood extent file %s cannot be opened.", path)
            continue
        for row in df_comp.itertuples():
            ours_ds = gdal.Open(
                str(
                    work_dir / "events" / row.event_id / f"{row.unit_id}_max_extent.tif"
                )
            )
            gt = ours_ds.GetGeoTransform()
            ours = ours_ds.ReadAsArray()
            b = (
                gt[0],
                gt[3] - ours.shape[0] * gt[1],
                gt[0] + ours.shape[1] * gt[1],
                gt[3],
            )
            try:
                official, _ = gis.Raster.read_window(
                    str(path),
                    b,
                    gt[1],
                    "near",
                    general["warp_memory_limit_mb"],
                    dst_wkt=ours_ds.GetProjection(),
                )
            except RuntimeError:
                continue
            covered = np.isfinite(official) & ~np.isin(
                official, sources["gfm_nodata_values"]
            )
            if not covered.any():
                continue
            off = covered & np.isin(official, sources["gfm_flood_values"])
            mine = (ours == 1) & covered
            cell = gt[1] ** 2 / 1e6
            list_cmp.append(
                {
                    "file": Path(path).name,
                    "event_id": row.event_id,
                    "unit_id": row.unit_id,
                    "both_km2": float((mine & off).sum() * cell),
                    "only_composite_km2": float((mine & ~off).sum() * cell),
                    "only_official_km2": float((off & ~mine).sum() * cell),
                    "official_on_cells_never_valid_km2": float(
                        (off & (ours == 255)).sum() * cell
                    ),
                }
            )
    if list_cmp:
        df_cmp = pd.DataFrame(list_cmp)
        df_cmp.to_csv(out_dir / "max_extent_comparison.csv", index=False)
        tot = df_cmp.groupby(["file", "event_id"])[
            [
                "both_km2",
                "only_composite_km2",
                "only_official_km2",
                "official_on_cells_never_valid_km2",
            ]
        ].sum()
        tot["iou"] = tot["both_km2"] / (
            tot["both_km2"] + tot["only_composite_km2"] + tot["only_official_km2"]
        ).clip(lower=1e-9)
        for (name, event_id), r in (
            tot.sort_values("iou", ascending=False).head(10).iterrows()
        ):
            log.info(
                "Official maximum flood extent %s against event %s: both %.2f km2, only composite %.2f, only official "
                "%.2f (%.2f of it on cells the composite never saw valid), IoU %.2f.",
                name,
                event_id,
                r["both_km2"],
                r["only_composite_km2"],
                r["only_official_km2"],
                r["official_on_cells_never_valid_km2"],
                r["iou"],
            )
    elif evidence["max_extent_files"]:
        log.warning("The maximum flood extent files overlap no composite.")

    # --- Maps: maximum extent, first wet day and number of valid observations, with the protected areas ---
    settings = cfg["implied"]["plots"]
    if settings["enabled"] and len(df_comp):
        plots = Plots(
            out_dir / "plots" / "events" / run["RUN_ID"],
            settings,
            run["SHOW_PLOTS"],
            log,
        )
        plt = plots.plt
        from matplotlib.colors import ListedColormap

        gdf_pa = (
            gpd.read_file(out_dir / "protected_areas.gpkg").to_crs(wkt)
            if (out_dir / "protected_areas.gpkg").exists()
            else None
        )
        df_comp["name"] = df_comp["event_id"] + "_" + df_comp["unit_id"]
        order = df_comp.sort_values("wet_km2", ascending=False)
        chosen = plots.select(
            "composite_maps",
            list(order["name"]),
            set(order["name"].head(settings["max_per_type"] // 2)),
        )
        for row in order[order["name"].isin(chosen)].itertuples():
            base = work_dir / "events" / row.event_id / row.unit_id
            arrays = {}
            for name in ("max_extent", "first_wet", "n_valid"):
                ds = gdal.Open(f"{base}_{name}.tif")
                arrays[name] = ds.ReadAsArray(
                    buf_xsize=min(ds.RasterXSize, 900),
                    buf_ysize=min(ds.RasterYSize, 900),
                )
                gt = ds.GetGeoTransform()
                extent = (
                    gt[0],
                    gt[0] + ds.RasterXSize * gt[1],
                    gt[3] - ds.RasterYSize * gt[1],
                    gt[3],
                )
            fig, axes = plt.subplots(1, 3, figsize=(15, 5.2), layout="constrained")
            ext = np.ma.masked_equal(arrays["max_extent"], 255)
            axes[0].imshow(
                ext,
                extent=extent,
                cmap=ListedColormap(["#ececec", CATEGORICAL[0]]),
                vmin=0,
                vmax=1,
                interpolation="nearest",
            )
            axes[0].set_title(
                "maximum extent (grey = observed dry, white = never valid)"
            )
            im = axes[1].imshow(
                np.ma.masked_less(arrays["first_wet"], 0),
                extent=extent,
                cmap="viridis",
                interpolation="nearest",
            )
            fig.colorbar(
                im, ax=axes[1], shrink=0.7, label="first wet day of the window"
            )
            axes[1].set_title("first wet day")
            im = axes[2].imshow(
                arrays["n_valid"], extent=extent, cmap="Greys", interpolation="nearest"
            )
            fig.colorbar(im, ax=axes[2], shrink=0.7, label="valid observations")
            axes[2].set_title("number of valid observations")
            for ax in axes:
                if gdf_pa is not None and len(gdf_pa):
                    gdf_pa.cx[
                        extent[0] : extent[1], extent[2] : extent[3]
                    ].boundary.plot(ax=ax, color=CATEGORICAL[1], linewidth=0.8)
                ax.set_xlim(extent[0], extent[1])
                ax.set_ylim(extent[2], extent[3])
                ax.grid(False)
                ax.tick_params(labelsize=6, colors=INK["muted"])
            fig.suptitle(
                f"{row.event_id}, unit {row.unit_id}: {row.scenes} scenes on {row.scene_days} days "
                f"({row.first_scene_day:+d} to {row.last_scene_day:+d} days from the peak), wet {row.wet_km2:.2f} km2; "
                "protected areas in orange",
                fontsize=9,
            )
            plots.save(
                fig,
                "composite_maps",
                row.name,
                {
                    "scenes": row.scenes,
                    "scene_days": row.scene_days,
                    "wet_km2": row.wet_km2,
                    "days_after_peak": row.scene_days_after_peak,
                },
            )
        log.info("Plots: %s", plots.finish())
    log.info("Done in %.1f s.", time.perf_counter() - start)
    if list_fail:
        raise SystemExit(1)
