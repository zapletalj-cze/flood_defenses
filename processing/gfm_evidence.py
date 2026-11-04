# Author: Jakub Zapletal
# Date: 2026-11-04
# Tests every compartment against GFM flood scenes: dry inside while the undefended floodplain at the same height
# around it is wet (held), or wet inside (flooded); optional hints of missing barriers from cells that stay dry.

import datetime
import logging
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from osgeo import gdal
from scipy import ndimage

import gis

# --- Run variables (no command-line parameters; DEFENSES_<NAME> environment variables override them) ---
CONFIG_PATH = "config/defenses.toml"
COUNTRY = "XX"
UNITS = []  # unit ids; empty = all units with GFM scenes
HINTS = (
    False  # also derive missing-barrier hints (reads the whole unit window per scene)
)
DATE_FROM = ""  # optional scene date filter, YYYY-MM-DD (e.g. one flood event)
DATE_TO = ""

KEYS = [
    "general.work_dir",
    "general.output_dir",
    "general.n_workers",
    "general.gdal_cachemax_mb",
    "general.warp_memory_limit_mb",
    "general.gdal_threads_per_worker",
    "general.log_level",
    "country.compartments.max_level_m",
    "country.compartments.min_area_m2",
    "evidence.comparison_distance_m",
    "evidence.min_valid_pixels",
    "evidence.f_wet_comparison",
    "evidence.f_dry_inside",
    "evidence.f_wet_inside",
    "evidence.river_level_percentile",
    "evidence.hint_margin_m",
    "evidence.hint_dry_share",
    "evidence.hint_min_scenes",
    "evidence.hint_exclusion_m",
    "sources.gfm_flood_values",
    "sources.gfm_exclusion_values",
    "sources.gfm_reference_water_values",
    "sources.gfm_nodata_values",
]
COLUMNS = [
    "compartment_id",
    "unit_id",
    "item_id",
    "datetime",
    "inside_valid_px",
    "inside_wet_frac",
    "comp_valid_px",
    "comp_wet_frac",
    "outcome",
]


def evidence_unit(task):
    start = time.perf_counter()
    p = task["params"]
    unit_dir = Path(task["unit_dir"])
    labels_ds = gdal.Open(str(unit_dir / "labels.tif"))
    gt, wkt = labels_ds.GetGeoTransform(), labels_ds.GetProjection()
    res = gt[1]
    labels = labels_ds.ReadAsArray()
    kept = dict(task["compartments"])  # unit_label -> compartment_id
    reach = int(np.ceil(p["comparison_distance_m"] / res))
    if task["hints"] or not kept:
        r0, r1, c0, c1 = 0, labels.shape[0], 0, labels.shape[1]
    else:
        mask = np.isin(labels, np.array(list(kept), dtype=labels.dtype))
        rows_any, cols_any = np.flatnonzero(mask.any(axis=1)), np.flatnonzero(
            mask.any(axis=0)
        )
        r0, r1 = max(rows_any[0] - reach, 0), min(
            rows_any[-1] + 1 + reach, labels.shape[0]
        )
        c0, c1 = max(cols_any[0] - reach, 0), min(
            cols_any[-1] + 1 + reach, labels.shape[1]
        )
    # Analysis window: the compartments and their comparison zones (whole unit window with hints).
    labels = labels[r0:r1, c0:c1]
    hand = gdal.Open(str(unit_dir / "hand_undef.tif")).ReadAsArray(
        int(c0), int(r0), int(c1 - c0), int(r1 - r0)
    )
    hand = hand.astype(np.float32)
    gt = (gt[0] + c0 * res, res, 0.0, gt[3] - r0 * res, 0.0, -res)
    bounds = (
        gt[0],
        gt[3] - labels.shape[0] * res,
        gt[0] + labels.shape[1] * res,
        gt[3],
    )

    # Comparison zone per compartment: outside all compartments, within the comparison distance, not higher than
    # the median HAND inside the compartment and not above max_level_m.
    zones = {}
    for label, compartment_id in kept.items():
        inside = labels == label
        if not inside.any():
            continue
        rows_any, cols_any = np.flatnonzero(inside.any(axis=1)), np.flatnonzero(
            inside.any(axis=0)
        )
        q0, q1 = max(rows_any[0] - reach, 0), min(
            rows_any[-1] + 1 + reach, labels.shape[0]
        )
        p0, p1 = max(cols_any[0] - reach, 0), min(
            cols_any[-1] + 1 + reach, labels.shape[1]
        )
        part = (slice(q0, q1), slice(p0, p1))
        near = gis.Hydro.distance_cells(inside[part], reach + 1) <= reach
        level = np.nanmedian(hand[part][inside[part]])
        zone = (
            near
            & (labels[part] == 0)
            & (hand[part] <= level)
            & (hand[part] <= p["max_level_m"])
        )
        zones[label] = (compartment_id, part, inside[part], zone)

    if task["hints"]:
        rivers = gpd.read_file(task["rivers_used"], bbox=bounds)
        stream = (
            gis.Raster.rasterize_geometries(
                rivers.geometry,
                np.ones(len(rivers)),
                labels.shape,
                gt,
                wkt,
                False,
                np.uint8,
            )
            > 0
        )
        # Levee bodies, their toes and the mixed GFM pixels along them stay dry whenever a levee holds: cells
        # close to a barrier line or a compartment are not tested.
        barriers = gpd.read_file(task["barrier_lines"], bbox=bounds)
        near = (labels > 0) | (
            gis.Raster.rasterize_geometries(
                barriers.geometry,
                np.ones(len(barriers)),
                labels.shape,
                gt,
                wkt,
                True,
                np.uint8,
            )
            > 0
        )
        exclusion = p["hint_exclusion_m"] / res
        testable = gis.Hydro.distance_cells(near, exclusion + 1) > exclusion
        del near
        dry_count = np.zeros(labels.shape, np.uint16)
        test_count = np.zeros(labels.shape, np.uint16)

    nodata = p["nodata_values"][0] if p["nodata_values"] else None
    rows = []
    for scene in task["scenes"]:
        layers = {}
        for role in ["flood", "exclusion", "reference_water"]:
            layers[role], _ = gis.Raster.read_window(
                scene[f"{role}_path"],
                bounds,
                res,
                "near",
                task["warp_memory_mb"],
                dst_wkt=wkt,
                src_nodata=nodata,
            )
        flood = layers["flood"]
        valid = (
            np.isfinite(flood)
            & ~np.isin(flood, p["nodata_values"])
            & ~np.isin(layers["exclusion"], p["exclusion_values"])
            & ~np.isin(layers["reference_water"], p["reference_water_values"])
            & np.isfinite(hand)
        )
        wet = valid & np.isin(flood, p["flood_values"])
        for label, (compartment_id, part, inside, zone) in zones.items():
            inside_valid = int((inside & valid[part]).sum())
            inside_wet = int((inside & wet[part]).sum())
            comp_valid = int((zone & valid[part]).sum())
            comp_wet = int((zone & wet[part]).sum())
            inside_frac = inside_wet / inside_valid if inside_valid else np.nan
            comp_frac = comp_wet / comp_valid if comp_valid else np.nan
            # A dry compartment counts only if the undefended floodplain around it is wet in the same scene.
            if (
                inside_valid < p["min_valid_pixels"]
                or comp_valid < p["min_valid_pixels"]
                or not comp_frac >= p["f_wet_comparison"]
            ):
                outcome = "no_test"
            elif inside_frac <= p["f_dry_inside"]:
                outcome = "held"
            elif inside_frac >= p["f_wet_inside"]:
                outcome = "flooded"
            else:
                outcome = "inconclusive"
            rows.append(
                [
                    compartment_id,
                    task["unit_id"],
                    scene["item_id"],
                    scene["datetime"],
                    inside_valid,
                    inside_frac,
                    comp_valid,
                    comp_frac,
                    outcome,
                ]
            )
        if task["hints"]:
            # Scene water level in HAND: percentile of HAND among wet valid cells connected to the river. The
            # river itself is permanent water (not valid, never "wet"), so connectivity runs through the wet
            # cells, the reference water and the river lines.
            reference = np.isin(layers["reference_water"], p["reference_water_values"])
            parts, _ = ndimage.label(wet | reference | stream)
            touching = np.unique(parts[stream])
            connected = np.isin(parts, touching[touching > 0]) & wet
            del parts, reference
            if connected.any():
                level = np.percentile(hand[connected], p["river_level_percentile"])
                tested = (
                    valid
                    & testable
                    & (labels == 0)
                    & (hand < level - p["hint_margin_m"])
                )
                test_count += tested
                dry_count += tested & ~wet
    pd.DataFrame(rows, columns=COLUMNS).to_csv(
        unit_dir / "gfm_evidence.csv", index=False
    )

    if task["hints"]:
        share = np.where(test_count > 0, dry_count / np.maximum(test_count, 1), 0.0)
        hint = (test_count >= p["hint_min_scenes"]) & (share >= p["hint_dry_share"])
        clusters, n = ndimage.label(hint)
        sizes = np.bincount(clusters.ravel())
        keep = np.flatnonzero(sizes * res * res >= p["min_area_m2"])
        keep = keep[keep > 0]
        relabel = np.zeros(n + 1, np.int32)
        relabel[keep] = np.arange(1, len(keep) + 1)
        clusters = relabel[clusters]
        hints = gis.Raster.polygonize_labels(clusters, gt, wkt)
        if len(hints):
            ids = hints["label"].to_numpy()
            hints["area_m2"] = np.bincount(clusters.ravel())[ids] * res * res
            hints["scenes_tested_mean"] = ndimage.mean(test_count, clusters, ids)
            hints["dry_share_mean"] = ndimage.mean(share, clusters, ids)
        hints["unit_id"] = task["unit_id"]
        gis.Vector.write_gpkg(hints, unit_dir / "missing_barrier_hints.gpkg")
    seconds = time.perf_counter() - start
    summary = f"unit {task['unit_id']}: {len(task['scenes'])} scenes x {len(kept)} compartments tested in {seconds:.1f} s"
    logging.getLogger("gfm_evidence").info(summary)
    return {"unit_id": task["unit_id"], "seconds": seconds, "summary": summary}


if __name__ == "__main__":
    run = gis.Run.env_overrides(
        CONFIG_PATH=CONFIG_PATH,
        COUNTRY=COUNTRY,
        UNITS=UNITS,
        HINTS=HINTS,
        DATE_FROM=DATE_FROM,
        DATE_TO=DATE_TO,
    )
    cfg = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], KEYS)
    general, evidence, sources = cfg["general"], cfg["evidence"], cfg["sources"]
    work_dir = Path(general["work_dir"]) / run["COUNTRY"]
    out_dir = Path(general["output_dir"]) / run["COUNTRY"]
    log = gis.Run.logger("gfm_evidence", out_dir / "logs", general["log_level"])
    gis.Config.warn_placeholders(cfg)
    gis.Run.configure_gdal(
        general["gdal_cachemax_mb"], general["gdal_threads_per_worker"]
    )
    if (
        not (out_dir / "gfm_index.csv").exists()
        or not (out_dir / "compartments.gpkg").exists()
    ):
        raise SystemExit(
            "gfm_index.csv or compartments.gpkg is missing; run compartments.py and download_gfm.py."
        )
    index = pd.read_csv(out_dir / "gfm_index.csv", dtype={"unit_id": str}).fillna("")
    if len(index):
        dates = pd.to_datetime(index["datetime"], utc=True, format="ISO8601").dt.date
        for bound, compare in [
            (run["DATE_FROM"], dates.__ge__),
            (run["DATE_TO"], dates.__le__),
        ]:
            if bound:
                index = index[compare(datetime.date.fromisoformat(bound)).to_numpy()]
    if run["UNITS"]:
        index = index[index["unit_id"].isin(run["UNITS"])]
    compartments = gpd.read_file(out_dir / "compartments.gpkg")
    params = {
        **{
            k: evidence[k]
            for k in [
                "comparison_distance_m",
                "min_valid_pixels",
                "f_wet_comparison",
                "f_dry_inside",
                "f_wet_inside",
                "river_level_percentile",
                "hint_margin_m",
                "hint_dry_share",
                "hint_min_scenes",
                "hint_exclusion_m",
            ]
        },
        "max_level_m": cfg["country"]["compartments"]["max_level_m"],
        "min_area_m2": cfg["country"]["compartments"]["min_area_m2"],
        "flood_values": sources["gfm_flood_values"],
        "exclusion_values": sources["gfm_exclusion_values"],
        "reference_water_values": sources["gfm_reference_water_values"],
        "nodata_values": sources["gfm_nodata_values"],
    }
    tasks = []
    for unit_id, scenes in index.groupby("unit_id"):
        members = compartments[compartments["unit_id"] == unit_id]
        if members.empty and not run["HINTS"]:
            continue
        tasks.append(
            {
                "unit_id": unit_id,
                "unit_dir": str(work_dir / "units" / unit_id),
                "compartments": list(
                    zip(members["unit_label"].astype(int), members["compartment_id"])
                ),
                "scenes": scenes.to_dict("records"),
                "hints": run["HINTS"],
                "rivers_used": str(work_dir / "rivers_used.gpkg"),
                "barrier_lines": str(work_dir / "barrier_lines.gpkg"),
                "warp_memory_mb": general["warp_memory_limit_mb"],
                "params": params,
            }
        )
    log.info(
        "Evidence: %d scene windows in %d units with compartments%s.",
        len(index),
        len(tasks),
        (
            f" (dates {run['DATE_FROM'] or '...'} to {run['DATE_TO'] or '...'})"
            if run["DATE_FROM"] or run["DATE_TO"]
            else ""
        ),
    )
    results, failed = gis.Run.run_units(
        evidence_unit,
        tasks,
        general["n_workers"],
        general["gdal_cachemax_mb"],
        general["gdal_threads_per_worker"],
        [len(t["scenes"]) for t in tasks],
        log,
    )
    if failed:
        raise SystemExit(f"Units failed: {failed}; rerun to retry them.")

    tables = [pd.read_csv(Path(t["unit_dir"]) / "gfm_evidence.csv") for t in tasks]
    table = (
        pd.concat(tables, ignore_index=True)
        if tables
        else pd.DataFrame(columns=COLUMNS)
    )
    table.to_csv(out_dir / "gfm_evidence.csv", index=False)
    # One feature per compartment and scene, for visual checks in QGIS.
    layer = compartments[["compartment_id", "geometry"]].merge(
        table, on="compartment_id", how="inner"
    )
    gis.Vector.write_gpkg(
        gpd.GeoDataFrame(layer, geometry="geometry", crs=compartments.crs),
        out_dir / "gfm_evidence.gpkg",
    )
    counts = table["outcome"].value_counts()
    log.info(
        "gfm_evidence.csv: %d tests (%s) for %d compartments.",
        len(table),
        ", ".join(f"{k} {v}" for k, v in counts.items()),
        table["compartment_id"].nunique() if len(table) else 0,
    )

    if run["HINTS"]:
        cores = gpd.read_file(out_dir / "units_core.gpkg").set_index("unit_id")
        frames = [
            gpd.read_file(Path(t["unit_dir"]) / "missing_barrier_hints.gpkg")
            for t in tasks
        ]
        frames = [f for f in frames if len(f)]
        if frames:
            hints = gpd.GeoDataFrame(
                pd.concat(frames, ignore_index=True),
                geometry="geometry",
                crs=frames[0].crs,
            )
            inside = shapely.covers(
                cores.geometry.reindex(hints["unit_id"]).values,
                shapely.point_on_surface(hints.geometry.values),
            )
            hints = hints[inside].reset_index(drop=True)
        else:
            hints = gpd.GeoDataFrame(
                {"unit_id": pd.Series(dtype=str)},
                geometry=gpd.GeoSeries([], crs=compartments.crs),
            )
        gis.Vector.write_gpkg(hints, out_dir / "missing_barrier_hints.gpkg")
        log.info(
            "missing_barrier_hints.gpkg: %d areas that stay dry below the scene water level.",
            len(hints),
        )
