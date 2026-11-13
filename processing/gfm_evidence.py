# Author: Jakub Zapletal
# Date: 2026-11-04
# Tests every protected area against GFM flood scenes: dry inside while the undefended floodplain at the same height
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
LAYER = "A"  # "A": protected areas of all water sources together, "B": per river group (no hints)

KEYS = [
    "general.work_dir",
    "general.output_dir",
    "general.n_workers",
    "general.gdal_cachemax_mb",
    "general.warp_memory_limit_mb",
    "general.gdal_threads_per_worker",
    "general.log_level",
    "country.protected_areas.max_level_m",
    "country.protected_areas.min_area_m2",
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
    "sources.gfm_advisory_values",
    "evidence.use_advisory_flags",
]
COLUMNS = [
    "protected_area_id",
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
    labels_ds = gdal.Open(task["labels_path"])
    gt, wkt = labels_ds.GetGeoTransform(), labels_ds.GetProjection()
    res = gt[1]
    labels = labels_ds.ReadAsArray()
    kept = dict(task["protected_areas"])  # unit_label -> protected_area_id
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
    # Analysis window: the protected areas and their comparison zones (whole unit window with hints).
    labels = labels[r0:r1, c0:c1]
    hand = gdal.Open(task["hand_path"]).ReadAsArray(
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

    # Comparison zone per protected area: outside all protected areas, within the comparison distance, not higher than
    # the median HAND inside the protected area and not above max_level_m.
    zones = gis.Evidence.comparison_zones(labels, hand, kept, reach, p["max_level_m"])

    if task["hints"]:
        gdf_rn = gpd.read_file(task["rivers_used"], bbox=bounds)
        stream = (
            gis.Raster.rasterize_geometries(
                gdf_rn.geometry,
                np.ones(len(gdf_rn)),
                labels.shape,
                gt,
                wkt,
                False,
                np.uint8,
            )
            > 0
        )
        # Levee bodies, their toes and the mixed GFM pixels along them stay dry whenever a levee holds: cells
        # close to a barrier line or a protected area are not tested.
        gdf_bl = gpd.read_file(task["barrier_lines"], bbox=bounds)
        near = (labels > 0) | (
            gis.Raster.rasterize_geometries(
                gdf_bl.geometry,
                np.ones(len(gdf_bl)),
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

    list_rows = []
    for scene in task["scenes"]:
        valid, wet, reference = gis.Evidence.scene_masks(
            scene, bounds, res, wkt, hand, p, task["warp_memory_mb"]
        )
        for t in gis.Evidence.test_regions(zones, valid, wet, p):
            list_rows.append(
                [
                    t["region"],
                    task["unit_id"],
                    scene["item_id"],
                    scene["datetime"],
                    t["inside_valid"],
                    t["inside_frac"],
                    t["comp_valid"],
                    t["comp_frac"],
                    t["outcome"],
                ]
            )
        if task["hints"]:
            # Scene water level in HAND: percentile of HAND among wet valid cells connected to the river. The
            # river itself is permanent water (not valid, never "wet"), so connectivity runs through the wet
            # cells, the reference water and the river lines.
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
    pd.DataFrame(list_rows, columns=COLUMNS).to_csv(task["csv_path"], index=False)

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
        gdf_hint = gis.Raster.polygonize_labels(clusters, gt, wkt)
        if len(gdf_hint):
            ids = gdf_hint["label"].to_numpy()
            gdf_hint["area_m2"] = np.bincount(clusters.ravel())[ids] * res * res
            gdf_hint["scenes_tested_mean"] = ndimage.mean(test_count, clusters, ids)
            gdf_hint["dry_share_mean"] = ndimage.mean(share, clusters, ids)
        gdf_hint["unit_id"] = task["unit_id"]
        gis.Vector.write_gpkg(gdf_hint, unit_dir / "missing_barrier_hints.gpkg")
    seconds = time.perf_counter() - start
    summary = f"unit {task['unit_id']}{task['group_text']}: {len(task['scenes'])} scenes x {len(kept)} protected areas tested in {seconds:.1f} s"
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
        LAYER=LAYER,
    )
    cfg = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], KEYS)
    files = gis.Run.layer_files(run["LAYER"])
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
        or not (out_dir / files["areas"]).exists()
    ):
        raise SystemExit(
            f"gfm_index.csv or {files['areas']} is missing; run protected_areas.py and download_gfm.py."
        )
    if run["HINTS"] and files["layer"] == "B":
        log.warning(
            "Missing-barrier hints come from file A; HINTS is ignored for LAYER = B."
        )
        run["HINTS"] = False
    df_idx = pd.read_csv(out_dir / "gfm_index.csv", dtype={"unit_id": str}).fillna("")
    if len(df_idx):
        dates = pd.to_datetime(df_idx["datetime"], utc=True, format="ISO8601").dt.date
        for bound, compare in [
            (run["DATE_FROM"], dates.__ge__),
            (run["DATE_TO"], dates.__le__),
        ]:
            if bound:
                df_idx = df_idx[compare(datetime.date.fromisoformat(bound)).to_numpy()]
    if run["UNITS"]:
        df_idx = df_idx[df_idx["unit_id"].isin(run["UNITS"])]
    gdf_pa = gpd.read_file(out_dir / files["areas"])
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
        "max_level_m": cfg["country"]["protected_areas"]["max_level_m"],
        "min_area_m2": cfg["country"]["protected_areas"]["min_area_m2"],
        "flood_values": sources["gfm_flood_values"],
        "exclusion_values": sources["gfm_exclusion_values"],
        "reference_water_values": sources["gfm_reference_water_values"],
        "nodata_values": sources["gfm_nodata_values"],
        # switch: cells with an advisory flag (frost, dry soil, wind) are not valid observations
        "advisory_values": (
            sources["gfm_advisory_values"] if evidence["use_advisory_flags"] else []
        ),
    }
    if evidence["use_advisory_flags"]:
        n_adv = (
            int((df_idx["advisory_path"] != "").sum())
            if "advisory_path" in df_idx.columns
            else 0
        )
        log.info(
            "use_advisory_flags: %d of %d scene windows have the advisory layer; the others are tested without it "
            "(run download_gfm.py with gfm_assets.advisory set to get it).",
            n_adv,
            len(df_idx),
        )
    list_task = []
    for unit_id, df_scn in df_idx.groupby("unit_id"):
        unit_dir = work_dir / "units" / unit_id
        gdf_unit = gdf_pa[gdf_pa["unit_id"] == unit_id]
        if gdf_unit.empty and not run["HINTS"]:
            continue
        # file A: one task per unit; file B: one per river group of the unit (its own labels and HAND)
        list_group = (
            sorted(gdf_unit["group_id"].unique()) if files["layer"] == "B" else [None]
        )
        for group_id in list_group:
            gdf_mem = (
                gdf_unit
                if group_id is None
                else gdf_unit[gdf_unit["group_id"] == group_id]
            )
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
                    "csv_path": gis.Run.unit_layer_path(
                        unit_dir, "gfm_evidence.csv", group_id
                    ),
                    "unit_dir": str(unit_dir),
                    "protected_areas": list(
                        zip(
                            gdf_mem["unit_label"].astype(int),
                            gdf_mem["protected_area_id"],
                        )
                    ),
                    "scenes": df_scn.to_dict("records"),
                    "hints": run["HINTS"],
                    "rivers_used": str(work_dir / "rivers_used.gpkg"),
                    "barrier_lines": str(work_dir / "barrier_lines.gpkg"),
                    "warp_memory_mb": general["warp_memory_limit_mb"],
                    "params": params,
                }
            )
    log.info(
        "Evidence: %d scene windows in %d units with protected areas%s.",
        len(df_idx),
        len(list_task),
        (
            f" (dates {run['DATE_FROM'] or '...'} to {run['DATE_TO'] or '...'})"
            if run["DATE_FROM"] or run["DATE_TO"]
            else ""
        ),
    )
    list_res, list_fail = gis.Run.run_units(
        evidence_unit,
        list_task,
        general["n_workers"],
        general["gdal_cachemax_mb"],
        general["gdal_threads_per_worker"],
        [len(t["scenes"]) for t in list_task],
        log,
    )
    if list_fail:
        raise SystemExit(f"Units failed: {list_fail}; rerun to retry them.")

    list_ev = [pd.read_csv(t["csv_path"]) for t in list_task]
    df_ev = (
        pd.concat(list_ev, ignore_index=True)
        if list_ev
        else pd.DataFrame(columns=COLUMNS)
    )
    df_ev.to_csv(out_dir / f"gfm_evidence{files['suffix']}.csv", index=False)
    # One feature per protected area and scene, for visual checks in QGIS.
    gdf_ev = gdf_pa[["protected_area_id", "geometry"]].merge(
        df_ev, on="protected_area_id", how="inner"
    )
    gis.Vector.write_gpkg(
        gpd.GeoDataFrame(gdf_ev, geometry="geometry", crs=gdf_pa.crs),
        out_dir / f"gfm_evidence{files['suffix']}.gpkg",
    )
    counts = df_ev["outcome"].value_counts()
    log.info(
        "gfm_evidence%s.csv: %d tests (%s) for %d protected areas.",
        files["suffix"],
        len(df_ev),
        ", ".join(f"{k} {v}" for k, v in counts.items()),
        df_ev["protected_area_id"].nunique() if len(df_ev) else 0,
    )

    if run["HINTS"]:
        gdf_core = gpd.read_file(out_dir / "units_core.gpkg").set_index("unit_id")
        list_hint = [
            gpd.read_file(Path(t["unit_dir"]) / "missing_barrier_hints.gpkg")
            for t in list_task
            if (Path(t["unit_dir"]) / "missing_barrier_hints.gpkg").exists()
        ]
        list_hint = [f for f in list_hint if len(f)]
        if list_hint:
            gdf_hint = gpd.GeoDataFrame(
                pd.concat(list_hint, ignore_index=True),
                geometry="geometry",
                crs=list_hint[0].crs,
            )
            inside = shapely.covers(
                gdf_core.geometry.reindex(gdf_hint["unit_id"]).values,
                shapely.point_on_surface(gdf_hint.geometry.values),
            )
            gdf_hint = gdf_hint[inside].reset_index(drop=True)
        else:
            gdf_hint = gpd.GeoDataFrame(
                {"unit_id": pd.Series(dtype=str)},
                geometry=gpd.GeoSeries([], crs=gdf_pa.crs),
            )
        gis.Vector.write_gpkg(gdf_hint, out_dir / "missing_barrier_hints.gpkg")
        log.info(
            "missing_barrier_hints.gpkg: %d areas that stay dry below the scene water level.",
            len(gdf_hint),
        )
