# Author: Jakub Zapletal
# Date: 2026-11-05
# Stage E, step E2: bank subsegments per unit. Envelope of low HAND connected to the snapped rivers, nearest river
# point of every envelope cell (EDT with indices), bank side and chainage bins; per-unit labels, merged polygons.

import hashlib
import json
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
FORCE = (
    False  # recompute units whose checkpoint matches the current inputs and parameters
)

KEYS = [
    "general.work_dir",
    "general.output_dir",
    "general.log_level",
    "general.n_workers",
    "general.gdal_cachemax_mb",
    "general.gdal_threads_per_worker",
    "country.working_resolution_m",
    "country.compartments.stream_burn_m",
    "country.compartments.stream_burn_buffer_m",
    "implied.units",
    "implied.recompute_hand",
    "implied.envelope_hand_m",
    "implied.subsegment_length_m",
]
SIDES = np.array(["left", "right"])  # looking downstream


def segment_unit(task):
    # Worker: subsegment labels of one unit (core cells only; the halo serves the nearest-point search).
    start = time.perf_counter()
    unit_dir = Path(task["unit_dir"])
    secs = {}
    t = time.perf_counter()
    ds = gdal.Open(str(unit_dir / "hand_undef.tif"))
    gt, wkt = ds.GetGeoTransform(), ds.GetProjection()
    shape = (ds.RasterYSize, ds.RasterXSize)
    res = gt[1]
    lines = shapely.from_wkb(task["rivers"])
    river = (
        gis.Raster.rasterize_geometries(
            lines, np.ones(len(lines)), shape, gt, wkt, True, np.uint8
        )
        > 0
        if len(lines)
        else np.zeros(shape, bool)
    )
    if task["recompute_hand"]:
        # HAND with the snapped rivers as streams, on the undefended DTM kept by stage A.
        ds_undef = gdal.Open(str(unit_dir / "dtm_undefended.tif"))
        dtm_undef = ds_undef.GetRasterBand(1).ReadAsArray().astype(np.float32)
        ds_undef = None
        hand, _, n_unreached = gis.Hydro.hand_above_streams(
            dtm_undef, river, task["burn_m"], task["burn_buffer_m"] / res
        )
        del dtm_undef
        gis.Raster.write_array(unit_dir / "hand_undef_snapped.tif", hand, gt, wkt)
    else:
        hand = ds.GetRasterBand(1).ReadAsArray().astype(np.float32)
        n_unreached = 0
    secs["hand"] = time.perf_counter() - t

    # --- Envelope: HAND at most envelope_hand_m, 4-connected to river cells ---
    t = time.perf_counter()
    low = np.isfinite(hand) & (hand <= task["envelope_hand_m"])
    lab, n_lab = ndimage.label(low | river)
    keep = np.zeros(n_lab + 1, bool)
    keep[np.unique(lab[river])] = True
    keep[0] = False
    envelope = keep[lab] & low
    del lab, low, hand
    secs["envelope"] = time.perf_counter() - t

    # --- Nearest river point of every cell (EDT with indices over the whole window) ---
    t = time.perf_counter()
    x, y = task["points"]["x"], task["points"]["y"]
    pr = np.floor((gt[3] - y) / res).astype(np.int64)
    pc = np.floor((x - gt[0]) / res).astype(np.int64)
    inside = (pr >= 0) & (pr < shape[0]) & (pc >= 0) & (pc < shape[1])
    core = np.zeros(shape, bool)
    r0 = int(round((gt[3] - task["core"][3]) / res))
    r1 = int(round((gt[3] - task["core"][1]) / res))
    c0 = int(round((task["core"][0] - gt[0]) / res))
    c1 = int(round((task["core"][2] - gt[0]) / res))
    core[max(r0, 0) : r1, max(c0, 0) : c1] = True
    cells = envelope & core
    labels = np.zeros(shape, np.int32)
    df_lut = pd.DataFrame(
        {
            "label": pd.Series(dtype=np.int32),
            "river_id": pd.Series(dtype=np.int64),
            "side": pd.Series(dtype=str),
            "k": pd.Series(dtype=np.int64),
            "n_cells": pd.Series(dtype=np.int64),
        }
    )
    if inside.any() and cells.any():
        idx = np.full(shape, -1, np.int32)
        idx[pr[inside], pc[inside]] = np.flatnonzero(inside)
        ir, ic = ndimage.distance_transform_edt(
            idx < 0, return_distances=False, return_indices=True
        )
        rr, cc = np.nonzero(cells)
        near = idx[ir[rr, cc], ic[rr, cc]]
        del ir, ic, idx
        p = task["points"]
        vx = gt[0] + (cc + 0.5) * res - x[near]
        vy = gt[3] - (rr + 0.5) * res - y[near]
        # Cross product of the downstream tangent and the vector to the cell: positive = left bank.
        side = (p["tx"][near] * vy - p["ty"][near] * vx < 0).astype(np.int64)
        k = np.floor(p["chainage"][near] / task["length_m"]).astype(np.int64)
        n_k = task["n_k"]
        key = (p["river_id"][near].astype(np.int64) * 2 + side) * n_k + k
        uniq, inv = np.unique(key, return_inverse=True)
        labels[rr, cc] = inv + 1
        df_lut = pd.DataFrame(
            {
                "label": np.arange(1, len(uniq) + 1, dtype=np.int32),
                "river_id": uniq // (2 * n_k),
                "side": SIDES[(uniq // n_k) % 2],
                "k": uniq % n_k,
                "n_cells": np.bincount(inv, minlength=len(uniq)),
            }
        )
    df_lut["chainage_from_m"] = df_lut["k"] * task["length_m"]
    df_lut["chainage_to_m"] = (df_lut["k"] + 1) * task["length_m"]
    df_lut["envelope_area_m2"] = df_lut["n_cells"] * res * res
    secs["allocation"] = time.perf_counter() - t

    # --- Outputs: labels, lookup table, polygons of the core cells ---
    t = time.perf_counter()
    gis.Raster.write_array(unit_dir / "subsegments.tif", labels, gt, wkt, nodata=0)
    df_lut.to_csv(unit_dir / "subsegments.csv", index=False)
    gdf_poly = gis.Raster.polygonize_labels(labels, gt, wkt)
    gdf_poly = gdf_poly.merge(df_lut, on="label", how="inner")
    gis.Vector.write_gpkg(gdf_poly, unit_dir / "subsegments_unit.gpkg")
    secs["write"] = time.perf_counter() - t
    (unit_dir / "subsegments_done.json").write_text(
        json.dumps({"hash": task["hash"], "secs": secs})
    )
    return {
        "unit_id": task["unit_id"],
        "summary": (
            f"unit {task['unit_id']}: {len(df_lut)} subsegments, envelope "
            f"{df_lut['envelope_area_m2'].sum() / 1e6:.2f} km2 in the core, "
            f"{time.perf_counter() - start:.1f} s"
        ),
        "n_unreached": n_unreached,
        "secs": secs,
        "peak_rss_mb": gis.Run.peak_rss_mb(),
    }


if __name__ == "__main__":
    run = gis.Run.env_overrides(CONFIG_PATH=CONFIG_PATH, COUNTRY=COUNTRY, FORCE=FORCE)
    cfg = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], KEYS)
    general, implied = cfg["general"], cfg["implied"]
    burn = cfg["country"]["compartments"]
    work_dir = Path(general["work_dir"]) / run["COUNTRY"]
    out_dir = Path(general["output_dir"]) / run["COUNTRY"]
    log = gis.Run.logger("implied_segments", out_dir / "logs", general["log_level"])
    gis.Config.warn_placeholders(cfg)
    start = time.perf_counter()
    list_need = [
        out_dir / "units.csv",
        out_dir / "rivers_snapped.gpkg",
        out_dir / "river_points.gpkg",
    ]
    list_miss = [str(p) for p in list_need if not p.exists()]
    if list_miss:
        raise SystemExit(
            f"Missing inputs {list_miss}; run build_units.py, compartments.py and snap_rivers.py first."
        )
    df_units = pd.read_csv(out_dir / "units.csv", dtype={"unit_id": str})
    if implied["units"]:
        list_unk = sorted(set(implied["units"]) - set(df_units["unit_id"]))
        if list_unk:
            raise SystemExit(f"implied.units not found in units.csv: {list_unk}")
        df_units = df_units[df_units["unit_id"].isin(implied["units"])]
    has_hand = df_units["unit_id"].map(
        lambda u: (work_dir / "units" / u / "hand_undef.tif").exists()
    )
    if not has_hand.all():
        log.warning(
            "No hand_undef.tif (compartments.py) for units %s; they are skipped.",
            sorted(df_units.loc[~has_hand, "unit_id"]),
        )
    df_units = df_units[has_hand]
    if implied["recompute_hand"]:
        list_miss = [
            u
            for u in df_units["unit_id"]
            if not (work_dir / "units" / u / "dtm_undefended.tif").exists()
        ]
        if list_miss:
            raise SystemExit(
                "recompute_hand needs dtm_undefended.tif of stage A; run compartments.py with "
                f"keep_intermediate = true (missing for units {list_miss})."
            )
    if df_units.empty:
        raise SystemExit("No unit to process.")

    gdf_rn = gpd.read_file(out_dir / "rivers_snapped.gpkg")
    gdf_pts = gpd.read_file(out_dir / "river_points.gpkg")
    xy = shapely.get_coordinates(gdf_pts.geometry.values)
    length_m = implied["subsegment_length_m"]
    n_k = (
        int(np.ceil(gdf_pts["chainage_m"].max() / length_m)) + 1 if len(gdf_pts) else 1
    )
    pts_tree = shapely.STRtree(gdf_pts.geometry.values)
    rn_tree = shapely.STRtree(gdf_rn.geometry.values)
    params = {
        "envelope_hand_m": implied["envelope_hand_m"],
        "subsegment_length_m": length_m,
        "recompute_hand": implied["recompute_hand"],
        "n_k": n_k,
        "burn": [burn["stream_burn_m"], burn["stream_burn_buffer_m"]],
        "points_mtime": (out_dir / "river_points.gpkg").stat().st_mtime,
    }
    list_task, list_w, list_skip = [], [], []
    for u in df_units.itertuples():
        unit_dir = work_dir / "units" / u.unit_id
        digest = hashlib.sha256(
            json.dumps(
                {**params, "hand_mtime": (unit_dir / "hand_undef.tif").stat().st_mtime},
                sort_keys=True,
            ).encode()
        ).hexdigest()
        done = unit_dir / "subsegments_done.json"
        if (
            not run["FORCE"]
            and done.exists()
            and json.loads(done.read_text()).get("hash") == digest
        ):
            list_skip.append(u.unit_id)
            continue
        window = shapely.box(u.win_xmin, u.win_ymin, u.win_xmax, u.win_ymax)
        # rivers clipped two cells beyond the window, so lines leaving it still touch its edge cells
        reach = window.buffer(
            2 * cfg["country"]["working_resolution_m"], join_style="mitre"
        )
        sel = pts_tree.query(window, predicate="intersects")
        list_task.append(
            {
                "unit_id": u.unit_id,
                "unit_dir": str(unit_dir),
                "core": (u.core_xmin, u.core_ymin, u.core_xmax, u.core_ymax),
                "rivers": list(
                    shapely.to_wkb(
                        shapely.intersection(
                            gdf_rn.geometry.values[
                                rn_tree.query(reach, predicate="intersects")
                            ],
                            reach,
                        )
                    )
                ),
                "points": {
                    "x": xy[sel, 0],
                    "y": xy[sel, 1],
                    "river_id": gdf_pts["river_id"].to_numpy()[sel],
                    "chainage": gdf_pts["chainage_m"].to_numpy()[sel],
                    "tx": gdf_pts["tangent_dx"].to_numpy()[sel],
                    "ty": gdf_pts["tangent_dy"].to_numpy()[sel],
                },
                "recompute_hand": implied["recompute_hand"],
                "burn_m": burn["stream_burn_m"],
                "burn_buffer_m": burn["stream_burn_buffer_m"],
                "envelope_hand_m": implied["envelope_hand_m"],
                "length_m": length_m,
                "n_k": n_k,
                "hash": digest,
            }
        )
        list_w.append(float(u.est_window_cells))
    if list_skip:
        log.info(
            "%d units unchanged since their checkpoint (FORCE = True recomputes them).",
            len(list_skip),
        )
    t = time.perf_counter()
    list_res, list_fail = gis.Run.run_units(
        segment_unit,
        list_task,
        general["n_workers"],
        general["gdal_cachemax_mb"],
        general["gdal_threads_per_worker"],
        list_w,
        log,
    )
    wall = time.perf_counter() - t
    for res_u in list_res:
        if res_u["n_unreached"]:
            log.warning(
                "Unit %s: %d cells reached no snapped river cell; the nearest one was used.",
                res_u["unit_id"],
                res_u["n_unreached"],
            )

    # --- Merged subsegments: per-unit core polygons dissolved by (river_id, side, k) ---
    t = time.perf_counter()
    list_gdf = [
        gpd.read_file(work_dir / "units" / u / "subsegments_unit.gpkg")
        for u in df_units["unit_id"]
        if (work_dir / "units" / u / "subsegments_unit.gpkg").exists()
        and u not in list_fail
    ]
    list_gdf = [g for g in list_gdf if len(g)]
    if list_gdf:
        gdf_seg = gpd.GeoDataFrame(
            pd.concat(list_gdf, ignore_index=True),
            geometry="geometry",
            crs=list_gdf[0].crs,
        )
        gdf_seg = gdf_seg.dissolve(
            by=["river_id", "side", "k"],
            aggfunc={
                "chainage_from_m": "first",
                "chainage_to_m": "first",
                "envelope_area_m2": "sum",
                "n_cells": "sum",
            },
        ).reset_index()
    else:
        gdf_seg = gpd.GeoDataFrame(
            {
                c: pd.Series(dtype=d)
                for c, d in [
                    ("river_id", np.int64),
                    ("side", str),
                    ("k", np.int64),
                    ("chainage_from_m", float),
                    ("chainage_to_m", float),
                    ("envelope_area_m2", float),
                    ("n_cells", np.int64),
                ]
            },
            geometry=gpd.GeoSeries([], crs=gdf_rn.crs),
        )
    gdf_seg.insert(
        0,
        "subsegment_id",
        [
            f"{r}_{s}_{k}"
            for r, s, k in zip(gdf_seg["river_id"], gdf_seg["side"], gdf_seg["k"])
        ],
    )
    gis.Vector.write_gpkg(gdf_seg, out_dir / "subsegments.gpkg")
    merge_s = time.perf_counter() - t
    secs = {
        k: round(sum(r["secs"].get(k, 0.0) for r in list_res), 2)
        for k in ("hand", "envelope", "allocation", "write")
    }
    perf = {
        "units": len(df_units),
        "processed": len(list_task),
        "skipped_checkpoint": len(list_skip),
        "failed": list_fail,
        **{f"{k}_s": v for k, v in secs.items()},
        "units_wall_s": round(wall, 2),
        "merge_s": round(merge_s, 2),
        "total_s": round(time.perf_counter() - start, 2),
        "peak_rss_mb_worker": round(max([r["peak_rss_mb"] for r in list_res] or [0])),
        "peak_rss_mb_main": round(gis.Run.peak_rss_mb()),
    }
    (out_dir / "implied_segments_performance.json").write_text(
        json.dumps(perf, indent=1)
    )
    log.info(
        "subsegments.gpkg: %d subsegments of %.0f m (%d left, %d right), envelope %.2f km2 in %d units. "
        "Timings: HAND %.1f s, envelope %.1f s, allocation %.1f s, writing %.1f s (sums over units); "
        "units wall %.1f s, merge %.1f s, total %.1f s; peak RAM worker %d MB.",
        len(gdf_seg),
        length_m,
        int((gdf_seg["side"] == "left").sum()),
        int((gdf_seg["side"] == "right").sum()),
        gdf_seg["envelope_area_m2"].sum() / 1e6,
        len(df_units),
        secs["hand"],
        secs["envelope"],
        secs["allocation"],
        secs["write"],
        wall,
        merge_s,
        perf["total_s"],
        perf["peak_rss_mb_worker"],
    )
    if list_fail:
        raise SystemExit(f"Failed units: {list_fail}")
