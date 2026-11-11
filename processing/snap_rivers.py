# Author: Jakub Zapletal
# Date: 2026-11-04
# Stage E, step E1: snaps river lines to the thalweg on the DTM (direction check, least-cost path per overlapping
# chunk, stitching, smoothing, chainage); writes rivers_snapped.gpkg, river_points.gpkg and the snap plots.

import json
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import shapely
from matplotlib.patches import Rectangle
from osgeo import gdal
from scipy import ndimage
from shapely.ops import substring

import gis
from plots import COLOURS, DASHED, INK, Plots

# --- Run variables (no command-line parameters; DEFENSES_<NAME> environment variables override them) ---
CONFIG_PATH = "config/defenses.toml"
COUNTRY = "XX"
RUN_ID = "base"  # plot folder <output_dir>/<country>/plots/implied/<RUN_ID>
SHOW_PLOTS = False  # show up to implied.plots.show_max figures per type while running (needs a display)

KEYS = [
    "general.work_dir",
    "general.output_dir",
    "general.log_level",
    "general.n_workers",
    "general.gdal_cachemax_mb",
    "general.warp_memory_limit_mb",
    "general.gdal_threads_per_worker",
    "country.rivers",
    "country.water_polygons",
    "implied.units",
    *[
        f"implied.snap.{k}"
        for k in (
            "step_m",
            "resolution_m",
            "max_offset_m",
            "min_window_m",
            "cost_epsilon",
            "chunk_length_m",
            "chunk_overlap_m",
            "max_length_ratio",
            "wide_river_m",
            "direction_check_m",
            "max_nodata_share",
            "chainage_origin",
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


def snap_batch(task):
    # Worker: snaps the rivers of one batch; returns geometries, chunk records, profiles and timings.
    start = time.perf_counter()
    p = task["snap"]
    res, eps, off = p["resolution_m"], p["cost_epsilon"], p["max_offset_m"]
    size, overlap, step = p["chunk_length_m"], p["chunk_overlap_m"], p["step_m"]
    win = max(3, int(round(p["min_window_m"] / res)))
    win += 1 - win % 2  # odd, centred window
    ds = gdal.Open(task["vrt"])
    wkt = ds.GetProjection()
    water = shapely.from_wkb(task["water"]) if task["water"] else np.empty(0, object)
    tree = shapely.STRtree(water)
    secs = {"read": 0.0, "cost": 0.0, "path": 0.0, "stitch": 0.0}
    list_out = []
    for river_id, wkb in task["rivers"]:
        line = shapely.from_wkb(wkb)
        length = line.length

        # --- Direction: mean DTM elevation over the first and last direction_check_m ---
        list_end = []
        for a, b in (
            (0.0, min(p["direction_check_m"], length)),
            (max(length - p["direction_check_m"], 0.0), length),
        ):
            part = substring(line, a, b)
            x0, y0, x1, y1 = part.buffer(res).bounds
            bounds = (
                np.floor(x0 / res) * res,
                np.floor(y0 / res) * res,
                np.ceil(x1 / res) * res,
                np.ceil(y1 / res) * res,
            )
            t = time.perf_counter()
            z, gt = gis.Raster.read_window(ds, bounds, res, "average", task["warp_mb"])
            secs["read"] += time.perf_counter() - t
            xy = gis.Vector.interpolate_xy(
                part, np.arange(0.0, part.length + 1e-9, res)
            )
            r = np.clip(((gt[3] - xy[:, 1]) // res).astype(int), 0, z.shape[0] - 1)
            c = np.clip(((xy[:, 0] - gt[0]) // res).astype(int), 0, z.shape[1] - 1)
            values = z[r, c]
            list_end.append(
                float(np.mean(values[np.isfinite(values)]))
                if np.isfinite(values).any()
                else np.nan
            )
        flipped = bool(np.isfinite(list_end).all() and list_end[0] < list_end[1])
        if flipped:
            line = shapely.reverse(line)
            list_end = list_end[::-1]

        # --- Chunks of chunk_length_m overlapping by chunk_overlap_m ---
        n = max(1, int(np.ceil(max(length - overlap, 1e-9) / size)))
        list_chunk = []
        for k in range(n):
            a = k * size
            b = length if k == n - 1 else min((k + 1) * size + overlap, length)
            piece = substring(line, a, b)
            corr = piece.buffer(off)
            x0, y0, x1, y1 = corr.bounds
            bounds = (
                np.floor(x0 / res) * res,
                np.floor(y0 / res) * res,
                np.ceil(x1 / res) * res,
                np.ceil(y1 / res) * res,
            )
            t = time.perf_counter()
            z, gt = gis.Raster.read_window(ds, bounds, res, "average", task["warp_mb"])
            secs["read"] += time.perf_counter() - t
            t = time.perf_counter()
            shape = z.shape
            inside = (
                gis.Raster.rasterize_geometries(
                    [corr], [1], shape, gt, wkt, False, np.uint8
                )
                > 0
            )
            hits = tree.query(corr, predicate="intersects")
            wet = (
                shapely.intersection(shapely.union_all(water[hits]), corr)
                if len(hits)
                else None
            )
            wide = bool(
                wet is not None
                and wet.area / max(piece.length, res) > p["wide_river_m"]
            )
            if wide:
                # Wide river: follow the centre of the water polygons. Distances to their edge come from a window
                # padded by 2 * max_offset_m with a dry border, so the window edge never acts as a bank.
                pad = int(np.ceil(2 * off / res))
                gt_w = (gt[0] - pad * res, res, 0.0, gt[3] + pad * res, 0.0, -res)
                shape_w = (shape[0] + 2 * pad, shape[1] + 2 * pad)
                near = tree.query(
                    shapely.box(
                        gt_w[0],
                        gt_w[3] - shape_w[0] * res,
                        gt_w[0] + shape_w[1] * res,
                        gt_w[3],
                    ),
                    predicate="intersects",
                )
                mask_w = (
                    gis.Raster.rasterize_geometries(
                        water[near],
                        np.ones(len(near)),
                        shape_w,
                        gt_w,
                        wkt,
                        False,
                        np.uint8,
                    )
                    > 0
                )
                dist = ndimage.distance_transform_edt(np.pad(mask_w, 1))[
                    1 + pad : 1 + pad + shape[0], 1 + pad : 1 + pad + shape[1]
                ]
                mask = mask_w[pad : pad + shape[0], pad : pad + shape[1]]
                cost = np.where(mask, 1.0 / (dist * res + eps), np.nan)
            else:
                # Narrow river: height above the local minimum over min_window_m.
                low = ndimage.minimum_filter(
                    np.where(np.isfinite(z), z, np.inf), size=win, mode="nearest"
                )
                cost = z - low + eps
                cost[~np.isfinite(z)] = (
                    eps  # nodata in the corridor: water surface removed from the DTM
                )
            cost[~inside] = np.nan
            nan_share = float(np.mean(~np.isfinite(z[inside]))) if inside.any() else 1.0
            secs["cost"] += time.perf_counter() - t

            # Start and end: cells of the chunk end points, else the nearest passable cell.
            list_cell = []
            ok = np.isfinite(cost)
            radius = int(np.ceil(2 * off / res))
            for x, y in (piece.coords[0], piece.coords[-1]):
                r = min(max(int((gt[3] - y) // res), 0), shape[0] - 1)
                c = min(max(int((x - gt[0]) // res), 0), shape[1] - 1)
                r0, c0 = max(r - radius, 0), max(c - radius, 0)
                rr, cc = np.nonzero(ok[r0 : r + radius + 1, c0 : c + radius + 1])
                if len(rr) == 0:
                    list_cell.append(None)
                    continue
                j = int(np.argmin((rr + r0 - r) ** 2 + (cc + c0 - c) ** 2))
                list_cell.append((int(rr[j] + r0), int(cc[j] + c0)))
            t = time.perf_counter()
            path = np.empty(0, np.int64)
            if None not in list_cell and (wide or nan_share <= p["max_nodata_share"]):
                path, _ = gis.Hydro.least_cost_path(cost, list_cell[0], list_cell[1])
            secs["path"] += time.perf_counter() - t

            # Quality check: length ratio of the simplified path; offsets from the original chunk.
            pr, pc = np.divmod(path, shape[1])
            px, py = gt[0] + (pc + 0.5) * res, gt[3] - (pr + 0.5) * res
            reason, ratio, offsets = "", np.nan, np.zeros(0)
            if not wide and nan_share > p["max_nodata_share"]:
                reason = f"DTM nodata in {100 * nan_share:.0f}% of the corridor"
            elif len(path) < 2:
                reason = "no path"
            else:
                snapped = shapely.LineString(np.c_[px, py])
                ratio = shapely.simplify(snapped, res / 2).length / max(
                    piece.length, res
                )
                offsets = shapely.distance(piece, shapely.points(px, py))
                if ratio > p["max_length_ratio"]:
                    reason = f"length ratio {ratio:.2f} > {p['max_length_ratio']:.2f}"
            # DTM samples for the profile: (distance along the chunk line, elevation) on the original chunk and
            # on the path cells.
            u = np.arange(0.0, piece.length + 1e-9, res)
            xy = gis.Vector.interpolate_xy(piece, u)
            r = np.clip(((gt[3] - xy[:, 1]) // res).astype(int), 0, shape[0] - 1)
            c = np.clip(((xy[:, 0] - gt[0]) // res).astype(int), 0, shape[1] - 1)
            sample_o = np.c_[u, z[r, c]]
            sample_s = (
                sample_o
                if reason
                else np.c_[
                    np.r_[0.0, np.cumsum(np.hypot(np.diff(px), np.diff(py)))], z[pr, pc]
                ]
            )
            list_chunk.append(
                {
                    "k": k,
                    "from_m": a,
                    "to_m": b,
                    "mode": "wide" if wide else "narrow",
                    "kept_original": bool(reason),
                    "reason": reason,
                    "length_ratio": float(ratio),
                    "mean_offset_m": float(offsets.mean()) if len(offsets) else 0.0,
                    "max_offset_m": float(offsets.max()) if len(offsets) else 0.0,
                    "n_cells": int(len(offsets)),
                    "nan_share": nan_share,
                    "bounds": [float(v) for v in bounds],
                    "piece": shapely.to_wkb(piece),
                    "path": None if reason else shapely.to_wkb(snapped),
                    "geom": piece if reason else snapped,
                    "sample_o": sample_o,
                    "sample_s": sample_s,
                }
            )

        # --- Stitch at the middle of each overlap, then simplify and smooth ---
        t = time.perf_counter()
        list_xy, list_span = [], []
        for ch in list_chunk:
            g, k = ch["geom"], ch["k"]
            u0 = k * size + overlap / 2 if k > 0 else 0.0
            u1 = (k + 1) * size + overlap / 2 if k < n - 1 else length
            t0 = g.project(line.interpolate(u0)) if k > 0 else 0.0
            t1 = g.project(line.interpolate(u1)) if k < n - 1 else g.length
            # kept part: [u0, u1] on the original line, [t0, t1] on the chunk geometry
            list_span.append((u0 - ch["from_m"], u1 - ch["from_m"], t0, t1))
            if t1 > t0:
                list_xy.append(shapely.get_coordinates(substring(g, t0, t1)))
        joints = [part[-1] for part in list_xy[:-1]]
        xy = np.concatenate(list_xy)
        xy = xy[np.r_[True, np.any(np.diff(xy, axis=0) != 0, axis=1)]]
        smooth = gis.Vector.chaikin(
            shapely.simplify(shapely.LineString(xy), res / 2), 2, res
        )
        total = smooth.length

        # --- Profile against the chainage of the snapped line, binned by step_m ---
        # Kept part of every chunk mapped linearly onto the stretch of the snapped line between its joints.
        ends = np.maximum.accumulate(
            np.r_[0.0, [smooth.project(shapely.Point(j)) for j in joints], total]
        )
        bins = int(total // step) + 1
        profile = {"chainage_m": ((np.arange(bins) + 0.5) * step).astype(np.float32)}
        for name, key, col in (
            ("z_original", "sample_o", 0),
            ("z_snapped", "sample_s", 2),
        ):
            list_d, list_z = [], []
            for k, ch in enumerate(list_chunk):
                lo, hi = list_span[k][col], list_span[k][col + 1]
                if key == "sample_s" and ch["kept_original"]:
                    lo, hi = list_span[k][0], list_span[k][1]
                sample = ch[key]
                keep = (sample[:, 0] >= lo) & (sample[:, 0] <= hi) & (hi > lo)
                list_d.append(
                    ends[k]
                    + (sample[keep, 0] - lo)
                    / max(hi - lo, 1e-9)
                    * (ends[k + 1] - ends[k])
                )
                list_z.append(sample[keep, 1])
            d, zs = np.concatenate(list_d), np.concatenate(list_z)
            chain = total - d if p["chainage_origin"] == "downstream" else d
            idx = np.clip((chain // step).astype(int), 0, bins - 1)
            f = np.isfinite(zs)
            counts = np.bincount(idx[f], minlength=bins)
            sums = np.bincount(idx[f], weights=zs[f], minlength=bins)
            profile[name] = np.divide(
                sums, counts, out=np.full(bins, np.nan), where=counts > 0
            ).astype(np.float32)
        cuts = [
            float(total - d if p["chainage_origin"] == "downstream" else d)
            for d in ends[1:-1]
        ]
        edges = [0.0 if p["chainage_origin"] == "upstream" else total] + cuts
        edges += [total if p["chainage_origin"] == "upstream" else 0.0]
        profile["stitch_m"] = [float(v) for v in cuts]
        profile["kept_spans"] = [
            [float(min(edges[k], edges[k + 1])), float(max(edges[k], edges[k + 1]))]
            for k, ch in enumerate(list_chunk)
            if ch["kept_original"]
        ]
        secs["stitch"] += time.perf_counter() - t
        n_cells = sum(ch["n_cells"] for ch in list_chunk)
        list_out.append(
            {
                "river_id": river_id,
                "geometry": shapely.to_wkb(smooth),
                "flipped": flipped,
                "n_chunks": n,
                "n_kept_original": sum(ch["kept_original"] for ch in list_chunk),
                "mean_offset_m": (
                    sum(ch["mean_offset_m"] * ch["n_cells"] for ch in list_chunk)
                    / n_cells
                    if n_cells
                    else 0.0
                ),
                "max_offset_m": max(ch["max_offset_m"] for ch in list_chunk),
                "wide_river_chunks": sum(ch["mode"] == "wide" for ch in list_chunk),
                "length_m": total,
                "length_original_m": length,
                "upstream_elev_m": list_end[0],
                "downstream_elev_m": list_end[1],
                "chunks": [
                    {
                        k: v
                        for k, v in ch.items()
                        if k not in ("geom", "sample_o", "sample_s")
                    }
                    for ch in list_chunk
                ],
                "profile": profile,
            }
        )
    km = sum(r["length_original_m"] for r in list_out) / 1000
    return {
        "summary": f"{task['unit_id']}: {len(list_out)} rivers, {km:.1f} km in {time.perf_counter() - start:.1f} s",
        "rivers": list_out,
        "secs": secs,
        "peak_rss_mb": gis.Run.peak_rss_mb(),
    }


if __name__ == "__main__":
    run = gis.Run.env_overrides(
        CONFIG_PATH=CONFIG_PATH, COUNTRY=COUNTRY, RUN_ID=RUN_ID, SHOW_PLOTS=SHOW_PLOTS
    )
    cfg = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], KEYS)
    general, country, implied = cfg["general"], cfg["country"], cfg["implied"]
    snap, settings = implied["snap"], implied["plots"]
    work_dir = Path(general["work_dir"]) / run["COUNTRY"]
    out_dir = Path(general["output_dir"]) / run["COUNTRY"]
    log = gis.Run.logger("snap_rivers", out_dir / "logs", general["log_level"])
    gis.Config.warn_placeholders(cfg)
    gis.Run.configure_gdal(
        general["gdal_cachemax_mb"], general["gdal_threads_per_worker"]
    )
    start = time.perf_counter()
    vrt_path = work_dir / "dtm" / "dtm.vrt"
    if not vrt_path.exists():
        raise SystemExit("dtm.vrt is missing; run prepare_dtm.py first.")
    wkt = gdal.Open(str(vrt_path)).GetProjection()
    res, off = snap["resolution_m"], snap["max_offset_m"]

    # --- Rivers: snapped into their channel by protected_areas.py (rivers_channel.gpkg), else the rivers input ---
    channel_path = out_dir / "rivers_channel.gpkg"
    if channel_path.exists():
        gdf_rn = gpd.read_file(channel_path).to_crs(wkt)
        log.info("Rivers: rivers_channel.gpkg (%d features).", len(gdf_rn))
    else:
        gdf_rn = gis.Vector.read_merged(country["rivers"], wkt, keep_columns=True)
        log.info(
            "Rivers: the rivers input (%d features); rivers_channel.gpkg not found.",
            len(gdf_rn),
        )
    gdf_rn = gdf_rn[gdf_rn.geometry.notna() & ~gdf_rn.geometry.is_empty].explode(
        ignore_index=True
    )
    if "river_id" not in gdf_rn.columns:
        gdf_rn["river_id"] = np.arange(len(gdf_rn))
    if "name" not in gdf_rn.columns:
        gdf_rn["name"] = ""
    # Features joined end to end become one line; confluences and crossings stay breaks.
    lines = shapely.get_parts(
        shapely.line_merge(shapely.MultiLineString(list(gdf_rn.geometry.values)))
    )
    # Source features of every merged line: the line passes exactly through the midpoint of their first segment
    # (prepared line geometries against small boxes, linear in the number of vertices).
    xy, ix = shapely.get_coordinates(gdf_rn.geometry.values, return_index=True)
    first = np.flatnonzero(np.r_[True, ix[1:] != ix[:-1]])
    mid = (xy[first] + xy[first + 1]) / 2
    e = res / 10
    l_idx, f_idx = shapely.STRtree(
        shapely.box(mid[:, 0] - e, mid[:, 1] - e, mid[:, 0] + e, mid[:, 1] + e)
    ).query(lines, predicate="intersects")
    list_src, list_name = [], []
    for i in range(len(lines)):
        rows = f_idx[l_idx == i]
        list_src.append(",".join(str(v) for v in gdf_rn["river_id"].values[rows]))
        names = [
            str(v)
            for v in gdf_rn["name"].values[rows]
            if str(v) not in ("", "None", "nan")
        ]
        list_name.append(max(set(names), key=names.count) if names else "")
    parent = np.arange(len(lines))

    # --- Pilot units: rivers clipped to the unit windows ---
    if implied["units"]:
        win_path = out_dir / "units_window.gpkg"
        if not win_path.exists():
            raise SystemExit("units_window.gpkg is missing; run build_units.py first.")
        gdf_win = gpd.read_file(win_path).to_crs(wkt)
        gdf_win["unit_id"] = gdf_win["unit_id"].astype(str)
        list_unk = sorted(set(implied["units"]) - set(gdf_win["unit_id"]))
        if list_unk:
            raise SystemExit(
                f"implied.units not found in units_window.gpkg: {list_unk}"
            )
        area = shapely.union_all(
            gdf_win[gdf_win["unit_id"].isin(implied["units"])].geometry.values
        )
        lines, parent = shapely.get_parts(
            shapely.intersection(lines, area), return_index=True
        )
        log.info(
            "Pilot units %s: rivers clipped to their windows; chainage starts at the clipped ends.",
            implied["units"],
        )
    keep = (shapely.get_type_id(lines) == 1) & (shapely.length(lines) >= 2 * res)
    lines, parent = lines[keep], parent[keep]
    if not len(lines):
        raise SystemExit("No river line to snap.")
    list_src = [list_src[i] for i in parent]
    list_name = [list_name[i] for i in parent]
    lengths = shapely.length(lines)
    log.info(
        "%d river lines to snap, %.1f km (from %d features).",
        len(lines),
        lengths.sum() / 1000,
        len(gdf_rn),
    )

    # --- Water polygons for wide rivers ---
    gdf_wb = None
    if country["water_polygons"]:
        gdf_wb = gis.Vector.read_merged(country["water_polygons"], wkt)
        log.info("Water polygons: water_polygons input (%d).", len(gdf_wb))
    else:
        log.warning("No water polygons: the wide-river cost is disabled.")
    water = np.empty(0, object)
    if gdf_wb is not None and len(gdf_wb):
        water = shapely.get_parts(shapely.make_valid(gdf_wb.geometry.values))
        water = water[shapely.get_type_id(water) == 3]
        # Polygons within 3 * max_offset_m: the corridor plus the padding of the wide-river distance window.
        near = shapely.STRtree(water).query(
            shapely.buffer(lines, 3 * off), predicate="intersects"
        )[1]
        water = water[np.unique(near)]
        log.info("%d water polygons within 3 * max_offset_m of the rivers.", len(water))

    # --- Batches of similar total length (longest first), run in the process pool ---
    n_batch = max(1, min(len(lines), 4 * general["n_workers"]))
    loads, members = np.zeros(n_batch), [[] for _ in range(n_batch)]
    for i in np.argsort(-lengths):
        b = int(np.argmin(loads))
        members[b].append(int(i))
        loads[b] += lengths[i]
    water_tree = shapely.STRtree(water)
    list_task, list_w = [], []
    for b, idx in enumerate(members):
        if not idx:
            continue
        near = (
            np.unique(
                water_tree.query(
                    shapely.buffer(lines[idx], 3 * off), predicate="intersects"
                )[1]
            )
            if len(water)
            else np.empty(0, int)
        )
        list_task.append(
            {
                "unit_id": f"batch_{b:03d}",
                "rivers": [(i, shapely.to_wkb(lines[i])) for i in idx],
                "water": [shapely.to_wkb(g) for g in water[near]],
                "vrt": str(vrt_path),
                "warp_mb": general["warp_memory_limit_mb"],
                "snap": snap,
            }
        )
        list_w.append(float(loads[b]))
    t = time.perf_counter()
    list_res, list_fail = gis.Run.run_units(
        snap_batch,
        list_task,
        general["n_workers"],
        general["gdal_cachemax_mb"],
        general["gdal_threads_per_worker"],
        list_w,
        log,
    )
    wall_snap = time.perf_counter() - t
    list_riv = sorted(
        (r for res_b in list_res for r in res_b["rivers"]), key=lambda r: r["river_id"]
    )
    if not list_riv:
        raise SystemExit(f"No river was snapped; failed batches: {list_fail}")

    # --- Outputs: snapped rivers and points every step_m with chainage and downstream tangent ---
    t = time.perf_counter()
    geoms = shapely.from_wkb([r["geometry"] for r in list_riv])
    gdf_snap = gpd.GeoDataFrame(
        {
            "river_id": [r["river_id"] for r in list_riv],
            "source_river_ids": [list_src[r["river_id"]] for r in list_riv],
            "name": [list_name[r["river_id"]] for r in list_riv],
            **{
                key: [r[key] for r in list_riv]
                for key in (
                    "flipped",
                    "n_chunks",
                    "n_kept_original",
                    "mean_offset_m",
                    "max_offset_m",
                    "wide_river_chunks",
                    "length_m",
                    "length_original_m",
                    "upstream_elev_m",
                    "downstream_elev_m",
                )
            },
        },
        geometry=geoms,
        crs=wkt,
    )
    gis.Vector.write_gpkg(gdf_snap, out_dir / "rivers_snapped.gpkg")
    step = snap["step_m"]
    list_xy, list_tan, list_ch = [], [], []
    for g in geoms:
        total = g.length
        chain = np.arange(0.0, total + 1e-9, step)
        d = total - chain if snap["chainage_origin"] == "downstream" else chain
        list_xy.append(gis.Vector.interpolate_xy(g, d))
        # Tangent in the flow direction (lines run downstream) over one step around the point.
        tangent = gis.Vector.interpolate_xy(
            g, d + step / 2
        ) - gis.Vector.interpolate_xy(g, d - step / 2)
        list_tan.append(
            tangent / np.maximum(np.hypot(tangent[:, 0], tangent[:, 1]), 1e-12)[:, None]
        )
        list_ch.append(chain)
    counts = [len(c) for c in list_ch]
    xy, tangent = np.concatenate(list_xy), np.concatenate(list_tan)
    chain = np.concatenate(list_ch)
    pts = shapely.points(xy[:, 0], xy[:, 1])
    gdf_pts = gpd.GeoDataFrame(
        {
            "river_id": np.repeat(gdf_snap["river_id"].to_numpy(), counts),
            "chainage_m": chain,
            "tangent_dx": tangent[:, 0],
            "tangent_dy": tangent[:, 1],
        },
        geometry=pts,
        crs=wkt,
    )
    gis.Vector.write_gpkg(gdf_pts, out_dir / "river_points.gpkg")
    write_s = time.perf_counter() - t
    list_ck = [ch for r in list_riv for ch in r["chunks"]]
    log.info(
        "rivers_snapped.gpkg: %d rivers, %.1f km; %d chunks, %d kept original, %d wide; %d flipped; "
        "offset mean %.1f m, max %.1f m. river_points.gpkg: %d points every %.1f m.",
        len(gdf_snap),
        gdf_snap["length_m"].sum() / 1000,
        len(list_ck),
        sum(ch["kept_original"] for ch in list_ck),
        sum(ch["mode"] == "wide" for ch in list_ck),
        int(gdf_snap["flipped"].sum()),
        float(np.average(gdf_snap["mean_offset_m"], weights=gdf_snap["length_m"])),
        float(gdf_snap["max_offset_m"].max()),
        len(gdf_pts),
        step,
    )
    for r in list_riv:
        for ch in r["chunks"]:
            if ch["kept_original"]:
                log.warning(
                    "River %d chunk %d kept original: %s.",
                    r["river_id"],
                    ch["k"],
                    ch["reason"],
                )

    # --- Plots: snap_profile per river, snap_map per chunk ---
    t = time.perf_counter()
    if settings["enabled"]:
        plots = Plots(
            out_dir / "plots" / "implied" / run["RUN_ID"],
            settings,
            run["SHOW_PLOTS"],
            log,
        )
        plt = plots.plt
        units_of, units_of_ck = {}, {}
        if settings["plot_all_in_units"]:
            gdf_core = gpd.read_file(out_dir / "units_core.gpkg").to_crs(wkt)
            core_tree = shapely.STRtree(gdf_core.geometry.values)
            core_ids = gdf_core["unit_id"].astype(str).to_numpy()
            for r, g in zip(list_riv, geoms):
                units_of[f"river_{r['river_id']}"] = core_ids[
                    core_tree.query(g, predicate="intersects")
                ]
                for ch in r["chunks"]:
                    units_of_ck[f"river_{r['river_id']}_chunk_{ch['k']:03d}"] = (
                        core_ids[
                            core_tree.query(
                                shapely.from_wkb(ch["piece"]), predicate="intersects"
                            )
                        ]
                    )
        flag_ck = {
            f"river_{r['river_id']}_chunk_{ch['k']:03d}"
            for r in list_riv
            for ch in r["chunks"]
            if ch["kept_original"]
            or ch["mode"] == "wide"
            or ch["nan_share"] > 0
            or ch["max_offset_m"] >= off - res
        }
        label = (
            "Chainage from the downstream end (m)"
            if snap["chainage_origin"] == "downstream"
            else "Chainage from the upstream end (m)"
        )
        by_id = {f"river_{r['river_id']}": r for r in list_riv}
        list_sel = plots.select(
            "snap_profile",
            list(by_id),
            [
                f"river_{r['river_id']}"
                for r in list_riv
                if r["flipped"]
                or r["n_kept_original"]
                or r["wide_river_chunks"]
                or any(
                    f"river_{r['river_id']}_chunk_{ch['k']:03d}" in flag_ck
                    for ch in r["chunks"]
                )
            ],
            units_of,
        )
        for name in list_sel:
            r = by_id[name]
            prof = r["profile"]
            fig, ax = plt.subplots(figsize=(10, 4.4), layout="constrained")
            ax.plot(
                prof["chainage_m"],
                prof["z_original"],
                color=COLOURS["original"],
                ls=DASHED,
                label="original line",
            )
            ax.plot(
                prof["chainage_m"],
                prof["z_snapped"],
                color=COLOURS["snapped"],
                label="snapped line (thalweg)",
            )
            for i, s_m in enumerate(prof["stitch_m"]):
                ax.axvline(
                    s_m,
                    color=INK["muted"],
                    lw=0.6,
                    ls=":",
                    zorder=1,
                    label="chunk boundary" if i == 0 else None,
                )
            for i, (a, b) in enumerate(prof["kept_spans"]):
                ax.axvspan(
                    a,
                    b,
                    color=INK["grid"],
                    lw=0,
                    zorder=0,
                    label="chunk kept original" if i == 0 else None,
                )
            ax.set_xlabel(label)
            ax.set_ylabel("DTM elevation (m)")
            ax.set_title(
                f"snap_profile · river {r['river_id']} · run {run['RUN_ID']}\n"
                f"{r['length_m'] / 1000:.2f} km, {r['n_chunks']} chunks, {r['n_kept_original']} kept original, "
                f"{'flipped' if r['flipped'] else 'not flipped'}, offset mean {r['mean_offset_m']:.1f} m, "
                f"max {r['max_offset_m']:.1f} m"
            )
            ax.legend(loc="upper left", ncols=2)
            plots.save(
                fig,
                "snap_profile",
                name,
                {
                    key: r[key]
                    for key in (
                        "flipped",
                        "n_chunks",
                        "n_kept_original",
                        "wide_river_chunks",
                        "mean_offset_m",
                        "max_offset_m",
                        "length_m",
                    )
                },
            )
        by_ck = {
            f"river_{r['river_id']}_chunk_{ch['k']:03d}": (r, ch)
            for r in list_riv
            for ch in r["chunks"]
        }
        list_sel = plots.select("snap_map", list(by_ck), flag_ck, units_of_ck)
        ds = gdal.Open(str(vrt_path))
        wb_tree = shapely.STRtree(water)
        half = 4 * off  # detail crop: 8 * max_offset_m wide
        for name in list_sel:
            r, ch = by_ck[name]
            piece = shapely.from_wkb(ch["piece"])
            path_g = shapely.from_wkb(ch["path"]) if ch["path"] is not None else None
            # Detail at the path point farthest from the original chunk, else at the chunk middle.
            if path_g is not None:
                pxy = shapely.get_coordinates(path_g)
                gaps = shapely.distance(piece, shapely.points(pxy))
                centre, far = pxy[int(np.argmax(gaps))], float(gaps.max())
            else:
                centre = shapely.get_coordinates(
                    piece.interpolate(0.5, normalized=True)
                )[0]
                far = 0.0
            fig, (ax_o, ax_d) = plt.subplots(
                1, 2, figsize=(12, 5.8), width_ratios=[1.6, 1], layout="constrained"
            )
            detail = (
                centre[0] - half,
                centre[1] - half,
                centre[0] + half,
                centre[1] + half,
            )
            for ax, window in ((ax_o, ch["bounds"]), (ax_d, detail)):
                bounds_p, cell = Plots.window(window, res)
                Plots.hillshade_layer(
                    ax, ds, bounds_p, cell, general["warp_memory_limit_mb"]
                )
                box = shapely.box(*bounds_p)
                wet = shapely.get_parts(
                    shapely.intersection(
                        water[wb_tree.query(box, predicate="intersects")], box
                    )
                )
                for i, part in enumerate(wet[shapely.get_type_id(wet) == 3]):
                    ax.fill(*part.exterior.xy, color=COLOURS["water"], alpha=0.1, lw=0)
                    ax.plot(
                        *part.exterior.xy,
                        color=COLOURS["water"],
                        lw=0.8,
                        label="water polygon" if i == 0 else None,
                    )
                ax.plot(
                    *piece.xy,
                    color=COLOURS["original"],
                    ls=DASHED,
                    path_effects=Plots.halo(),
                    label="original line",
                )
                if path_g is not None:
                    ax.plot(
                        *path_g.xy,
                        color=COLOURS["snapped"],
                        path_effects=Plots.halo(),
                        label="least-cost path (cells)",
                    )
                ax.set_xlabel("Easting (m)")
            ax_o.set_ylabel("Northing (m)")
            ax_o.add_patch(
                Rectangle(
                    (centre[0] - half, centre[1] - half),
                    2 * half,
                    2 * half,
                    fill=False,
                    ec=INK["primary"],
                    lw=0.8,
                )
            )
            ax_o.set_title("Chunk overview (box: detail)")
            ax_d.set_title(
                f"Detail at the largest offset ({far:.1f} m)"
                if path_g is not None
                else "Detail at the chunk middle"
            )
            fig.suptitle(
                f"snap_map · river {r['river_id']} chunk {ch['k']} · run {run['RUN_ID']}\n"
                f"{ch['mode']} river, "
                f"{'kept original: ' + ch['reason'] if ch['kept_original'] else 'snapped'}, "
                f"length ratio {ch['length_ratio']:.2f}, offset mean {ch['mean_offset_m']:.1f} m, "
                f"max {ch['max_offset_m']:.1f} m, DTM nodata {100 * ch['nan_share']:.0f} % of the corridor",
                x=0.01,
                ha="left",
                fontsize=10,
            )
            handles, labels = ax_o.get_legend_handles_labels()
            fig.legend(handles, labels, loc="outside lower center", ncols=len(labels))
            plots.save(
                fig,
                "snap_map",
                name,
                {
                    "mode": ch["mode"],
                    "kept_original": ch["kept_original"],
                    "reason": ch["reason"],
                    "length_ratio": ch["length_ratio"],
                    "mean_offset_m": ch["mean_offset_m"],
                    "max_offset_m": ch["max_offset_m"],
                    "nan_share": ch["nan_share"],
                },
            )
        index = plots.finish()
        log.info("Plots: %s", index)
    plots_s = time.perf_counter() - t

    # --- Measured timings and memory ---
    secs = {
        k: sum(res_b["secs"][k] for res_b in list_res)
        for k in ("read", "cost", "path", "stitch")
    }
    km = gdf_snap["length_original_m"].sum() / 1000
    perf = {
        "rivers": len(gdf_snap),
        "river_km": round(km, 3),
        "chunks": len(list_ck),
        "batches": len(list_task),
        "failed_batches": list_fail,
        "n_workers": general["n_workers"],
        **{f"{k}_s": round(v, 3) for k, v in secs.items()},
        "worker_s_per_km": round(sum(secs.values()) / max(km, 1e-9), 4),
        "snap_wall_s": round(wall_snap, 2),
        "write_s": round(write_s, 2),
        "plots_s": round(plots_s, 2),
        "total_s": round(time.perf_counter() - start, 2),
        "peak_rss_mb_main": round(gis.Run.peak_rss_mb()),
        "peak_rss_mb_worker": round(max(res_b["peak_rss_mb"] for res_b in list_res)),
    }
    (out_dir / "snap_rivers_performance.json").write_text(json.dumps(perf, indent=1))
    log.info(
        "Timings: read %.2f s, cost %.2f s, path %.2f s, stitch %.2f s (%.3f s per km in the workers); "
        "snapping wall %.1f s, outputs %.1f s, plots %.1f s, total %.1f s; peak RAM main %d MB, worker %d MB.",
        secs["read"],
        secs["cost"],
        secs["path"],
        secs["stitch"],
        perf["worker_s_per_km"],
        wall_snap,
        write_s,
        plots_s,
        perf["total_s"],
        perf["peak_rss_mb_main"],
        perf["peak_rss_mb_worker"],
    )
    if list_fail:
        raise SystemExit(
            f"Failed batches: {list_fail}; their rivers are missing from the outputs."
        )
