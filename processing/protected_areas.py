# Author: Jakub Zapletal
# Date: 2026-11-03
# Derives areas protected by flood levees (protected areas) per processing unit and merges the units: priority-flood
# spill levels above the river (HAND on the undefended DTM) with and without barriers; margin = difference.
# File A: all rivers are water sources at once; file B: the same per river group (connected reaches of similar
# upstream area), the other rivers then being terrain.

import json
import logging
import math
import os
import shutil
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from osgeo import gdal
from scipy import ndimage
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from shapely.geometry import LineString, box

import gis

# --- Run variables (no command-line parameters; DEFENSES_<NAME> environment variables override them) ---
CONFIG_PATH = "config/defenses.toml"
COUNTRY = "XX"
FORCE = True  # recompute units that already have a done marker
UNITS = []  # unit ids to process, e.g. ["012_034", "012_035"]; empty = all units

KEYS = [
    "general.work_dir",
    "general.output_dir",
    "general.n_workers",
    "general.max_ram_gb_per_worker",
    "general.gdal_cachemax_mb",
    "general.warp_memory_limit_mb",
    "general.gdal_threads_per_worker",
    "general.log_level",
    "country.levees",
    "country.rivers",
    "country.river_guides",
    "country.river_upstream_area_field",
    "river_groups.max_join_ratio",
    "river_groups.max_length_m",
    "river_groups.min_upstream_area_km2",
    "levee_snap.densify_m",
    "levee_snap.max_shift_m",
    "levee_snap.side_m",
    "levee_snap.min_height_m",
    "levee_snap.level_tolerance_m",
    "levee_snap.smooth_vertices",
    "levee_snap.min_ridge_share",
    "levee_gaps.max_gap_m",
    "levee_gaps.crest_tolerance_m",
    "leaks.sample_spacing_m",
    "evidence.comparison_distance_m",
    "river_snap.densify_m",
    "river_snap.guide_search_m",
    "river_snap.guide_max_shift_m",
    "river_snap.window_m",
    "river_snap.max_shift_m",
    "river_snap.max_shift_without_guide_m",
    "river_snap.min_angle_deg",
    "river_snap.angle_passes",
    "river_snap.spike_angle_deg",
    "river_snap.spike_length_m",
    "river_snap.straight_angle_deg",
    "river_snap.straight_length_factor",
    "river_snap.straight_max_deviation_m",
    "river_snap.spike_passes",
    "river_snap.band_rows",
    "country.working_resolution_m",
    "country.bytes_per_cell",
    "country.crest_box_size_m",
    "country.merge_overlap_ratio",
    "country.drop_contradicted_truncated",
    "country.contradiction_ratio",
    "country.id_rounding_m",
    "country.keep_intermediate",
    "country.protected_areas.gap_tolerance_m",
    "country.protected_areas.touch_tolerance_m",
    "country.protected_areas.self_close_min_length_factor",
    "country.protected_areas.crest_search_radius_m",
    "country.protected_areas.barrier_removal_width_m",
    "country.protected_areas.min_margin_m",
    "country.protected_areas.max_level_m",
    "country.protected_areas.min_area_m2",
    "country.protected_areas.gated_crossings",
    "country.protected_areas.fill_nodata_holes",
    "country.protected_areas.stream_burn_m",
    "country.protected_areas.stream_burn_buffer_m",
    "country.protected_areas.preserve_river_cells",
    "country.protected_areas.cap_river_holes",
    "country.protected_areas.max_fill_distance_m",
    "country.rp.levee_adjacent_distance_m",
]
UNIT_COLUMNS = {
    "unit_label": "int64",
    "area_m2": "float64",
    "protection_level_m": "float64",
    "margin_median_m": "float64",
    "margin_max_m": "float64",
    "river_ref_median_m": "float64",
    "entry_x": "float64",
    "entry_y": "float64",
    "source_x": "float64",
    "source_y": "float64",
    "truncated": "bool",
    "touches_corridor_edge": "bool",
    "touches_dtm_edge": "bool",
}
LEAK_COLUMNS = {
    "cause": "str",
    "x": "float64",
    "y": "float64",
    "levee_length_m": "float64",
    "n_samples": "int64",
    "line_ids": "str",
    "crest_hand_median_m": "float64",
    "level_behind_median_m": "float64",
    "lost_margin_median_m": "float64",
}
OUTPUT_COLUMNS = [
    "protected_area_id",
    "area_m2",
    "protection_level_m",
    "margin_median_m",
    "margin_max_m",
    "river_ref_median_m",
    "entry_x",
    "entry_y",
    "source_x",
    "source_y",
    "truncated",
    "touches_corridor_edge",
    "touches_dtm_edge",
    "no_levee",
    "unit_id",
    "unit_label",
    "geometry",
]


def protected_regions(f, barrier, domain, dtm, gt, wkt, p):
    """
    Protected regions of one flood pair and their statistics (shared by all water sources together and by each
    river group)
    :param f: dict margin, spill_undef, spill_def, river_ref, order_def, seed_def (2D arrays of the unit window)
    :param barrier: barrier cells
    :param domain: processing domain cells
    :param dtm: terrain (NaN = missing coverage)
    :param gt: geotransform of the window
    :param wkt: CRS of the window
    :param p: parameters (min_margin_m, max_level_m, min_area_m2)
    :return: (label array, GeoDataFrame of the regions with UNIT_COLUMNS)
    """
    margin, spill_undef, spill_def = f["margin"], f["spill_undef"], f["spill_def"]
    river_ref, order_def, seed_def = f["river_ref"], f["order_def"], f["seed_def"]
    res = gt[1]
    cell_area = res * res
    shape = margin.shape
    protected = (
        (margin >= p["min_margin_m"]) & (spill_undef <= p["max_level_m"]) & ~barrier
    )
    labels, n_candidates = ndimage.label(protected)
    sizes = np.bincount(labels.ravel())
    keep = np.flatnonzero(sizes * cell_area >= p["min_area_m2"])
    keep = keep[keep > 0]
    relabel = np.zeros(n_candidates + 1, np.int32)
    relabel[keep] = np.arange(1, len(keep) + 1)
    labels = relabel[labels]
    ids = np.arange(1, len(keep) + 1)
    if len(ids):
        # Statistics over protected area cells only (scipy.ndimage label statistics scan the whole window).
        cells = np.flatnonzero(labels.ravel())
        df_cells = pd.DataFrame(
            {
                "label": labels.ravel()[cells],
                "spill": spill_def.ravel()[cells],
                "margin": margin.ravel()[cells],
                "river_ref": river_ref.ravel()[cells],
                "order": order_def.ravel()[cells],
                "cell": cells,
            }
        )
        grouped = df_cells.groupby("label")
        # Entry = first protected area cell reached by the defended flood; source = seed river cell of its path.
        entry_rows, entry_cols = np.divmod(
            df_cells["cell"].to_numpy()[
                grouped["order"].idxmin().reindex(ids).to_numpy()
            ],
            shape[1],
        )
        seed_rows, seed_cols = np.divmod(seed_def[entry_rows, entry_cols], shape[1])
        del df_cells["order"], df_cells["cell"]
        border = np.unique(
            np.concatenate([labels[0], labels[-1], labels[:, 0], labels[:, -1]])
        )
        touching = np.unique(labels[ndimage.binary_dilation(~domain) & (labels > 0)])
        # Missing DTM coverage acts as a wall: such protected areas may be closed artificially.
        at_dtm_edge = np.unique(
            labels[ndimage.binary_dilation(np.isnan(dtm) & domain) & (labels > 0)]
        )
        df_stats = pd.DataFrame(
            {
                "unit_label": ids,
                "area_m2": sizes[keep] * cell_area,
                "protection_level_m": grouped["spill"].min().reindex(ids).to_numpy(),
                "margin_median_m": grouped["margin"].median().reindex(ids).to_numpy(),
                "margin_max_m": grouped["margin"].max().reindex(ids).to_numpy(),
                "river_ref_median_m": grouped["river_ref"]
                .median()
                .reindex(ids)
                .to_numpy(),
                "entry_x": gt[0] + (entry_cols + 0.5) * res,
                "entry_y": gt[3] - (entry_rows + 0.5) * res,
                "source_x": gt[0] + (seed_cols + 0.5) * res,
                "source_y": gt[3] - (seed_rows + 0.5) * res,
                "truncated": np.isin(ids, border),
                "touches_corridor_edge": np.isin(ids, touching),
                "touches_dtm_edge": np.isin(ids, at_dtm_edge),
            }
        ).astype(UNIT_COLUMNS)
        del df_cells, grouped
    else:
        df_stats = pd.DataFrame(
            {c: pd.Series(dtype=d) for c, d in UNIT_COLUMNS.items()}
        )
    gdf_poly = gis.Raster.polygonize_labels(labels, gt, wkt).rename(
        columns={"label": "unit_label"}
    )
    return labels, gdf_poly.merge(df_stats, on="unit_label", how="inner")


def merge_units(gdf_cand, gdf_core, gdf_win, list_done, country, code, wkt, gdf_bl):
    """
    Merges the protected area candidates of all finished units: one area seen by several windows is kept once,
    areas truncated everywhere are united or dropped when neighbouring windows contradict them; ids and no_levee
    :param gdf_cand: candidates of all units (unit columns and unit_id)
    :param gdf_core: unit cores indexed by unit_id
    :param gdf_win: unit windows indexed by unit_id
    :param list_done: finished unit ids
    :param country: country configuration
    :param code: country code (id prefix)
    :param wkt: CRS
    :param gdf_bl: barrier lines
    :return: (merged GeoDataFrame, number dropped as contradicted)
    """
    if len(gdf_cand):
        geoms = gdf_cand.geometry.values
        areas = shapely.area(geoms)
        rep_pt = shapely.point_on_surface(geoms)
        gdf_cand["in_core"] = shapely.covers(
            gdf_core.geometry.reindex(gdf_cand.unit_id).values, rep_pt
        )
        gdf_cand["edge_distance"] = shapely.distance(
            rep_pt,
            shapely.boundary(gdf_win.geometry.reindex(gdf_cand.unit_id).values),
        )
        # The same protected area seen from several unit windows: candidates of different units overlapping by at
        # least merge_overlap_ratio of the smaller one form one group.
        left, right = gdf_cand.sindex.query(geoms, predicate="intersects")
        pair = (left < right) & (
            gdf_cand.unit_id.values[left] != gdf_cand.unit_id.values[right]
        )
        left, right = left[pair], right[pair]
        overlap = shapely.area(
            shapely.intersection(geoms[left], geoms[right])
        ) / np.minimum(areas[left], areas[right])
        same = overlap >= country["merge_overlap_ratio"]
        graph = coo_matrix(
            (np.ones(same.sum()), (left[same], right[same])),
            shape=(len(gdf_cand),) * 2,
        )
        _, gdf_cand["group"] = connected_components(graph, directed=False)
        group_size = gdf_cand.groupby("group")["group"].transform("size")
        list_keep = [gdf_cand[group_size == 1]]
        for _, gdf_grp in gdf_cand[group_size > 1].groupby("group"):
            gdf_full = gdf_grp[~gdf_grp["truncated"].astype(bool)]
            if len(gdf_full):
                # Prefer a complete candidate whose representative point lies in its own core, then the most
                # central one (largest distance from its window edge).
                list_keep.append(
                    gdf_full.sort_values(
                        ["in_core", "edge_distance"], ascending=False
                    ).iloc[[0]]
                )
            else:
                # Truncated in every window: union of the pieces (recall), lowest protection level among them.
                gdf_mrg = gdf_grp.loc[[gdf_grp.geometry.area.idxmax()]].copy()
                lowest = gdf_grp.loc[gdf_grp["protection_level_m"].idxmin()]
                gdf_mrg = gdf_mrg.set_geometry(
                    gpd.GeoSeries(
                        [shapely.union_all(gdf_grp.geometry.values)],
                        index=gdf_mrg.index,
                        crs=wkt,
                    )
                )
                gdf_mrg["area_m2"] = gdf_mrg.geometry.area
                for column in [
                    "protection_level_m",
                    "entry_x",
                    "entry_y",
                    "source_x",
                    "source_y",
                ]:
                    gdf_mrg[column] = lowest[column]
                for column in ["touches_corridor_edge", "touches_dtm_edge"]:
                    gdf_mrg[column] = gdf_grp[column].astype(bool).any()
                list_keep.append(gdf_mrg)
        gdf_final = gpd.GeoDataFrame(
            pd.concat(list_keep, ignore_index=True), geometry="geometry", crs=wkt
        )
        # A protected area truncated in every window may exist only because window edges
        # closed it. Units whose windows see the same area and found nothing there have
        # the better view (an edge can only add walls); drop it if they cover enough.
        n_contra = 0
        if country["drop_contradicted_truncated"]:
            gdf_win_done = gdf_win.loc[list_done]
            keep_rows = np.ones(len(gdf_final), bool)
            for k in np.flatnonzero(gdf_final["truncated"].astype(bool).to_numpy()):
                geometry = gdf_final.geometry.iloc[k]
                members = set(
                    gdf_cand.unit_id[gdf_cand.group == gdf_final.group.iloc[k]]
                )
                gdf_oth = gdf_win_done[
                    ~gdf_win_done.index.isin(members)
                    & gdf_win_done.intersects(geometry)
                ]
                hits = gdf_cand.unit_id.values[
                    gdf_cand.sindex.query(geometry, predicate="intersects")
                ]
                gdf_sil = gdf_oth[~gdf_oth.index.isin(hits)]
                if len(gdf_sil):
                    covered = geometry.intersection(
                        shapely.union_all(gdf_sil.geometry.values)
                    ).area / max(geometry.area, 1e-9)
                    keep_rows[k] = covered < country["contradiction_ratio"]
            n_contra = int((~keep_rows).sum())
            gdf_final = gdf_final[keep_rows].reset_index(drop=True)
    else:
        gdf_final = gdf_cand.copy()
        n_contra = 0
    rep_pt = shapely.point_on_surface(gdf_final.geometry.values)
    step = country["id_rounding_m"]
    ids = pd.Series(
        [
            f"{code}_{int(round(x / step) * step)}_{int(round(y / step) * step)}"
            for x, y in shapely.get_coordinates(rep_pt)
        ],
        dtype=str,
    )
    dup_rank = ids.groupby(ids).cumcount()
    gdf_final["protected_area_id"] = np.where(
        dup_rank > 0, ids + "_" + (dup_rank + 1).astype(str), ids
    )
    distance = country["rp"]["levee_adjacent_distance_m"]
    gdf_final["no_levee"] = (
        gis.Vector.length_near_boundary(gdf_final, gdf_bl, distance) == 0
    )
    return gdf_final, n_contra


def process_unit(task):
    start = time.perf_counter()
    times, counts = {}, {}
    unit_id, p = task["unit_id"], task["params"]
    unit_dir = Path(task["unit_dir"])
    unit_dir.mkdir(parents=True, exist_ok=True)
    (unit_dir / "done.json").unlink(missing_ok=True)
    log = logging.getLogger(f"unit.{unit_id}")
    log.setLevel(logging.INFO)
    log.propagate = False
    handler = logging.FileHandler(unit_dir / "unit.log", mode="w", encoding="utf-8")
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    log.addHandler(handler)
    res, wkt = task["resolution"], task["wkt"]
    cell_area = res * res
    bounds = (task["win_xmin"], task["win_ymin"], task["win_xmax"], task["win_ymax"])
    core = box(
        task["core_xmin"], task["core_ymin"], task["core_xmax"], task["core_ymax"]
    )

    # --- Terrain at working resolution and processing domain ---
    t = time.perf_counter()
    dtm, gt = gis.Raster.read_window(
        task["dtm_vrt"], bounds, res, "average", task["warp_memory_mb"]
    )
    times["read_terrain"] = time.perf_counter() - t
    shape = dtm.shape
    t = time.perf_counter()
    gdf_dom = gpd.read_file(task["domain_path"], where=f"unit_id = '{unit_id}'")
    domain = (
        gis.Raster.rasterize_geometries(
            gdf_dom.geometry,
            np.ones(len(gdf_dom)),
            shape,
            gt,
            wkt,
            True,
            np.uint8,
        )
        > 0
    )
    data = ~np.isnan(dtm)
    dtm[~domain] = np.nan
    counts["domain_cells"] = int(domain.sum())
    gdf_rn = gpd.read_file(task["rivers_used"], bbox=bounds)
    river_cells = (
        gis.Raster.rasterize_geometries(
            gdf_rn.geometry, np.ones(len(gdf_rn)), shape, gt, wkt, False, np.uint8
        )
        > 0
    )
    # river group + 1 per river cell (0 = no river), for the protected areas per group
    river_group = (
        gis.Raster.rasterize_geometries(
            gdf_rn.geometry, gdf_rn["group_id"] + 1, shape, gt, wkt, False, np.int32
        )
        if "group_id" in gdf_rn.columns
        else None
    )
    if p["fill_nodata_holes"]:
        # Nodata inside the DTM coverage (lidar water surfaces, data gaps) would act as a wall and
        # leave rivers without cells; it is interpolated, also where it reaches the window edge.
        # Nodata outside the coverage stays a wall (flagged by touches_dtm_edge).
        # Coverage = morphological closing of the data cells (two chamfer distance passes): gaps
        # narrower than twice max_fill_distance_m are inside it, a straight coverage edge stays.
        # Cells beyond the window are unknown and count as no data (padding), so that a coverage
        # edge near the window border is not closed against the border.
        radius = max(p["max_fill_distance_m"] / res, 1.0)
        pad = int(np.ceil(radius)) + 1
        grown = gis.Hydro.distance_cells(np.pad(data, pad), radius + 1) <= radius
        covered = (gis.Hydro.distance_cells(~grown, radius + 1) > radius)[
            pad:-pad, pad:-pad
        ]
        del grown
        gaps = np.isnan(dtm) & domain & covered
        del covered
        if gaps.any():
            rows_any = np.flatnonzero(gaps.any(axis=1))
            cols_any = np.flatnonzero(gaps.any(axis=0))
            q0, q1 = max(rows_any[0] - 2, 0), min(rows_any[-1] + 3, shape[0])
            p0, p1 = max(cols_any[0] - 2, 0), min(cols_any[-1] + 3, shape[1])
            limit = max(p["max_fill_distance_m"] / res, 1.0)
            # Distance from gap cells to the nearest valid cell: half-width of the gaps.
            to_valid = gis.Hydro.distance_cells(~np.isnan(dtm[q0:q1, p0:p1]), limit)
            reach = int(np.ceil(to_valid[gaps[q0:q1, p0:p1]].max())) + 2
            del to_valid
            before = dtm
            dtm, unfilled = gis.Raster.fill_nodata_cells(before, gaps, gt, wkt, reach)
            counts["filled_gap_cells"] = int(gaps.sum()) - unfilled
            river_gap = gaps & river_cells
            if p["cap_river_holes"] and river_gap.any():
                # Water surfaces: interpolation from the banks lies above the water. River cells are
                # interpolated along each river line between valid river cells (the river profile).
                along = np.zeros(shape, bool)
                for geometry in gdf_rn.geometry:
                    for part in shapely.get_parts(geometry):
                        steps = np.append(
                            np.arange(0.0, part.length, res / 2), part.length
                        )
                        xy = shapely.get_coordinates(
                            shapely.line_interpolate_point(part, steps)
                        )
                        cols = np.floor((xy[:, 0] - gt[0]) / res).astype(np.int64)
                        rows = np.floor((gt[3] - xy[:, 1]) / res).astype(np.int64)
                        inside = (
                            (rows >= 0)
                            & (rows < shape[0])
                            & (cols >= 0)
                            & (cols < shape[1])
                        )
                        steps, rows, cols = steps[inside], rows[inside], cols[inside]
                        known = ~np.isnan(before[rows, cols])
                        missing = river_gap[rows, cols]
                        if known.any() and missing.any():
                            dtm[rows[missing], cols[missing]] = np.interp(
                                steps[missing],
                                steps[known],
                                before[rows[known], cols[known]],
                            )
                            along[rows[missing], cols[missing]] = True
                counts["river_cells_interpolated_along"] = int(along.sum())
                # Without valid river cells (whole river nodata in the window), and for the other cells
                # of gaps crossed by a river, the fill is capped at the lowest valid cell within reach:
                # the waterline at the edge of the gap.
                gap_parts, _ = ndimage.label(gaps)
                river_parts = np.setdiff1d(np.unique(gap_parts[river_gap]), [0])
                target = np.isin(gap_parts, river_parts) & ~along
                del gap_parts
                if target.any():
                    rows_any = np.flatnonzero(target.any(axis=1))
                    cols_any = np.flatnonzero(target.any(axis=0))
                    r0, r1 = max(rows_any[0] - reach, 0), min(
                        rows_any[-1] + 1 + reach, shape[0]
                    )
                    c0, c1 = max(cols_any[0] - reach, 0), min(
                        cols_any[-1] + 1 + reach, shape[1]
                    )
                    original = before[r0:r1, c0:c1]
                    lowest = ndimage.minimum_filter(
                        np.where(np.isnan(original), np.inf, original),
                        size=2 * reach + 1,
                        mode="constant",
                        cval=np.inf,
                    )
                    cap = target[r0:r1, c0:c1] & np.isfinite(lowest)
                    window = dtm[r0:r1, c0:c1]
                    window[cap] = np.minimum(window[cap], lowest[cap])
                    counts["capped_gap_cells"] = int(cap.sum())
            del before, river_gap
        del gaps
    del data
    valid = ~np.isnan(dtm)
    times["domain"] = time.perf_counter() - t

    # --- Barrier cells and crest elevations ---
    t = time.perf_counter()
    gdf_bl = gpd.read_file(task["barrier_lines"], bbox=bounds)
    bridged = (
        gdf_bl["bridged"].astype(bool).to_numpy() if len(gdf_bl) else np.zeros(0, bool)
    )
    gdf_lv, gdf_br = gdf_bl[~bridged], gdf_bl[bridged]
    levee_cells = np.zeros(shape, bool)
    crest_full = np.full(shape, np.nan, np.float32)
    counts["crest_boxes"] = 0
    if len(gdf_lv):
        levee_cells = (
            gis.Raster.rasterize_geometries(
                gdf_lv.geometry,
                np.ones(len(gdf_lv)),
                shape,
                gt,
                wkt,
                True,
                np.uint8,
            )
            > 0
        )
        # Crest = maximum of the native DTM near the line. Native pixels are read only in boxes around levees:
        # gdal.Warp with "max" per working cell (overviews disabled), then a maximum filter of the search radius.
        box_cells = max(1, int(round(p["crest_box_size_m"] / res)))
        radius_cells = max(1, int(math.ceil(p["crest_search_radius_m"] / res)))
        yy, xx = np.mgrid[
            -radius_cells : radius_cells + 1, -radius_cells : radius_cells + 1
        ]
        footprint = np.hypot(yy, xx) * res <= p["crest_search_radius_m"] + 1e-9
        footprint[radius_cells, radius_cells] = True
        box_shape = (-(-shape[0] // box_cells), -(-shape[1] // box_cells))
        box_gt = (gt[0], gt[1] * box_cells, 0.0, gt[3], 0.0, gt[5] * box_cells)
        needed = (
            gis.Raster.rasterize_geometries(
                shapely.buffer(
                    gdf_lv.geometry.values, p["crest_search_radius_m"] + res
                ),
                np.ones(len(gdf_lv)),
                box_shape,
                box_gt,
                wkt,
                True,
                np.uint8,
            )
            > 0
        )
        for box_row in range(box_shape[0]):
            on = np.flatnonzero(needed[box_row])
            # Boxes in a row separated by at most one empty box are read together.
            for run_cols in (
                np.split(on, np.flatnonzero(np.diff(on) > 2) + 1) if len(on) else []
            ):
                r0, r1 = box_row * box_cells, min((box_row + 1) * box_cells, shape[0])
                c0, c1 = run_cols[0] * box_cells, min(
                    (run_cols[-1] + 1) * box_cells, shape[1]
                )
                read_bounds = (
                    gt[0] + (c0 - radius_cells) * res,
                    gt[3] - (r1 + radius_cells) * res,
                    gt[0] + (c1 + radius_cells) * res,
                    gt[3] - (r0 - radius_cells) * res,
                )
                cell_max, _ = gis.Raster.read_window(
                    task["dtm_vrt"],
                    read_bounds,
                    res,
                    "max",
                    task["warp_memory_mb"],
                    overview="NONE",
                )
                cell_max = ndimage.maximum_filter(
                    np.where(np.isnan(cell_max), -np.inf, cell_max),
                    footprint=footprint,
                    mode="constant",
                    cval=-np.inf,
                )
                block = cell_max[
                    radius_cells : radius_cells + r1 - r0,
                    radius_cells : radius_cells + c1 - c0,
                ]
                crest_full[r0:r1, c0:c1] = np.where(np.isfinite(block), block, np.nan)
                counts["crest_boxes"] += 1
    times["crest_read"] = time.perf_counter() - t

    t = time.perf_counter()
    crest = np.where(levee_cells, crest_full, np.nan).astype(np.float32)
    bridge_cells = np.zeros(shape, bool)
    if len(gdf_br):
        bridge_id = gis.Raster.rasterize_geometries(
            gdf_br.geometry,
            np.arange(1, len(gdf_br) + 1),
            shape,
            gt,
            wkt,
            True,
            np.int32,
        )
        bridge_cells = (bridge_id > 0) & ~levee_cells
        ends = np.array([[g.coords[0], g.coords[-1]] for g in gdf_br.geometry])
        end_cols = np.floor((ends[..., 0] - gt[0]) / res).astype(np.int64)
        end_rows = np.floor((gt[3] - ends[..., 1]) / res).astype(np.int64)
        inside = (
            (end_rows >= 0)
            & (end_rows < shape[0])
            & (end_cols >= 0)
            & (end_cols < shape[1])
        )
        end_crest = np.full(end_cols.shape, np.nan, np.float32)
        end_crest[inside] = crest_full[end_rows[inside], end_cols[inside]]
        # A missing crest at one end takes the value of the other end.
        end_crest = np.where(np.isnan(end_crest), end_crest[:, ::-1], end_crest)
        rows, cols = np.nonzero(bridge_cells)
        seg = bridge_id[rows, cols] - 1
        a = ends[seg, 0]
        ab = ends[seg, 1] - a
        along = np.clip(
            (
                (gt[0] + (cols + 0.5) * res - a[:, 0]) * ab[:, 0]
                + (gt[3] - (rows + 0.5) * res - a[:, 1]) * ab[:, 1]
            )
            / np.maximum((ab**2).sum(axis=1), 1e-12),
            0,
            1,
        )
        crest[rows, cols] = end_crest[seg, 0] + along * (
            end_crest[seg, 1] - end_crest[seg, 0]
        )
    barrier = levee_cells | bridge_cells
    # Close diagonal-only contacts so that water cannot pass between two barrier cells touching at a corner.
    height = np.where(valid, dtm, -np.inf)
    counts["added_4conn"] = 0
    while True:
        diag = (
            barrier[:-1, :-1] & barrier[1:, 1:] & ~barrier[:-1, 1:] & ~barrier[1:, :-1]
        )
        anti = (
            barrier[:-1, 1:] & barrier[1:, :-1] & ~barrier[:-1, :-1] & ~barrier[1:, 1:]
        )
        if not diag.any() and not anti.any():
            break
        i, j = np.nonzero(diag)
        first = height[i, j + 1] >= height[i + 1, j]
        barrier[np.where(first, i, i + 1), np.where(first, j + 1, j)] = True
        i, j = np.nonzero(anti)
        first = height[i, j] >= height[i + 1, j + 1]
        barrier[np.where(first, i, i + 1), np.where(first, j, j + 1)] = True
        counts["added_4conn"] += int(diag.sum() + anti.sum())
    del height
    gaps = barrier & np.isnan(crest)
    crest[gaps] = crest_full[gaps]
    gaps = barrier & np.isnan(crest)
    if gaps.any() and (barrier & ~np.isnan(crest)).any():
        nearest = ndimage.distance_transform_edt(
            np.isnan(crest), return_distances=False, return_indices=True
        )
        crest[gaps] = crest[nearest[0][gaps], nearest[1][gaps]]
        del nearest
    gaps = barrier & np.isnan(crest)
    crest[gaps] = dtm[gaps]
    times["barriers"] = time.perf_counter() - t

    # --- Rivers, enclosed streams, water sources ---
    t = time.perf_counter()
    stream = river_cells & valid
    del river_cells
    barrier_full = barrier.copy()
    counts.update(
        {
            "stream_cells": int(stream.sum()),
            "barrier_cells": int(barrier.sum()),
            "crossing_cells": int((barrier & stream).sum()),
        }
    )
    # Zones separated by barriers (crossings counted as barrier). Streams entering the window edge belong to the
    # main network; streams in other zones are enclosed by levees (e.g. tributaries with flap gates).
    zones, _ = ndimage.label(~barrier_full & valid)
    # Rivers on the window edge or next to nodata (missing DTM coverage) enter from outside.
    stream_edge = stream & ndimage.binary_dilation(~valid)
    stream_edge[[0, -1], :] |= stream[[0, -1], :]
    stream_edge[:, [0, -1]] |= stream[:, [0, -1]]
    outer_zones = np.setdiff1d(np.unique(zones[stream_edge]), [0])
    stream_zones = np.setdiff1d(np.unique(zones[stream]), [0])
    if not len(outer_zones) and len(stream_zones):
        outer_zones = stream_zones[
            [
                np.argmax(
                    np.bincount(zones[stream], minlength=zones.max() + 1)[stream_zones]
                )
            ]
        ]
    inner_zones = np.setdiff1d(stream_zones, outer_zones)
    if p["gated_crossings"]:
        seeds = stream & np.isin(zones, outer_zones)
    else:
        barrier &= ~stream
        seeds = stream
    inner_labels = np.zeros(shape, np.int32)
    for k, zone in enumerate(inner_zones, start=1):
        inner_labels[zones == zone] = k
    gdf_encl = gis.Raster.polygonize_labels(inner_labels, gt, wkt).rename(
        columns={"label": "zone"}
    )
    gdf_encl["area_m2"] = (
        np.bincount(inner_labels.ravel())[gdf_encl["zone"]] * cell_area
    )
    gdf_encl["stream_cells"] = np.bincount(
        inner_labels[stream], minlength=len(inner_zones) + 1
    )[gdf_encl["zone"]]
    gdf_encl["gated"] = p["gated_crossings"]
    gdf_encl["unit_id"] = unit_id
    counts["enclosed_stream_zones"] = len(inner_zones)
    if len(inner_zones):
        log.warning(
            "%d levee-enclosed areas contain river cells; treated as %s.",
            len(inner_zones),
            (
                "gated (not water sources)"
                if p["gated_crossings"]
                else "unprotected (river wins)"
            ),
        )
    times["rivers"] = time.perf_counter() - t

    if stream.any():
        # --- Undefended DTM: barrier cells and their surroundings removed and interpolated ---
        t = time.perf_counter()
        dtm_undef = dtm
        if barrier_full.any():
            reach = int(math.ceil(p["barrier_removal_width_m"] / res)) + 1
            rows_any = np.flatnonzero(barrier_full.any(axis=1))
            cols_any = np.flatnonzero(barrier_full.any(axis=0))
            r0, r1 = max(rows_any[0] - reach, 0), min(
                rows_any[-1] + 1 + reach, shape[0]
            )
            c0, c1 = max(cols_any[0] - reach, 0), min(
                cols_any[-1] + 1 + reach, shape[1]
            )
            removal = np.zeros(shape, bool)
            removal[r0:r1, c0:c1] = (
                ndimage.distance_transform_edt(
                    ~barrier_full[r0:r1, c0:c1], sampling=res
                )
                <= p["barrier_removal_width_m"]
            )
            removal &= valid
            if p["preserve_river_cells"]:
                # River cells next to a levee keep their elevation (the bed); only crossing cells,
                # where the DTM shows the levee body, are interpolated.
                removal &= ~(stream & ~barrier_full)
            dtm_undef, unfilled = gis.Raster.fill_nodata_cells(
                dtm,
                removal,
                gt,
                wkt,
                math.ceil(4 * p["barrier_removal_width_m"] / res) + 10,
            )
            counts["removal_cells"], counts["removal_unfilled"] = (
                int(removal.sum()),
                unfilled,
            )
            del removal
        times["undefended"] = time.perf_counter() - t

        # --- HAND on the undefended DTM ---
        t = time.perf_counter()
        # outlets of the drainage: the water sources (closed streams in gated mode drain nowhere of their own)
        hand_undef, river_ref, n_unreached = gis.Hydro.hand_above_streams(
            dtm_undef,
            seeds,
            p["stream_burn_m"],
            p["stream_burn_buffer_m"] / res,
            valid,
        )
        if n_unreached:
            log.warning(
                "%d cells are not connected to any stream cell; using the nearest stream cell.",
                n_unreached,
            )
        counts["unreached_cells"] = n_unreached
        times["hand"] = time.perf_counter() - t

        # --- Spill levels with and without barriers ---
        t = time.perf_counter()
        surface_def = hand_undef.copy()
        surface_def[barrier] = np.maximum(
            crest[barrier] - river_ref[barrier], hand_undef[barrier]
        )
        spill_undef, _, _ = gis.Hydro.priority_flood(hand_undef, seeds)
        spill_def, seed_def, order_def = gis.Hydro.priority_flood(surface_def, seeds)
        del surface_def
        margin = spill_def - spill_undef
        times["flood"] = time.perf_counter() - t
    else:
        log.warning("No river cells in the window; no protected areas.")
        dtm_undef = dtm
        river_ref = hand_undef = spill_undef = spill_def = margin = np.full(
            shape, np.nan, np.float32
        )
        seed_def = order_def = np.full(shape, -1, np.int32)

    # --- Protected areas: all water sources together ---
    t = time.perf_counter()
    flood = {
        "margin": margin,
        "spill_undef": spill_undef,
        "spill_def": spill_def,
        "river_ref": river_ref,
        "order_def": order_def,
        "seed_def": seed_def,
    }
    labels, gdf_pa = protected_regions(flood, barrier, domain, dtm, gt, wkt, p)
    gdf_pa["unit_id"] = unit_id
    counts["protected_areas"] = len(gdf_pa)
    times["protected_areas"] = time.perf_counter() - t

    # --- Leak diagnostics (all water sources together): levees checked every leak_spacing_m inside the core. A
    # levee holds water back where one side is protected. Otherwise, when the terrain behind it lies above its
    # crest, it is not needed; else the flood path to the cell behind it is followed back to where the water gets
    # past the levee: an open river crossing, a low crest, the free end of a line, or a water source behind it. ---
    t = time.perf_counter()
    list_rows = []
    status = {"protected": 0, "terrain_above_crest": 0, "leak": 0, "unclear": 0}
    gdf_bl_core = gdf_bl[gdf_bl.intersects(core)] if len(gdf_bl) else gdf_bl
    if stream.any() and len(gdf_bl_core):
        owner, s_xy, s_norm, weight = gis.Vector.line_samples(
            shapely.intersection(gdf_bl_core.geometry.values, core),
            p["leak_spacing_m"],
            res,
        )
        side = p["barrier_removal_width_m"] + 2 * res
        list_cell = []
        for xy in (s_xy, s_xy + side * s_norm, s_xy - side * s_norm):
            cols = np.floor((xy[:, 0] - gt[0]) / res).astype(np.int64)
            rows = np.floor((gt[3] - xy[:, 1]) / res).astype(np.int64)
            inside = (rows >= 0) & (rows < shape[0]) & (cols >= 0) & (cols < shape[1])
            list_cell.append(np.where(inside, rows * shape[1] + cols, -1))
        i0, ia, ib = list_cell
        flat_zone = zones.ravel()
        use = (i0 >= 0) & (ia >= 0) & (ib >= 0)
        use[use] = (
            barrier.ravel()[i0[use]]
            & (flat_zone[ia[use]] > 0)
            & (flat_zone[ib[use]] > 0)
        )
        owner, s_xy, weight = owner[use], s_xy[use], weight[use]
        i0, ia, ib = i0[use], ia[use], ib[use]
        crest_h = (crest - river_ref).ravel()[i0]
        flat_margin = margin.ravel()
        protected = (flat_margin[ia] >= p["min_margin_m"]) | (
            flat_margin[ib] >= p["min_margin_m"]
        )
        # behind = the side not reached by rivers entering the window; else the side farther from them
        is_outer = np.zeros(int(flat_zone.max()) + 1, bool)
        is_outer[outer_zones] = True
        outer_seeds = seeds & is_outer[zones]
        to_river = gis.Hydro.distance_cells(
            outer_seeds if outer_seeds.any() else seeds, float(sum(shape))
        ).ravel()
        out_a, out_b = is_outer[flat_zone[ia]], is_outer[flat_zone[ib]]
        behind = np.where(
            out_a != out_b,
            np.where(out_b, ia, ib),
            np.where(to_river[ia] >= to_river[ib], ia, ib),
        )
        del to_river, outer_seeds
        level = spill_def.ravel()[behind]
        terrain = ~protected & (
            hand_undef.ravel()[behind] >= crest_h - p["min_margin_m"]
        )
        leak = ~protected & ~terrain & np.isfinite(level)
        status["protected"] = int(protected.sum())
        status["terrain_above_crest"] = int(terrain.sum())
        # events on the flood path: 1 = open river crossing, 2 = raised barrier cell, 3 = free end of a line
        event = np.zeros(shape, np.int8)
        event[barrier_full & ~barrier] = 1
        event[barrier] = 2
        list_g = list(gdf_bl.geometry)
        ends = shapely.points([g.coords[k] for g in list_g for k in (0, -1)])
        end_owner = np.repeat(np.arange(len(list_g)), 2)
        e_idx, l_idx = gdf_bl.sindex.query(
            ends, predicate="dwithin", distance=p["touch_tolerance_m"]
        )
        tied = np.zeros(len(ends), bool)
        tied[e_idx[l_idx != end_owner[e_idx]]] = True
        closed = shapely.distance(ends[0::2], ends[1::2]) <= p["touch_tolerance_m"]
        tied |= np.repeat(closed, 2)
        free_xy = shapely.get_coordinates(ends[~tied])
        end_id = np.zeros(shape, np.int32)
        if len(free_xy):
            end_id = gis.Raster.rasterize_geometries(
                shapely.buffer(shapely.points(free_xy), side),
                np.arange(1, len(free_xy) + 1),
                shape,
                gt,
                wkt,
                True,
                np.int32,
            )
            event[(end_id > 0) & (event == 0)] = 3
        code, cell = gis.Hydro.first_event_on_path(
            order_def, seeds, event, behind[leak]
        )
        del event
        # open crossings next to each zone, for paths that start at a water source inside a closed zone
        cross = np.flatnonzero((barrier_full & stream).ravel())
        c_rows, c_cols = np.divmod(cross, shape[1])
        list_cz = []
        for dr, dc in ((-1, 0), (1, 0), (0, -1), (0, 1)):
            r, c = c_rows + dr, c_cols + dc
            ok = (r >= 0) & (r < shape[0]) & (c >= 0) & (c < shape[1])
            list_cz.append(
                pd.DataFrame(
                    {"cell": cross[ok], "zone": flat_zone[r[ok] * shape[1] + c[ok]]}
                )
            )
        df_cz = pd.concat(list_cz).query("zone > 0").drop_duplicates()
        cause = np.full(len(code), "", dtype=object)
        loc = np.full(len(code), -1, np.int64)
        loc_xy = np.full((len(code), 2), np.nan)
        cause[code == 1], loc[code == 1] = "river_crossing", cell[code == 1]
        cause[code == 2], loc[code == 2] = "low_crest", cell[code == 2]
        at_end = code == 3
        cause[at_end] = "line_end"
        loc_xy[at_end] = free_xy[end_id.ravel()[cell[at_end]] - 1]
        zone_b = flat_zone[behind[leak]]
        for k in np.flatnonzero(code == -1):
            if is_outer[zone_b[k]]:
                continue  # the side choice or an opening not found on the path: unclear
            near = df_cz.loc[df_cz.zone == zone_b[k], "cell"].to_numpy()
            if len(near):
                sr, sc = divmod(int(cell[k]), shape[1])
                nr, nc = np.divmod(near, shape[1])
                cause[k] = "river_crossing"
                loc[k] = near[np.argmin((nr - sr) ** 2 + (nc - sc) ** 2)]
            else:
                cause[k], loc[k] = "water_source_behind", cell[k]
        has_cell = loc >= 0
        lr, lc = np.divmod(loc[has_cell], shape[1])
        loc_xy[has_cell, 0] = gt[0] + (lc + 0.5) * res
        loc_xy[has_cell, 1] = gt[3] - (lr + 0.5) * res
        status["leak"] = int((cause != "").sum())
        status["unclear"] = int(leak.sum()) - status["leak"]
        line_ids = gdf_bl_core["line_id"].to_numpy()[owner]
        known = cause != ""
        df_leak = pd.DataFrame(
            {
                "cause": cause[known],
                "x": loc_xy[known, 0],
                "y": loc_xy[known, 1],
                "line_id": line_ids[leak][known],
                "weight": weight[leak][known],
                "crest_h": crest_h[leak][known],
                "level": level[leak][known],
            }
        )
        df_terr = pd.DataFrame(
            {
                "cause": "terrain_above_crest",
                "x": s_xy[terrain, 0],
                "y": s_xy[terrain, 1],
                "line_id": line_ids[terrain],
                "weight": weight[terrain],
                "crest_h": crest_h[terrain],
                "level": hand_undef.ravel()[behind[terrain]],
            }
        )
        # terrain above the crest: one record per line, at its first such sample
        df_terr[["x", "y"]] = df_terr.groupby("line_id")[["x", "y"]].transform("first")
        df_s = pd.concat([df_leak, df_terr], ignore_index=True)
        if len(df_s):
            step = p["id_rounding_m"]
            df_s["kx"] = np.round(df_s["x"] / step)
            df_s["ky"] = np.round(df_s["y"] / step)
            df_s["lost"] = df_s["crest_h"] - df_s["level"]
            list_rows.append(
                df_s.groupby(["cause", "kx", "ky"], sort=False)
                .agg(
                    x=("x", "first"),
                    y=("y", "first"),
                    levee_length_m=("weight", "sum"),
                    n_samples=("weight", "size"),
                    line_ids=("line_id", lambda v: ",".join(map(str, sorted(set(v))))),
                    crest_hand_median_m=("crest_h", "median"),
                    level_behind_median_m=("level", "median"),
                    lost_margin_median_m=("lost", "median"),
                )
                .reset_index()
                .drop(columns=["kx", "ky"])
            )
    df_lk = (
        pd.concat(list_rows, ignore_index=True)
        if list_rows
        else pd.DataFrame({c: pd.Series(dtype=d) for c, d in LEAK_COLUMNS.items()})
    )
    gdf_leaks = gpd.GeoDataFrame(
        df_lk.astype(LEAK_COLUMNS),
        geometry=gpd.points_from_xy(df_lk["x"], df_lk["y"]),
        crs=wkt,
    )
    gdf_leaks["unit_id"] = unit_id
    counts.update({f"levee_samples_{k}": v for k, v in status.items()})
    counts["leaks"] = len(gdf_leaks)
    del zones
    times["leaks"] = time.perf_counter() - t

    # --- Protected areas per river group: water sources and drainage outlets are the cells of one group; the
    # other rivers are terrain (a flood does not come from all rivers at once) ---
    t = time.perf_counter()
    list_grp = []
    pad_g = int(math.ceil(p["comparison_distance_m"] / res)) + 1
    shutil.rmtree(unit_dir / "groups", ignore_errors=True)
    if river_group is not None and stream.any():
        for value in np.setdiff1d(np.unique(river_group[seeds]), [0]):
            seeds_g = seeds & (river_group == value)
            hand_g, ref_g, _ = gis.Hydro.hand_above_streams(
                dtm_undef,
                seeds_g,
                p["stream_burn_m"],
                p["stream_burn_buffer_m"] / res,
                valid,
            )
            surface = hand_g.copy()
            surface[barrier] = np.maximum(
                crest[barrier] - ref_g[barrier], hand_g[barrier]
            )
            spill_u, _, _ = gis.Hydro.priority_flood(hand_g, seeds_g)
            spill_d, seed_g, order_g = gis.Hydro.priority_flood(surface, seeds_g)
            del surface
            labels_g, gdf_g = protected_regions(
                {
                    "margin": spill_d - spill_u,
                    "spill_undef": spill_u,
                    "spill_def": spill_d,
                    "river_ref": ref_g,
                    "order_def": order_g,
                    "seed_def": seed_g,
                },
                barrier,
                domain,
                dtm,
                gt,
                wkt,
                p,
            )
            del spill_u, spill_d, ref_g, seed_g, order_g
            gdf_g["group_id"] = int(value) - 1
            if len(gdf_g):
                # labels and HAND of the group for the later steps (file B), cropped to its protected areas and
                # their comparison zones
                rows_any = np.flatnonzero((labels_g > 0).any(axis=1))
                cols_any = np.flatnonzero((labels_g > 0).any(axis=0))
                r0, r1 = max(rows_any[0] - pad_g, 0), min(
                    rows_any[-1] + 1 + pad_g, shape[0]
                )
                c0, c1 = max(cols_any[0] - pad_g, 0), min(
                    cols_any[-1] + 1 + pad_g, shape[1]
                )
                gt_g = (gt[0] + c0 * res, res, 0.0, gt[3] - r0 * res, 0.0, -res)
                gid = int(value) - 1
                gis.Raster.write_array(
                    gis.Run.unit_layer_path(unit_dir, "labels.tif", gid),
                    labels_g[r0:r1, c0:c1].astype(np.int32),
                    gt_g,
                    wkt,
                    nodata=0,
                )
                gis.Raster.write_array(
                    gis.Run.unit_layer_path(unit_dir, "hand_undef.tif", gid),
                    hand_g[r0:r1, c0:c1],
                    gt_g,
                    wkt,
                )
            del hand_g, labels_g
            list_grp.append(gdf_g)
    gdf_grp = (
        gpd.GeoDataFrame(pd.concat(list_grp, ignore_index=True), crs=wkt)
        if list_grp
        else gpd.GeoDataFrame(
            {c: pd.Series(dtype=d) for c, d in UNIT_COLUMNS.items()}
            | {"group_id": pd.Series(dtype="int64")},
            geometry=gpd.GeoSeries([], crs=wkt),
        )
    )
    gdf_grp["unit_id"] = unit_id
    counts["river_groups"] = len(list_grp)
    counts["protected_areas_by_group"] = len(gdf_grp)
    times["groups"] = time.perf_counter() - t

    # --- Crest statistics of barrier lines inside the core ---
    t = time.perf_counter()
    crest_attr = np.where(barrier, crest, np.nan)
    crest_hand = crest_attr - river_ref
    list_crest = []
    gdf_bl_core = gdf_bl[gdf_bl.intersects(core)] if len(gdf_bl) else gdf_bl
    for line_id, geometry in (
        zip(gdf_bl_core["line_id"], gdf_bl_core.geometry) if len(gdf_bl_core) else []
    ):
        xy = shapely.get_coordinates(
            shapely.segmentize(geometry.intersection(core), res / 2)
        )
        cols = np.clip(
            np.floor((xy[:, 0] - gt[0]) / res).astype(np.int64), 0, shape[1] - 1
        )
        rows = np.clip(
            np.floor((gt[3] - xy[:, 1]) / res).astype(np.int64), 0, shape[0] - 1
        )
        values, hand_values = crest_attr[rows, cols], crest_hand[rows, cols]
        list_crest.append(
            (
                line_id,
                float(np.nanmin(values)) if np.isfinite(values).any() else np.nan,
                (
                    float(np.nanmin(hand_values))
                    if np.isfinite(hand_values).any()
                    else np.nan
                ),
            )
        )
    times["barrier_crest"] = time.perf_counter() - t

    # --- Outputs ---
    t = time.perf_counter()
    gis.Raster.write_array(unit_dir / "hand_undef.tif", hand_undef, gt, wkt)
    gis.Raster.write_array(
        unit_dir / "labels.tif", labels.astype(np.int32), gt, wkt, nodata=0
    )
    if p["keep_intermediate"]:
        for name, array in [
            ("dtm", dtm),
            ("dtm_undefended", dtm_undef),
            ("river_ref", river_ref),
            ("spill_undef", spill_undef),
            ("spill_def", spill_def),
            ("margin", margin),
            ("barrier", barrier.astype(np.uint8)),
            ("crest", np.where(barrier, crest, np.nan)),
        ]:
            gis.Raster.write_array(unit_dir / f"{name}.tif", array, gt, wkt)
    gis.Vector.write_gpkg(gdf_pa, unit_dir / "protected_areas.gpkg")
    if river_group is not None:
        gis.Vector.write_gpkg(gdf_grp, unit_dir / "protected_areas_by_group.gpkg")
    gis.Vector.write_gpkg(gdf_encl, unit_dir / "enclosed_streams.gpkg")
    gis.Vector.write_gpkg(gdf_leaks, unit_dir / "leaks.gpkg")
    pd.DataFrame(
        list_crest, columns=["line_id", "crest_min_m", "crest_hand_min_m"]
    ).to_csv(unit_dir / "barrier_crest.csv", index=False)
    times["write"] = time.perf_counter() - t

    seconds = time.perf_counter() - start
    peak_mb = gis.Run.peak_rss_mb()
    window_cells = shape[0] * shape[1]
    done = {
        "unit_id": unit_id,
        "seconds": seconds,
        "times": {k: round(v, 3) for k, v in times.items()},
        "peak_rss_mb": peak_mb,
        "window_cells": window_cells,
        "bytes_per_cell": peak_mb * 2**20 / window_cells,
        **counts,
    }
    done["summary"] = (
        f"unit {unit_id}: {len(gdf_pa)} protected areas in {seconds:.0f} s (read {times['read_terrain']:.1f} s, "
        f"crest {times['crest_read']:.1f} s in {counts['crest_boxes']} boxes, hand {times.get('hand', 0):.1f} s, "
        f"flood {times.get('flood', 0):.1f} s, {counts['leaks']} leaks {times['leaks']:.1f} s, "
        f"{counts['river_groups']} river groups {times['groups']:.1f} s, "
        f"write {times['write']:.1f} s), peak RSS {peak_mb / 1024:.2f} GB "
        f"({done['bytes_per_cell']:.0f} B per window cell)"
    )
    marker = unit_dir / "done.json.tmp"
    marker.write_text(json.dumps(done, indent=2))
    os.replace(marker, unit_dir / "done.json")
    log.info(done["summary"])
    log.removeHandler(handler)
    handler.close()
    return done


if __name__ == "__main__":
    run = gis.Run.env_overrides(
        CONFIG_PATH=CONFIG_PATH, COUNTRY=COUNTRY, FORCE=FORCE, UNITS=UNITS
    )
    cfg = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], KEYS)
    general, country = cfg["general"], cfg["country"]
    params = dict(country["protected_areas"])
    params.update(
        {
            "crest_box_size_m": country["crest_box_size_m"],
            "keep_intermediate": country["keep_intermediate"],
            "leak_spacing_m": cfg["leaks"]["sample_spacing_m"],
            "comparison_distance_m": cfg["evidence"]["comparison_distance_m"],
            "id_rounding_m": country["id_rounding_m"],
        }
    )
    work_dir = Path(general["work_dir"]) / run["COUNTRY"]
    out_dir = Path(general["output_dir"]) / run["COUNTRY"]
    log = gis.Run.logger("protected_areas", out_dir / "logs", general["log_level"])
    gis.Config.warn_placeholders(cfg)
    if params["stream_burn_m"] <= 0:
        log.warning(
            "stream_burn_m = 0 reproduces the previous HAND routing, in which drainage "
            "paths can run along channel cells parallel to the rasterized river and reach "
            "it only far downstream (river_ref too low, HAND too high). See README, bug B1."
        )
    gis.Run.configure_gdal(
        general["gdal_cachemax_mb"], general["gdal_threads_per_worker"]
    )
    vrt_path = work_dir / "dtm" / "dtm.vrt"
    list_need = [
        vrt_path,
        out_dir / "units.csv",
        out_dir / "units_core.gpkg",
        out_dir / "units_window.gpkg",
        work_dir / "units_domain.gpkg",
    ]
    if not all(path.exists() for path in list_need):
        raise SystemExit(
            f"Missing inputs {[str(p) for p in list_need if not p.exists()]}; "
            "run prepare_dtm.py and build_units.py first."
        )
    wkt = gdal.Open(str(vrt_path)).GetProjection()
    resolution = country["working_resolution_m"]
    df_units = pd.read_csv(out_dir / "units.csv", dtype={"unit_id": str})
    all_windows = (
        df_units.win_xmin.min(),
        df_units.win_ymin.min(),
        df_units.win_xmax.max(),
        df_units.win_ymax.max(),
    )

    # --- Levee lines moved onto the crest seen in the DTM (once per country, shared by all units) ---
    t = time.perf_counter()
    gap, touch = params["gap_tolerance_m"], params["touch_tolerance_m"]
    snap, lg = cfg["levee_snap"], cfg["levee_gaps"]
    crest_path = out_dir / "levees_crest.gpkg"
    if crest_path.exists() and not run["FORCE"]:
        gdf_lv = gpd.read_file(crest_path).to_crs(wkt)
        log.info(
            "Levee lines: %d from levees_crest.gpkg (set FORCE to snap again).",
            len(gdf_lv),
        )
    else:
        gdf_lv = gis.Vector.read_merged(
            country["levees"], wkt, bbox=all_windows
        ).explode(index_parts=False)
        gdf_lv = gdf_lv[
            (gdf_lv.geom_type == "LineString") & (gdf_lv.length > 0)
        ].reset_index(drop=True)
        geoms, df_ridge, crest_stats = gis.Hydro.snap_to_crest(
            gdf_lv.geometry.values, str(vrt_path), snap, touch=touch
        )
        df_ridge["no_ridge"] = df_ridge["ridge_share"] < snap["min_ridge_share"]
        gdf_lv = gpd.GeoDataFrame(
            pd.concat([gdf_lv[["source"]], df_ridge], axis=1),
            geometry=gpd.GeoSeries(geoms, crs=wkt),
        )
        gis.Vector.write_gpkg(gdf_lv, crest_path)
        length_km = gdf_lv.geometry.length / 1000
        log.info(
            "Levee lines onto the DTM crest (levees_crest.gpkg): %d lines, %d vertices, %.0f %% of them on a ridge "
            "at least %.1f m high; shift mean %.1f m, p95 %.1f m (at most %.0f m); %d lines moved; %d lines "
            "(%.1f km) with a ridge on less than %.0f %% of their length (no_ridge), %d of them (%.1f km) without "
            "any; %.1f s.",
            crest_stats["lines"],
            crest_stats["vertices"],
            100 * crest_stats.get("ridge_share", 0.0),
            snap["min_height_m"],
            crest_stats.get("shift_mean_m", 0.0),
            crest_stats.get("shift_p95_m", 0.0),
            snap["max_shift_m"],
            crest_stats.get("lines_moved", 0),
            int(gdf_lv["no_ridge"].sum()),
            length_km[gdf_lv["no_ridge"]].sum(),
            100 * snap["min_ridge_share"],
            crest_stats.get("lines_without_ridge", 0),
            length_km[gdf_lv["ridge_share"] == 0].sum(),
            time.perf_counter() - t,
        )

    # --- Barrier lines for the whole country: bridge gaps (deterministic, shared by all units) ---
    t = time.perf_counter()
    geoms = gdf_lv.geometry.values
    starts, ends = shapely.get_point(geoms, 0), shapely.get_point(geoms, -1)
    bridges = {}
    own_gap = shapely.distance(starts, ends)
    # A nearly closed ring digitised as one line; the length condition keeps short lines from closing on themselves.
    self_close = (
        (own_gap > touch)
        & (own_gap <= gap)
        & (shapely.length(geoms) > params["self_close_min_length_factor"] * gap)
    )
    for i in np.flatnonzero(self_close):
        key = tuple(
            sorted(
                [
                    (round(starts[i].x, 3), round(starts[i].y, 3)),
                    (round(ends[i].x, 3), round(ends[i].y, 3)),
                ]
            )
        )
        bridges[key] = (
            gdf_lv.at[i, "source"],
            LineString([starts[i], ends[i]]),
            "short",
        )
    open_lines = np.flatnonzero(~self_close)
    points = np.concatenate([starts[open_lines], ends[open_lines]])
    owners = np.concatenate([open_lines, open_lines])
    n_long_cand = n_long = 0
    if len(points):
        point_idx, line_idx = gdf_lv.sindex.query(
            points, predicate="dwithin", distance=max(lg["max_gap_m"], gap)
        )
        other = line_idx != owners[point_idx]
        df_pairs = pd.DataFrame({"point": point_idx[other], "line": line_idx[other]})
        df_pairs["distance"] = shapely.distance(
            points[df_pairs.point], geoms[df_pairs.line]
        )
        connected = df_pairs.groupby("point")["distance"].min() <= touch
        df_pairs = df_pairs[~df_pairs.point.map(connected).astype(bool).to_numpy()]
        df_near = df_pairs.loc[df_pairs.groupby("point")["distance"].idxmin()]
        df_near = df_near[df_near["distance"] <= gap]
        for point, line in zip(df_near.point, df_near.line):
            endpoint, other_line = points[point], geoms[line]
            list_ends = [starts[line], ends[line]]
            list_dist = [endpoint.distance(e) for e in list_ends]
            # Prefer endpoint-to-endpoint, otherwise connect to the nearest point on the other line.
            target = (
                list_ends[int(np.argmin(list_dist))]
                if min(list_dist) <= gap
                else shapely.line_interpolate_point(
                    other_line, shapely.line_locate_point(other_line, endpoint)
                )
            )
            key = tuple(
                sorted(
                    [
                        (round(endpoint.x, 3), round(endpoint.y, 3)),
                        (round(target.x, 3), round(target.y, 3)),
                    ]
                )
            )
            list_src = sorted(
                {gdf_lv.at[owners[point], "source"], gdf_lv.at[line, "source"]}
            )
            bridges[key] = (
                "|".join(list_src),
                LineString([endpoint, target]),
                "short",
            )
        # Longer gaps up to max_gap_m are bridged where the DTM holds the crest along the gap (the levee exists,
        # its line is missing); otherwise the removal of the levee next to the line ends would open the levee.
        # Candidates of a free end in order of distance; the first that passes is taken.
        df_long = df_pairs[~df_pairs.point.isin(df_near.point)].sort_values(
            ["point", "distance"], kind="stable"
        )
        if len(df_long):
            pts = points[df_long.point.to_numpy()]
            lines_c = geoms[df_long.line.to_numpy()]
            targets = shapely.line_interpolate_point(
                lines_c, shapely.line_locate_point(lines_c, pts)
            )
            segs = shapely.linestrings(
                np.stack(
                    [shapely.get_coordinates(pts), shapely.get_coordinates(targets)],
                    axis=1,
                )
            )
            end_low, along_min = gis.Hydro.bridge_profile(
                segs, str(vrt_path), snap["max_shift_m"]
            )
            passed = np.isfinite(end_low) & (
                along_min >= end_low - lg["crest_tolerance_m"]
            )
            n_long_cand = int(df_long.point.nunique())
            first_ok = (
                df_long.assign(k=np.arange(len(df_long)))[passed]
                .groupby("point")["k"]
                .first()
                .to_numpy()
            )
            # the bridge follows the crest between the two fixed ends
            snapped, _, _ = gis.Hydro.snap_to_crest(
                segs[first_ok], str(vrt_path), snap, fix_ends=True
            )
            for k, geometry in zip(first_ok, snapped):
                a, b = shapely.get_coordinates(segs[k])
                key = tuple(sorted([tuple(np.round(a, 3)), tuple(np.round(b, 3))]))
                if key in bridges:
                    continue
                list_src = sorted(
                    {
                        gdf_lv.at[owners[df_long.point.iat[k]], "source"],
                        gdf_lv.at[df_long.line.iat[k], "source"],
                    }
                )
                bridges[key] = ("|".join(list_src), geometry, "dtm")
                n_long += 1
    gdf_bl = gpd.GeoDataFrame(
        pd.concat(
            [
                gdf_lv.drop(columns="geometry").assign(bridged=False, bridge_type=""),
                pd.DataFrame(
                    {
                        "source": [b[0] for b in bridges.values()],
                        "bridged": True,
                        "bridge_type": [b[2] for b in bridges.values()],
                        "no_ridge": False,
                    }
                ),
            ],
            ignore_index=True,
        ),
        geometry=list(geoms) + [b[1] for b in bridges.values()],
        crs=wkt,
    )
    gdf_bl["no_ridge"] = gdf_bl["no_ridge"].astype(bool)
    gdf_bl.insert(0, "line_id", np.arange(len(gdf_bl)))
    gis.Vector.write_gpkg(gdf_bl, work_dir / "barrier_lines.gpkg")
    log.info(
        "Barrier lines: %d levee lines, %d bridged gaps (%d self-closed rings; %d gaps longer than gap_tolerance_m, "
        "up to %.0f m, where the DTM holds the crest, of %d free ends with such a candidate) in %.1f s.",
        len(gdf_lv),
        len(bridges),
        int(self_close.sum()),
        n_long,
        lg["max_gap_m"],
        n_long_cand,
        time.perf_counter() - t,
    )

    # --- Rivers snapped into their channel (always): once per country, shared by all units ---
    t = time.perf_counter()
    channel_path = out_dir / "rivers_channel.gpkg"
    if channel_path.exists() and not run["FORCE"]:
        gdf_rn = gpd.read_file(channel_path).to_crs(wkt)
        log.info(
            "Rivers: %d lines from rivers_channel.gpkg (set FORCE to snap again).",
            len(gdf_rn),
        )
    else:
        gdf_rn = gis.Vector.read_merged(
            country["rivers"], wkt, bbox=all_windows, keep_columns=True
        )
        gdf_rn = gdf_rn[
            gdf_rn.geometry.notna() & ~gdf_rn.geometry.is_empty
        ].reset_index(drop=True)
        guides = None
        if country["river_guides"]:
            guides = gis.Vector.read_merged(
                country["river_guides"], wkt, bbox=all_windows
            ).geometry.values
        geoms, snap_stats = gis.Hydro.snap_to_channel(
            gdf_rn.geometry.values, str(vrt_path), guides, cfg["river_snap"]
        )
        gdf_rn = gdf_rn.set_geometry(gpd.GeoSeries(geoms, crs=wkt))
        gis.Vector.write_gpkg(gdf_rn, channel_path)
        log.info(
            "Rivers snapped into their channel (rivers_channel.gpkg): %d lines, %d vertices (%.0f %% near a guide "
            "line%s), shift mean %.1f m, p95 %.1f m, %d parts unchanged; %.1f s.",
            snap_stats["lines"],
            snap_stats["vertices"],
            100 * snap_stats.get("guided_share", 0.0),
            "" if guides is not None else ", no river_guides given",
            snap_stats.get("shift_mean_m", 0.0),
            snap_stats.get("shift_p95_m", 0.0),
            snap_stats.get("parts_unchanged", 0),
            time.perf_counter() - t,
        )
    log.info(
        "Water sources: all %d river lines (rivers crossing a barrier: gated_crossings = %s).",
        len(gdf_rn),
        params["gated_crossings"],
    )
    # --- River groups: connected reaches of similar upstream area, for the protected areas per group ---
    field, rg = country["river_upstream_area_field"], cfg["river_groups"]
    groups_on = bool(field) and field in gdf_rn.columns
    if groups_on:
        gdf_rn["group_id"] = gis.Vector.river_groups(
            gdf_rn.geometry.values,
            pd.to_numeric(gdf_rn[field], errors="coerce"),
            rg["max_join_ratio"],
            rg["max_length_m"],
            rg["min_upstream_area_km2"],
            params["touch_tolerance_m"],
        )
        gdf_rn["length_m"] = gdf_rn.geometry.length
        gis.Vector.write_gpkg(
            gdf_rn[["group_id", field, "length_m", "geometry"]],
            out_dir / "river_groups.gpkg",
        )
        g_len = (
            gdf_rn[gdf_rn["group_id"] > 0].groupby("group_id")["length_m"].sum() / 1000
        )
        log.info(
            "River groups (river_groups.gpkg): %d groups of reaches from %.0f km2 (join ratio %.1f, at most %.0f km), "
            "length per group %s km; %d minor reaches (group 0) are handled together.",
            len(g_len),
            rg["min_upstream_area_km2"],
            rg["max_join_ratio"],
            rg["max_length_m"] / 1000,
            (
                ", ".join(
                    f"p{q} {v:.1f}"
                    for q, v in zip((5, 50, 95), np.percentile(g_len, (5, 50, 95)))
                )
                if len(g_len)
                else "-"
            ),
            int((gdf_rn["group_id"] == 0).sum()),
        )
    else:
        log.warning(
            "No river_upstream_area_field in the rivers: protected areas per river group are not computed."
        )
    gis.Vector.write_gpkg(
        gdf_rn[(["group_id"] if groups_on else []) + ["geometry"]],
        work_dir / "rivers_used.gpkg",
    )

    # --- Units ---
    list_unk = sorted(set(run["UNITS"]) - set(df_units.unit_id))
    if list_unk:
        raise SystemExit(f"Unknown unit ids in UNITS: {list_unk}")
    df_sel = df_units[df_units.unit_id.isin(run["UNITS"])] if run["UNITS"] else df_units
    list_task, list_w, list_skip = [], [], []
    for row in df_sel.to_dict("records"):
        unit_dir = work_dir / "units" / row["unit_id"]
        if (unit_dir / "done.json").exists() and not run["FORCE"]:
            list_skip.append(row["unit_id"])
            continue
        row.update(
            {
                "unit_dir": str(unit_dir),
                "dtm_vrt": str(vrt_path),
                "wkt": wkt,
                "resolution": resolution,
                "domain_path": str(work_dir / "units_domain.gpkg"),
                "barrier_lines": str(work_dir / "barrier_lines.gpkg"),
                "rivers_used": str(work_dir / "rivers_used.gpkg"),
                "warp_memory_mb": general["warp_memory_limit_mb"],
                "params": params,
            }
        )
        list_task.append(row)
        list_w.append(max(row["est_domain_cells"], 1))
    log.info(
        "Units: %d to process, %d already finished and skipped (set FORCE to recompute).",
        len(list_task),
        len(list_skip),
    )
    if list_task:
        # Compile the numba kernels once so that worker processes load them from the cache.
        gis.Hydro.breach_to_streams(np.zeros((3, 3), np.float32), np.eye(3, dtype=bool))
        gis.Hydro.priority_flood(np.zeros((3, 3), np.float32), np.eye(3, dtype=bool))
    list_res, list_fail = gis.Run.run_units(
        process_unit,
        list_task,
        general["n_workers"],
        general["gdal_cachemax_mb"],
        general["gdal_threads_per_worker"],
        list_w,
        log,
    )
    if list_fail:
        log.error(
            "%d units failed: %s. Rerun to retry them; finished units are skipped.",
            len(list_fail),
            list_fail,
        )

    # --- Merge all finished units ---
    t = time.perf_counter()
    gdf_core = gpd.read_file(out_dir / "units_core.gpkg").set_index("unit_id")
    gdf_win = gpd.read_file(out_dir / "units_window.gpkg").set_index("unit_id")
    list_done = [
        u for u in df_units.unit_id if (work_dir / "units" / u / "done.json").exists()
    ]
    list_pa = [
        gpd.read_file(work_dir / "units" / u / "protected_areas.gpkg")
        for u in list_done
    ]
    list_pa = [f for f in list_pa if len(f)]
    if list_pa:
        gdf_cand = gpd.GeoDataFrame(
            pd.concat(list_pa, ignore_index=True), geometry="geometry", crs=wkt
        )
    else:
        gdf_cand = gpd.GeoDataFrame(
            {c: pd.Series(dtype=d) for c, d in UNIT_COLUMNS.items()}
            | {"unit_id": pd.Series(dtype=str)},
            geometry=gpd.GeoSeries([], crs=wkt),
        )
    gdf_final, n_contra = merge_units(
        gdf_cand, gdf_core, gdf_win, list_done, country, run["COUNTRY"], wkt, gdf_bl
    )
    gdf_final = gdf_final[OUTPUT_COLUMNS].astype(
        {
            "truncated": bool,
            "touches_corridor_edge": bool,
            "touches_dtm_edge": bool,
            "no_levee": bool,
        }
    )
    gis.Vector.write_gpkg(gdf_final, out_dir / "protected_areas.gpkg")

    # --- Protected areas per river group: merged group by group (areas of different groups may overlap) ---
    n_by_group = 0
    if groups_on:
        list_g = [
            gpd.read_file(work_dir / "units" / u / "protected_areas_by_group.gpkg")
            for u in list_done
            if (work_dir / "units" / u / "protected_areas_by_group.gpkg").exists()
        ]
        list_g = [f for f in list_g if len(f)]
        list_merged = []
        if list_g:
            gdf_cg = pd.concat(list_g, ignore_index=True)
            for group_id, gdf_one in gdf_cg.groupby("group_id"):
                gdf_m, _ = merge_units(
                    gpd.GeoDataFrame(
                        gdf_one.reset_index(drop=True), geometry="geometry", crs=wkt
                    ),
                    gdf_core,
                    gdf_win,
                    list_done,
                    country,
                    run["COUNTRY"],
                    wkt,
                    gdf_bl,
                )
                gdf_m["protected_area_id"] = (
                    gdf_m["protected_area_id"] + f"_g{int(group_id)}"
                )
                gdf_m["group_id"] = int(group_id)
                list_merged.append(gdf_m)
        columns = OUTPUT_COLUMNS[:-1] + ["group_id", "geometry"]
        gdf_byg = (
            gpd.GeoDataFrame(
                pd.concat(list_merged, ignore_index=True), geometry="geometry", crs=wkt
            )[columns]
            if list_merged
            else gpd.GeoDataFrame(columns=columns, geometry="geometry", crs=wkt)
        )
        gis.Vector.write_gpkg(gdf_byg, out_dir / "protected_areas_by_group.gpkg")
        n_by_group = len(gdf_byg)
        log.info(
            "Protected areas per river group (protected_areas_by_group.gpkg): %d areas in %d groups, %.1f km2 in "
            "total (overlaps between groups counted per group), against %.1f km2 with all water sources together.",
            n_by_group,
            gdf_byg["group_id"].nunique() if n_by_group else 0,
            gdf_byg.geometry.area.sum() / 1e6 if n_by_group else 0.0,
            gdf_final.geometry.area.sum() / 1e6,
        )
    n_truncated = int(gdf_final["truncated"].sum())
    log.info(
        "Protected areas: %d from %d unit candidates in %d finished units; truncated %d (%.1f %%; many truncated "
        "protected areas mean the halo is too small), dropped as contradicted by neighbouring windows %d, touching "
        "the corridor edge %d, no levee %d. Merge took %.1f s.",
        len(gdf_final),
        len(gdf_cand),
        len(list_done),
        n_truncated,
        100 * n_truncated / max(len(gdf_final), 1),
        n_contra,
        int(gdf_final["touches_corridor_edge"].sum()),
        int(gdf_final["no_levee"].sum()),
        time.perf_counter() - t,
    )

    # --- Barrier lines with crest attributes (minimum over the units' cores) ---
    list_crest = [
        pd.read_csv(work_dir / "units" / u / "barrier_crest.csv") for u in list_done
    ]
    list_crest = [table for table in list_crest if len(table)]
    df_crest = (
        pd.concat(list_crest)
        .groupby("line_id")[["crest_min_m", "crest_hand_min_m"]]
        .min()
        if list_crest
        else pd.DataFrame(columns=["crest_min_m", "crest_hand_min_m"])
    )
    gdf_bar = gdf_bl.merge(df_crest, left_on="line_id", right_index=True, how="left")
    gis.Vector.write_gpkg(gdf_bar, out_dir / "barriers.gpkg")

    # --- Leaks: the same place found from several units counted once; largest first ---
    list_lk = [
        gpd.read_file(work_dir / "units" / u / "leaks.gpkg")
        for u in list_done
        if (work_dir / "units" / u / "leaks.gpkg").exists()
    ]
    list_lk = [f for f in list_lk if len(f)]
    if list_lk:
        df_lk = pd.concat(list_lk, ignore_index=True)
        step = country["id_rounding_m"]
        df_lk["kx"], df_lk["ky"] = np.round(df_lk["x"] / step), np.round(
            df_lk["y"] / step
        )
        grouped = df_lk.sort_values("levee_length_m", ascending=False).groupby(
            ["cause", "kx", "ky"], sort=False
        )
        df_lk = grouped.first()
        df_lk["levee_length_m"] = grouped["levee_length_m"].sum()
        df_lk["n_samples"] = grouped["n_samples"].sum()
        df_lk["line_ids"] = grouped["line_ids"].agg(
            lambda v: ",".join(
                sorted({i for s in v for i in str(s).split(",") if i}, key=int)
            )
        )
        df_lk = df_lk.reset_index().sort_values(
            "levee_length_m", ascending=False, ignore_index=True
        )
    else:
        df_lk = pd.DataFrame({c: pd.Series(dtype=d) for c, d in LEAK_COLUMNS.items()})
        df_lk["unit_id"] = pd.Series(dtype=str)
    df_lk.insert(
        0, "leak_id", [f"{run['COUNTRY']}_L{k + 1}" for k in range(len(df_lk))]
    )
    gdf_leaks = gpd.GeoDataFrame(
        df_lk[["leak_id"] + list(LEAK_COLUMNS) + ["unit_id"]],
        geometry=gpd.points_from_xy(df_lk["x"], df_lk["y"]),
        crs=wkt,
    )
    gis.Vector.write_gpkg(gdf_leaks, out_dir / "leaks.gpkg")
    by_cause = gdf_leaks.groupby("cause")["levee_length_m"].agg(["size", "sum"])
    log.info(
        "Levees that protect nothing (leaks.gpkg), by cause: %s. Largest: %s.",
        ", ".join(
            f"{c} {int(r['size'])} places, {r['sum'] / 1000:.1f} km of levee"
            for c, r in by_cause.sort_values("sum", ascending=False).iterrows()
        )
        or "none",
        "; ".join(
            f"{r.leak_id} {r.cause} at ({r.x:.0f}, {r.y:.0f}), {r.levee_length_m:.0f} m"
            for r in gdf_leaks.head(5).itertuples()
        )
        or "-",
    )

    list_encl = [
        gpd.read_file(work_dir / "units" / u / "enclosed_streams.gpkg")
        for u in list_done
    ]
    list_encl = [f for f in list_encl if len(f)]
    if list_encl:
        gdf_encl = gpd.GeoDataFrame(
            pd.concat(list_encl, ignore_index=True), geometry="geometry", crs=wkt
        )
        gdf_encl = gdf_encl[
            shapely.covers(
                gdf_core.geometry.reindex(gdf_encl.unit_id).values,
                shapely.point_on_surface(gdf_encl.geometry.values),
            )
        ]
        gis.Vector.write_gpkg(
            gdf_encl.reset_index(drop=True), out_dir / "enclosed_streams.gpkg"
        )
        log.warning(
            "%d levee-enclosed areas contain river cells (see enclosed_streams.gpkg).",
            len(gdf_encl),
        )

    # --- Measured performance and extrapolation to the whole country ---
    df_perf = pd.DataFrame(
        [
            json.loads((work_dir / "units" / u / "done.json").read_text())
            for u in list_done
        ]
    )
    if len(df_perf):
        s_cell = df_perf["seconds"].sum() / max(df_perf["domain_cells"].sum(), 1)
        b_cell = df_perf["bytes_per_cell"].max()
        df_rem = df_units[~df_units.unit_id.isin(list_done)]
        step_share = (
            pd.DataFrame(list(df_perf["times"])).sum() / df_perf["seconds"].sum()
        )
        total_s = df_units["est_domain_cells"].sum() * s_cell / general["n_workers"]
        remaining_s = df_rem["est_domain_cells"].sum() * s_cell / general["n_workers"]
        log.info(
            "Measured on %d units: %.2f microseconds per domain cell, up to %.0f bytes per window cell (configured "
            "%.0f), peak RSS up to %.2f GB. Time shares: %s.",
            len(df_perf),
            1e6 * s_cell,
            b_cell,
            country["bytes_per_cell"],
            df_perf["peak_rss_mb"].max() / 1024,
            ", ".join(
                f"{k} {100 * v:.0f} %"
                for k, v in step_share.sort_values(ascending=False).items()
            ),
        )
        log.info(
            "Extrapolation (measured s/cell x estimated domain cells / %d workers): whole country %s, remaining "
            "%d units %s; largest unit window needs about %.1f GB.",
            general["n_workers"],
            gis.Run.format_seconds(total_s),
            len(df_rem),
            gis.Run.format_seconds(remaining_s),
            df_units["est_window_cells"].max() * b_cell / 1e9,
        )
        (out_dir / "performance.json").write_text(
            json.dumps(
                {
                    "units_measured": len(df_perf),
                    "seconds_per_domain_cell": s_cell,
                    "bytes_per_window_cell": b_cell,
                    "country_seconds": total_s,
                    "remaining_seconds": remaining_s,
                    "time_shares": step_share.to_dict(),
                },
                indent=2,
            )
        )
