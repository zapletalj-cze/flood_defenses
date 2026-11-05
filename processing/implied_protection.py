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
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
from osgeo import gdal
from scipy import ndimage

import gis
from plots import CATEGORICAL, COLOURS, INK, Plots

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
            "component_built_min_m2",
            "merge_max_dlevel_m",
            "merge_max_gap",
            "corroboration_distance_m",
            "corroboration_min",
            "explicit_overlap_max",
            "regional_sop_rp",
            "implied_corroborated_max_rp",
            "historic_low_core_uncorroborated_cap",
            "epoch_cutoff",
            "growth_low_share_min",
            "min_cells_temporal",
            "rp_distribution_classes",
            "rp_distribution_shares",
        )
    ],
    "country.rivers",
    "country.river_upstream_area_field",
    "country.rp.classes",
    "country.rp.storey_height_m",
    "country.rp.max_damage_csv",
    *[
        f"country.context.params.{c}.{k}"
        for c in ("industrial", "dense_urban", "urban", "rural")
        for k in (
            "score_weight_value",
            "score_weight_population",
            "score_weight_upstream",
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
# Stage C keys, loaded only when GFM scenes exist (gfm_index.csv)
GFM_KEYS = [
    "country.compartments.max_level_m",
    "evidence.comparison_distance_m",
    "evidence.min_valid_pixels",
    "evidence.f_wet_comparison",
    "evidence.f_dry_inside",
    "evidence.f_wet_inside",
    "sources.gfm_flood_values",
    "sources.gfm_exclusion_values",
    "sources.gfm_reference_water_values",
    "sources.gfm_nodata_values",
]


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


def implied_unit(task):
    # Worker, pass 2: implied cells, 4-connected components per group, explicit overlap, edge cells with crest and
    # closure support, value and population (no depth weighting), GFM tests; polygons without explicit cells.
    start = time.perf_counter()
    unit_dir = Path(task["unit_dir"])
    arrays = {}
    for name in ("subsegments.tif", task["hand_name"], "implied_built.tif"):
        ds = gdal.Open(str(unit_dir / name))
        arrays[name] = ds.GetRasterBand(1).ReadAsArray()
        gt, wkt = ds.GetGeoTransform(), ds.GetProjection()
        ds = None
    labels = arrays.pop("subsegments.tif")
    hand = arrays.pop(task["hand_name"]).astype(np.float32, copy=False)
    built = arrays.pop("implied_built.tif").astype(np.float32, copy=False)
    shape, res = labels.shape, gt[1]
    cell = res * res
    bounds = (gt[0], gt[3] - shape[0] * res, gt[0] + shape[1] * res, gt[3])
    lev = np.asarray(task["level_lab"], np.float32)[labels]
    grp = np.asarray(task["group_lab"], np.int32)[labels]
    implied = (
        (grp >= 0)
        & np.isfinite(hand)
        & (hand >= lev)
        & (hand <= task["envelope_hand_m"])
    )
    explicit = (
        gis.Raster.rasterize_geometries(
            shapely.from_wkb(task["explicit"]),
            np.ones(len(task["explicit"])),
            shape,
            gt,
            wkt,
            False,
            np.uint8,
        )
        > 0
        if task["explicit"]
        else np.zeros(shape, bool)
    )
    # Components per group (a component never mixes groups).
    comp = np.zeros(shape, np.int32)
    list_group = [-1]
    for g in np.unique(grp[implied]):
        m = implied & (grp == g)
        rr, cc = np.flatnonzero(m.any(axis=1)), np.flatnonzero(m.any(axis=0))
        sub = (slice(rr[0], rr[-1] + 1), slice(cc[0], cc[-1] + 1))
        lab_g, n = ndimage.label(m[sub])
        view = comp[sub]
        view[lab_g > 0] = lab_g[lab_g > 0] + len(list_group) - 1
        list_group += [int(g)] * n
    n_comp = len(list_group)
    del grp, implied, labels
    comp_ne = np.where(explicit, 0, comp)
    # Counts over implied cells only (bincount weights are float64: whole-window weights would cost 8 bytes a cell).
    cells = np.flatnonzero(comp.ravel())
    flat = comp.ravel()[cells]
    flat_ne = comp_ne.ravel()[cells]
    built_c = built.ravel()[cells]
    stats = {
        "n_cells": np.bincount(flat, minlength=n_comp).astype(float),
        "n_explicit": np.bincount(
            flat[explicit.ravel()[cells]], minlength=n_comp
        ).astype(float),
        "built_all_m2": np.bincount(flat, weights=built_c, minlength=n_comp) * cell,
        "built_m2": np.bincount(flat_ne, weights=built_c, minlength=n_comp) * cell,
        "value": np.zeros(n_comp),
        "population": np.zeros(n_comp),
    }
    del built, built_c, flat
    volume = {}
    for product, path in task["ghsl"].items():
        source = gdal.Open(path)
        sgt = source.GetGeoTransform()
        warped, _ = gis.Raster.read_window(
            source, bounds, res, "average", task["warp_mb"], dst_wkt=wkt
        )
        volume[product] = np.nan_to_num(warped, nan=0.0) * cell / abs(sgt[1] * sgt[5])
        source = None
    if "total" in volume and "non_residential" in volume:
        residential = np.clip(volume["total"] - volume["non_residential"], 0.0, None)
        cell_value = (
            residential * task["value_m2"]["residential"]
            + volume["non_residential"] * task["value_m2"]["non_residential"]
        ) / task["storey_height_m"]
        stats["value"] = np.bincount(
            flat_ne, weights=cell_value.ravel()[cells], minlength=n_comp
        )
    if "pop" in volume:
        stats["population"] = np.bincount(
            flat_ne, weights=volume["pop"].ravel()[cells], minlength=n_comp
        )
    # Edge cells: implied cells with a 4-neighbour outside the implied area whose HAND is below their level.
    inner = comp_ne > 0
    edge = np.zeros(shape, bool)
    for a, b in ((slice(1, None), slice(None, -1)), (slice(None, -1), slice(1, None))):
        for axis in (0, 1):
            here = (a, slice(None)) if axis == 0 else (slice(None), a)
            there = (b, slice(None)) if axis == 0 else (slice(None), b)
            out_nb = ~inner[there] & (hand[there] < lev[here])
            edge[here] |= inner[here] & out_nb
    edge_ids = comp_ne[edge]
    stats["edge_n"] = np.bincount(edge_ids, minlength=n_comp).astype(float)
    stats["edge_supported"] = np.zeros(n_comp)
    reach_m = task["corroboration_distance_m"]
    near = np.zeros(shape, bool)
    closures = [[] for _ in range(n_comp)]
    if edge.any() and (task["lines"] or task["closures"]):
        # Distances are needed only near edge cells: window around them, padded by the corroboration distance.
        pad = int(np.ceil(reach_m / res)) + 1
        er_, ec_ = np.nonzero(edge)
        win = (
            slice(max(er_.min() - pad, 0), min(er_.max() + 1 + pad, shape[0])),
            slice(max(ec_.min() - pad, 0), min(ec_.max() + 1 + pad, shape[1])),
        )
        wgt = (
            gt[0] + win[1].start * res,
            res,
            0.0,
            gt[3] - win[0].start * res,
            0.0,
            -res,
        )
        wshape = (win[0].stop - win[0].start, win[1].stop - win[1].start)
        if task["lines"]:
            lines = (
                gis.Raster.rasterize_geometries(
                    shapely.from_wkb(task["lines"]),
                    np.ones(len(task["lines"])),
                    wshape,
                    wgt,
                    wkt,
                    True,
                    np.uint8,
                )
                > 0
            )
            if lines.any():
                near[win] = gis.Hydro.distance_cells(lines, pad + 1) * res <= reach_m
                stats["edge_supported"] = np.bincount(
                    comp_ne[edge & near], minlength=n_comp
                ).astype(float)
        if task["closures"]:
            ids = np.asarray(task["closure_ids"])
            burnt = gis.Raster.rasterize_geometries(
                shapely.from_wkb(task["closures"]),
                np.arange(1, len(ids) + 1),
                wshape,
                wgt,
                wkt,
                True,
                np.int32,
            )
            if (burnt > 0).any():
                dist, (ir, ic) = ndimage.distance_transform_edt(
                    burnt == 0, return_indices=True
                )
                hit = edge[win] & (dist * res <= reach_m)
                pairs = set(
                    zip(comp_ne[win][hit].tolist(), burnt[ir[hit], ic[hit]].tolist())
                )
                for c, j in pairs:
                    closures[c].append(str(ids[j - 1]))
                del dist, ir, ic
    # GFM tests: comparison zone outside both implied and explicit compartments.
    list_gfm = []
    if task["scenes"] and n_comp > 1:
        p = task["gfm"]
        zones = gis.Evidence.comparison_zones(
            comp_ne,
            hand,
            {c: c for c in range(1, n_comp)},
            int(np.ceil(p["comparison_distance_m"] / res)),
            p["max_level_m"],
            outside=~explicit,
        )
        for scene in task["scenes"]:
            valid, wet, _ = gis.Evidence.scene_masks(
                scene, bounds, res, wkt, hand, p, task["warp_mb"]
            )
            for t in gis.Evidence.test_regions(zones, valid, wet, p):
                list_gfm.append(
                    {
                        "comp": t["region"],
                        "item_id": scene["item_id"],
                        "datetime": scene["datetime"],
                        **{
                            k: t[k]
                            for k in (
                                "inside_valid",
                                "inside_wet",
                                "comp_valid",
                                "comp_wet",
                            )
                        },
                    }
                )
    # Polygons of the components without explicit cells; edge cells for the maps.
    gdf_poly = (
        gis.Raster.polygonize_labels(comp_ne, gt, wkt)
        .dissolve(by="label")
        .reset_index()
    )
    geom = {
        int(l): shapely.to_wkb(g) for l, g in zip(gdf_poly["label"], gdf_poly.geometry)
    }
    er, ec = np.nonzero(edge)
    sup = near[er, ec]
    comps = [
        {
            "unit_id": task["unit_id"],
            "comp": c,
            "group": list_group[c],
            **{k: float(v[c]) for k, v in stats.items()},
            "closures": closures[c],
            "geometry": geom.get(c),
        }
        for c in range(1, n_comp)
    ]
    return {
        "unit_id": task["unit_id"],
        "summary": (
            f"unit {task['unit_id']}: {n_comp - 1} implied components, {int(inner.sum())} cells outside explicit "
            f"compartments, {len(list_gfm)} GFM tests, {time.perf_counter() - start:.1f} s"
        ),
        "comps": comps,
        "gfm": list_gfm,
        "edges": {
            "comp": comp_ne[er, ec],
            "x": gt[0] + (ec + 0.5) * res,
            "y": gt[3] - (er + 0.5) * res,
            "supported": sup,
        },
        "secs": {"pass2": time.perf_counter() - start},
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

    # --- C. Merging along each river side; gaps of at most merge_max_gap bridged ---
    t = time.perf_counter()
    group = np.full(n_seg, -1, np.int64)
    level_used = np.full(n_seg, np.nan)
    bridged = np.zeros(n_seg, bool)
    n_group = 0
    for _, idx in gdf_seg.groupby(["river_id", "side"]).indices.items():
        g, used, br = gis.Stats.merge_levels(
            gdf_seg["k"].to_numpy()[idx],
            eligible[idx],
            primary[idx],
            imp["merge_max_dlevel_m"],
            imp["merge_max_gap"],
        )
        group[idx] = np.where(g >= 0, g + n_group, -1)
        level_used[idx], bridged[idx] = used, br
        n_group += int(g.max()) + 1 if (g >= 0).any() else 0
    gdf_seg["group"], gdf_seg["level_used_m"], gdf_seg["bridged"] = (
        group,
        level_used,
        bridged,
    )
    gis.Vector.write_gpkg(gdf_seg, run_dir / "implied_subsegments.gpkg")
    log.info(
        "Merging: %d groups from %d eligible subsegments, %d bridged.",
        n_group,
        int(eligible.sum()),
        int(bridged.sum()),
    )

    # --- Inputs of pass 2; optional ones switch features off ---
    explicit_path = out_dir / "compartments_rp.gpkg"
    if not explicit_path.exists():
        explicit_path = out_dir / "compartments.gpkg"
        log.warning(
            "compartments_rp.gpkg missing: precedence uses compartments.gpkg of stage A."
        )
    gdf_expl = (
        gpd.read_file(explicit_path).to_crs(wkt) if explicit_path.exists() else None
    )
    if gdf_expl is None:
        log.warning("No explicit compartments: no precedence.")
    list_line = []
    for name in ("barriers_all.gpkg", "crest_lines.gpkg"):
        if (out_dir / name).exists():
            list_line.append(gpd.read_file(out_dir / name).to_crs(wkt))
    if not list_line and (out_dir / "barriers.gpkg").exists():
        log.warning(
            "barriers_all.gpkg and crest_lines.gpkg (Stage D, D2) missing: crest support uses the stage A "
            "barrier lines (barriers.gpkg); closure support is off."
        )
        list_line.append(gpd.read_file(out_dir / "barriers.gpkg").to_crs(wkt))
    lines = (
        np.concatenate([g.geometry.values for g in list_line])
        if list_line
        else np.empty(0, object)
    )
    closures, closure_ids = np.empty(0, object), np.empty(0, object)
    if (out_dir / "barriers_all.gpkg").exists():
        gdf_bar = list_line[0]
        if "closure_gap" in gdf_bar.columns:
            sel = gdf_bar["closure_gap"].astype(bool).to_numpy()
            closures = gdf_bar.geometry.values[sel]
            id_col = "barrier_id" if "barrier_id" in gdf_bar.columns else None
            closure_ids = (
                gdf_bar[id_col].to_numpy() if id_col else np.arange(len(gdf_bar))
            )[sel]
    ghsl = {
        product: str(out_dir / f"ghsl_{name}.vrt")
        for product, name in (
            ("total", "total"),
            ("non_residential", "non_residential"),
            ("pop", "pop"),
        )
        if (out_dir / f"ghsl_{name}.vrt").exists()
    }
    df_maxd = pd.read_csv(country["rp"]["max_damage_csv"])
    df_row = df_maxd[
        df_maxd["country"].astype(str).str.upper() == run["COUNTRY"].upper()
    ]
    value_m2 = (
        {k: float(df_row[k].iloc[0]) for k in ("residential", "non_residential")}
        if len(df_row) == 1
        else None
    )
    if value_m2 is None or not {"total", "non_residential"} <= set(ghsl):
        log.warning(
            "GHSL building volume or the max damage row missing: V is not computed (score without it)."
        )
        ghsl = {k: v for k, v in ghsl.items() if k == "pop"}
    if "pop" not in ghsl:
        log.warning("ghsl_pop.vrt missing: P is not computed (score without it).")
    scenes_by_unit, gfm = {}, None
    if (out_dir / "gfm_index.csv").exists():
        cfg_gfm = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], GFM_KEYS)
        ev, src = cfg_gfm["evidence"], cfg_gfm["sources"]
        gfm = {
            **{
                k: ev[k]
                for k in (
                    "comparison_distance_m",
                    "min_valid_pixels",
                    "f_wet_comparison",
                    "f_dry_inside",
                    "f_wet_inside",
                )
            },
            "max_level_m": cfg_gfm["country"]["compartments"]["max_level_m"],
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
        log.warning("gfm_index.csv missing: GFM corroboration is off.")

    # --- Pass 2 per unit ---
    expl_tree = (
        shapely.STRtree(gdf_expl.geometry.values) if gdf_expl is not None else None
    )
    line_tree, closure_tree = shapely.STRtree(lines), shapely.STRtree(closures)
    list_task2 = []
    for u in df_units.itertuples():
        unit_dir = work_dir / "units" / u.unit_id
        df_lut = pd.read_csv(unit_dir / "subsegments.csv")
        glob_idx = np.array(
            [
                pos[f"{r}_{s}_{k}"]
                for r, s, k in zip(df_lut["river_id"], df_lut["side"], df_lut["k"])
            ],
            dtype=np.int64,
        )
        window = shapely.box(u.win_xmin, u.win_ymin, u.win_xmax, u.win_ymax)
        sel_c = closure_tree.query(window, predicate="intersects")
        list_task2.append(
            {
                "unit_id": u.unit_id,
                "unit_dir": str(unit_dir),
                "hand_name": hand_name,
                "level_lab": [np.nan] + list(level_used[glob_idx]),
                "group_lab": [-1] + list(group[glob_idx]),
                "envelope_hand_m": imp["envelope_hand_m"],
                "explicit": (
                    list(
                        shapely.to_wkb(
                            gdf_expl.geometry.values[
                                expl_tree.query(window, predicate="intersects")
                            ]
                        )
                    )
                    if expl_tree is not None
                    else []
                ),
                "lines": list(
                    shapely.to_wkb(
                        lines[line_tree.query(window, predicate="intersects")]
                    )
                ),
                "closures": list(shapely.to_wkb(closures[sel_c])),
                "closure_ids": list(closure_ids[sel_c]),
                "corroboration_distance_m": imp["corroboration_distance_m"],
                "ghsl": ghsl,
                "value_m2": value_m2,
                "storey_height_m": country["rp"]["storey_height_m"],
                "scenes": scenes_by_unit.get(u.unit_id, []) if gfm else [],
                "gfm": gfm,
                "warp_mb": general["warp_memory_limit_mb"],
            }
        )
    t2 = time.perf_counter()
    list_res2, list_fail = gis.Run.run_units(
        implied_unit,
        list_task2,
        general["n_workers"],
        general["gdal_cachemax_mb"],
        general["gdal_threads_per_worker"],
        list_w,
        log,
    )
    pass2_s = time.perf_counter() - t2
    if list_fail:
        raise SystemExit(f"Failed units in pass 2: {list_fail}")

    # --- Implied compartments: components joined across units, built-area filter, precedence ---
    df_comp = pd.DataFrame([c for r in list_res2 for c in r["comps"]])
    df_gfm = pd.DataFrame(
        [{**g, "unit_id": r["unit_id"]} for r in list_res2 for g in r["gfm"]]
    )
    list_row, list_drop = [], []
    seg_of = pd.Series(np.arange(n_seg)).groupby(group).apply(list).to_dict()
    for g in range(n_group):
        cg = df_comp[df_comp["group"] == g] if len(df_comp) else df_comp
        if cg.empty:
            list_drop.append({"group": g, "reason": "no implied cells"})
            continue
        with_geom = cg[cg["geometry"].notna()]
        geoms = shapely.from_wkb(with_geom["geometry"].to_numpy())
        parts = (
            shapely.get_parts(shapely.union_all(geoms))
            if len(geoms)
            else np.empty(0, object)
        )
        owner = (
            shapely.STRtree(parts).query_nearest(
                shapely.point_on_surface(geoms), all_matches=False
            )[1]
            if len(parts)
            else np.empty(0, int)
        )
        part_built = np.bincount(
            owner, weights=with_geom["built_all_m2"].to_numpy(), minlength=len(parts)
        )
        keep_part = part_built >= imp["component_built_min_m2"]
        kept = with_geom[keep_part[owner]] if len(parts) else with_geom.iloc[:0]
        only_expl = cg[
            cg["geometry"].isna()
            & (cg["built_all_m2"] >= imp["component_built_min_m2"])
        ]
        kept_all = pd.concat([kept, only_expl])
        if kept_all.empty:
            list_drop.append(
                {
                    "group": g,
                    "reason": "no component with component_built_min_m2 of built area",
                }
            )
            continue
        overlap = kept_all["n_explicit"].sum() / max(kept_all["n_cells"].sum(), 1.0)
        if overlap > imp["explicit_overlap_max"]:
            list_drop.append(
                {
                    "group": g,
                    "reason": f"explicit_overlap_share {overlap:.2f} > explicit_overlap_max",
                }
            )
            continue
        geometry = shapely.union_all(parts[keep_part]) if keep_part.any() else None
        if geometry is None or geometry.is_empty:
            list_drop.append({"group": g, "reason": "covered by explicit compartments"})
            continue
        members = np.asarray(seg_of[g])
        list_row.append(
            {
                "group": g,
                "geometry": geometry,
                "explicit_overlap_share": float(overlap),
                **{
                    k: float(kept_all[k].sum())
                    for k in (
                        "built_m2",
                        "value",
                        "population",
                        "edge_n",
                        "edge_supported",
                    )
                },
                "closure_support": len({c for cl in kept_all["closures"] for c in cl}),
                "units": ",".join(sorted(set(kept_all["unit_id"]))),
                "members": members,
                "keys": set(zip(kept_all["unit_id"], kept_all["comp"])),
            }
        )
    log.info(
        "Implied compartments: %d kept, %d groups dropped (%s).",
        len(list_row),
        len(list_drop),
        (
            "; ".join(
                f"{k} {v}"
                for k, v in pd.Series([d["reason"].split(" ")[0] for d in list_drop])
                .value_counts()
                .items()
            )
            if list_drop
            else "none"
        ),
    )
    # --- Group attributes from their subsegments ---
    gdf_imp = (
        gpd.GeoDataFrame(list_row, geometry="geometry", crs=wkt) if list_row else None
    )
    n_imp = len(list_row)
    attrs = []
    for row in list_row:
        m = row["members"]
        el = m[eligible[m]]
        c_sum = ctx[m].sum(axis=0)
        attrs.append(
            {
                "river_id": int(gdf_seg["river_id"].iloc[m[0]]),
                "side": gdf_seg["side"].iloc[m[0]],
                "k_from": int(gdf_seg["k"].iloc[m].min()),
                "k_to": int(gdf_seg["k"].iloc[m].max()),
                "chainage_from_m": float(gdf_seg["chainage_from_m"].iloc[m].min()),
                "chainage_to_m": float(gdf_seg["chainage_to_m"].iloc[m].max()),
                "n_subsegments": len(m),
                "n_bridged": int(bridged[m].sum()),
                "level_min_m": float(np.nanmin(level_used[m])),
                "level_median_m": float(np.nanmedian(level_used[m])),
                "level_max_m": float(np.nanmax(level_used[m])),
                "sharpness": float(np.nanmedian(sharp[el])) if len(el) else np.nan,
                "context": (
                    contexts[int(np.argmax(c_sum))] if c_sum.sum() > 0 else "none"
                ),
            }
        )

    # --- D. Temporal signal from the GHSL epochs ---
    epochs = sorted(imp["ghsl_epochs"]) if n_ep else []
    if n_ep:
        sgt = gdal.Open(ghsl_epochs[0]).GetGeoTransform()
        ghsl_cells = env_total / abs(sgt[1] * sgt[5])
        cut = [e for e in epochs if e <= imp["epoch_cutoff"]]
    for a, row in zip(attrs, list_row):
        a.update(
            {
                "onset_epochs": "",
                "delta_onset_m": np.nan,
                "growth_low_share": np.nan,
                "growth_into_floodplain": False,
                "historic_low_core": False,
            }
        )
        if not n_ep:
            continue
        m = row["members"][ghsl_cells[row["members"]] >= imp["min_cells_temporal"]]
        if not len(m) or not cut:
            continue
        env_g = env[m].sum(axis=0)
        ep_g = ep[:, m].sum(axis=1)
        share_e = np.where(
            env_g >= imp["min_bin_area_m2"],
            np.divide(ep_g, env_g, out=np.zeros_like(ep_g), where=env_g > 0),
            np.nan,
        )
        ob, _ = gis.Stats.onset(
            share_e,
            imp["smooth_bins"],
            imp["onset_fraction"],
            imp["onset_persist_bins"],
        )
        onset_e = np.where(ob >= 0, edges[np.maximum(ob, 0)], np.nan)
        i_cut = epochs.index(cut[-1])
        added = np.clip(ep_g[-1] - ep_g[i_cut], 0.0, None)
        below = (
            edges[1:] <= onset_e[i_cut]
            if np.isfinite(onset_e[i_cut])
            else np.zeros(n_bins, bool)
        )
        growth = (
            float(added[below].sum() / ep_g[-1].sum()) if ep_g[-1].sum() > 0 else np.nan
        )
        into = (
            bool(growth >= imp["growth_low_share_min"])
            if np.isfinite(growth)
            else False
        )
        a.update(
            {
                "onset_epochs": ";".join(
                    f"{e}:{v:.2f}" for e, v in zip(epochs, onset_e)
                ),
                "delta_onset_m": float(onset_e[0] - onset_e[-1]),
                "growth_low_share": growth,
                "growth_into_floodplain": into,
                "historic_low_core": bool(
                    np.isfinite(onset_e[0])
                    and np.isfinite(onset_e[-1])
                    and onset_e[0] <= onset_e[-1] + imp["merge_max_dlevel_m"]
                    and not into
                ),
            }
        )

    # --- E. Corroboration: crest lines and barriers, closure barriers, GFM ---
    list_test = []
    for a, row in zip(attrs, list_row):
        a["crest_support"] = (
            row["edge_supported"] / row["edge_n"] if row["edge_n"] else 0.0
        )
        a["closure_support"] = row["closure_support"]
        a["implied_id"] = f"implied_{a['river_id']}_{a['side']}_{a['k_from']}"
        counts = {"held": 0, "flooded": 0, "no_test": 0, "inconclusive": 0}
        if len(df_gfm):
            mine = df_gfm[
                [k in row["keys"] for k in zip(df_gfm["unit_id"], df_gfm["comp"])]
            ]
            for (item, when), df_t in mine.groupby(["item_id", "datetime"]):
                sums = df_t[
                    ["inside_valid", "inside_wet", "comp_valid", "comp_wet"]
                ].sum()
                outcome = gis.Evidence.outcome(*sums.astype(int).tolist(), gfm)
                counts[outcome] += 1
                list_test.append(
                    {
                        "compartment_id": a["implied_id"],
                        "item_id": item,
                        "datetime": when,
                        "inside_valid_px": int(sums["inside_valid"]),
                        "inside_wet_frac": (
                            sums["inside_wet"] / sums["inside_valid"]
                            if sums["inside_valid"]
                            else np.nan
                        ),
                        "comp_valid_px": int(sums["comp_valid"]),
                        "comp_wet_frac": (
                            sums["comp_wet"] / sums["comp_valid"]
                            if sums["comp_valid"]
                            else np.nan
                        ),
                        "outcome": outcome,
                    }
                )
        a.update({f"gfm_{k}": v for k, v in counts.items()})
        a["corroborated"] = bool(
            a["crest_support"] >= imp["corroboration_min"]
            or a["closure_support"] > 0
            or a["gfm_held"] > 0
        )
    pd.DataFrame(
        list_test,
        columns=[
            "compartment_id",
            "item_id",
            "datetime",
            "inside_valid_px",
            "inside_wet_frac",
            "comp_valid_px",
            "comp_wet_frac",
            "outcome",
        ],
    ).to_csv(run_dir / "implied_gfm_evidence.csv", index=False)
    # Event return periods of flooded tests, from event_rp.csv of rp_bounds.py (Stage F) when present.
    event_rp = {}
    if (out_dir / "event_rp.csv").exists():
        df_ev = pd.read_csv(out_dir / "event_rp.csv")
        df_ev = df_ev[(df_ev["outcome"] == "flooded") & df_ev["rp_event"].notna()]
        event_rp = df_ev.groupby("compartment_id")["rp_event"].min().to_dict()
    else:
        log.info(
            "event_rp.csv missing: no event caps for flooded implied compartments."
        )

    # --- F. RP of the stratum "implied" ---
    classes = np.asarray(country["rp"]["classes"], float)
    sop = float(imp["regional_sop_rp"])
    above_sop = float(classes[classes > sop].min()) if (classes > sop).any() else sop
    upstream = {}
    field = country["river_upstream_area_field"]
    if field:
        src_path = out_dir / "rivers_sources.gpkg"
        gdf_src = (
            gpd.read_file(src_path)
            if src_path.exists()
            else gis.Vector.read_merged(
                country["rivers"], wkt, keep_columns=[field, "river_id"]
            )
        )
        if field in gdf_src.columns and "river_id" in gdf_src.columns:
            area_of = dict(
                zip(
                    gdf_src["river_id"].astype(str),
                    pd.to_numeric(gdf_src[field], errors="coerce"),
                )
            )
            gdf_sn = gpd.read_file(out_dir / "rivers_snapped.gpkg")
            for rid, ids in zip(
                gdf_sn["river_id"], gdf_sn["source_river_ids"].astype(str)
            ):
                vals = [
                    area_of.get(i) for i in ids.split(",") if area_of.get(i) is not None
                ]
                upstream[int(rid)] = float(np.nanmax(vals)) if vals else np.nan
    if not upstream:
        log.info(
            "No upstream area per river (river_upstream_area_field): the score has no upstream term."
        )
    if n_imp:
        front = np.array([a["chainage_to_m"] - a["chainage_from_m"] for a in attrs])
        terms = {
            "value": np.array([r["value"] for r in list_row]) / front,
            "population": np.array([r["population"] for r in list_row]) / front,
            "upstream": np.array([upstream.get(a["river_id"], np.nan) for a in attrs]),
        }
        with np.errstate(divide="ignore", invalid="ignore"):
            ranks = {
                k: gis.Stats.percentile_rank(np.where(v > 0, np.log(v), np.nan))
                for k, v in terms.items()
            }
        params = country["context"]["params"]
        num, den = np.zeros(n_imp), np.zeros(n_imp)
        for k, rank in ranks.items():
            w = np.array(
                [
                    params[a["context"] if a["context"] != "none" else "rural"][
                        f"score_weight_{k}"
                    ]
                    for a in attrs
                ]
            )
            has = np.isfinite(rank)
            num += np.where(has, w * np.nan_to_num(rank), 0.0)
            den += np.where(has, w, 0.0)
        score = np.divide(num, den, out=np.full(n_imp, np.nan), where=den > 0)
        rp_model = gis.Stats.quantile_classes(
            score, imp["rp_distribution_classes"], imp["rp_distribution_shares"]
        )
        for i, a in enumerate(attrs):
            a["front_length_m"] = float(front[i])
            a["value"] = float(terms["value"][i] * front[i])
            a["population"] = float(terms["population"][i] * front[i])
            a["upstream_area_km2"] = float(terms["upstream"][i])
            a["score"] = float(score[i])
            a["rp_model"] = (
                float(rp_model[i])
                if np.isfinite(rp_model[i])
                else float(min(imp["rp_distribution_classes"]))
            )
            cap, reasons = np.inf, []
            if not a["corroborated"]:
                cap, reasons = min(cap, above_sop), reasons + ["uncorroborated"]
            if a["historic_low_core"] and imp["historic_low_core_uncorroborated_cap"]:
                cap, reasons = min(cap, above_sop), reasons + ["historic_low_core"]
            if a["corroborated"]:
                cap, reasons = min(
                    cap, float(imp["implied_corroborated_max_rp"])
                ), reasons + ["corroborated_max"]
            rp_val = min(a["rp_model"], cap)
            high_cap = cap
            ev = event_rp.get(a["implied_id"])
            if ev is not None and np.isfinite(ev):
                high_cap = min(high_cap, float(ev))
                rp_val = min(rp_val, high_cap)
                reasons.append("flooded_event")
            above = classes[classes > rp_val]
            below = classes[classes < rp_val]
            a["rp"] = float(rp_val)
            a["rp_low"] = float(max(sop, below.max() if len(below) else sop))
            a["rp_high"] = float(min(above.min() if len(above) else rp_val, high_cap))
            a["rp_cap"] = float(cap) if np.isfinite(cap) else np.nan
            a["cap_reason"] = ",".join(
                r for r in reasons if r != "corroborated_max" or a["rp_model"] > cap
            )
            a["confidence"] = (
                "low"
                if (not a["corroborated"] or a["historic_low_core"])
                else (
                    "high"
                    if np.isfinite(a["sharpness"])
                    and a["sharpness"] >= imp["sharp_high"]
                    else "medium"
                )
            )
            a["rp_source"] = "implied"
        for k, v in attrs[0].items():
            gdf_imp[k] = [a[k] for a in attrs]
        gdf_imp["area_m2"] = gdf_imp.area
        gdf_imp = gdf_imp.drop(
            columns=["members", "keys", "edge_n", "edge_supported", "group"]
        ).rename(columns={"built_m2": "built_area_m2"})
        lead = ["implied_id", "river_id", "side", "k_from", "k_to"]
        gdf_imp = gdf_imp[lead + [c for c in gdf_imp.columns if c not in lead]]
    else:
        gdf_imp = gpd.GeoDataFrame(
            {"implied_id": pd.Series(dtype=str)}, geometry=gpd.GeoSeries([], crs=wkt)
        )
    gis.Vector.write_gpkg(gdf_imp, run_dir / "implied_compartments.gpkg")
    pd.DataFrame(list_drop, columns=["group", "reason"]).to_csv(
        run_dir / "implied_dropped.csv", index=False
    )
    parts_s = time.perf_counter() - t
    if n_imp:
        log.info(
            "implied_compartments.gpkg: %d compartments, %.2f km2; confidence %s; corroborated %d; rp %s.",
            n_imp,
            gdf_imp.area.sum() / 1e6,
            gdf_imp["confidence"].value_counts().to_dict(),
            int(gdf_imp["corroborated"].sum()),
            gdf_imp["rp"].value_counts().sort_index().to_dict(),
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
        if n_imp:
            conf = (
                pd.crosstab(gdf_imp["context"], gdf_imp["confidence"])
                .reindex(order)
                .fillna(0)
            )
            stacks = (
                ("high", INK["primary"]),
                ("medium", INK["muted"]),
                ("low", INK["grid"]),
            )
            panel_title = "Implied compartments per context and confidence"
        else:
            conf = (
                pd.crosstab(
                    pd.Series(context, name="context"),
                    pd.Series(
                        np.where(eligible, "eligible", "not eligible"), name="status"
                    ),
                )
                .reindex(order)
                .fillna(0)
            )
            stacks = (("eligible", INK["primary"]), ("not eligible", INK["grid"]))
            panel_title = "Subsegments per context"
        x = np.arange(len(order))
        bottom = np.zeros(len(order))
        for status, colour in stacks:
            vals = conf[status].to_numpy() if status in conf else np.zeros(len(order))
            axes[1, 1].bar(
                x, vals, bottom=bottom, color=colour, width=0.6, label=status
            )
            bottom += vals
        axes[1, 1].set_xticks(x, order)
        axes[1, 1].set_title(panel_title)
        axes[1, 1].set_ylabel("Count")
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
        # longitudinal: levels along each river, both banks, groups shaded, bridged subsegments marked
        mid = (gdf_seg["chainage_from_m"] + gdf_seg["chainage_to_m"]).to_numpy() / 2
        rivers = sorted(gdf_seg["river_id"].unique())
        flag_r = sorted({f"river_{r}" for r in gdf_seg["river_id"][bridged]})
        for name in plots.select(
            "longitudinal", [f"river_{r}" for r in rivers], flag_r
        ):
            rid = int(name.split("_")[1])
            fig, axes2 = plt.subplots(
                2, 1, figsize=(11, 6.5), sharex=True, layout="constrained"
            )
            for ax, side in zip(axes2, ("left", "right")):
                sel = (gdf_seg["river_id"] == rid).to_numpy() & (
                    gdf_seg["side"] == side
                ).to_numpy()
                for g in np.unique(group[sel & (group >= 0)]):
                    gm = sel & (group == g)
                    ax.axvspan(
                        gdf_seg["chainage_from_m"][gm].min(),
                        gdf_seg["chainage_to_m"][gm].max(),
                        color=INK["grid"],
                        lw=0,
                        zorder=0,
                    )
                ax.plot(
                    mid[sel & eligible],
                    primary[sel & eligible],
                    "o",
                    color=COLOURS["water"],
                    ms=5,
                    label="eligible",
                )
                ax.plot(
                    mid[sel & ~eligible],
                    primary[sel & ~eligible],
                    "o",
                    mfc="none",
                    color=INK["muted"],
                    ms=5,
                    label="not eligible",
                )
                ax.plot(
                    mid[sel & bridged],
                    level_used[sel & bridged],
                    "^",
                    color=COLOURS["snapped"],
                    ms=7,
                    label="bridged (level used)",
                )
                ax.set_ylabel(
                    f"{side.capitalize()} bank\n{imp['primary_level_method']} level (m)"
                )
                if not np.isfinite(primary[sel]).any():
                    ax.text(
                        0.5,
                        0.5,
                        "no levels on this bank",
                        transform=ax.transAxes,
                        ha="center",
                        color=INK["muted"],
                    )
                    ax.set_yticks([])
            axes2[1].set_xlabel("Chainage from the downstream end (m)")
            handles, labels_ = axes2[0].get_legend_handles_labels()
            handles.append(Patch(color=INK["grid"]))
            labels_.append("merged group")
            fig.legend(handles, labels_, loc="outside lower center", ncols=4)
            fig.suptitle(
                f"longitudinal · river {rid} · run {run['RUN_ID']}", x=0.01, ha="left"
            )
            sel_r = (gdf_seg["river_id"] == rid).to_numpy()
            plots.save(
                fig,
                "longitudinal",
                name,
                {
                    "subsegments": int(sel_r.sum()),
                    "eligible": int((sel_r & eligible).sum()),
                    "bridged": int((sel_r & bridged).sum()),
                },
            )
        if n_imp:
            ids_c = list(gdf_imp["implied_id"])
            by_c = gdf_imp.set_index("implied_id")
            if n_ep:
                flag_t = list(
                    by_c.index[
                        by_c["historic_low_core"] | by_c["growth_into_floodplain"]
                    ]
                )
                for name in plots.select("temporal", ids_c, flag_t):
                    r = by_c.loc[name]
                    pairs = [v.split(":") for v in r["onset_epochs"].split(";") if v]
                    fig, ax = plt.subplots(figsize=(7, 4.2), layout="constrained")
                    if pairs:
                        ax.plot(
                            [int(e) for e, _ in pairs],
                            [float(v) for _, v in pairs],
                            "o-",
                            color=COLOURS["water"],
                            label="onset level",
                        )
                    ax.axvline(
                        imp["epoch_cutoff"],
                        color=INK["muted"],
                        ls=":",
                        lw=0.8,
                        label="epoch_cutoff",
                    )
                    ax.axhline(
                        r["level_median_m"],
                        color=COLOURS["snapped"],
                        ls="--",
                        lw=1.0,
                        label="median level used",
                    )
                    ax.set_xlabel("GHSL epoch")
                    ax.set_ylabel("Onset level (HAND, m)")
                    ax.legend(loc="best")
                    ax.set_title(
                        f"temporal · {name} · run {run['RUN_ID']}\n"
                        f"delta onset {r['delta_onset_m']:.2f} m, growth below {r['growth_low_share']:.2f}, "
                        f"historic low core {'yes' if r['historic_low_core'] else 'no'}"
                    )
                    plots.save(
                        fig,
                        "temporal",
                        name,
                        {
                            k: r[k]
                            for k in (
                                "delta_onset_m",
                                "growth_low_share",
                                "growth_into_floodplain",
                                "historic_low_core",
                            )
                        },
                    )
            # compartment_map: hillshade, built fraction, outlines, barriers, rivers, edge cells
            built_vrt = work_dir / "implied_built.vrt"
            gdal.BuildVRT(
                str(built_vrt),
                [
                    str(work_dir / "units" / u / "implied_built.tif")
                    for u in df_units["unit_id"]
                ],
            )
            ds_dtm = gdal.Open(str(work_dir / "dtm" / "dtm.vrt"))
            gdf_sn = gpd.read_file(out_dir / "rivers_snapped.gpkg")
            edges_all = {
                k: np.concatenate([r["edges"][k] for r in list_res2])
                for k in ("comp", "x", "y", "supported")
            }
            edge_unit = np.concatenate(
                [
                    np.full(len(r["edges"]["comp"]), r["unit_id"], dtype=object)
                    for r in list_res2
                ]
            )
            keys_of = {a["implied_id"]: row["keys"] for a, row in zip(attrs, list_row)}
            built_cmap = LinearSegmentedColormap.from_list(
                "built", [(0, 0, 0, 0.0), (0.04, 0.04, 0.04, 0.35)]
            )
            flag_m = list(
                by_c.index[
                    (by_c["confidence"] == "low") | (by_c["explicit_overlap_share"] > 0)
                ]
            )
            for name in plots.select("compartment_map", ids_c, flag_m):
                r = by_c.loc[name]
                geom = r.geometry
                bounds_p, cell = Plots.window(
                    geom.buffer(2 * imp["corroboration_distance_m"]).bounds,
                    country["working_resolution_m"],
                )
                aspect = (bounds_p[3] - bounds_p[1]) / (bounds_p[2] - bounds_p[0])
                fig, ax = plt.subplots(
                    figsize=(10, min(10.0, max(4.0, 10 * aspect + 2.0))),
                    layout="constrained",
                )
                extent = Plots.hillshade_layer(
                    ax, ds_dtm, bounds_p, cell, general["warp_memory_limit_mb"]
                )
                frac, _ = gis.Raster.read_window(
                    str(built_vrt),
                    bounds_p,
                    cell,
                    "average",
                    general["warp_memory_limit_mb"],
                )
                ax.imshow(
                    np.ma.masked_less_equal(np.nan_to_num(frac), 0.01),
                    cmap=built_cmap,
                    vmin=0,
                    vmax=1,
                    extent=extent,
                    interpolation="nearest",
                )
                box = shapely.box(*bounds_p)
                for g in shapely.get_parts(
                    shapely.intersection(gdf_sn.geometry.values, box)
                ):
                    if not g.is_empty and shapely.get_type_id(g) == 1:
                        ax.plot(
                            *g.xy,
                            color=COLOURS["water"],
                            lw=1.2,
                            label="river (snapped)",
                        )
                for g in (
                    shapely.get_parts(shapely.intersection(lines, box))
                    if len(lines)
                    else []
                ):
                    if not g.is_empty and shapely.get_type_id(g) == 1:
                        ax.plot(
                            *g.xy,
                            color=INK["secondary"],
                            lw=0.9,
                            label="barrier / crest line",
                        )
                if gdf_expl is not None:
                    for g in shapely.get_parts(
                        shapely.intersection(gdf_expl.geometry.values, box)
                    ):
                        if shapely.get_type_id(g) == 3:
                            ax.plot(
                                *g.exterior.xy,
                                color=INK["primary"],
                                lw=1.0,
                                ls=(0, (4, 2)),
                                label="explicit compartment",
                            )
                for g in shapely.get_parts(geom):
                    ax.plot(
                        *g.exterior.xy,
                        color=COLOURS["snapped"],
                        lw=1.6,
                        label="implied compartment",
                    )
                sel = np.array(
                    [
                        (u, c) in keys_of[name]
                        for u, c in zip(edge_unit, edges_all["comp"])
                    ],
                    bool,
                )
                sup = sel & edges_all["supported"]
                ax.plot(
                    edges_all["x"][sup],
                    edges_all["y"][sup],
                    "o",
                    ms=2.5,
                    color=CATEGORICAL[5],
                    label="edge cell, supported",
                )
                ax.plot(
                    edges_all["x"][sel & ~sup],
                    edges_all["y"][sel & ~sup],
                    "x",
                    ms=3,
                    color=CATEGORICAL[7],
                    label="edge cell, not supported",
                )
                ax.set_xlim(extent[0], extent[1])
                ax.set_ylim(extent[2], extent[3])
                ax.set_xlabel("Easting (m)")
                ax.set_ylabel("Northing (m)")
                handles, labels_ = ax.get_legend_handles_labels()
                uniq = dict(zip(labels_, handles))
                fig.legend(
                    uniq.values(),
                    uniq.keys(),
                    loc="outside lower center",
                    ncols=3,
                    fontsize=8,
                )
                ax.set_title(
                    f"compartment_map · {name} · run {run['RUN_ID']}\n"
                    f"level {r['level_min_m']:.2f}-{r['level_max_m']:.2f} m, crest support {r['crest_support']:.2f}, "
                    f"closures {r['closure_support']}, GFM held {r['gfm_held']}, rp {r['rp']:.0f} ({r['confidence']})"
                )
                plots.save(
                    fig,
                    "compartment_map",
                    name,
                    {
                        k: r[k]
                        for k in (
                            "level_median_m",
                            "explicit_overlap_share",
                            "crest_support",
                            "closure_support",
                            "gfm_held",
                            "gfm_flooded",
                            "corroborated",
                            "rp",
                            "rp_low",
                            "rp_high",
                            "confidence",
                        )
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
        "pass2_s": round(sum(r["secs"]["pass2"] for r in list_res2), 2),
        "pass2_wall_s": round(pass2_s, 2),
        "parts_c_f_s": round(parts_s, 2),
        "implied_compartments": n_imp,
        "levels_s": round(curves_s, 2),
        "plots_s": round(plots_s, 2),
        "total_s": round(time.perf_counter() - start, 2),
        "peak_rss_mb_worker": round(
            max(r["peak_rss_mb"] for r in list_res + list_res2)
        ),
        "peak_rss_mb_main": round(gis.Run.peak_rss_mb()),
    }
    (run_dir / "implied_protection_performance.json").write_text(
        json.dumps(perf, indent=1)
    )
    log.info("Timings: %s", ", ".join(f"{k} {v}" for k, v in perf.items()))
