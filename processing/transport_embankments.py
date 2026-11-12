# Author: Jakub Zapletal
# Date: 2026-11-12
# Stage G step 3: road and railway embankments from the DTM. The lines of country.transport_lines (a file you prepare)
# only show where to look; every decision comes from the DTM. Cross-sections every station_spacing_m give the crest,
# the terrain on both sides and a class per station (embankment, one-sided, cutting, at grade); only stations in the
# floodplain (HAND of the units at most max_level_m) count. Openings: bridges and underpasses (breaks in the crest)
# and probable culverts (valleys crossing an embankment, or a mapped stream crossing it). Outputs
# transport_embankments.gpkg, openings.gpkg and plots in plots/transport/<RUN_ID> with index.html.
# Run: pixi run transport-embankments

import math
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from osgeo import gdal
from shapely.ops import substring

import gis
from plots import CATEGORICAL, INK, Plots

# --- Run variables (no command-line parameters; DEFENSES_<NAME> environment variables override them) ---
CONFIG_PATH = "config/defenses.toml"
COUNTRY = "XX"
RUN_ID = "base"  # plots in <output_dir>/<country>/plots/transport/<RUN_ID>
SHOW_PLOTS = False  # show up to implied.plots.show_max figures per type while running (needs a display)

KEYS = [
    "general.work_dir",
    "general.output_dir",
    "general.log_level",
    "country.transport_lines",
    "country.rivers",
    "country.protected_areas.max_level_m",
    "transport.station_spacing_m",
    "transport.half_width_m",
    "transport.crest_search_m",
    "transport.side_near_m",
    "transport.side_far_m",
    "transport.min_height_m",
    "transport.footprint_fraction",
    "transport.min_segment_m",
    "transport.break_window_m",
    "transport.break_depth_m",
    "transport.max_opening_m",
    "transport.culvert_window_m",
    "transport.culvert_min_depth_m",
    "transport.stream_distance_m",
    "transport.cross_sections_per_class",
    "implied.plots.enabled",
    "implied.plots.dpi",
    "implied.plots.max_per_type",
    "implied.plots.sample_seed",
    "implied.plots.show_max",
    "implied.plots.plot_all_in_units",
]
CLASSES = {
    1: "embankment",
    2: "one_sided",
    3: "cutting",
    0: "at_grade",
    -1: "no_data",
    -2: "outside_floodplain",
}
COLOURS = {
    "embankment": CATEGORICAL[0],
    "one_sided": CATEGORICAL[3],
    "cutting": CATEGORICAL[1],
    "at_grade": INK["muted"],
    "no_data": INK["grid"],
    "outside_floodplain": INK["axis"],
    "bridge": CATEGORICAL[2],
    "underpass": CATEGORICAL[6],
    "culvert": CATEGORICAL[4],
}

if __name__ == "__main__":
    run = gis.Run.env_overrides(
        CONFIG_PATH=CONFIG_PATH, COUNTRY=COUNTRY, RUN_ID=RUN_ID, SHOW_PLOTS=SHOW_PLOTS
    )
    cfg = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], KEYS)
    general, country, p = cfg["general"], cfg["country"], cfg["transport"]
    work_dir = Path(general["work_dir"]) / run["COUNTRY"]
    out_dir = Path(general["output_dir"]) / run["COUNTRY"]
    log = gis.Run.logger(
        "transport_embankments", out_dir / "logs", general["log_level"]
    )
    gis.Config.warn_placeholders(cfg)
    start = time.perf_counter()
    times = {}
    if not country["transport_lines"]:
        log.warning(
            "countries.<code>.transport_lines is empty: no road or railway lines to examine."
        )
        raise SystemExit(0)
    vrt_path = work_dir / "dtm" / "dtm.vrt"
    if not vrt_path.exists():
        raise SystemExit(f"{vrt_path} is missing; run prepare_dtm.py first.")
    wkt = gdal.Open(str(vrt_path)).GetProjection()
    spacing = p["station_spacing_m"]

    # --- Candidate lines: the file you prepared, inside the corridor, single lines ---
    t = time.perf_counter()
    gdf_ln = gis.Vector.read_merged(country["transport_lines"], wkt)
    n_input, km_input = len(gdf_ln), gdf_ln.geometry.length.sum() / 1000
    if (out_dir / "corridor.gpkg").exists():
        corridor = shapely.union_all(
            gpd.read_file(out_dir / "corridor.gpkg").to_crs(wkt).geometry.values
        )
        gdf_ln = gdf_ln.set_geometry(gdf_ln.geometry.intersection(corridor))
    gdf_ln = gdf_ln.explode(index_parts=False)
    gdf_ln = gdf_ln[
        (gdf_ln.geom_type == "LineString") & (gdf_ln.length >= spacing)
    ].reset_index(drop=True)
    gdf_ln.insert(0, "line_id", np.arange(len(gdf_ln)))
    geoms = gdf_ln.geometry.values
    log.info(
        "Transport lines: %d lines (%.1f km) in the input, %d single lines (%.1f km) inside the corridor.",
        n_input,
        km_input,
        len(gdf_ln),
        gdf_ln.geometry.length.sum() / 1000,
    )
    if not len(gdf_ln):
        raise SystemExit(0)

    # --- Stations and their cross-sections (native DTM, read block by block) ---
    owner, xy, normal, weight = gis.Vector.line_samples(geoms, spacing, spacing / 2)
    first = np.r_[True, owner[1:] != owner[:-1]]
    start_idx = np.maximum.accumulate(np.where(first, np.arange(len(owner)), 0))
    dist = (np.arange(len(owner)) - start_idx + 0.5) * weight
    n = len(owner)
    st = {
        k: np.full(n, np.nan)
        for k in [
            "crest_z",
            "line_z",
            "terrain_left",
            "terrain_right",
            "height_left",
            "height_right",
            "footprint_m",
        ]
    }
    st["cls"] = np.full(n, -1, np.int8)
    order = gis.Hydro._block_order(xy, vrt_path)
    params = {
        k: p[k]
        for k in [
            "crest_search_m",
            "side_near_m",
            "side_far_m",
            "min_height_m",
            "footprint_fraction",
        ]
    }
    for s in range(0, n, 100000):
        sel = order[s : s + 100000]
        offsets, z = gis.Hydro.crest_profiles(
            xy[sel], normal[sel], str(vrt_path), p["half_width_m"]
        )
        for k, v in gis.Hydro.embankment_stations(offsets, z, params).items():
            st[k][sel] = v
    times["profiles"] = time.perf_counter() - t

    # --- Floodplain: HAND of the units at most max_level_m ---
    t = time.perf_counter()
    list_hand = sorted(str(x) for x in (work_dir / "units").glob("*/hand_undef.tif"))
    if list_hand:
        (work_dir / "transport").mkdir(parents=True, exist_ok=True)
        hand_vrt = gis.Raster.build_vrt(list_hand, work_dir / "transport" / "hand.vrt")
        hand = gis.Raster.sample_points(hand_vrt, xy[:, 0], xy[:, 1])
        floodplain = hand <= country["protected_areas"]["max_level_m"]
    else:
        log.warning(
            "No unit HAND rasters (protected_areas.py): every station counts as floodplain."
        )
        hand = np.full(n, np.nan)
        floodplain = np.ones(n, bool)
    cls = np.where(floodplain, st["cls"], -2).astype(np.int8)
    times["floodplain"] = time.perf_counter() - t

    # --- Segments of one class along each line ---
    t = time.perf_counter()
    min_run = max(1, int(math.ceil(p["min_segment_m"] / spacing)))
    list_seg = []
    seg_of = np.full(n, -1, np.int64)
    cls_s = (
        cls.copy()
    )  # class of the segment each station belongs to (short changes absorbed)
    bounds = np.flatnonzero(first).tolist() + [n]
    for a, b in zip(bounds[:-1], bounds[1:]):
        line = geoms[owner[a]]
        for r0, r1, value in gis.Stats.runs(cls[a:b], min_run):
            i0, i1 = a + r0, a + r1
            w = weight[i0]
            d0 = 0.0 if r0 == 0 else dist[i0] - w / 2
            d1 = line.length if i1 == b else dist[i1 - 1] + w / 2
            seg_of[i0:i1] = len(list_seg)
            cls_s[i0:i1] = value
            list_seg.append(
                {
                    "segment_id": len(list_seg),
                    "line_id": int(owner[a]),
                    "source": gdf_ln.at[int(owner[a]), "source"],
                    "class": CLASSES[int(value)],
                    "length_m": d1 - d0,
                    "crest_z_min": (
                        float(np.nanmin(st["crest_z"][i0:i1]))
                        if np.isfinite(st["crest_z"][i0:i1]).any()
                        else np.nan
                    ),
                    "crest_z_median": (
                        float(np.nanmedian(st["crest_z"][i0:i1]))
                        if np.isfinite(st["crest_z"][i0:i1]).any()
                        else np.nan
                    ),
                    "height_left_m": (
                        float(np.nanmedian(st["height_left"][i0:i1]))
                        if np.isfinite(st["height_left"][i0:i1]).any()
                        else np.nan
                    ),
                    "height_right_m": (
                        float(np.nanmedian(st["height_right"][i0:i1]))
                        if np.isfinite(st["height_right"][i0:i1]).any()
                        else np.nan
                    ),
                    "footprint_m": (
                        float(np.nanmedian(st["footprint_m"][i0:i1]))
                        if np.isfinite(st["footprint_m"][i0:i1]).any()
                        else np.nan
                    ),
                    "hand_median_m": (
                        float(np.nanmedian(hand[i0:i1]))
                        if np.isfinite(hand[i0:i1]).any()
                        else np.nan
                    ),
                    "stations": i1 - i0,
                    "geometry": substring(line, max(d0, 0.0), min(d1, line.length)),
                }
            )
    gdf_seg = gpd.GeoDataFrame(list_seg, geometry="geometry", crs=wkt)
    gdf_seg["height_min_m"] = gdf_seg[["height_left_m", "height_right_m"]].min(axis=1)
    times["segments"] = time.perf_counter() - t

    # --- Openings: breaks in the crest (bridges, underpasses), probable culverts ---
    t = time.perf_counter()
    for name in ("rivers_channel.gpkg", ""):
        if not name or (out_dir / name).exists():
            break
    rivers = (
        gpd.read_file(out_dir / name).to_crs(wkt).geometry.values
        if name
        else gis.Vector.read_merged(country["rivers"], wkt).geometry.values
    )
    tree = shapely.STRtree(rivers)
    near_river = np.zeros(n, bool)
    hit_pt, _ = tree.query(
        shapely.points(xy), predicate="dwithin", distance=p["stream_distance_m"]
    )
    near_river[hit_pt] = True
    cw = max(1, int(round(p["culvert_window_m"] / spacing)))
    # dense crest along the lines with embankments, every native pixel: openings narrower than the station spacing
    # (an underpass of a few metres) lie between two cross-sections
    step_n = abs(gdal.Open(str(vrt_path)).GetGeoTransform()[1])
    # openings belong to embankment segments only (a single station on a crossing levee is not an embankment)
    lines_bank = np.unique(owner[(cls_s == 1) | (cls_s == 2)])
    d_own, d_xy, d_norm, d_w = gis.Vector.line_samples(
        geoms[lines_bank], step_n, step_n / 2
    )
    d_own = lines_bank[d_own] if len(d_own) else d_own
    d_first = np.r_[True, d_own[1:] != d_own[:-1]] if len(d_own) else d_own.astype(bool)
    d_dist = (
        np.arange(len(d_own))
        - np.maximum.accumulate(np.where(d_first, np.arange(len(d_own)), 0))
        + 0.5
    ) * d_w
    d_crest, d_line = np.full(len(d_own), np.nan), np.full(len(d_own), np.nan)
    d_order = gis.Hydro._block_order(d_xy, vrt_path) if len(d_own) else d_own
    for s in range(0, len(d_own), 200000):
        sel = d_order[s : s + 200000]
        offsets, zd = gis.Hydro.crest_profiles(
            d_xy[sel], d_norm[sel], str(vrt_path), p["crest_search_m"]
        )
        d_crest[sel] = np.where(np.isnan(zd), -np.inf, zd).max(axis=1)
        d_line[sel] = zd[:, len(offsets) // 2]
    d_crest[~np.isfinite(d_crest)] = np.nan
    d_bounds = {
        int(d_own[i]): (i, j)
        for i, j in zip(
            np.flatnonzero(d_first), np.r_[np.flatnonzero(d_first)[1:], len(d_own)]
        )
    }
    window_d = max(1, int(round(p["break_window_m"] / step_n)))
    flank_d = max(1, int(round(3 * spacing / step_n)))
    list_op = []
    for a, b in zip(bounds[:-1], bounds[1:]):
        c = cls_s[a:b]
        bank = (c == 1) | (c == 2)
        if not bank.any() or int(owner[a]) not in d_bounds:
            continue
        da, db = d_bounds[int(owner[a])]
        # station of every dense sample (stations lie at (k + 0.5) x spacing)
        idx = np.clip(
            np.round(d_dist[da:db] / weight[a] - 0.5).astype(np.int64), 0, b - a - 1
        )
        bank_d = bank[idx]
        crest = pd.Series(d_crest[da:db])
        ref = (
            crest.where(pd.Series(bank_d))
            .rolling(2 * window_d + 1, center=True, min_periods=1)
            .median()
        )
        brk = ((ref - crest).to_numpy() >= p["break_depth_m"]) & (c[idx] != -2)
        taken = np.zeros(b - a, bool)
        for r0, r1, value in gis.Stats.runs(brk, 1):
            if not value:
                continue
            width = (r1 - r0) * d_w[da]
            # an opening has embankment on both sides of it
            left_ok = bank_d[max(r0 - flank_d, 0) : r0].any()
            right_ok = bank_d[r1 : r1 + flank_d].any()
            if not (left_ok and right_ok) or width > p["max_opening_m"]:
                continue
            part = d_line[da + r0 : da + r1]
            kk = r0 + (
                int(np.nanargmin(part)) if np.isfinite(part).any() else (r1 - r0) // 2
            )
            k = a + int(idx[kk])
            hit, _ = tree.query(
                shapely.points(d_xy[da + r0 : da + r1]),
                predicate="dwithin",
                distance=p["stream_distance_m"],
            )
            river = len(hit) > 0
            list_op.append(
                {
                    "type": "bridge" if river else "underpass",
                    "x": d_xy[da + kk, 0],
                    "y": d_xy[da + kk, 1],
                    "width_m": width,
                    "invert_z": float(d_line[da + kk]),
                    "crest_z": float(ref.iloc[kk]),
                    "confidence": "high",
                    "evidence": "break in the crest" + (" at a river" if river else ""),
                    "line_id": int(owner[a]),
                    "segment_id": int(seg_of[k]),
                    "station": k,
                }
            )
            taken[max(int(idx[r0]) - cw, 0) : int(idx[r1 - 1]) + 1 + cw] = True
        # probable culverts: a mapped stream crossing the embankment, or a valley under it (terrain of both sides
        # lowest along the line, by at least culvert_min_depth_m)
        side = np.fmax(st["terrain_left"][a:b], st["terrain_right"][a:b])
        stream_runs = [r for r in gis.Stats.runs(near_river[a:b] & bank, 1) if r[2]]
        for r0, r1, _ in stream_runs:
            mid = (r0 + r1) // 2
            if taken[mid]:
                continue
            k = a + mid
            list_op.append(
                {
                    "type": "culvert",
                    "x": xy[k, 0],
                    "y": xy[k, 1],
                    "width_m": np.nan,
                    "invert_z": float(
                        np.fmin(st["terrain_left"][k], st["terrain_right"][k])
                    ),
                    "crest_z": float(st["crest_z"][k]),
                    "confidence": "high",
                    "evidence": "mapped stream crosses the embankment",
                    "line_id": int(owner[a]),
                    "segment_id": int(seg_of[k]),
                    "station": k,
                }
            )
            taken[max(mid - cw, 0) : mid + cw + 1] = True
        for i in np.flatnonzero(bank & ~taken & np.isfinite(side)):
            lo, hi = max(i - cw, 0), min(i + cw + 1, b - a)
            if side[i] > np.nanmin(side[lo:hi]) or i - lo < 1 or hi - i < 2:
                continue
            depth = min(
                np.nanmax(side[lo:i]) - side[i], np.nanmax(side[i + 1 : hi]) - side[i]
            )
            if depth < p["culvert_min_depth_m"] or taken[i]:
                continue
            k = a + i
            list_op.append(
                {
                    "type": "culvert",
                    "x": xy[k, 0],
                    "y": xy[k, 1],
                    "width_m": np.nan,
                    "invert_z": float(
                        np.fmin(st["terrain_left"][k], st["terrain_right"][k])
                    ),
                    "crest_z": float(st["crest_z"][k]),
                    "confidence": (
                        "medium" if depth >= 2 * p["culvert_min_depth_m"] else "low"
                    ),
                    "evidence": f"valley under the embankment, {depth:.2f} m deep along the line",
                    "line_id": int(owner[a]),
                    "segment_id": int(seg_of[k]),
                    "station": k,
                }
            )
            taken[lo:hi] = True
    columns = [
        "opening_id",
        "type",
        "x",
        "y",
        "width_m",
        "invert_z",
        "crest_z",
        "confidence",
        "evidence",
        "line_id",
        "segment_id",
        "station",
    ]
    df_op = pd.DataFrame(list_op, columns=columns[1:])
    df_op.insert(
        0, "opening_id", [f"{run['COUNTRY']}_O{k + 1}" for k in range(len(df_op))]
    )
    gdf_op = gpd.GeoDataFrame(
        df_op, geometry=gpd.points_from_xy(df_op["x"], df_op["y"]), crs=wkt
    )
    gdf_seg["openings"] = (
        gdf_seg["segment_id"]
        .map(gdf_op["segment_id"].value_counts())
        .fillna(0)
        .astype(int)
    )
    gis.Vector.write_gpkg(gdf_seg, out_dir / "transport_embankments.gpkg")
    gis.Vector.write_gpkg(gdf_op.drop(columns="station"), out_dir / "openings.gpkg")
    times["openings"] = time.perf_counter() - t
    by_class = gdf_seg.groupby("class")["length_m"].sum() / 1000
    log.info(
        "transport_embankments.gpkg: %d segments; %s km. Embankments: median height %.2f m (lower side), footprint "
        "%.0f m.",
        len(gdf_seg),
        ", ".join(
            f"{k} {v:.1f}" for k, v in by_class.sort_values(ascending=False).items()
        ),
        (
            gdf_seg.loc[gdf_seg["class"] == "embankment", "height_min_m"].median()
            if (gdf_seg["class"] == "embankment").any()
            else float("nan")
        ),
        (
            gdf_seg.loc[gdf_seg["class"] == "embankment", "footprint_m"].median()
            if (gdf_seg["class"] == "embankment").any()
            else float("nan")
        ),
    )
    log.info(
        "openings.gpkg: %d openings (%s).",
        len(gdf_op),
        ", ".join(
            f"{k[0]} {k[1]} {v}"
            for k, v in gdf_op.groupby(["type", "confidence"]).size().items()
        )
        or "none",
    )

    # --- Plots ---
    t = time.perf_counter()
    settings = cfg["implied"]["plots"]
    if settings["enabled"]:
        plots = Plots(
            out_dir / "plots" / "transport" / run["RUN_ID"],
            settings,
            run["SHOW_PLOTS"],
            log,
        )
        plt = plots.plt
        # cross-sections: a sample of stations of every class
        rng = np.random.default_rng(settings["sample_seed"])
        per = max(1, int(p["cross_sections_per_class"]))
        list_pick = []
        for value in (1, 2, 3, 0):
            idx = np.flatnonzero(cls == value)
            list_pick += (
                list(rng.choice(idx, min(per, len(idx)), replace=False))
                if len(idx)
                else []
            )
        names = [f"{CLASSES[int(cls[i])]}_{i}" for i in list_pick]
        chosen = set(plots.select("cross_sections", names))
        list_pick = [i for i, name in zip(list_pick, names) if name in chosen]
        if list_pick:
            pick = np.array(list_pick)
            offsets, z = gis.Hydro.crest_profiles(
                xy[pick], normal[pick], str(vrt_path), p["half_width_m"]
            )
            for row, i in enumerate(pick):
                name = f"{CLASSES[int(cls[i])]}_{i}"
                fig, ax = plt.subplots(figsize=(7, 3.2), layout="constrained")
                ax.plot(offsets, z[row], color=COLOURS[CLASSES[int(cls[i])]])
                for sgn, terrain in (
                    (-1, st["terrain_left"][i]),
                    (1, st["terrain_right"][i]),
                ):
                    ax.plot(
                        sorted([sgn * p["side_near_m"], sgn * p["side_far_m"]]),
                        [terrain] * 2,
                        color=INK["secondary"],
                        linestyle=(0, (4, 2)),
                    )
                ax.axhline(st["crest_z"][i], color=INK["muted"], linewidth=0.8)
                ax.set_xlabel("offset from the line [m] (left negative)")
                ax.set_ylabel("elevation [m]")
                ax.set_title(
                    f"{CLASSES[int(cls[i])].replace('_', ' ')}: line {owner[i]} at {dist[i]:.0f} m, height left "
                    f"{st['height_left'][i]:.2f} m, right {st['height_right'][i]:.2f} m"
                )
                plots.save(
                    fig,
                    "cross_sections",
                    name,
                    {
                        "class": CLASSES[int(cls[i])],
                        "height_left_m": st["height_left"][i],
                        "height_right_m": st["height_right"][i],
                        "footprint_m": st["footprint_m"][i],
                    },
                )
        # crest profiles along lines with embankments, openings marked
        list_line = sorted(
            {int(owner[i]) for i in np.flatnonzero((cls_s == 1) | (cls_s == 2))}
        )
        flagged = {f"line_{x}" for x in gdf_op["line_id"].unique()}
        for name in plots.select(
            "crest_profiles", [f"line_{x}" for x in list_line], flagged
        ):
            x_id = int(name.split("_")[1])
            idx = np.flatnonzero(owner == x_id)
            a, b = idx[0], idx[-1] + 1
            fig, ax = plt.subplots(figsize=(10, 3.6), layout="constrained")
            d = dist[a:b]
            for r0, r1, value in gis.Stats.runs(cls[a:b], min_run):
                ax.axvspan(
                    d[r0] - weight[a] / 2,
                    d[r1 - 1] + weight[a] / 2,
                    color=COLOURS[CLASSES[int(value)]],
                    alpha=0.12,
                    linewidth=0,
                )
            if x_id in d_bounds:
                da, db = d_bounds[x_id]
                ax.plot(
                    d_dist[da:db],
                    d_crest[da:db],
                    color=INK["primary"],
                    linewidth=0.9,
                    label="crest (every pixel)",
                )
            else:
                ax.plot(d, st["crest_z"][a:b], color=INK["primary"], label="crest")
            ax.plot(
                d,
                st["terrain_left"][a:b],
                color=CATEGORICAL[0],
                label="terrain left",
                linestyle=(0, (4, 2)),
            )
            ax.plot(
                d,
                st["terrain_right"][a:b],
                color=CATEGORICAL[1],
                label="terrain right",
                linestyle=(0, (1, 1.5)),
            )
            for r in gdf_op[gdf_op["line_id"] == x_id].itertuples():
                ax.axvline(dist[r.station], color=COLOURS[r.type], linewidth=1.0)
                ax.annotate(
                    f"{r.type} ({r.confidence})",
                    (dist[r.station], ax.get_ylim()[1]),
                    rotation=90,
                    va="top",
                    ha="right",
                    fontsize=7,
                    color=INK["secondary"],
                )
            ax.set_xlabel("distance along the line [m]")
            ax.set_ylabel("elevation [m]")
            ax.set_title(
                f"line {x_id}: crest and terrain of both sides; shading = class"
            )
            ax.legend(loc="lower right")
            plots.save(
                fig,
                "crest_profiles",
                name,
                {
                    "length_m": float(geoms[x_id].length),
                    "openings": int((gdf_op["line_id"] == x_id).sum()),
                },
            )
        # map of classes and openings, summary
        plots.select("openings_map", ["map"])
        fig, ax = plt.subplots(figsize=(8, 8), layout="constrained")
        for name_c, gdf_c in gdf_seg.groupby("class"):
            gdf_c.plot(
                ax=ax,
                color=COLOURS[name_c],
                linewidth=1.6 if name_c in ("embankment", "one_sided") else 0.8,
                label=name_c.replace("_", " "),
            )
        for name_t, gdf_t in gdf_op.groupby("type"):
            gdf_t.plot(
                ax=ax,
                color=COLOURS[name_t],
                markersize=22,
                marker="o",
                label=name_t,
                zorder=3,
            )
        ax.set_aspect("equal")
        ax.legend(loc="upper left", fontsize=7)
        ax.set_title("Transport lines by class and their openings")
        plots.save(
            fig,
            "openings_map",
            "map",
            {"segments": len(gdf_seg), "openings": len(gdf_op)},
        )
        plots.select("summary", ["summary"])
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.6), layout="constrained")
        km = by_class.reindex([c for c in CLASSES.values() if c in by_class.index])
        ax1.barh(
            km.index.str.replace("_", " "),
            km.to_numpy(),
            color=[COLOURS[c] for c in km.index],
        )
        ax1.set_xlabel("length [km]")
        ax1.set_title("Length per class")
        counts = gdf_op.groupby(["type", "confidence"]).size()
        labels = [f"{a} ({b})" for a, b in counts.index]
        ax2.barh(labels, counts.to_numpy(), color=[COLOURS[a] for a, _ in counts.index])
        ax2.set_xlabel("openings")
        ax2.set_title("Openings by type and confidence")
        plots.save(
            fig,
            "summary",
            "summary",
            {
                **{f"km_{k}": float(v) for k, v in by_class.items()},
                "openings": len(gdf_op),
            },
        )
        index = plots.finish()
        log.info("Plots: %s", index)
    times["plots"] = time.perf_counter() - t
    log.info(
        "Done in %.1f s (%s) for %d stations.",
        time.perf_counter() - start,
        ", ".join(f"{k} {v:.1f} s" for k, v in times.items()),
        n,
    )
