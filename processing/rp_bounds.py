# Author: Jakub Zapletal
# Date: 2026-11-05
# Stage F: matches protected areas (and Stage E implied protected areas) to GloFAS and EFAS river cells, turns the GFM
# tests into event return periods per source, combines the sources and bounds the model RP interval; plots.
# LAYER "A": protected_areas_rp.gpkg, "B": protected_areas_by_group_rp.gpkg (outputs with the suffix _b, no implied areas).

import json
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from osgeo import gdal, osr

import gis
from plots import (
    COLOURS,
    INK,
    LINE_WIDTH,
    OUTCOME_MARKERS,
    RATIO_COLOURS,
    SOURCE_NAMES,
    Plots,
)

gdal.UseExceptions()

# --- Run variables (no command-line parameters; DEFENSES_<NAME> environment variables override them) ---
CONFIG_PATH = "config/defenses.toml"
COUNTRY = "XX"
RUN_ID = "base"  # plots in <output_dir>/<country>/plots/discharge/<RUN_ID>
IMPLIED_RUN = "base"  # Stage E run whose implied protected areas get event return periods too (when it exists)
SHOW_PLOTS = False  # show up to implied.plots.show_max figures per type while running (needs a display)
LAYER = "A"  # "A": protected areas of all water sources together, "B": per river group

KEYS = [
    "general.work_dir",
    "general.output_dir",
    "general.log_level",
    "country.rp.classes",
    "country.rivers",
    "country.river_upstream_area_field",
    "country.upstream_area_raster",
    "country.working_resolution_m",
    "discharge.sources",
    "discharge.source_mode",
    "discharge.max_rp_extrapolation",
    "discharge.flooded_window_days",
    "discharge.rp_disagreement_ratio",
    *[
        f"discharge.{s}.{k}"
        for s in ("glofas", "efas")
        for k in (
            "search_radius_km",
            "upstream_ratio_min",
            "upstream_ratio_max",
            "match_river_distance_m",
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
MODEL_COLUMNS = ["rp", "rp_low", "rp_high", "rp_source"]
RP_LINES = [2, 10, 50, 100]  # return periods drawn in the event plots (as specified)
EVENT_DAYS = (
    15  # days before and after the acquisition shown in the event plots (as specified)
)

if __name__ == "__main__":
    run = gis.Run.env_overrides(
        CONFIG_PATH=CONFIG_PATH,
        COUNTRY=COUNTRY,
        RUN_ID=RUN_ID,
        IMPLIED_RUN=IMPLIED_RUN,
        SHOW_PLOTS=SHOW_PLOTS,
        LAYER=LAYER,
    )
    cfg = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], KEYS)
    general, country, dis = cfg["general"], cfg["country"], cfg["discharge"]
    files = gis.Run.layer_files(run["LAYER"])
    sfx = files["suffix"]
    out_dir = Path(general["output_dir"]) / run["COUNTRY"]
    log = gis.Run.logger("rp_bounds", out_dir / "logs", general["log_level"])
    gis.Config.warn_placeholders(cfg)
    start = time.perf_counter()
    list_need = [out_dir / files["rp"], out_dir / f"gfm_evidence{sfx}.csv"]
    if not all(p.exists() for p in list_need):
        raise SystemExit(
            f"Missing inputs {[str(p) for p in list_need if not p.exists()]}; run rp_index.py and gfm_evidence.py first."
        )
    gdf_pa = gpd.read_file(list_need[0])
    # The model values of rp_index.py are kept as *_prior, so that reruns start from them again.
    for column in MODEL_COLUMNS:
        if f"{column}_prior" not in gdf_pa.columns:
            gdf_pa[f"{column}_prior"] = gdf_pa[column]
        gdf_pa[column] = gdf_pa[f"{column}_prior"]
    wkt = gdf_pa.crs.to_wkt()
    classes = np.asarray(country["rp"]["classes"], float)
    sources = [
        s
        for s in dis["sources"]
        if (out_dir / f"discharge_{s}.npz").exists()
        and (out_dir / f"discharge_cells_{s}.gpkg").exists()
    ]
    for s in sorted(set(dis["sources"]) - set(sources)):
        log.warning(
            "%s: discharge_%s.npz missing (download_discharge.py, discharge_prepare.py); source not used.",
            s,
            s,
        )
    if not sources:
        log.warning("No discharge source prepared: no observed bounds.")
        raise SystemExit(0)

    # --- Sites: explicit protected areas at their source river cell, implied protected areas at their river ---
    field = country["river_upstream_area_field"]
    # river lines of the same-river test: the Stage E thalweg, else the channel lines of protected_areas.py
    for name in ("rivers_snapped.gpkg", "rivers_channel.gpkg", ""):
        if not name or (out_dir / name).exists():
            break
    gdf_line = (
        gpd.read_file(out_dir / name).to_crs(wkt)
        if name
        else gis.Vector.read_merged(country["rivers"], wkt, keep_columns=True)
    )
    gdf_area = None
    if field:
        gdf_area = gis.Vector.read_merged(country["rivers"], wkt, keep_columns=[field])
        if field not in gdf_area.columns:
            raise SystemExit(
                f"river_upstream_area_field {field!r} is not a column of the rivers."
            )
    sites = pd.DataFrame(
        {
            "protected_area_id": gdf_pa["protected_area_id"],
            "kind": "explicit",
            "x": gdf_pa["source_x"].astype(float),
            "y": gdf_pa["source_y"].astype(float),
        }
    )
    # implied protected areas (Stage E) belong to file A
    imp_dir = out_dir / "implied" / run["IMPLIED_RUN"]
    if files["layer"] == "B":
        imp_dir = out_dir / "implied" / "_none_"
    if (imp_dir / "implied_protected_areas.gpkg").exists() and (
        out_dir / "river_points.gpkg"
    ).exists():
        gdf_imp = gpd.read_file(imp_dir / "implied_protected_areas.gpkg").to_crs(wkt)
        if len(gdf_imp):
            gdf_pts = gpd.read_file(out_dir / "river_points.gpkg").to_crs(wkt)
            near = shapely.STRtree(gdf_pts.geometry.values).query_nearest(
                shapely.point_on_surface(gdf_imp.geometry.values), all_matches=False
            )[1]
            xy = shapely.get_coordinates(gdf_pts.geometry.values[near])
            sites = pd.concat(
                [
                    sites,
                    pd.DataFrame(
                        {
                            "protected_area_id": gdf_imp["implied_id"],
                            "kind": "implied",
                            "x": xy[:, 0],
                            "y": xy[:, 1],
                        }
                    ),
                ],
                ignore_index=True,
            )
            log.info(
                "Implied protected areas of run %s: %d sites at their nearest river point.",
                run["IMPLIED_RUN"],
                len(gdf_imp),
            )
    pts = shapely.points(sites["x"].to_numpy(), sites["y"].to_numpy())
    # at a confluence a river cell can lie nearest to the tributary: within one working cell of the nearest
    # line, the line with the largest upstream area is the river of the site
    tol = country["working_resolution_m"]
    sites["river"] = gis.Vector.dominant_nearby(
        pts,
        gdf_line.geometry.values,
        (
            pd.to_numeric(gdf_line[field], errors="coerce")
            if field and field in gdf_line.columns
            else np.zeros(len(gdf_line))
        ),
        tol,
    )
    if gdf_area is not None:
        sites["area_km2"] = pd.to_numeric(
            gdf_area[field].to_numpy()[
                gis.Vector.dominant_nearby(
                    pts,
                    gdf_area.geometry.values,
                    pd.to_numeric(gdf_area[field], errors="coerce"),
                    tol,
                )
            ],
            errors="coerce",
        )
    elif country["upstream_area_raster"]:
        ds_ua = gdal.Open(country["upstream_area_raster"])
        to_r = osr.CoordinateTransformation(
            osr.SpatialReference(wkt=wkt),
            osr.SpatialReference(wkt=ds_ua.GetProjection()),
        )
        to_r.GetSourceCS().SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        values = []
        ugt = ds_ua.GetGeoTransform()
        band = ds_ua.GetRasterBand(1)
        for x, y in zip(sites["x"], sites["y"]):
            rx, ry, _ = to_r.TransformPoint(float(x), float(y))
            c, r = int((rx - ugt[0]) / ugt[1]), int((ry - ugt[3]) / ugt[5])
            ok = 0 <= c < ds_ua.RasterXSize and 0 <= r < ds_ua.RasterYSize
            values.append(float(band.ReadAsArray(c, r, 1, 1)[0, 0]) if ok else np.nan)
        sites["area_km2"] = values
    else:
        log.warning(
            "No river_upstream_area_field and no upstream_area_raster: matching skipped, no observed bounds."
        )
        raise SystemExit(0)
    wgs84 = osr.SpatialReference()
    wgs84.ImportFromEPSG(4326)
    wgs84.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    local = osr.SpatialReference(wkt=wkt)
    local.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    to_geo = osr.CoordinateTransformation(local, wgs84)
    to_local = osr.CoordinateTransformation(wgs84, local)
    lonlat = np.array(to_geo.TransformPoints(np.c_[sites["x"], sites["y"]].tolist()))[
        :, :2
    ]
    sites["lon"], sites["lat"] = lonlat[:, 0], lonlat[:, 1]

    # --- 1. Matching per site and source ---
    t = time.perf_counter()
    data, list_match = {}, []
    for src in sources:
        # every array decompressed once (indexing an open npz file decompresses the array on each access)
        with np.load(out_dir / f"discharge_{src}.npz") as f:
            npz = {k: f[k] for k in f.files}
        gdf_cells = (
            gpd.read_file(out_dir / f"discharge_cells_{src}.gpkg")
            .set_index("cell_id")
            .loc[npz["cell_id"]]
            .reset_index()
        )
        cxy = np.array(
            to_local.TransformPoints(np.c_[npz["lon"], npz["lat"]].tolist())
        )[:, :2]
        data[src] = {
            "npz": npz,
            "cells": gdf_cells,
            "xy": cxy,
            "dates": npz["dates"].astype("datetime64[D]"),
        }
        c = dis[src]
        cell_pts = shapely.points(cxy[:, 0], cxy[:, 1])
        for col in ("cell", "distance_km", "ratio", "same_river", "match_quality"):
            sites[f"{src}_{col}"] = None
        for i, site in sites.iterrows():
            same = (
                shapely.distance(gdf_line.geometry.values[site["river"]], cell_pts)
                <= c["match_river_distance_m"]
            )
            m = gis.Discharge.match(
                site["lon"],
                site["lat"],
                site["area_km2"] if np.isfinite(site["area_km2"]) else -1.0,
                npz["lon"],
                npz["lat"],
                npz["upstream_area_km2"],
                c["search_radius_km"],
                c["upstream_ratio_min"],
                c["upstream_ratio_max"],
                same,
            )
            sites.at[i, f"{src}_cell"] = m["index"]
            sites.at[i, f"{src}_distance_km"] = m["distance_km"]
            sites.at[i, f"{src}_ratio"] = m["ratio"]
            sites.at[i, f"{src}_same_river"] = m["same_river"]
            sites.at[i, f"{src}_match_quality"] = m["quality"]
            if m["index"] >= 0:
                list_match.append(
                    {
                        "protected_area_id": site["protected_area_id"],
                        "kind": site["kind"],
                        "source": src,
                        "cell_id": int(npz["cell_id"][m["index"]]),
                        "distance_km": m["distance_km"],
                        "ratio": m["ratio"],
                        "same_river": m["same_river"],
                        "match_quality": m["quality"],
                        "geometry": shapely.LineString(
                            [(site["x"], site["y"]), tuple(cxy[m["index"]])]
                        ),
                    }
                )
        q = sites[f"{src}_match_quality"].value_counts().to_dict()
        log.info("%s matches: %s.", src, q)
    gdf_match = gpd.GeoDataFrame(
        list_match
        or [
            {
                "protected_area_id": "",
                "kind": "",
                "source": "",
                "cell_id": -1,
                "distance_km": np.nan,
                "ratio": np.nan,
                "same_river": False,
                "match_quality": "",
                "geometry": shapely.LineString([(0, 0), (0, 0)]),
            }
        ],
        geometry="geometry",
        crs=wkt,
    )
    gis.Vector.write_gpkg(
        gdf_match.iloc[: len(list_match)], out_dir / f"discharge_matches{sfx}.gpkg"
    )
    match_s = time.perf_counter() - t

    # --- 2. Event discharge and return period per test and matched source ---
    df_test = pd.read_csv(list_need[1])
    if (imp_dir / "implied_gfm_evidence.csv").exists():
        df_test = pd.concat(
            [df_test, pd.read_csv(imp_dir / "implied_gfm_evidence.csv")],
            ignore_index=True,
        )
    df_test = df_test[df_test["outcome"].isin(["held", "flooded"])].drop_duplicates(
        ["protected_area_id", "item_id"]
    )
    site_of = {cid: i for i, cid in enumerate(sites["protected_area_id"])}
    acq = (
        pd.to_datetime(df_test["datetime"], utc=True, format="ISO8601")
        .dt.tz_localize(None)
        .dt.floor("D")
    )
    window = dis["flooded_window_days"]
    list_event = []
    for (_, row), day in zip(df_test.iterrows(), acq):
        i = site_of.get(row["protected_area_id"])
        if i is None:
            continue
        for src in sources:
            cell = sites.at[i, f"{src}_cell"]
            if cell is None or cell < 0:
                continue
            d = data[src]
            k = int((np.datetime64(day.date(), "D") - d["dates"][0]).astype(int))
            if not 0 <= k < len(d["dates"]):
                continue
            series = d["npz"]["values"][:, cell]
            q = (
                series[k]
                if row["outcome"] == "held"
                else np.nanmax(series[max(k - window, 0) : k + 1])
            )
            xi, alpha = d["cells"].loc[cell, ["gumbel_xi", "gumbel_alpha"]]
            if not np.isfinite(q) or not np.isfinite(alpha):
                continue
            rp_event = float(gis.Stats.gumbel_rp(q, xi, alpha))
            extrapolated = rp_event > dis["max_rp_extrapolation"]
            list_event.append(
                {
                    "protected_area_id": row["protected_area_id"],
                    "item_id": row["item_id"],
                    "datetime": row["datetime"],
                    "outcome": row["outcome"],
                    "source": src,
                    "cell_id": int(d["npz"]["cell_id"][cell]),
                    "q_event": float(q),
                    "rp_event": min(rp_event, dis["max_rp_extrapolation"]),
                    "extrapolated": bool(extrapolated),
                    "match_quality": sites.at[i, f"{src}_match_quality"],
                }
            )
    df_ev = pd.DataFrame(
        list_event,
        columns=[
            "protected_area_id",
            "item_id",
            "datetime",
            "outcome",
            "source",
            "cell_id",
            "q_event",
            "rp_event",
            "extrapolated",
            "match_quality",
        ],
    )
    df_ev.drop(columns="match_quality").to_csv(
        out_dir / f"event_rp{sfx}.csv", index=False
    )

    # --- 3. One return period per test from the sources (source_mode) ---
    list_comb = []
    for (cid, item), df_t in df_ev.groupby(["protected_area_id", "item_id"]):
        by = df_t.set_index("source")
        outcome = df_t["outcome"].iloc[0]
        mode = dis["source_mode"]
        if mode in ("prefer_efas", "prefer_glofas"):
            first, other = (
                ("efas", "glofas") if mode == "prefer_efas" else ("glofas", "efas")
            )
            # the preferred source, unless only its match is weak and the other source matches well
            pick = (
                first
                if first in by.index
                and (
                    by.at[first, "match_quality"] == "good"
                    or other not in by.index
                    or by.at[other, "match_quality"] != "good"
                )
                else other
            )
            rp_c, used = float(by.at[pick, "rp_event"]), pick
        else:
            # conservative: held takes the lower RP, flooded the higher; the bounds stay valid when the sources disagree
            rp_c = float(
                by["rp_event"].min() if outcome == "held" else by["rp_event"].max()
            )
            used = "both" if len(by) == 2 else by.index[0]
        rps = by["rp_event"].to_numpy()
        disagree = (
            len(rps) == 2
            and max(rps) / max(min(rps), 1e-9) > dis["rp_disagreement_ratio"]
        )
        list_comb.append(
            {
                "protected_area_id": cid,
                "item_id": item,
                "outcome": outcome,
                "rp_combined": rp_c,
                "used": used,
                "source_disagreement": disagree,
                **{
                    f"rp_{s}": float(by.at[s, "rp_event"]) if s in by.index else np.nan
                    for s in ("glofas", "efas")
                },
            }
        )
    df_comb = pd.DataFrame(
        list_comb,
        columns=[
            "protected_area_id",
            "item_id",
            "outcome",
            "rp_combined",
            "used",
            "source_disagreement",
            "rp_glofas",
            "rp_efas",
        ],
    )
    df_comb.to_csv(out_dir / f"event_rp_combined{sfx}.csv", index=False)
    log.info(
        "Event return periods: %d per source and test, %d combined tests (%s), %d source disagreements.",
        len(df_ev),
        len(df_comb),
        dis["source_mode"],
        int(df_comb["source_disagreement"].sum()),
    )

    # --- 4. Bounds per explicit protected area and combination with the model interval (logic of Stage C) ---
    dc = df_comb[df_comb["protected_area_id"].isin(gdf_pa["protected_area_id"])]
    held = dc[dc["outcome"] == "held"].groupby("protected_area_id")
    flooded = dc[dc["outcome"] == "flooded"].groupby("protected_area_id")
    gdf_pa["n_held"] = (
        gdf_pa["protected_area_id"]
        .map(held["rp_combined"].count())
        .fillna(0)
        .astype(int)
    )
    gdf_pa["n_flooded"] = (
        gdf_pa["protected_area_id"]
        .map(flooded["rp_combined"].count())
        .fillna(0)
        .astype(int)
    )
    gdf_pa["obs_rp_low"] = gdf_pa["protected_area_id"].map(held["rp_combined"].max())
    gdf_pa["obs_rp_high"] = gdf_pa["protected_area_id"].map(
        flooded["rp_combined"].min()
    )
    for s in ("glofas", "efas"):
        gdf_pa[f"rp_obs_low_{s}"] = gdf_pa["protected_area_id"].map(
            held[f"rp_{s}"].max()
        )
        gdf_pa[f"rp_obs_high_{s}"] = gdf_pa["protected_area_id"].map(
            flooded[f"rp_{s}"].min()
        )
    used = dc.groupby("protected_area_id")["used"].agg(
        lambda u: "both" if "both" in set(u) or len(set(u)) > 1 else next(iter(set(u)))
    )
    gdf_pa["obs_source"] = gdf_pa["protected_area_id"].map(used).fillna("")
    gdf_pa["source_disagreement"] = (
        gdf_pa["protected_area_id"]
        .map(dc.groupby("protected_area_id")["source_disagreement"].any())
        .fillna(False)
        .astype(bool)
    )
    gdf_pa["obs_conflict"] = (gdf_pa["obs_rp_low"] > gdf_pa["obs_rp_high"]).fillna(
        False
    )
    gdf_pa["model_obs_conflict"] = False
    for s in sources:
        part = sites[sites["kind"] == "explicit"].set_index("protected_area_id")
        for col in ("cell", "distance_km", "ratio", "same_river", "match_quality"):
            values = gdf_pa["protected_area_id"].map(part[f"{s}_{col}"])
            if col == "cell":
                gdf_pa[f"{s}_cell_id"] = [
                    (
                        int(data[s]["npz"]["cell_id"][int(v)])
                        if v is not None and int(v) >= 0
                        else -1
                    )
                    for v in values
                ]
            else:
                gdf_pa[f"{s}_{col}"] = values.to_numpy()
    cases = {
        "combined": 0,
        "model_obs_conflict": 0,
        "obs_conflict": 0,
        "override_kept": 0,
    }
    log_classes = np.log(classes)
    for i in gdf_pa.index[(gdf_pa["n_held"] + gdf_pa["n_flooded"]) > 0]:
        row = gdf_pa.loc[i]
        if row["rp_source_prior"] == "override":
            cases["override_kept"] += 1
            continue
        if row["obs_conflict"]:
            cases["obs_conflict"] += 1
            continue
        obs_low = row["obs_rp_low"] if np.isfinite(row["obs_rp_low"]) else -np.inf
        obs_high = row["obs_rp_high"] if np.isfinite(row["obs_rp_high"]) else np.inf
        low = max(row["rp_low_prior"], obs_low)
        high = min(row["rp_high_prior"], obs_high)
        if low > high:
            # The observations contradict the model interval: only the observed interval is used, a side without
            # observations stays open.
            gdf_pa.loc[i, "model_obs_conflict"] = True
            cases["model_obs_conflict"] += 1
            low, high = obs_low, obs_high
        else:
            cases["combined"] += 1
        value = float(np.clip(row["rp_prior"], low, high))
        # Snap to the class inside [low, high] closest to the value (log scale), else to the closest class.
        inside = (classes >= low) & (classes <= high)
        pool = np.flatnonzero(inside) if inside.any() else np.arange(len(classes))
        k = pool[np.argmin(np.abs(log_classes[pool] - np.log(value)))]
        snapped_rp = classes[k]
        # An open side of the interval becomes the neighbouring class of the result.
        if not np.isfinite(low):
            low = classes[max(k - 1, 0)]
        if not np.isfinite(high):
            high = classes[min(k + 1, len(classes) - 1)]
        gdf_pa.loc[i, "rp"] = snapped_rp
        gdf_pa.loc[i, "rp_low"] = min(low, snapped_rp)
        gdf_pa.loc[i, "rp_high"] = max(high, snapped_rp)
        gdf_pa.loc[i, "rp_source"] = f"{row['rp_source_prior']}+obs"
    for column in ["rp", "rp_low", "rp_high"]:
        gdf_pa[column] = gdf_pa[column].astype(float)
    gis.Vector.write_gpkg(gdf_pa, out_dir / files["rp"])
    narrowed = int(
        (
            (gdf_pa["rp_low"] > gdf_pa["rp_low_prior"])
            | (gdf_pa["rp_high"] < gdf_pa["rp_high_prior"])
        ).sum()
    )
    log.info(
        "Observed bounds: %s; %d intervals narrowed. "
        + files["rp"]
        + " updated (%d protected areas).",
        ", ".join(f"{k} {v}" for k, v in cases.items()),
        narrowed,
        len(gdf_pa),
    )

    # --- Plots: match_map, fit, event, source_comparison, summary ---
    t = time.perf_counter()
    settings = cfg["implied"]["plots"]
    if settings["enabled"]:
        plots = Plots(
            out_dir / "plots" / "discharge" / (run["RUN_ID"] + sfx),
            settings,
            run["SHOW_PLOTS"],
            log,
        )
        plt = plots.plt
        from matplotlib.cm import ScalarMappable
        from matplotlib.colors import LinearSegmentedColormap, LogNorm
        from matplotlib.dates import AutoDateLocator, ConciseDateFormatter
        from matplotlib.patches import Circle
        from matplotlib.ticker import MaxNLocator

        cmap_ratio = LinearSegmentedColormap.from_list("ratio", RATIO_COLOURS)
        norm_ratio = LogNorm(0.25, 4.0)
        matched = sites[
            [
                any(
                    sites.at[i, f"{s}_match_quality"] in ("good", "weak")
                    for s in sources
                )
                for i in sites.index
            ]
        ]
        flag_m = [
            c
            for c in matched["protected_area_id"]
            if any(
                sites.loc[site_of[c], f"{s}_match_quality"] == "weak" for s in sources
            )
        ]
        # match_map: one panel per source in offsets from the source river cell, cells coloured by the upstream
        # area ratio, the search radius and the band that counts as the same river
        for name in plots.select(
            "match_map", list(matched["protected_area_id"]), flag_m
        ):
            site = sites.loc[site_of[name]]
            reach_km = max(dis[s]["search_radius_km"] for s in sources)
            box = shapely.box(
                site["x"] - 1000 * reach_km,
                site["y"] - 1000 * reach_km,
                site["x"] + 1000 * reach_km,
                site["y"] + 1000 * reach_km,
            )
            origin = np.array([site["x"], site["y"]])
            own = shapely.intersection(gdf_line.geometry.values[site["river"]], box)
            fig, axes = plt.subplots(
                1,
                len(sources),
                figsize=(4.4 * len(sources) + 1.4, 5.6),
                layout="constrained",
                squeeze=False,
            )
            for ax, s in zip(axes[0], sources):
                band = shapely.intersection(
                    shapely.buffer(own, dis[s]["match_river_distance_m"]), box
                )
                for g in shapely.get_parts(band):
                    if shapely.get_type_id(g) == 3 and not g.is_empty:
                        ax.fill(
                            *((np.asarray(g.exterior.coords) - origin) / 1000).T,
                            color=INK["grid"],
                            lw=0,
                            label="same-river band",
                        )
                for g in shapely.get_parts(
                    shapely.intersection(gdf_line.geometry.values, box)
                ):
                    if shapely.get_type_id(g) == 1 and not g.is_empty:
                        ax.plot(
                            *((np.asarray(g.coords) - origin) / 1000).T,
                            color=INK["muted"],
                            lw=0.8,
                            label="rivers",
                        )
                for g in shapely.get_parts(own):
                    if shapely.get_type_id(g) == 1 and not g.is_empty:
                        ax.plot(
                            *((np.asarray(g.coords) - origin) / 1000).T,
                            color=INK["secondary"],
                            lw=1.6,
                            label="river of the protected area",
                        )
                c = dis[s]
                ax.add_patch(
                    Circle(
                        (0.0, 0.0),
                        c["search_radius_km"],
                        fill=False,
                        ec=INK["secondary"],
                        lw=0.8,
                        ls=(0, (4, 2)),
                        label="search radius",
                    )
                )
                d = data[s]
                off = (d["xy"] - origin) / 1000
                in_box = np.all(np.abs(off) <= reach_km, axis=1)
                near = in_box & (
                    np.hypot(off[:, 0], off[:, 1]) <= c["search_radius_km"]
                )
                ax.scatter(
                    *off[in_box & ~near].T,
                    s=22,
                    facecolors="none",
                    edgecolors=INK["muted"],
                    linewidths=0.8,
                    label="cells beyond the radius",
                )
                ax.scatter(
                    *off[near].T,
                    c=d["npz"]["upstream_area_km2"][near] / max(site["area_km2"], 1e-9),
                    cmap=cmap_ratio,
                    norm=norm_ratio,
                    s=56,
                    edgecolors=INK["secondary"],
                    linewidths=0.6,
                    label="cells within the radius",
                )
                k = site[f"{s}_cell"]
                chosen = k is not None and k >= 0
                if chosen:
                    ax.scatter(
                        *off[k],
                        s=200,
                        facecolors="none",
                        edgecolors=INK["primary"],
                        linewidths=1.6,
                        label="chosen cell",
                    )
                ax.plot(
                    0.0,
                    0.0,
                    "*",
                    color=INK["primary"],
                    ms=13,
                    mec=INK["surface"],
                    label=(
                        "site: source river cell"
                        if site["kind"] == "explicit"
                        else "site: nearest river point"
                    ),
                )
                ax.set_xlim(-reach_km, reach_km)
                ax.set_ylim(-reach_km, reach_km)
                ax.set_aspect("equal")
                ax.set_xlabel("Easting from the site (km)")
                ax.set_title(
                    f"{SOURCE_NAMES[s]} · radius {c['search_radius_km']:g} km · ratio "
                    f"{c['upstream_ratio_min']:g}–{c['upstream_ratio_max']:g} · same river ≤ {c['match_river_distance_m']:g} m\n"
                    + (
                        f"cell {int(d['npz']['cell_id'][k])}: {site[f'{s}_match_quality']}, ratio "
                        f"{site[f'{s}_ratio']:.2f}, {site[f'{s}_distance_km']:.2f} km"
                        if chosen
                        else "no cell matched"
                    ),
                    loc="left",
                    fontsize=9,
                )
            axes[0][0].set_ylabel("Northing from the site (km)")
            bar = fig.colorbar(
                ScalarMappable(norm_ratio, cmap_ratio),
                ax=axes[0].tolist(),
                shrink=0.8,
                label="Upstream area ratio (cell / river)",
            )
            bar.set_ticks([0.25, 0.5, 1, 2, 4], labels=["0.25", "0.5", "1", "2", "4"])
            bar.ax.minorticks_off()
            handles = {}
            for ax in axes[0]:
                for h, lbl in zip(*ax.get_legend_handles_labels()):
                    handles.setdefault(lbl, h)
            fig.legend(
                handles.values(),
                handles.keys(),
                loc="outside lower center",
                ncols=4,
                fontsize=8,
            )
            fig.suptitle(
                f"match_map · {name} · run {run['RUN_ID']} · river upstream area {site['area_km2']:.0f} km²",
                x=0.01,
                ha="left",
            )
            plots.save(
                fig,
                "match_map",
                name,
                {
                    f"{s}_{k}": site[f"{s}_{k}"]
                    for s in sources
                    for k in ("match_quality", "ratio", "distance_km")
                },
            )
        # fit: annual maxima of the matched cells at Gringorten positions with their Gumbel curves
        rp_line = np.logspace(
            np.log10(1.01), np.log10(dis["max_rp_extrapolation"]), 200
        )
        for name in plots.select("fit", list(matched["protected_area_id"]), flag_m):
            site = sites.loc[site_of[name]]
            fig, ax = plt.subplots(figsize=(7.5, 4.8), layout="constrained")
            parts = []
            for s in sources:
                k = site[f"{s}_cell"]
                if k is None or k < 0:
                    continue
                cell = data[s]["cells"].loc[k]
                x = np.sort(
                    data[s]["npz"]["maxima"][:, k][
                        np.isfinite(data[s]["npz"]["maxima"][:, k])
                    ]
                )
                n = len(x)
                ax.plot(
                    rp_line,
                    gis.Stats.gumbel_q(
                        rp_line, cell["gumbel_xi"], cell["gumbel_alpha"]
                    ),
                    color=COLOURS["source"][s],
                    lw=LINE_WIDTH,
                    label=f"{SOURCE_NAMES[s]} Gumbel fit",
                )
                ax.plot(
                    1.0 / (1.0 - (np.arange(1, n + 1) - 0.44) / (n + 0.12)),
                    x,
                    "o",
                    mfc="none",
                    color=COLOURS["source"][s],
                    ms=5,
                    label=f"{SOURCE_NAMES[s]} annual maxima",
                )
                parts.append(
                    f"{SOURCE_NAMES[s]} cell {int(cell['cell_id'])}: {int(cell['n_years'])} years, "
                    f"fit_rmse {cell['fit_rmse']:.1f} m³/s"
                )
            Plots.rp_axis(ax, "x", dis["max_rp_extrapolation"])
            ax.set_xlabel("Return period (years)")
            ax.set_ylabel("Discharge (m³/s)")
            ax.legend(loc="upper left", fontsize=8)
            ax.set_title(f"fit · {name} · run {run['RUN_ID']}\n" + "; ".join(parts))
            plots.save(
                fig,
                "fit",
                name,
                {
                    "sources": ",".join(
                        s
                        for s in sources
                        if site[f"{s}_cell"] is not None and site[f"{s}_cell"] >= 0
                    )
                },
            )
        # event: one panel per matched source (own y-axis) with the RP thresholds of its fit and the discharge the
        # event return period comes from
        tests_ids = [
            f"{c}__{i}"
            for c, i in zip(df_comb["protected_area_id"], df_comb["item_id"])
        ]
        flag_e = [t_ for t_, d_ in zip(tests_ids, df_comb["source_disagreement"]) if d_]
        for name in plots.select("event", tests_ids, flag_e):
            cid, item = name.split("__")
            row = df_comb[
                (df_comb["protected_area_id"] == cid) & (df_comb["item_id"] == item)
            ].iloc[0]
            df_e = df_ev[
                (df_ev["protected_area_id"] == cid) & (df_ev["item_id"] == item)
            ]
            day = np.datetime64(
                pd.Timestamp(df_e["datetime"].iloc[0]).tz_localize(None).date(), "D"
            )
            list_src = [s for s in sources if (df_e["source"] == s).any()]
            fig, axes = plt.subplots(
                len(list_src),
                1,
                figsize=(9, 1.4 + 2.4 * len(list_src)),
                sharex=True,
                layout="constrained",
                squeeze=False,
            )
            for ax, s in zip(axes[:, 0], list_src):
                e = df_e[df_e["source"] == s].iloc[0]
                d = data[s]
                k = int(sites.at[site_of[cid], f"{s}_cell"])
                k0 = int((day - d["dates"][0]).astype(int))
                lo, hi = max(k0 - EVENT_DAYS, 0), min(
                    k0 + EVENT_DAYS + 1, len(d["dates"])
                )
                if row["outcome"] == "flooded" and window > 0:
                    ax.axvspan(
                        day - np.timedelta64(window, "D"),
                        day,
                        color=INK["grid"],
                        lw=0,
                        label=f"flooded window ({window} days before the acquisition)",
                    )
                ax.plot(
                    d["dates"][lo:hi],
                    d["npz"]["values"][lo:hi, k],
                    color=COLOURS["source"][s],
                    lw=LINE_WIDTH,
                    label="daily discharge",
                )
                cell = d["cells"].loc[k]
                for rp_l in RP_LINES:
                    q_l = gis.Stats.gumbel_q(
                        rp_l, cell["gumbel_xi"], cell["gumbel_alpha"]
                    )
                    ax.axhline(
                        q_l,
                        color=COLOURS["source"][s],
                        lw=0.7,
                        ls=(0, (1, 1.5)),
                        label="RP thresholds of the Gumbel fit",
                    )
                    ax.text(
                        1.006,
                        q_l,
                        f"RP {rp_l}",
                        transform=ax.get_yaxis_transform(),
                        fontsize=7,
                        color=INK["secondary"],
                        va="center",
                        ha="left",
                    )
                ax.axvline(day, color=INK["primary"], lw=1.0, label="acquisition")
                # the day whose discharge gave the event RP: the acquisition day (held) or the maximum of the window
                k_ev = (
                    k0
                    if row["outcome"] == "held"
                    else max(k0 - window, 0)
                    + int(
                        np.nanargmax(
                            d["npz"]["values"][max(k0 - window, 0) : k0 + 1, k]
                        )
                    )
                )
                ax.plot(
                    d["dates"][k_ev],
                    e["q_event"],
                    "o",
                    ms=8,
                    color=COLOURS["source"][s],
                    mec=INK["surface"],
                    mew=1.5,
                    label="event discharge",
                )
                ax.annotate(
                    f"{e['q_event']:.0f} m³/s · RP "
                    + (
                        f"≥ {e['rp_event']:g} (capped)"
                        if e["extrapolated"]
                        else f"{e['rp_event']:.1f}"
                    ),
                    (d["dates"][k_ev], e["q_event"]),
                    xytext=(7, 0),
                    textcoords="offset points",
                    va="center",
                    fontsize=8,
                    color=INK["primary"],
                    bbox={
                        "boxstyle": "square,pad=0.15",
                        "fc": INK["surface"],
                        "ec": "none",
                        "alpha": 0.85,
                    },
                )
                ax.set_ylabel("Discharge (m³/s)")
                ax.set_title(
                    f"{SOURCE_NAMES[s]} cell {int(e['cell_id'])} ({sites.at[site_of[cid], f'{s}_match_quality']} match)",
                    loc="left",
                    fontsize=9,
                )
            locator = AutoDateLocator()
            axes[-1, 0].xaxis.set_major_locator(locator)
            axes[-1, 0].xaxis.set_major_formatter(ConciseDateFormatter(locator))
            handles = {}
            for ax in axes[:, 0]:
                for h, lbl in zip(*ax.get_legend_handles_labels()):
                    if lbl != "daily discharge":
                        handles.setdefault(lbl, h)
            fig.legend(
                handles.values(),
                handles.keys(),
                loc="outside lower center",
                ncols=4,
                fontsize=8,
            )
            fig.suptitle(
                f"event · {cid} · {item} · run {run['RUN_ID']}\n{row['outcome']}; event RP "
                + ", ".join(f"{SOURCE_NAMES[s]} {row[f'rp_{s}']:.1f}" for s in list_src)
                + f" → {row['rp_combined']:.1f} ({dis['source_mode']})"
                + (", sources disagree" if row["source_disagreement"] else ""),
                x=0.01,
                ha="left",
            )
            plots.save(
                fig,
                "event",
                name,
                {
                    k: row[k]
                    for k in (
                        "outcome",
                        "rp_glofas",
                        "rp_efas",
                        "rp_combined",
                        "source_disagreement",
                    )
                },
            )
        # source_comparison: event RP of GloFAS against EFAS for the tests with both sources
        both = df_comb.dropna(subset=["rp_glofas", "rp_efas"])
        plots.select("source_comparison", ["source_comparison"])
        fig, ax = plt.subplots(figsize=(6.5, 6), layout="constrained")
        lim = np.array([1.0, dis["max_rp_extrapolation"]])
        ax.plot(lim, lim, color=INK["muted"], lw=0.8, label="equal")
        for j, f_ in enumerate(
            (dis["rp_disagreement_ratio"], 1 / dis["rp_disagreement_ratio"])
        ):
            ax.plot(
                lim,
                lim * f_,
                color=INK["muted"],
                lw=0.8,
                ls=(0, (4, 2)),
                label=(
                    f"disagreement ratio {dis['rp_disagreement_ratio']:g}"
                    if j == 0
                    else None
                ),
            )
        for outcome in ("held", "flooded"):
            sel = both[both["outcome"] == outcome]
            ax.scatter(
                sel["rp_glofas"],
                sel["rp_efas"],
                marker=OUTCOME_MARKERS[outcome],
                color=COLOURS["outcome"][outcome],
                s=48,
                edgecolors=INK["surface"],
                linewidths=1.0,
                label=f"{outcome} ({len(sel)})",
            )
        Plots.rp_axis(ax, "x", dis["max_rp_extrapolation"])
        Plots.rp_axis(ax, "y", dis["max_rp_extrapolation"])
        ax.set_xlim(lim * [0.9, 1.1])
        ax.set_ylim(lim * [0.9, 1.1])
        ax.set_aspect("equal")
        ax.set_xlabel("Event RP, GloFAS (years)")
        ax.set_ylabel("Event RP, EFAS (years)")
        ax.legend(loc="upper left", fontsize=8)
        ax.set_title(
            f"source_comparison · run {run['RUN_ID']}\n{len(both)} tests with both sources, "
            f"{int(both['source_disagreement'].sum())} beyond the disagreement ratio"
        )
        plots.save(
            fig, "source_comparison", "source_comparison", {"tests_both": len(both)}
        )
        # summary: matches by quality, event RP distribution per source, protected area flags
        plots.select("summary", ["summary"])
        fig, axes = plt.subplots(1, 3, figsize=(14, 4.6), layout="constrained")
        qual = (
            pd.DataFrame(
                {s: sites[f"{s}_match_quality"].value_counts() for s in sources}
            )
            .reindex(["good", "weak", "none"])
            .fillna(0)
        )
        x = np.arange(len(qual))
        width = 0.16  # about 24 px per bar at 120 dpi
        for j, s in enumerate(sources):
            bars = axes[0].bar(
                x + (j - (len(sources) - 1) / 2) * (width + 0.03),
                qual[s],
                width=width,
                color=COLOURS["source"][s],
                label=SOURCE_NAMES[s],
            )
            axes[0].bar_label(
                bars, fmt="%d", fontsize=8, color=INK["secondary"], padding=2
            )
        axes[0].set_xticks(x, qual.index)
        axes[0].set_ylabel("Sites (protected areas)")
        axes[0].set_title(f"Matches by quality ({len(sites)} sites)", loc="left")
        axes[0].legend(fontsize=8)
        bins = np.logspace(0, np.log10(dis["max_rp_extrapolation"]), 13)
        for s in sources:
            v = df_ev[df_ev["source"] == s]["rp_event"]
            if len(v):
                axes[1].hist(
                    np.clip(v, 1.0, dis["max_rp_extrapolation"]),
                    bins=bins,
                    histtype="step",
                    lw=LINE_WIDTH + 0.4,
                    color=COLOURS["source"][s],
                    label=f"{SOURCE_NAMES[s]} ({len(v)})",
                )
        Plots.rp_axis(axes[1], "x", dis["max_rp_extrapolation"])
        axes[1].set_xlabel("Event RP (years)")
        axes[1].set_ylabel("Tests")
        axes[1].set_title("Event return periods per source", loc="left")
        axes[1].legend(fontsize=8)
        counts = {
            "narrowed": narrowed,
            "obs_conflict": cases["obs_conflict"],
            "model_obs_conflict": cases["model_obs_conflict"],
            "source_disagreement": int(gdf_pa["source_disagreement"].sum()),
        }
        bars = axes[2].bar(
            range(len(counts)), list(counts.values()), color=INK["secondary"], width=0.3
        )
        axes[2].bar_label(bars, fmt="%d", fontsize=8, color=INK["secondary"], padding=2)
        axes[2].set_xticks(range(len(counts)), list(counts.keys()), rotation=15)
        axes[2].set_ylabel("Protected areas")
        axes[2].set_title(
            f"Explicit protected areas ({len(gdf_pa)}, {int(((gdf_pa['n_held'] + gdf_pa['n_flooded']) > 0).sum())} with tests)",
            loc="left",
        )
        for ax in axes:
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
            ax.margins(y=0.12)
        for ax in (axes[0], axes[2]):
            ax.grid(axis="x", visible=False)
        fig.suptitle(f"summary · run {run['RUN_ID']}", x=0.01, ha="left")
        plots.save(fig, "summary", "summary", counts)
        log.info("Plots: %s", plots.finish())
    perf = {
        "sites": len(sites),
        "sources": sources,
        "matching_s": round(match_s, 2),
        "events": len(df_ev),
        "plots_s": round(time.perf_counter() - t, 2),
        "total_s": round(time.perf_counter() - start, 2),
        "peak_rss_mb": round(gis.Run.peak_rss_mb()),
    }
    (out_dir / f"rp_bounds_performance{sfx}.json").write_text(
        json.dumps(perf, indent=1)
    )
    log.info("Timings: %s", json.dumps(perf))
