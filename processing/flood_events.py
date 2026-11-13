# Author: Jakub Zapletal
# Date: 2026-11-13
# Stage H step 2: catalogue of flood events from the prepared discharge series (GloFAS, EFAS). Every river group
# (river_groups.gpkg of protected_areas.py) is matched to a discharge cell at its largest reach; peaks of at least
# events.min_peak_rp start an event with a window from days_before before the peak to days_after after it;
# overlapping windows are merged. Output: events.gpkg (one feature per event with the lines of its group),
# event_groups.csv (the match of every group) and hydrograph plots in plots/events/<RUN_ID> with index.html.
# Run: pixi run flood-events

import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely

import gis
from plots import CATEGORICAL, INK, Plots

# --- Run variables (no command-line parameters; DEFENSES_<NAME> environment variables override them) ---
CONFIG_PATH = "config/defenses.toml"
COUNTRY = "XX"
RUN_ID = "base"  # plots in <output_dir>/<country>/plots/events/<RUN_ID>
SHOW_PLOTS = False  # show up to implied.plots.show_max figures per type while running (needs a display)

KEYS = [
    "general.output_dir",
    "general.log_level",
    "country.river_upstream_area_field",
    "events.min_peak_rp",
    "events.days_before",
    "events.days_after",
    "discharge.sources",
    "discharge.max_rp_extrapolation",
    "discharge.glofas.search_radius_km",
    "discharge.glofas.upstream_ratio_min",
    "discharge.glofas.upstream_ratio_max",
    "discharge.glofas.match_river_distance_m",
    "discharge.efas.search_radius_km",
    "discharge.efas.upstream_ratio_min",
    "discharge.efas.upstream_ratio_max",
    "discharge.efas.match_river_distance_m",
    "implied.plots.enabled",
    "implied.plots.dpi",
    "implied.plots.max_per_type",
    "implied.plots.sample_seed",
    "implied.plots.show_max",
    "implied.plots.plot_all_in_units",
]
EVENT_COLUMNS = [
    "event_id",
    "group_id",
    "source",
    "cell_id",
    "peak_date",
    "peak_discharge",
    "peak_rp",
    "extrapolated",
    "window_start",
    "window_end",
    "n_peaks",
    "match_quality",
    "upstream_area_km2",
    "geometry",
]

if __name__ == "__main__":
    run = gis.Run.env_overrides(
        CONFIG_PATH=CONFIG_PATH, COUNTRY=COUNTRY, RUN_ID=RUN_ID, SHOW_PLOTS=SHOW_PLOTS
    )
    cfg = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], KEYS)
    general, ev, dis = cfg["general"], cfg["events"], cfg["discharge"]
    out_dir = Path(general["output_dir"]) / run["COUNTRY"]
    log = gis.Run.logger("flood_events", out_dir / "logs", general["log_level"])
    gis.Config.warn_placeholders(cfg)
    start = time.perf_counter()
    field = cfg["country"]["river_upstream_area_field"]
    if not (out_dir / "river_groups.gpkg").exists():
        raise SystemExit(
            "river_groups.gpkg is missing; run protected_areas.py with river_upstream_area_field set."
        )
    # sources in the order of discharge.sources: the first one with a match is the source of a group
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
        raise SystemExit("No discharge source prepared: no events.")

    # --- River groups: lines, the largest reach and its mid point ---
    gdf_rg = gpd.read_file(out_dir / "river_groups.gpkg")
    wkt = gdf_rg.crs.to_wkt()
    gdf_rg = gdf_rg[gdf_rg["group_id"] > 0]
    n_minor = int(
        (
            gpd.read_file(out_dir / "river_groups.gpkg", columns=["group_id"])[
                "group_id"
            ]
            == 0
        ).sum()
    )
    area = pd.to_numeric(gdf_rg[field], errors="coerce")
    top = gdf_rg.loc[area.groupby(gdf_rg["group_id"]).idxmax()]
    groups = pd.DataFrame(
        {
            "group_id": top["group_id"].to_numpy(),
            "area_km2": area.loc[top.index].to_numpy(),
            "length_km": gdf_rg.groupby("group_id")["length_m"]
            .sum()
            .loc[top["group_id"]]
            .to_numpy()
            / 1000,
        }
    )
    mid = shapely.line_interpolate_point(
        [shapely.get_parts(g)[0] for g in top.geometry.values], 0.5, normalized=True
    )
    lonlat = shapely.get_coordinates(gpd.GeoSeries(mid, crs=wkt).to_crs(4326).values)
    groups["lon"], groups["lat"] = lonlat[:, 0], lonlat[:, 1]
    lines = gdf_rg.dissolve("group_id").geometry
    log.info(
        "River groups: %d groups (group 0 with %d minor reaches has no single river and gets no events).",
        len(groups),
        n_minor,
    )

    # --- Matching: per group the first source with a discharge cell; events of its series ---
    t = time.perf_counter()
    data = {}
    for src in sources:
        with np.load(out_dir / f"discharge_{src}.npz") as f:
            npz = {k: f[k] for k in f.files}
        cells = (
            gpd.read_file(out_dir / f"discharge_cells_{src}.gpkg")
            .set_index("cell_id")
            .loc[npz["cell_id"]]
            .reset_index()
        )
        data[src] = {
            "npz": npz,
            "cells": cells,
            "points": cells.geometry.to_crs(wkt).values,
            "dates": npz["dates"].astype("datetime64[D]"),
        }
    list_event, list_match, series = [], [], {}
    for g in groups.itertuples():
        row = {
            "group_id": g.group_id,
            "area_km2": g.area_km2,
            "source": "",
            "cell_id": -1,
            "match_quality": "",
            "events": 0,
        }
        for src in sources:
            d, c = data[src], dis[src]
            same = (
                shapely.distance(lines.loc[g.group_id], d["points"])
                <= c["match_river_distance_m"]
            )
            m = gis.Discharge.match(
                g.lon,
                g.lat,
                g.area_km2 if np.isfinite(g.area_km2) else -1.0,
                d["npz"]["lon"],
                d["npz"]["lat"],
                d["npz"]["upstream_area_km2"],
                c["search_radius_km"],
                c["upstream_ratio_min"],
                c["upstream_ratio_max"],
                same,
            )
            if m["index"] < 0:
                continue
            k = int(m["index"])
            xi, alpha = d["cells"].loc[k, ["gumbel_xi", "gumbel_alpha"]]
            values = d["npz"]["values"][:, k].astype(float)
            rp = gis.Stats.gumbel_rp(values, xi, alpha)
            found = gis.Discharge.detect_events(
                d["dates"],
                values,
                rp,
                ev["min_peak_rp"],
                ev["days_before"],
                ev["days_after"],
            )
            row.update(
                source=src,
                cell_id=int(d["npz"]["cell_id"][k]),
                match_quality=m["quality"],
                events=len(found),
            )
            series[g.group_id] = (src, k, values, rp)
            for e in found:
                list_event.append(
                    {
                        "event_id": f"{run['COUNTRY']}_E{g.group_id}_{str(e['peak_date']).replace('-', '')}",
                        "group_id": g.group_id,
                        "source": src,
                        "cell_id": row["cell_id"],
                        "peak_date": str(e["peak_date"]),
                        "peak_discharge": e["peak_value"],
                        "peak_rp": min(e["peak_rp"], dis["max_rp_extrapolation"]),
                        "extrapolated": bool(
                            e["peak_rp"] > dis["max_rp_extrapolation"]
                        ),
                        "window_start": str(e["start"]),
                        "window_end": str(e["end"]),
                        "n_peaks": e["n_peaks"],
                        "match_quality": m["quality"],
                        "upstream_area_km2": g.area_km2,
                        "geometry": lines.loc[g.group_id],
                    }
                )
            break
        list_match.append(row)
    df_match = pd.DataFrame(list_match)
    df_match.to_csv(out_dir / "event_groups.csv", index=False)
    gdf_ev = gpd.GeoDataFrame(
        pd.DataFrame(list_event, columns=EVENT_COLUMNS), geometry="geometry", crs=wkt
    )
    gis.Vector.write_gpkg(gdf_ev, out_dir / "events.gpkg")
    matched = df_match[df_match["source"] != ""]
    log.info(
        "Matching: %d of %d groups have a discharge cell (%s); %d groups without one get no events. %.1f s.",
        len(matched),
        len(df_match),
        ", ".join(f"{k} {v}" for k, v in matched["source"].value_counts().items())
        or "-",
        len(df_match) - len(matched),
        time.perf_counter() - t,
    )
    if len(gdf_ev):
        years = pd.to_datetime(gdf_ev["peak_date"]).dt.year
        span = (data[sources[0]]["dates"][-1] - data[sources[0]]["dates"][0]).astype(
            int
        ) / 365.25
        log.info(
            "events.gpkg: %d events in %d groups (peaks of at least RP %.1f, window -%d / +%d days), %.2f events per "
            "matched group and year; peak RP p50 %.1f, p95 %.1f, max %.0f; %d merged from several peaks; %d "
            "extrapolated beyond RP %.0f. Years with most events: %s.",
            len(gdf_ev),
            gdf_ev["group_id"].nunique(),
            ev["min_peak_rp"],
            ev["days_before"],
            ev["days_after"],
            len(gdf_ev) / max(len(matched), 1) / max(span, 1e-9),
            gdf_ev["peak_rp"].median(),
            gdf_ev["peak_rp"].quantile(0.95),
            gdf_ev["peak_rp"].max(),
            int((gdf_ev["n_peaks"] > 1).sum()),
            int(gdf_ev["extrapolated"].sum()),
            dis["max_rp_extrapolation"],
            ", ".join(f"{y} ({n})" for y, n in years.value_counts().head(5).items()),
        )
    else:
        log.warning(
            "events.gpkg: no events (no peak reaches RP %.1f).", ev["min_peak_rp"]
        )

    # --- Plots: hydrograph of a group with its events (largest peaks first), events per year ---
    t = time.perf_counter()
    settings = cfg["implied"]["plots"]
    if settings["enabled"] and len(gdf_ev):
        plots = Plots(
            out_dir / "plots" / "events" / run["RUN_ID"],
            settings,
            run["SHOW_PLOTS"],
            log,
        )
        plt = plots.plt
        top_rp = (
            gdf_ev.groupby("group_id")["peak_rp"].max().sort_values(ascending=False)
        )
        names = [f"group_{g}" for g in top_rp.index]
        for name in plots.select(
            "hydrographs", names, set(names[: settings["max_per_type"] // 2])
        ):
            gid = int(name.split("_")[1])
            src, k, values, rp = series[gid]
            dates = data[src]["dates"]
            own = gdf_ev[gdf_ev["group_id"] == gid]
            fig, ax = plt.subplots(figsize=(11, 3.6), layout="constrained")
            ax.plot(
                dates.astype("datetime64[ms]").astype(object),
                values,
                color=CATEGORICAL[0],
                linewidth=0.7,
            )
            xi, alpha = data[src]["cells"].loc[k, ["gumbel_xi", "gumbel_alpha"]]
            ax.axhline(
                gis.Stats.gumbel_q(ev["min_peak_rp"], xi, alpha),
                color=INK["muted"],
                linestyle=(0, (4, 2)),
                linewidth=0.8,
            )
            for e in own.itertuples():
                ax.axvspan(
                    pd.Timestamp(e.window_start),
                    pd.Timestamp(e.window_end),
                    color=CATEGORICAL[1],
                    alpha=0.18,
                    linewidth=0,
                )
                ax.plot(
                    [pd.Timestamp(e.peak_date)],
                    [e.peak_discharge],
                    "o",
                    color=CATEGORICAL[1],
                    markersize=4,
                )
            ax.set_ylabel("discharge [m3/s]")
            ax.set_title(
                f"group {gid}: {len(own)} events at {src} cell {int(own['cell_id'].iloc[0])} "
                f"(upstream {own['upstream_area_km2'].iloc[0]:.0f} km2); dashed = RP {ev['min_peak_rp']:g}"
            )
            plots.save(
                fig,
                "hydrographs",
                name,
                {
                    "events": len(own),
                    "max_peak_rp": float(own["peak_rp"].max()),
                    "source": src,
                },
            )
        plots.select("summary", ["events_per_year"])
        fig, ax = plt.subplots(figsize=(9, 3.4), layout="constrained")
        per_year = (
            pd.to_datetime(gdf_ev["peak_date"]).dt.year.value_counts().sort_index()
        )
        ax.bar(per_year.index, per_year.to_numpy(), color=CATEGORICAL[0])
        ax.set_ylabel("events")
        ax.set_title(
            f"Events per year ({len(gdf_ev)} events in {gdf_ev['group_id'].nunique()} river groups)"
        )
        plots.save(
            fig,
            "summary",
            "events_per_year",
            {"events": len(gdf_ev), "groups": int(gdf_ev["group_id"].nunique())},
        )
        log.info("Plots: %s (%.1f s)", plots.finish(), time.perf_counter() - t)
    log.info("Done in %.1f s.", time.perf_counter() - start)
