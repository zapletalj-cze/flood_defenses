# Author: Jakub Zapletal
# Date: 2026-11-13
# Stage H on synthetic data. Step 2, event catalogue (flood_events.py): river groups with discharge cells of two
# sources and series with known peaks: one peak, two overlapping peaks (merged), two separate events, a group below
# the threshold, a group without a discharge cell, a group matched only by the second source.
# Run: pixi run python -m tests.test_stage_h  (DEFENSES_TEST_DIR sets the working folder, default: temporary)

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import LineString

import gis

EPSG = 32633  # synthetic test data only
root = Path(
    os.environ.get("DEFENSES_TEST_DIR") or tempfile.mkdtemp(prefix="defenses_stage_h_")
)
out = root / "output" / "T8"
out.mkdir(parents=True, exist_ok=True)
print(f"Test folder: {root}")

# --- River groups: six groups of one 10 km reach each, 20 km apart, and a minor reach (group 0) ---
X0, Y0 = 500000.0, 5500000.0
reaches = {
    g: LineString([(X0, Y0 + 20000 * g), (X0 + 10000, Y0 + 20000 * g)])
    for g in range(7)
}
areas = {0: 40.0, 1: 5000.0, 2: 3000.0, 3: 2000.0, 4: 1500.0, 5: 1200.0, 6: 900.0}
gdf_rg = gpd.GeoDataFrame(
    {
        "group_id": list(reaches),
        "uparea": [areas[g] for g in reaches],
        "length_m": 10000.0,
    },
    geometry=list(reaches.values()),
    crs=EPSG,
)
gdf_rg.to_file(out / "river_groups.gpkg")
mid = (
    gdf_rg.set_index("group_id").geometry.interpolate(0.5, normalized=True).to_crs(4326)
)

# --- Discharge: daily series 2000-2023, base flow with a yearly cycle; Gumbel fit from the annual maxima ---
days = np.arange("2000-01-01", "2024-01-01", dtype="datetime64[D]")
doy = (days - days.astype("datetime64[Y]")).astype(int)
rng = np.random.default_rng(3)


def base(scale):
    return scale * (1.0 + 0.4 * np.sin(2 * np.pi * doy / 365.25)) + rng.normal(
        0, 0.02 * scale, len(days)
    )


def add_peak(series, date, height, width=3):
    k = int(np.flatnonzero(days == np.datetime64(date))[0])
    series[k - 3 * width : k + 3 * width + 1] += height * np.exp(
        -0.5 * (np.arange(-3 * width, 3 * width + 1) / width) ** 2
    )


series = {g: base(100.0) for g in range(1, 7)}
add_peak(series[1], "2010-06-10", 400.0)  # one peak
add_peak(
    series[2], "2013-05-14", 300.0
)  # two peaks 18 days apart: windows overlap, the second is higher
add_peak(series[2], "2013-06-01", 450.0)
add_peak(series[3], "2005-03-15", 350.0)  # two events three years apart
add_peak(series[3], "2008-08-02", 500.0)
add_peak(series[5], "2018-07-07", 420.0)  # group 5: only the second source has a cell
add_peak(
    series[6], "2020-10-10", 380.0
)  # group 6: no cell within the upstream area ratio
# group 4: no added peak; its ordinary yearly maxima exceed RP 5 about once in five years


def write_source(name, list_group, cell_area):
    values = np.column_stack([series[g] for g in list_group]).astype(np.float32)
    years, maxima = gis.Discharge.annual_maxima(days, values, 1)
    n_years, xi, alpha, rmse = gis.Discharge.gumbel_fits(maxima, 20)
    cell_id = np.array([100 + g for g in list_group])
    lon = np.array([mid.loc[g].x for g in list_group])
    lat = np.array([mid.loc[g].y for g in list_group])
    area = np.array([cell_area.get(g, areas[g]) for g in list_group])
    np.savez_compressed(
        out / f"discharge_{name}.npz",
        cell_id=cell_id,
        lon=lon,
        lat=lat,
        upstream_area_km2=area,
        dates=days.astype(np.int64),
        values=values,
        years=years,
        maxima=maxima,
    )
    gpd.GeoDataFrame(
        {
            "cell_id": cell_id,
            "upstream_area_km2": area,
            "n_years": n_years,
            "gumbel_xi": xi,
            "gumbel_alpha": alpha,
            "fit_rmse": rmse,
            "lon": lon,
            "lat": lat,
        },
        geometry=gpd.points_from_xy(lon, lat),
        crs=4326,
    ).to_file(out / f"discharge_cells_{name}.gpkg")
    return dict(zip(list_group, zip(xi, alpha)))


fit_efas = write_source(
    "efas", [1, 2, 3, 4, 6], {6: 5000.0}
)  # cell of group 6 has a wrong upstream area
fit_glofas = write_source("glofas", [1, 2, 3, 4, 5], {})

config = f"""advanced = ""
[general]
work_dir = {json.dumps(str(root / 'work'))}
output_dir = {json.dumps(str(root / 'output'))}
log_level = "INFO"

[countries.T8]
river_upstream_area_field = "uparea"

[events]
min_peak_rp = 5.0
days_before = 5
days_after = 21

[discharge]
sources = ["efas", "glofas"]
max_rp_extrapolation = 1000

[discharge.glofas]
search_radius_km = 5.0
upstream_ratio_min = 0.7
upstream_ratio_max = 1.3
match_river_distance_m = 3000.0

[discharge.efas]
search_radius_km = 5.0
upstream_ratio_min = 0.7
upstream_ratio_max = 1.3
match_river_distance_m = 3000.0

[implied.plots]
enabled = true
dpi = 80
max_per_type = 6
sample_seed = 1
show_max = 0
plot_all_in_units = []
"""
(root / "test.toml").write_text(config)

env = dict(
    os.environ, DEFENSES_CONFIG_PATH=str(root / "test.toml"), DEFENSES_COUNTRY="T8"
)
proc = subprocess.run(
    [sys.executable, "-m", "processing.flood_events"],
    env=env,
    capture_output=True,
    text=True,
    cwd=Path(__file__).resolve().parents[1],
)
(root / "flood_events.log").write_text(proc.stdout + proc.stderr)
if proc.returncode != 0:
    print(proc.stderr[-3000:])

list_res = [("flood_events runs", proc.returncode == 0)]
gdf_ev = (
    gpd.read_file(out / "events.gpkg")
    if proc.returncode == 0
    else gpd.GeoDataFrame(columns=["group_id"])
)
df_match = (
    pd.read_csv(out / "event_groups.csv") if proc.returncode == 0 else pd.DataFrame()
)


def events_of(group):
    return (
        gdf_ev[gdf_ev["group_id"] == group]
        .sort_values("peak_date")
        .reset_index(drop=True)
    )


if len(gdf_ev):
    e1, e2, e3, e5 = events_of(1), events_of(2), events_of(3), events_of(5)
    k = int(np.flatnonzero(days == np.datetime64("2010-06-10"))[0])
    rp_expected = float(
        gis.Stats.gumbel_rp(series[1][k - 2 : k + 3].max(), *fit_efas[1])
    )
    list_res.append(
        (
            f"one peak: one event on {e1['peak_date'].iloc[0] if len(e1) else None}, window "
            f"{e1['window_start'].iloc[0] if len(e1) else None} to {e1['window_end'].iloc[0] if len(e1) else None}, "
            f"peak RP {e1['peak_rp'].iloc[0] if len(e1) else float('nan'):.0f} = Gumbel RP of the peak discharge "
            f"({min(rp_expected, 1000):.0f})",
            len(e1) == 1
            and abs(
                (
                    np.datetime64(e1["peak_date"].iloc[0]) - np.datetime64("2010-06-10")
                ).astype(int)
            )
            <= 2
            and (
                np.datetime64(e1["peak_date"].iloc[0])
                - np.datetime64(e1["window_start"].iloc[0])
            ).astype(int)
            == 5
            and (
                np.datetime64(e1["window_end"].iloc[0])
                - np.datetime64(e1["peak_date"].iloc[0])
            ).astype(int)
            == 21
            and abs(e1["peak_rp"].iloc[0] - min(rp_expected, 1000))
            <= 1e-6 * min(rp_expected, 1000)
            and e1["source"].iloc[0] == "efas",
        )
    )
    list_res.append(
        (
            f"overlapping peaks: merged into one event ({len(e2)}) with the higher peak "
            f"({e2['peak_date'].iloc[0] if len(e2) else None}, {int(e2['n_peaks'].iloc[0]) if len(e2) else 0} peaks), "
            "window from the first peak's start to the last peak's end",
            len(e2) == 1
            and abs(
                (
                    np.datetime64(e2["peak_date"].iloc[0]) - np.datetime64("2013-06-01")
                ).astype(int)
            )
            <= 2
            and e2["n_peaks"].iloc[0] == 2
            and np.datetime64(e2["window_start"].iloc[0]) <= np.datetime64("2013-05-11")
            and np.datetime64(e2["window_end"].iloc[0]) >= np.datetime64("2013-06-20"),
        )
    )
    list_res.append(
        (
            f"separate peaks: two events ({list(e3['peak_date'])})",
            len(e3) == 2
            and e3["peak_date"].iloc[0].startswith("2005-03")
            and e3["peak_date"].iloc[1].startswith("2008-0"),
        )
    )
    e4 = events_of(4)
    list_res.append(
        (
            f"ordinary years: a series without a flood has {len(e4)} events in 24 years at RP 5 (about one in five "
            f"years, peak RP at most {e4['peak_rp'].max() if len(e4) else 0:.0f}); the minor group 0 gets none",
            2 <= len(e4) <= 9
            and e4["peak_rp"].min() >= 5.0
            and e4["peak_rp"].max() < 200
            and len(events_of(0)) == 0
            and 0 not in set(df_match["group_id"]),
        )
    )
    list_res.append(
        (
            f"source by the order of discharge.sources: group 5 has no efas cell and takes glofas "
            f"({list(e5['source'])}); group 6 has no cell within the upstream area ratio and no events "
            f"({len(events_of(6))})",
            len(e5) == 1
            and e5["source"].iloc[0] == "glofas"
            and len(events_of(6)) == 0
            and df_match.set_index("group_id").loc[6, "cell_id"] == -1,
        )
    )
    list_res.append(
        (
            "event ids are unique; geometry is the line of the group",
            gdf_ev["event_id"].is_unique
            and all(
                g.equals(reaches[int(i)]) or g.equals_exact(reaches[int(i)], 1e-6)
                for i, g in zip(gdf_ev["group_id"], gdf_ev.geometry)
            ),
        )
    )
index = out / "plots" / "events" / "base" / "index.html"
list_res.append(
    (
        "plots with index.html",
        index.exists() and (index.parent / "hydrographs").is_dir(),
    )
)

# detect_events alone: two peaks around a missing day are merged by their windows; a later peak is its own event
d = np.arange("2020-01-01", "2020-03-01", dtype="datetime64[D]")
v = np.full(len(d), 10.0)
v[10:13] = [50, np.nan, 60]
v[40] = 70
rp = np.where(np.nan_to_num(v) >= 50, 10.0, 1.0)
found = gis.Discharge.detect_events(d, v, rp, 5.0, 2, 3)
list_res.append(
    (
        f"detect_events: {len(found)} events, peaks {[str(e['peak_date']) for e in found]}, peaks merged in the "
        f"first {found[0]['n_peaks'] if found else 0}",
        len(found) == 2
        and str(found[0]["peak_date"]) == "2020-01-13"
        and found[0]["n_peaks"] == 2
        and str(found[1]["peak_date"]) == "2020-02-10",
    )
)

for label, ok in list_res:
    print(f"{'PASS' if ok else 'FAIL'}  {label}")
sys.exit(0 if all(ok for _, ok in list_res) else 1)
