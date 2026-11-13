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

# --- Step 3: dry scene selection and event composites ---
df_scr = pd.DataFrame(
    {
        "unit_id": ["U1"] * 7 + ["U2"] * 2,
        "item_id": ["a", "b", "c", "d", "e", "f", "g", "h", "i"],
        "datetime": [
            f"2020-05-{d:02d}T05:00:00Z" for d in (8, 10, 12, 14, 16, 18, 20, 12, 14)
        ],
        "kept": [False, True, False, False, False, False, False, False, False],
    }
)
df_eu = pd.DataFrame(
    {
        "event_id": ["E1", "E1"],
        "unit_id": ["U1", "U2"],
        "peak_date": ["2020-05-11", "2020-05-11"],
        "window_start": ["2020-05-06", "2020-05-06"],
        "window_end": ["2020-06-01", "2020-06-01"],
    }
)
dry = gis.Evidence.select_dry_scenes(df_scr, df_eu, 3)
list_res.append(
    (
        f"select_dry_scenes: the first three dry days after the peak in the unit with a wet scene "
        f"({sorted(dry['item_id'])}); none before the peak, none in a unit without a wet scene",
        sorted(dry["item_id"]) == ["c", "d", "e"] and set(dry["unit_id"]) == {"U1"},
    )
)

from osgeo import gdal, osr  # noqa: E402
from pyproj import Transformer  # noqa: E402

# unit grid 100 x 100 cells of 20 m in EPSG:32633; scenes on a 20 m grid in EPSG:3035 (another CRS)
C = root / "comp"
unit_dir = C / "work" / "T8" / "units" / "U1"
unit_dir.mkdir(parents=True, exist_ok=True)
(C / "output" / "T8").mkdir(parents=True, exist_ok=True)
ux0, uy1 = 500000.0, 5502000.0
srs_u = osr.SpatialReference()
srs_u.ImportFromEPSG(EPSG)
hand = np.ones((100, 100), np.float32)
gis.Raster.write_array(
    unit_dir / "hand_undef.tif",
    hand,
    (ux0, 20.0, 0, uy1, 0, -20.0),
    srs_u.ExportToWkt(),
)
core = gpd.GeoDataFrame(
    {"unit_id": ["U1"]},
    geometry=[__import__("shapely").box(ux0, uy1 - 2000, ux0 + 2000, uy1)],
    crs=EPSG,
)
core.to_file(C / "output" / "T8" / "units_core.gpkg")
gpd.GeoDataFrame(
    {
        "event_id": ["T8_E1_20200511", "T8_E1_20210101"],
        "group_id": [1, 1],
        "peak_date": ["2020-05-11", "2021-01-01"],
        "window_start": ["2020-05-06", "2020-12-27"],
        "window_end": ["2020-06-01", "2021-01-22"],
    },
    geometry=[LineString([(ux0 - 500, uy1 - 1000), (ux0 + 2500, uy1 - 1000)])] * 2,
    crs=EPSG,
).to_file(C / "output" / "T8" / "events.gpkg")
to_s = Transformer.from_crs(EPSG, 3035, always_xy=True)
sx, sy = to_s.transform(
    [ux0, ux0 + 2000, ux0, ux0 + 2000], [uy1, uy1, uy1 - 2000, uy1 - 2000]
)
gx0, gy1 = np.floor(min(sx) / 20) * 20 - 200, np.ceil(max(sy) / 20) * 20 + 200
n_s = int((max(sx) - min(sx)) / 20) + 30
srs_s = osr.SpatialReference()
srs_s.ImportFromEPSG(3035)
cx, cy = np.meshgrid(
    gx0 + (np.arange(n_s) + 0.5) * 20, gy1 - (np.arange(n_s) + 0.5) * 20
)
back = Transformer.from_crs(3035, EPSG, always_xy=True)
lx, ly = back.transform(cx.ravel(), cy.ravel())
lx, ly = lx.reshape(cx.shape) - ux0, uy1 - ly.reshape(
    cy.shape
)  # local metres from the unit's upper left


def scene_files(name, wet, excluded=None, water=None):
    paths = {}
    for role, mask in (
        ("flood", wet),
        ("exclusion", excluded),
        ("reference_water", water),
    ):
        array = np.zeros(cx.shape, np.uint8) if mask is None else mask.astype(np.uint8)
        path = C / f"{name}_{role}.tif"
        gis.Raster.write_array(
            path, array, (gx0, 20.0, 0, gy1, 0, -20.0), srs_s.ExportToWkt(), nodata=255
        )
        paths[f"{role}_path"] = str(path)
    return paths


water = (lx > 1800) & (lx < 2000)  # reference water on the right edge
rows = [
    {
        "unit_id": "U1",
        "item_id": "S1",
        "datetime": "2020-05-08T05:00:00Z",
        "event_dry": False,
        **scene_files("s1", (ly > 200) & (ly < 600) & ~water, None, water),
    },
    {
        "unit_id": "U1",
        "item_id": "S2",
        "datetime": "2020-05-12T05:00:00Z",
        "event_dry": False,
        **scene_files(
            "s2", (ly > 400) & (ly < 1000) & (lx > 400) & ~water, lx < 400, water
        ),
    },
    {
        "unit_id": "U1",
        "item_id": "S3",
        "datetime": "2020-05-16T05:00:00Z",
        "event_dry": True,
        **scene_files("s3", np.zeros(cx.shape, bool), None, water),
    },
    {
        "unit_id": "U1",
        "item_id": "S0",
        "datetime": "2020-03-01T05:00:00Z",
        "event_dry": False,
        **scene_files("s0", ly < 2000, None, None),
    },  # outside every event window
]
pd.DataFrame(rows).to_csv(C / "output" / "T8" / "gfm_index.csv", index=False)
# an "official" maximum extent: the union of the wet scenes without the rows 800-1000
official = (((ly > 200) & (ly < 600)) | ((ly > 400) & (ly < 800) & (lx > 400))) & ~water
gis.Raster.write_array(
    C / "official.tif",
    official.astype(np.uint8),
    (gx0, 20.0, 0, gy1, 0, -20.0),
    srs_s.ExportToWkt(),
    nodata=255,
)
(C / "test.toml").write_text(
    config.replace(json.dumps(str(root / "work")), json.dumps(str(C / "work")))
    .replace(json.dumps(str(root / "output")), json.dumps(str(C / "output")))
    .replace(
        'log_level = "INFO"',
        'log_level = "INFO"\nn_workers = 1\ngdal_cachemax_mb = 128\nwarp_memory_limit_mb = 128\ngdal_threads_per_worker = 1',
    )
    .replace(
        'river_upstream_area_field = "uparea"',
        'river_upstream_area_field = "uparea"\ncorridor_buffer_m = 500.0',
    )
    + f"""
[evidence]
use_advisory_flags = false
max_extent_files = [{json.dumps(str(C / 'official.tif'))}]

[sources]
gfm_flood_values = [1]
gfm_exclusion_values = [1]
gfm_reference_water_values = [1, 2]
gfm_nodata_values = [255]
gfm_advisory_values = [1, 2, 3]
"""
)
proc = subprocess.run(
    [sys.executable, "-m", "processing.event_composites"],
    env=dict(
        os.environ, DEFENSES_CONFIG_PATH=str(C / "test.toml"), DEFENSES_COUNTRY="T8"
    ),
    capture_output=True,
    text=True,
    cwd=Path(__file__).resolve().parents[1],
)
(root / "event_composites.log").write_text(proc.stdout + proc.stderr)
if proc.returncode != 0:
    print(proc.stderr[-3000:])
list_res.append(("event_composites runs", proc.returncode == 0))
ev_dir = C / "work" / "T8" / "events" / "T8_E1_20200511"
if proc.returncode == 0 and (ev_dir / "U1_n_wet.tif").exists():
    a = {
        name: gdal.Open(str(ev_dir / f"U1_{name}.tif")).ReadAsArray()
        for name in ["max_extent", "n_valid", "n_wet", "first_wet", "last_wet"]
    }
    yy = (np.arange(100) + 0.5)[:, None] * 20 + np.zeros((1, 100))
    xx = (np.arange(100) + 0.5)[None, :] * 20 + np.zeros((100, 1))
    inner = (xx > 460) & (
        xx < 1740
    )  # away from the mask edges (nearest-neighbour warp of another CRS)
    band = lambda y0, y1: inner & (yy > y0 + 40) & (yy < y1 - 40)  # noqa: E731
    df_comp = pd.read_csv(C / "output" / "T8" / "event_composites.csv")
    list_res.append(
        (
            "composite on the unit grid from scenes in another CRS: three valid observations, two where a scene is "
            f"excluded, none on reference water (n_valid values {sorted(np.unique(a['n_valid']).tolist())})",
            (a["n_valid"][inner] == 3).all()
            and (a["n_valid"][(xx < 360) & (xx > 40)] == 2).all()
            and (a["n_valid"][(xx > 1840) & (xx < 1960)] == 0).all()
            and (a["max_extent"][(xx > 1840) & (xx < 1960)] == 255).all(),
        )
    )
    list_res.append(
        (
            "wet counts, first and last wet day: rows 200-400 m wet once on day 2, 400-600 m twice (days 2 and 6), "
            "600-1000 m once on day 6, the rest dry; the scene outside the window is not used",
            (a["n_wet"][band(200, 400)] == 1).all()
            and (a["first_wet"][band(200, 400)] == 2).all()
            and (a["n_wet"][band(400, 600)] == 2).all()
            and (a["first_wet"][band(400, 600)] == 2).all()
            and (a["last_wet"][band(400, 600)] == 6).all()
            and (a["first_wet"][band(600, 1000)] == 6).all()
            and (a["n_wet"][band(1000, 2000)] == 0).all()
            and (a["first_wet"][band(1000, 2000)] == -1).all()
            and (a["max_extent"][band(1000, 2000)] == 0).all()
            and (a["max_extent"][band(200, 1000)] == 1).all(),
        )
    )
    r = df_comp.iloc[0] if len(df_comp) else None
    list_res.append(
        (
            f"event_composites.csv: one composite ({len(df_comp)}); 3 scenes on 3 days, from "
            f"{r['first_scene_day'] if r is not None else None} to {r['last_scene_day'] if r is not None else None} days "
            "from the peak, two days after it; the event without scenes has none",
            len(df_comp) == 1
            and r["scenes"] == 3
            and r["scene_days"] == 3
            and r["first_scene_day"] == -3
            and r["last_scene_day"] == 5
            and r["scene_days_after_peak"] == 2
            and not (C / "work" / "T8" / "events" / "T8_E1_20210101").exists(),
        )
    )
    df_cmp = (
        pd.read_csv(C / "output" / "T8" / "max_extent_comparison.csv")
        if (C / "output" / "T8" / "max_extent_comparison.csv").exists()
        else pd.DataFrame()
    )
    expected_only = (
        ((yy > 800) & (yy < 1000) & (xx > 400) & (xx < 1800)).sum() * 400 / 1e6
    )
    list_res.append(
        (
            f"comparison with the official maximum extent: only in the composite "
            f"{df_cmp['only_composite_km2'].iloc[0] if len(df_cmp) else float('nan'):.3f} km2 (expected about "
            f"{expected_only:.3f}: the rows left out of the official file), only official "
            f"{df_cmp['only_official_km2'].iloc[0] if len(df_cmp) else float('nan'):.3f} km2",
            len(df_cmp) == 1
            and abs(df_cmp["only_composite_km2"].iloc[0] - expected_only) < 0.03
            and df_cmp["only_official_km2"].iloc[0] < 0.02
            and df_cmp["both_km2"].iloc[0] > 0.5,
        )
    )
    maps = C / "output" / "T8" / "plots" / "events" / "base" / "composite_maps"
    list_res.append(
        (
            "composite maps with index.html",
            maps.is_dir() and len(list(maps.glob("*.png"))) == 1,
        )
    )

for label, ok in list_res:
    print(f"{'PASS' if ok else 'FAIL'}  {label}")
sys.exit(0 if all(ok for _, ok in list_res) else 1)
