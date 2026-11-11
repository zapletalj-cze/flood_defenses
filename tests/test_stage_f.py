# Author: Jakub Zapletal
# Date: 2026-11-05
# Checks of Stage F (F.7): Gumbel fit, time convention, EFAS aggregation, matching at a confluence, static check; a
# synthetic country T6 with GloFAS (daily) and EFAS (6-hourly) served from a local EWDS folder: resume of the
# download, conservative combination of the sources, bounds and plots.
# Run: pixi run python -m tests.test_stage_f  (DEFENSES_TEST_DIR sets the working folder, default: temporary)

import ast
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from pyproj import Transformer
from osgeo import gdal

import gis

gdal.UseExceptions()
repo = Path(__file__).resolve().parents[1]
root = Path(
    os.environ.get("DEFENSES_TEST_DIR") or tempfile.mkdtemp(prefix="defenses_stage_f_")
)
data, out = root / "data_f", root / "output" / "T6"
for folder in (data, out):
    folder.mkdir(parents=True, exist_ok=True)
print(f"Test folder: {root}")
list_res, timings = [], {}
EPSG = 32633  # synthetic test data only

# --- Unit checks ---
sample = np.random.default_rng(3).gumbel(1000.0, 300.0, 5000)
xi, alpha = gis.Stats.gumbel_lmoments(sample)
rp_median = float(gis.Stats.gumbel_rp(np.median(sample), xi, alpha))
list_res.append(
    (
        f"Gumbel by L-moments: xi {xi:.1f} (1000), alpha {alpha:.1f} (300), RP of the median {rp_median:.2f}",
        abs(xi - 1000) < 15 and abs(alpha - 300) < 12 and abs(rp_median - 2.0) < 0.1,
    )
)
ts = [
    pd.Timestamp(v, tz="UTC").timestamp()
    for v in (
        "2020-01-02 00:00",
        "2020-01-01 06:00",
        "2020-01-01 18:00",
        "2020-01-02 00:00",
        "2020-01-02 06:00",
    )
]
glofas_day = gis.Discharge.period_days(ts[:1], "end")
efas_days = gis.Discharge.period_days(ts[1:], "end")
list_res.append(
    (
        f"time convention (end of period): GloFAS 24 h stamped 2020-01-02 00 UTC -> {glofas_day[0]}; EFAS 6 h stamped "
        f"06, 18, 24 h and 30 h -> {[str(d) for d in efas_days]}",
        str(glofas_day[0]) == "2020-01-01"
        and [str(d) for d in efas_days]
        == ["2020-01-01", "2020-01-01", "2020-01-01", "2020-01-02"]
        and str(gis.Discharge.period_days(ts[:1], "start")[0]) == "2020-01-02",
    )
)
d4 = np.array(["2020-01-01"] * 4, "datetime64[D]")
v4 = np.array([[10.0], [14.0], [12.0], [20.0]])
mean = gis.Discharge.aggregate_daily(d4, v4, "daily_mean", 4)[1][0, 0]
peak = gis.Discharge.aggregate_daily(d4, v4, "daily_max", 4)[1][0, 0]
short = gis.Discharge.aggregate_daily(d4[:3], v4[:3], "daily_mean", 4)[1][0, 0]
list_res.append(
    (
        f"EFAS aggregation: daily_mean {mean} (mean 14), daily_max {peak} (20), a day with 3 of 4 steps {short}",
        mean == 14.0 and peak == 20.0 and np.isnan(short),
    )
)
# Confluence: the source lies on the main river just upstream of a tributary of similar upstream area; the tributary
# cell is closer and has the better ratio, the main river cell is on the same river.
main_line = shapely.LineString([(0.0, 0.0), (6000.0, 0.0)])
trib_line = shapely.LineString([(5000.0, 6000.0), (5000.0, 0.0)])
cells = np.array(
    [[1500.0, 100.0], [5000.0, 1400.0], [3000.0, 7000.0]]
)  # main, tributary, beyond the radius
area = np.array([4700.0, 4790.0, 4800.0])
same = shapely.distance(main_line, shapely.points(cells)) <= 500.0
lon0, lat0 = 15.0, 49.0
deg = 1 / 111320.0
lon = lon0 + cells[:, 0] * deg / np.cos(np.radians(lat0))
lat = lat0 + cells[:, 1] * deg
src_lon, src_lat = lon0 + 4000.0 * deg / np.cos(np.radians(lat0)), lat0
m_pref = gis.Discharge.match(
    src_lon, src_lat, 4800.0, lon, lat, area, 5.0, 0.7, 1.3, same
)
m_none = gis.Discharge.match(
    src_lon, src_lat, 4800.0, lon, lat, area, 5.0, 0.7, 1.3, None
)
list_res.append(
    (
        f"matching at a confluence: the main river cell is chosen ({m_pref['quality']}, ratio {m_pref['ratio']:.3f}, "
        f"{m_pref['distance_km']:.2f} km); without the river preference the tributary cell is chosen (index "
        f"{m_none['index']}, {m_none['quality']}); the cell with ratio 1.0 at 7.1 km lies beyond the 5 km radius",
        m_pref["index"] == 0
        and m_pref["quality"] == "good"
        and m_none["index"] == 1
        and m_none["quality"] == "weak",
    )
)
# Configuration split: defaults in advanced.toml, user values on top (tables of the schema type "table" replaced
# whole), country defaults for every country, the country of the run taken from [general] country.
cdir = root / "config_split"
cdir.mkdir(exist_ok=True)
(cdir / "advanced.toml").write_text(
    'placeholders = ["general.log_level"]\n[general]\nlog_level = "INFO"\nn_workers = 4\n'
    "[country_defaults]\ncorridor_buffer_m = 3000.0\n[country_defaults.rp]\nmax_rp = 200\nstorey_height_m = 3.0\n"
    '[discharge.glofas]\nrequest = { a = ["1"], b = ["2"] }\n'
)
(cdir / "defenses.toml").write_text(
    '[general]\ncountry = "ZZ"\nn_workers = 8\n[countries.ZZ]\ncorridor_buffer_m = 500.0\n'
    "[countries.ZZ.rp]\nmax_rp = 100\n[countries.ZZ.compartments]\nmin_area_m2 = 5.0\n"
    '[discharge.glofas]\nrequest = { c = ["3"] }\n'
)
os.environ.pop("DEFENSES_COUNTRY", None)
run_split = gis.Run.env_overrides(CONFIG_PATH=str(cdir / "defenses.toml"), COUNTRY="XX")
cfg_split = gis.Config.load(
    run_split["CONFIG_PATH"],
    run_split["COUNTRY"],
    [
        "general.n_workers",
        "general.log_level",
        "country.corridor_buffer_m",
        "country.rp.max_rp",
        "country.rp.storey_height_m",
        "country.protected_areas.min_area_m2",
        "discharge.glofas.request",
    ],
)
list_res.append(
    (
        f"configuration split (old table compartments read as protected_areas): country {run_split['COUNTRY']}, n_workers {cfg_split['general']['n_workers']}, "
        f"buffer {cfg_split['country']['corridor_buffer_m']}, rp {cfg_split['country']['rp']}, request "
        f"{cfg_split['discharge']['glofas']['request']}",
        run_split["COUNTRY"] == "ZZ"
        and cfg_split["general"]["n_workers"] == 8
        and cfg_split["general"]["log_level"] == "INFO"
        and cfg_split["country"]["corridor_buffer_m"] == 500.0
        and cfg_split["country"]["rp"] == {"max_rp": 100, "storey_height_m": 3.0}
        and cfg_split["discharge"]["glofas"]["request"] == {"c": ["3"]}
        and cfg_split["placeholders"] == ["general.log_level"]
        and cfg_split["country"]["protected_areas"]["min_area_m2"] == 5.0,
    )
)
list_bad = []
for path in sorted(repo.rglob("*.py")):
    if any(part.startswith(".") or part in ("__pycache__",) for part in path.parts):
        continue
    tree = ast.parse(path.read_text())
    for node in ast.walk(tree):
        if isinstance(node, ast.Import) and any(
            a.name.split(".")[0] == "rasterio" for a in node.names
        ):
            list_bad.append(f"{path.name}: import rasterio")
        if (
            isinstance(node, ast.ImportFrom)
            and (node.module or "").split(".")[0] == "rasterio"
        ):
            list_bad.append(f"{path.name}: from rasterio")
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr in ("read_file", "to_file")
        ):
            if any(k.arg == "layer" for k in node.keywords):
                list_bad.append(
                    f"{path.name}:{node.lineno} {node.func.attr}(layer=...)"
                )
list_res.append(
    (
        "static check: no rasterio import, no layer argument in read_file / to_file "
        + str(list_bad),
        not list_bad,
    )
)

# --- Synthetic country T6: footprint, rivers, protected areas, GFM tests ---
X0, Y0 = 500000.0, 5500000.0
crs = f"EPSG:{EPSG}"
gis.Vector.write_gpkg(
    gpd.GeoDataFrame(
        {"path": ["dtm"]},
        geometry=[shapely.box(X0, Y0, X0 + 20000.0, Y0 + 20000.0)],
        crs=crs,
    ),
    out / "dtm_footprint.gpkg",
)
river = shapely.LineString([(X0, Y0 + 10000.0), (X0 + 20000.0, Y0 + 10000.0)])
gis.Vector.write_gpkg(
    gpd.GeoDataFrame(
        {
            "river_id": [0],
            "name": ["main"],
            "is_source": [True],
            "upstream_km2": [5000.0],
        },
        geometry=[river],
        crs=crs,
    ),
    data / "rivers.gpkg",
)
src_pts = [(X0 + 8000.0, Y0 + 10000.0), (X0 + 12000.0, Y0 + 10000.0)]
gis.Vector.write_gpkg(
    gpd.GeoDataFrame(
        {
            "protected_area_id": ["C1", "C2"],
            "source_x": [p[0] for p in src_pts],
            "source_y": [p[1] for p in src_pts],
            "rp": [50.0, 50.0],
            "rp_low": [20.0, 20.0],
            "rp_high": [100.0, 100.0],
            "rp_source": ["model", "model"],
        },
        geometry=[
            shapely.box(p[0] - 500, p[1] + 200, p[0] + 500, p[1] + 1200)
            for p in src_pts
        ],
        crs=crs,
    ),
    out / "protected_areas_rp.gpkg",
)
# Discharge: annual maxima with Gumbel(1000, 300) quantiles, EFAS 30 % higher on event days; within a day the four
# EFAS steps are q - 3, q - 1, q + 1, q + 3 (mean q, max q + 3).
years = list(range(2001, 2021))
days = np.arange(np.datetime64("2001-01-01"), np.datetime64("2021-01-01"))
yr = days.astype("datetime64[Y]").astype(int) + 1970
doy = (days - days.astype("datetime64[Y]")).astype(int)
pos = (np.arange(1, 21) - 0.44) / 20.12
planted = dict(
    zip(
        years,
        np.random.default_rng(9).permutation(1000.0 - 300.0 * np.log(-np.log(pos))),
    )
)
q = 300.0 + 100.0 * np.sin(2 * np.pi * doy / 365.25)
q[doy == 100] = [planted[y] for y in yr[doy == 100]]
events = {"2019-06-10": ("held", 1400.0), "2020-03-05": ("flooded", 2000.0)}
for day, (_, value) in events.items():
    q[days == np.datetime64(day)] = value
pd.DataFrame(
    {
        "protected_area_id": ["C1", "C1", "C2"],
        "unit_id": ["000_000"] * 3,
        "item_id": ["S1", "S2", "S1"],
        "datetime": [
            "2019-06-10T10:00:00Z",
            "2020-03-05T10:00:00Z",
            "2019-06-10T10:00:00Z",
        ],
        "inside_valid_px": [1000] * 3,
        "inside_wet_frac": [0.0, 0.9, 0.0],
        "comp_valid_px": [1000] * 3,
        "comp_wet_frac": [0.9] * 3,
        "outcome": ["held", "flooded", "held"],
    }
).to_csv(out / "gfm_evidence.csv", index=False)


def write_nc(path, values, lat_c, lon_c, time_hours=None, name="dis"):
    ds = gdal.GetDriverByName("netCDF").CreateMultiDimensional(str(path))
    group = ds.GetRootGroup()
    f64 = gdal.ExtendedDataType.Create(gdal.GDT_Float64)
    text = gdal.ExtendedDataType.CreateString()
    axes = [
        ("lat", lat_c, "degrees_north", gdal.DIM_TYPE_HORIZONTAL_Y),
        ("lon", lon_c, "degrees_east", gdal.DIM_TYPE_HORIZONTAL_X),
    ]
    if time_hours is not None:
        axes.insert(
            0,
            (
                "time",
                time_hours,
                "hours since 1970-01-01 00:00:00",
                gdal.DIM_TYPE_TEMPORAL,
            ),
        )
    dims = []
    for nm, coord, units, kind in axes:
        dim = group.CreateDimension(nm, kind, None, len(coord))
        arr = group.CreateMDArray(nm, [dim], f64)
        arr.Write(np.asarray(coord, np.float64))
        arr.CreateAttribute("units", [], text).Write(units)
        dims.append(dim)
    var = group.CreateMDArray(
        name, dims, gdal.ExtendedDataType.Create(gdal.GDT_Float32)
    )
    var.Write(np.ascontiguousarray(values, np.float32))
    ds = None


ewds = data / "ewds"
centre, grid = {}, {}
for src, step, res_deg in (("glofas", 24, 0.05), ("efas", 6, 1 / 60)):
    folder = ewds / src
    folder.mkdir(parents=True, exist_ok=True)
    lon_c = np.arange(14.9, 15.4, res_deg) + res_deg / 2
    lat_c = np.arange(49.75, 49.55, -res_deg) - res_deg / 2

    tr = Transformer.from_crs(EPSG, 4326, always_xy=True)
    s_lon, s_lat = tr.transform(*src_pts[0])
    c0, r0 = int(np.argmin(np.abs(lon_c - s_lon))), int(
        np.argmin(np.abs(lat_c - s_lat))
    )
    centre[src] = Transformer.from_crs(4326, EPSG, always_xy=True).transform(
        lon_c[c0], lat_c[r0]
    )
    grid[src] = (lat_c, lon_c, r0, c0)
    up = np.full((len(lat_c), len(lon_c)), 1e7)
    up[r0, c0] = 5e9
    write_nc(folder / "upstream_area_src.nc", up, lat_c, lon_c, name="upArea")
    base = np.full((len(lat_c), len(lon_c)), 5.0, np.float32)
    if src == "glofas":
        for y in years:
            sel = yr == y
            vals = np.repeat(base[None], sel.sum(), axis=0)
            vals[:, r0, c0] = q[sel]
            hours = ((days[sel] - np.datetime64("1970-01-01")).astype(int) + 1) * 24.0
            write_nc(folder / f"{y}.nc", vals, lat_c, lon_c, hours)
    else:
        offsets = np.array([-3.0, -1.0, 1.0, 3.0])
        scale = np.where(
            np.isin(days, np.array(list(events), "datetime64[D]")), 1.3, 1.0
        )
        for y in years:
            for m in range(1, 13):
                sel = (yr == y) & (
                    days.astype("datetime64[M]") == np.datetime64(f"{y}-{m:02d}")
                )
                n = int(sel.sum())
                vals = np.repeat(base[None], n * 4, axis=0)
                vals[:, r0, c0] = (
                    (q[sel] * scale[sel])[:, None] + offsets[None, :]
                ).ravel()
                hours = (
                    (days[sel] - np.datetime64("1970-01-01")).astype(int)[:, None] * 24
                    + np.array([6, 12, 18, 24])[None, :]
                ).ravel()
                write_nc(
                    folder / f"{y}_{m:02d}.nc", vals, lat_c, lon_c, hours.astype(float)
                )


def write_config(
    name,
    aggregation="daily_mean",
    mode="conservative_both",
    sources='["efas", "glofas"]',
    raw_dir="",
    min_years=20,
    output="output",
):
    def source(src, unit, step, radius, first):
        return f"""
[discharge.{src}]
dataset_id = "{'cems-glofas-historical' if src == 'glofas' else 'efas-historical'}"
request = {{ product_type = ["consolidated"] }}
data_format = "netcdf"
download_format = ""
first_year = {first}
last_year = 2020
request_unit = "{unit}"
time_keys = {{ year = "year", month = "month", day = "day" }}
area_key = "area"
step_hours = {step}
valid_time_marks = "end"
netcdf_variable = "dis"
upstream_area_file = "{(ewds / src / 'upstream_area_src.nc').as_posix()}"
upstream_area_variable = "upArea"
upstream_area_to_km2 = 1e-6
min_upstream_area_km2 = 100.0
search_radius_km = {radius}
upstream_ratio_min = 0.7
upstream_ratio_max = 1.3
match_river_distance_m = 500.0
"""

    (root / name).write_text(
        f"""[general]
work_dir = "{(root / 'work').as_posix()}"
output_dir = "{(root / output).as_posix()}"
log_level = "INFO"

[countries.T6]
working_resolution_m = 5.0
rivers = "{(data / 'rivers.gpkg').as_posix()}"
river_upstream_area_field = "upstream_km2"
upstream_area_raster = ""
[countries.T6.rp]
classes = [10, 20, 50, 100, 200, 500]

[discharge]
sources = {sources}
source_mode = "{mode}"
ewds_api_url = "file://{ewds.as_posix()}"
raw_dir = "{raw_dir}"
aggregation = "{aggregation}"
hydrological_year_start_month = 1
min_years = {min_years}
max_rp_extrapolation = 1000
flooded_window_days = 0
max_parallel_requests = 4
bbox_buffer_deg = 0.1
rp_disagreement_ratio = 3.0
"""
        + source("glofas", "year", 24, 5.0, 2001)
        + source("efas", "month", 6, 2.0, 2001)
        + """
[implied.plots]
enabled = true
dpi = 70
max_per_type = 50
sample_seed = 42
show_max = 2
plot_all_in_units = []
"""
    )


def run(label, module, config="test_f.toml", **variables):
    env = dict(
        os.environ,
        DEFENSES_CONFIG_PATH=str(root / config),
        DEFENSES_COUNTRY="T6",
        **{f"DEFENSES_{k}": str(v) for k, v in variables.items()},
    )
    begin = time.perf_counter()
    proc = subprocess.run(
        [sys.executable, "-m", module],
        env=env,
        capture_output=True,
        text=True,
        cwd=repo,
    )
    timings[label] = time.perf_counter() - begin
    (root / f"f_{label}.log").write_text(proc.stdout + proc.stderr)
    if proc.returncode != 0:
        print(proc.stderr[-3000:])
    return proc


write_config("test_f.toml")
raw = root / "work" / "T6" / "discharge" / "raw"
proc = run("download_1", "download.download_discharge")
proc2 = run("download_2", "download.download_discharge")
(raw / "efas" / "2005_03.nc").unlink(missing_ok=True)
(raw / "glofas" / "2010.nc").write_bytes(b"interrupted")
proc3 = run("download_3", "download.download_discharge")
log_csv = (
    pd.read_csv(raw / "download_log.csv")
    if (raw / "download_log.csv").exists()
    else pd.DataFrame()
)
if len(log_csv) == 3 * 260:
    r1, r2, r3 = log_csv.iloc[:260], log_csv.iloc[260:520], log_csv.iloc[520:]
    list_res.append(
        (
            f"download: 260 requests (20 GloFAS years, 240 EFAS months); rerun skips {int((r2['status'] == 'complete, skipped').sum())}; "
            f"after one deleted and one truncated file, {int((r3['status'] == 'downloaded').sum())} downloaded again",
            (r1["status"] == "downloaded").all()
            and (r2["status"] == "complete, skipped").all()
            and set(r3.loc[r3["status"] == "downloaded", "file"])
            == {"2005_03.nc", "2010.nc"},
        )
    )
else:
    list_res.append(
        (f"download_discharge runs three times ({len(log_csv)} log rows)", False)
    )
proc = run("prepare", "processing.discharge_prepare")
if proc.returncode == 0:
    g, e = np.load(out / "discharge_glofas.npz"), np.load(out / "discharge_efas.npz")
    eday = int(
        np.flatnonzero(
            e["dates"].astype("datetime64[D]") == np.datetime64("2019-06-10")
        )[0]
    )
    gday = int(
        np.flatnonzero(
            g["dates"].astype("datetime64[D]") == np.datetime64("2019-06-10")
        )[0]
    )
    cells_g = gpd.read_file(out / "discharge_cells_glofas.gpkg")
    list_res.append(
        (
            f"discharge_prepare: GloFAS {len(g['cell_id'])} cell, EFAS {len(e['cell_id'])} cell, 20 complete years; "
            f"GloFAS on 2019-06-10 {g['values'][gday, 0]:.0f} (1400, its own day), EFAS daily mean {e['values'][eday, 0]:.0f} (1820)",
            len(g["cell_id"]) == 1
            and len(e["cell_id"]) == 1
            and int(cells_g["n_years"].iloc[0]) == 20
            and abs(g["values"][gday, 0] - 1400.0) < 1e-3
            and abs(e["values"][eday, 0] - 1820.0) < 1e-2,
        )
    )
    write_config("test_f_max.toml", aggregation="daily_max")
    proc_max = run(
        "prepare_max", "processing.discharge_prepare", config="test_f_max.toml"
    )
    if proc_max.returncode == 0:
        e_max = np.load(out / "discharge_efas.npz")
        list_res.append(
            (
                f"daily_max for EFAS: {e_max['values'][eday, 0]:.0f} (1823); GloFAS stays daily_mean",
                abs(e_max["values"][eday, 0] - 1823.0) < 1e-2
                and "daily_max is not available" in proc_max.stderr,
            )
        )
    run("prepare_again", "processing.discharge_prepare")
else:
    list_res.append(("discharge_prepare runs", False))
# A Stage E implied protected area on the same river with a flooded test (implied/base and river_points.gpkg)
imp = out / "implied" / "base"
imp.mkdir(parents=True, exist_ok=True)
IMP_ID = "implied_0_L_0001"
gis.Vector.write_gpkg(
    gpd.GeoDataFrame(
        {"implied_id": [IMP_ID], "river_id": [0]},
        geometry=[shapely.box(X0 + 8000.0, Y0 + 10200.0, X0 + 9000.0, Y0 + 11000.0)],
        crs=crs,
    ),
    imp / "implied_protected_areas.gpkg",
)
chain = np.arange(0.0, 20001.0, 100.0)
gis.Vector.write_gpkg(
    gpd.GeoDataFrame(
        {"river_id": np.zeros(len(chain), int), "chainage_m": chain},
        geometry=shapely.points(X0 + chain, np.full(len(chain), Y0 + 10000.0)),
        crs=crs,
    ),
    out / "river_points.gpkg",
)
pd.DataFrame(
    {
        "protected_area_id": [IMP_ID],
        "unit_id": ["000_000"],
        "item_id": ["S2"],
        "datetime": ["2020-03-05T10:00:00Z"],
        "outcome": ["flooded"],
    }
).to_csv(imp / "implied_gfm_evidence.csv", index=False)
proc = run("bounds", "processing.rp_bounds")
if proc.returncode == 0:
    ev = pd.read_csv(out / "event_rp.csv")
    comb = pd.read_csv(out / "event_rp_combined.csv").set_index(
        ["protected_area_id", "item_id"]
    )
    held = comb.loc[("C1", "S1")]
    flooded = comb.loc[("C1", "S2")]
    list_res.append(
        (
            f"conservative_both: held takes the lower RP ({held['rp_glofas']:.1f} GloFAS, {held['rp_efas']:.1f} EFAS -> "
            f"{held['rp_combined']:.1f}), flooded the higher ({flooded['rp_glofas']:.1f}, {flooded['rp_efas']:.1f} -> "
            f"{flooded['rp_combined']:.1f})",
            held["rp_combined"] == min(held["rp_glofas"], held["rp_efas"])
            and flooded["rp_combined"] == max(flooded["rp_glofas"], flooded["rp_efas"])
            and held["rp_glofas"] != held["rp_efas"],
        )
    )
    gdf_c = gpd.read_file(out / "protected_areas_rp.gpkg").set_index(
        "protected_area_id"
    )
    c1, c2 = gdf_c.loc["C1"], gdf_c.loc["C2"]
    # Same river: the cell centre lies within match_river_distance_m (500 m) of the river line.
    expect = {
        s: "good" if river.distance(shapely.Point(centre[s])) <= 500.0 else "weak"
        for s in centre
    }
    c1_ev = comb.loc["C1"]
    disagree = bool(
        (
            (
                c1_ev[["rp_glofas", "rp_efas"]].max(axis=1)
                / c1_ev[["rp_glofas", "rp_efas"]].min(axis=1)
            )
            > 3.0
        ).any()
    )
    list_res.append(
        (
            f"bounds: C1 obs {c1['obs_rp_low']:.1f}-{c1['obs_rp_high']:.1f} from {c1['obs_source']} (GloFAS "
            f"{c1['rp_obs_low_glofas']:.1f}-{c1['rp_obs_high_glofas']:.1f}, EFAS {c1['rp_obs_low_efas']:.1f}-"
            f"{c1['rp_obs_high_efas']:.1f}), source_disagreement {c1['source_disagreement']} ({disagree} expected, ratio 3), obs_conflict "
            f"{c1['obs_conflict']}; rp {c1['rp']:.0f} ({c1['rp_low']:.0f}-{c1['rp_high']:.0f}); matches GloFAS "
            f"{c1['glofas_match_quality']} ({expect['glofas']} expected), EFAS {c1['efas_match_quality']} ({expect['efas']})",
            c1["obs_source"] == "both"
            and c1["n_held"] == 1
            and c1["n_flooded"] == 1
            and c1["obs_rp_low"] == min(c1["rp_obs_low_glofas"], c1["rp_obs_low_efas"])
            and c1["obs_rp_high"]
            == max(c1["rp_obs_high_glofas"], c1["rp_obs_high_efas"])
            and c1["obs_rp_low"] < c1["obs_rp_high"]
            and not c1["obs_conflict"]
            and bool(c1["source_disagreement"]) == disagree
            and c1["glofas_match_quality"] == expect["glofas"]
            and c1["efas_match_quality"] == expect["efas"],
        )
    )
    list_res.append(
        (
            f"one source only: C2 has no EFAS cell within 2 km (EFAS {c2['efas_match_quality']}), obs "
            f"{c2['obs_rp_low']:.1f}-{c2['obs_rp_high']:.1f} from {c2['obs_source']}, n_held {c2['n_held']}, n_flooded {c2['n_flooded']}",
            c2["obs_source"] == "glofas"
            and c2["efas_match_quality"] == "none"
            and c2["obs_rp_low"] == c2["rp_obs_low_glofas"]
            and c2["n_held"] == 1
            and c2["n_flooded"] == 0,
        )
    )
    imp_row = comb.loc[(IMP_ID, "S2")] if (IMP_ID, "S2") in comb.index else None
    list_res.append(
        (
            "implied protected area (Stage E): its flooded test on S2 gets "
            + (
                f"GloFAS {imp_row['rp_glofas']:.1f}, EFAS {imp_row['rp_efas']:.1f} -> {imp_row['rp_combined']:.1f} "
                "(the higher) in event_rp_combined.csv"
                if imp_row is not None
                else "no row in event_rp_combined.csv"
            )
            + "; protected_areas_rp.gpkg keeps the explicit protected areas only",
            imp_row is not None
            and imp_row["rp_combined"] == max(imp_row["rp_glofas"], imp_row["rp_efas"])
            and IMP_ID not in set(gdf_c.index)
            and len(gdf_c) == 2,
        )
    )
    matches = gpd.read_file(out / "discharge_matches.gpkg")
    m_exp = matches[matches["kind"] == "explicit"]
    list_res.append(
        (
            f"discharge_matches.gpkg: {len(m_exp)} lines from explicit sites (C1 both, C2 GloFAS) and "
            f"{len(matches) - len(m_exp)} from the implied site to their matched cells",
            len(m_exp) == 3
            and sorted(zip(m_exp["protected_area_id"], m_exp["source"]))
            == [("C1", "efas"), ("C1", "glofas"), ("C2", "glofas")]
            and set(matches.loc[matches["kind"] == "implied", "source"])
            == {"efas", "glofas"},
        )
    )
    pdir = out / "plots" / "discharge" / "base"
    kinds = {
        k: len(list((pdir / k).glob("*.png")))
        for k in ("match_map", "fit", "event", "source_comparison", "summary")
    }
    pp = out / "plots" / "discharge" / "prepare"
    kinds_p = {k: len(list((pp / k).glob("*.png"))) for k in ("fit_glofas", "fit_efas")}
    list_res.append(
        (
            f"plots: {kinds}, prepare {kinds_p}; index.html lists them",
            all(v > 0 for v in kinds.values())
            and all(v > 0 for v in kinds_p.values())
            and "match_map" in (pdir / "index.html").read_text(),
        )
    )
else:
    list_res.append(("rp_bounds runs", False))

# --- Days spanning two files: EFAS steps stamped 00, 06, 12, 18 UTC of each requested day (the layout when the time
# keys of the form are valid times). With end-marks the 00 step closes the day before, so the last day of every
# month is completed by the first step of the next file; one extra file closes 2001-12-31. ---
lay = data / "layout"
(lay / "efas").mkdir(parents=True, exist_ok=True)
lat_e, lon_e, r_e, c_e = grid["efas"]
up = np.full((len(lat_e), len(lon_e)), 1e7)
up[r_e, c_e] = 5e9
write_nc(lay / "efas" / "upstream_area.nc", up, lat_e, lon_e, name="upArea")
d2001 = np.arange(np.datetime64("2001-01-01"), np.datetime64("2002-01-01"))
base_e = 500.0 + np.arange(len(d2001), dtype=float)  # daily mean of period day i
for m in list(range(1, 13)) + [13]:
    if m == 13:
        stamps_h = np.array(
            [
                (np.datetime64("2002-01-01") - np.datetime64("1970-01-01")).astype(int)
                * 24.0
            ]
        )
        name = "2002_01.nc"
    else:
        sel = d2001.astype("datetime64[M]") == np.datetime64(f"2001-{m:02d}")
        stamps_h = (
            (
                (d2001[sel] - np.datetime64("1970-01-01")).astype(int)[:, None] * 24
                + np.array([0, 6, 12, 18])[None, :]
            )
            .ravel()
            .astype(float)
        )
        name = f"2001_{m:02d}.nc"
    # step value: mean of its period day plus -3, -1, +1, +3 for the steps ending 06, 12, 18, 24 h
    period = (np.floor((stamps_h * 3600 - 1) / 86400)).astype(int) - int(
        (d2001[0] - np.datetime64("1970-01-01")).astype(int)
    )
    offset = np.array([3.0, -3.0, -1.0, 1.0])[(np.round(stamps_h) % 24 / 6).astype(int)]
    vals = np.repeat(
        np.full((len(lat_e), len(lon_e)), 5.0, np.float32)[None], len(stamps_h), axis=0
    )
    vals[:, r_e, c_e] = (
        np.where(period >= 0, base_e[np.clip(period, 0, len(base_e) - 1)], 400.0)
        + offset
    )
    write_nc(lay / "efas" / name, vals, lat_e, lon_e, stamps_h)
(root / "output_layout" / "T6").mkdir(parents=True, exist_ok=True)
(root / "output_layout" / "T6" / "dtm_footprint.gpkg").write_bytes(
    (out / "dtm_footprint.gpkg").read_bytes()
)
write_config(
    "test_f_layout.toml",
    sources='["efas"]',
    raw_dir=lay.as_posix(),
    min_years=1,
    output="output_layout",
)
proc = run(
    "prepare_layout", "processing.discharge_prepare", config="test_f_layout.toml"
)
if proc.returncode == 0:
    lz = np.load(root / "output_layout" / "T6" / "discharge_efas.npz")
    ld = lz["dates"].astype("datetime64[D]")
    lv = lz["values"][:, 0]
    in_2001 = (ld >= d2001[0]) & (ld <= d2001[-1])
    boundary = [
        np.datetime64(f"2001-{m:02d}-01") - np.timedelta64(1, "D") for m in range(2, 13)
    ] + [d2001[-1]]
    err = float(np.nanmax(np.abs(lv[in_2001] - base_e)))
    list_res.append(
        (
            f"days spanning two files: 2001 has {int(np.isfinite(lv[in_2001]).sum())} complete days (365), largest error "
            f"{err:.4f}; month ends {[float(lv[ld == b][0]) for b in boundary[:2]]}... equal their means; 2000-12-31 "
            f"(one step) is {lv[ld == np.datetime64('2000-12-31')]}; maximum of 2001 {float(lz['maxima'][list(lz['years']).index(2001), 0]):.0f} (864)",
            int(np.isfinite(lv[in_2001]).sum()) == 365
            and err < 1e-3
            and all(
                abs(float(lv[ld == b][0]) - base_e[(b - d2001[0]).astype(int)]) < 1e-3
                for b in boundary
            )
            and np.isnan(lv[ld == np.datetime64("2000-12-31")]).all()
            and abs(float(lz["maxima"][list(lz["years"]).index(2001), 0]) - 864.0)
            < 1e-3,
        )
    )
else:
    list_res.append(
        ("discharge_prepare runs on files whose days span two files", False)
    )

print("Timings (s): " + ", ".join(f"{k} {v:.1f}" for k, v in timings.items()))
for label, ok in list_res:
    print(f"{'PASS' if ok else 'FAIL'}  {label}")
sys.exit(0 if all(ok for _, ok in list_res) else 1)
