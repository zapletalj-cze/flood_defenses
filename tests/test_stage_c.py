# Author: Jakub Zapletal
# Date: 2026-11-04
# End-to-end check of Stage C on the Stage B test country: STAC search with token and paging, GFM screening and
# clipping, evidence outcomes, missing-barrier hints, GloFAS NetCDF and GRIB reading, Gumbel fit, observed bounds,
# their combination with the model interval, and validation with an RP confusion matrix.
# Run: pixi run python -m tests.test_stage_c  (DEFENSES_TEST_DIR sets the working folder, default: temporary;
# the Stage B test runs first in the same folder when its outputs are missing)

import json
import os
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import urllib.parse
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from osgeo import gdal, osr
from pyproj import Transformer
from shapely import affinity
from shapely.geometry import Point, box

import gis

gdal.UseExceptions()

# The Stage B test country (tests/test_stage_b.py): synthetic test data only.
EPSG = 32633
GFM_EPSG = (
    3035  # the synthetic GFM scenes use another CRS than the DTM, like the real product
)
X0, Y_TOP = 500000.0, 5501600.0
WIDTH_M, HEIGHT_M = 2400, 1600
RIVER_Y = 800.0
probes = {"R1": (500, 1100), "R2": (1200, 1050), "R3": (2000, 1100)}
interiors = {
    "R1": (302, 952, 698, 1248),
    "R2": (1102, 952, 1498, 1248),
    "R3": (1802, 952, 2198, 1248),
}
TOKEN = "test-token"
GFM_RES = 20.0
WET_HAND_M = 3.8  # synthetic flood: undefended floodplain up to this HAND
POCKET = (780, 812, 1020, 885)  # dry in every scene, no levee: missing-barrier hint
CLASSES = [10, 20, 50, 100, 200, 500]
repo = Path(__file__).resolve().parents[1]

root = Path(
    os.environ.get("DEFENSES_TEST_DIR") or tempfile.mkdtemp(prefix="defenses_stage_c_")
)
out = root / "output" / "T4"
data = root / "data_c"
data.mkdir(parents=True, exist_ok=True)
print(f"Test folder: {root}")
if not (out / "compartments_rp.gpkg").exists():
    print("Stage B outputs missing: running tests.test_stage_b in the same folder.")
    proc = subprocess.run(
        [sys.executable, "-m", "tests.test_stage_b"],
        env=dict(os.environ, DEFENSES_TEST_DIR=str(root)),
        capture_output=True,
        text=True,
        cwd=repo,
    )
    (root / "test_stage_b.log").write_text(proc.stdout + proc.stderr)
    if proc.returncode != 0:
        print(proc.stdout[-3000:])
        sys.exit("The Stage B test failed; Stage C needs its outputs.")
# rp_bounds.py rewrites compartments_rp.gpkg: every run starts from the Stage B result.
stage_b = root / "compartments_rp_stage_b.gpkg"
if not stage_b.exists():
    shutil.copy(out / "compartments_rp.gpkg", stage_b)
shutil.copy(stage_b, out / "compartments_rp.gpkg")
model = gpd.read_file(stage_b)


def world(geometry):
    return affinity.translate(geometry, X0, Y_TOP - HEIGHT_M)


def covering(frame, label):
    return frame[frame.covers(world(Point(*probes[label])))]


ids = {k: covering(model, k).compartment_id.iloc[0] for k in probes}
before = {k: covering(model, k).iloc[0] for k in probes}

# --- GFM scenes: 20 m grid in EPSG:3035, flood extent from the undefended HAND of the Stage A units ---
to_gfm = Transformer.from_crs(EPSG, GFM_EPSG, always_xy=True)
from_gfm = Transformer.from_crs(GFM_EPSG, EPSG, always_xy=True)
corners = np.array(
    [
        to_gfm.transform(X0 + x, Y_TOP - HEIGHT_M + y)
        for x in (-200, WIDTH_M + 200)
        for y in (-200, HEIGHT_M + 200)
    ]
)
gx0, gy0 = np.floor(corners.min(axis=0) / GFM_RES) * GFM_RES
gx1, gy1 = np.ceil(corners.max(axis=0) / GFM_RES) * GFM_RES
gcols, grows = int(round((gx1 - gx0) / GFM_RES)), int(round((gy1 - gy0) / GFM_RES))
cx, cy = np.meshgrid(
    gx0 + (np.arange(gcols) + 0.5) * GFM_RES, gy1 - (np.arange(grows) + 0.5) * GFM_RES
)
ux, uy = from_gfm.transform(cx, cy)
lx, ly = ux - X0, uy - (Y_TOP - HEIGHT_M)
hand_tifs = sorted(
    str(p) for p in (root / "work" / "T4" / "units").glob("*/hand_undef.tif")
)
vrt = gdal.BuildVRT(str(data / "hand.vrt"), hand_tifs)
vrt = None
warped = gdal.Warp(
    "",
    str(data / "hand.vrt"),
    format="MEM",
    dstSRS=f"EPSG:{GFM_EPSG}",
    outputBounds=[gx0, gy0, gx1, gy1],
    xRes=GFM_RES,
    yRes=GFM_RES,
    resampleAlg="near",
    outputType=gdal.GDT_Float32,
    dstNodata=float("nan"),
)
hand = warped.GetRasterBand(1).ReadAsArray()
warped = None


def in_box(bounds, margin=0.0):
    return (
        (lx > bounds[0] - margin)
        & (lx < bounds[2] + margin)
        & (ly > bounds[1] - margin)
        & (ly < bounds[3] + margin)
    )


channel = np.abs(ly - RIVER_Y) <= 6
# Dry around the rings up to 30 m beyond the interior (levee body, its toe and mixed pixels): the bands lie below
# the scene water level and become hints unless hint_exclusion_m removes them.
ring_dry = in_box(interiors["R1"], 30) | in_box(interiors["R2"], 30)
ring_dry |= in_box(interiors["R3"], 30)
base = np.isfinite(hand) & (hand <= WET_HAND_M) & ~ring_dry & ~in_box(POCKET) & ~channel
nodata_block = in_box(
    (1260, 940, 1520, 1260)
)  # scene B: GFM nodata over most of ring 2
scenes = {
    # name: (datetime, flood, all excluded, with likelihood)
    "A": ("2016-11-15T05:12:00Z", base, False, True),
    "B": ("2021-06-10T05:12:00Z", base | in_box(interiors["R2"], 30), False, False),
    "C": ("2019-05-05T05:12:00Z", in_box((1500, 300, 1520, 320)), False, False),
    "D": ("2018-03-03T05:12:00Z", base, True, False),
    "E": ("2016-02-01T05:12:00Z", base & (lx < 1000), False, False),
}
gfm_srs = osr.SpatialReference()
gfm_srs.ImportFromEPSG(GFM_EPSG)
lonlat = Transformer.from_crs(EPSG, 4326, always_xy=True)
lon_w, lat_s = lonlat.transform(X0 - 300, Y_TOP - HEIGHT_M - 300)
lon_e, lat_n = lonlat.transform(X0 + WIDTH_M + 300, Y_TOP + 300)
items = []
for name, (stamp, flood, excluded, likelihood) in scenes.items():
    layers = {
        "ensemble_flood_extent": np.where(
            name == "B" and nodata_block, 255, flood.astype(np.uint8)
        ),
        "exclusion_mask": np.full(flood.shape, 1 if excluded else 0, np.uint8),
        "reference_water_mask": channel.astype(np.uint8),
    }
    if likelihood:
        layers["ensemble_likelihood"] = (flood * 90).astype(np.uint8)
    assets = {}
    for key, values in layers.items():
        path = data / f"GFM_{name}_{key}.tif"
        ds = gdal.GetDriverByName("GTiff").Create(
            str(path), gcols, grows, 1, gdal.GDT_Byte, ["COMPRESS=DEFLATE"]
        )
        ds.SetGeoTransform([gx0, GFM_RES, 0, gy1, 0, -GFM_RES])
        ds.SetProjection(gfm_srs.ExportToWkt())
        ds.GetRasterBand(1).WriteArray(values.astype(np.uint8))
        ds = None  # no band nodata: the flood layer's 255 must come from gfm_nodata_values
        assets[key] = {"href": f"file://{path}", "type": "image/tiff"}
    items.append(
        {
            "type": "Feature",
            "stac_version": "1.0.0",
            "id": f"GFM_TEST_{name}",
            "collection": "GFM_TEST",
            "bbox": [lon_w, lat_s, lon_e, lat_n],
            "geometry": shapely.geometry.mapping(box(lon_w, lat_s, lon_e, lat_n)),
            "properties": {"datetime": stamp},
            "assets": assets,
        }
    )
requests = []


class StacMock(BaseHTTPRequestHandler):
    # Item search with paging: page 1 links to page 2 by GET, page 2 to page 3 by POST (both occur in STAC APIs).
    def search(self, query, method):
        requests.append((method, int(query.get("page", 1))))
        if self.headers.get("Authorization") != f"Bearer {TOKEN}":
            self.send_response(401)
            self.end_headers()
            return
        start, end = query["datetime"].split("/")
        west, south, east, north = [float(v) for v in query["bbox"]]
        hits = sorted(
            (
                item
                for item in items
                if item["collection"] in query["collections"]
                and start <= item["properties"]["datetime"] <= end
                and item["bbox"][0] <= east
                and west <= item["bbox"][2]
                and item["bbox"][1] <= north
                and south <= item["bbox"][3]
            ),
            key=lambda item: item["id"],
        )
        page, limit = int(query.get("page", 1)), int(query["limit"])
        links = []
        if page * limit < len(hits):
            following = {**query, "page": page + 1}
            href = f"http://127.0.0.1:{self.server.server_address[1]}/stac/search"
            if page % 2:
                params = {
                    k: ",".join(map(str, v)) if isinstance(v, list) else v
                    for k, v in following.items()
                }
                links.append(
                    {"rel": "next", "href": f"{href}?{urllib.parse.urlencode(params)}"}
                )
            else:
                links.append(
                    {"rel": "next", "href": href, "method": "POST", "body": following}
                )
        body = json.dumps(
            {
                "type": "FeatureCollection",
                "features": hits[(page - 1) * limit : page * limit],
                "links": links,
                "numberMatched": len(hits),
            }
        ).encode()
        self.send_response(200)
        self.send_header("Content-Type", "application/geo+json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        url = urllib.parse.urlsplit(self.path)
        query = dict(urllib.parse.parse_qsl(url.query))
        query["collections"] = query["collections"].split(",")
        query["bbox"] = query["bbox"].split(",")
        self.search(query, "GET")

    def do_POST(self):
        length = int(self.headers.get("Content-Length", 0))
        self.search(json.loads(self.rfile.read(length)), "POST")

    def log_message(self, *args):
        pass


server = ThreadingHTTPServer(("127.0.0.1", 0), StacMock)
threading.Thread(target=server.serve_forever, daemon=True).start()
port = server.server_address[1]

# --- GloFAS: 0.05 degree grid around the country; the main river cell is the one closest to ring 2's source ---
lon_c, lat_c = lonlat.transform(X0 + WIDTH_M / 2, Y_TOP - HEIGHT_M / 2)
grid_west = np.floor(lon_c / 0.05) * 0.05 - 0.2
grid_north = np.ceil(lat_c / 0.05) * 0.05 + 0.15
ncols, nrows = 9, 7
lon_centers = grid_west + (np.arange(ncols) + 0.5) * 0.05
lat_centers = grid_north - (np.arange(nrows) + 0.5) * 0.05
glon, glat = np.meshgrid(lon_centers, lat_centers)
gx, gy = Transformer.from_crs(4326, EPSG, always_xy=True).transform(glon, glat)
source = before["R2"]
order = np.argsort(np.hypot(gx - source.source_x, gy - source.source_y), axis=None)
main_rc, decoy_rc = np.unravel_index(order[0], glon.shape), np.unravel_index(
    order[1], glon.shape
)
uparea = np.full((nrows, ncols), 1e7)  # m2 (10 km2)
uparea[(main_rc[0] + 1) % nrows, :] = 1e6  # 1 km2: below glofas_min_upstream_area_km2
uparea[main_rc] = 5e9  # 5000 km2 = upstream area of the test river
uparea[decoy_rc] = 4.6e9  # within the ratio range, but farther from ratio 1
days = np.arange(np.datetime64("1981-01-01"), np.datetime64("2022-01-01"))
years = days.astype("datetime64[Y]").astype(int) + 1970
day_of_year = (days - days.astype("datetime64[Y]")).astype(int)
# Annual maxima: Gumbel(1000, 300) quantiles at Gringorten plotting positions. The largest is replaced by the
# discharge of scene E, the 2016 maximum (RP above the model interval of ring 1); 2021 gets the second largest,
# above the days of scene B. Scene A (2016, held) gets the discharge of RP 70 of the resulting fit.
n_years = 41
positions = (np.arange(1, n_years + 1) - 0.44) / (n_years + 0.12)
planted = 1000.0 - 300.0 * np.log(-np.log(positions))
others = np.random.default_rng(11).permutation(planted[:-2])
peak = dict(zip([y for y in range(1981, 2022) if y not in (2016, 2021)], others))
peak[2021] = planted[-2]
q_main = 300.0 + 150.0 * np.sin(2 * np.pi * (day_of_year - 60) / 365.25)
peak_days = (day_of_year == 120) & (years != 2016)
q_main[peak_days] = [peak[y] for y in years[peak_days]]
q_main[days == np.datetime64("2016-02-01")] = 2481.0  # scene E
design = gis.Stats.gumbel_lmoments(
    pd.Series(q_main.astype(np.float32)).groupby(years).max().to_numpy()
)
special = {
    "2016-11-15": design[0] - design[1] * np.log(-np.log(1 - 1 / 70)),  # scene A
    "2021-06-08": 1300.0,
    "2021-06-09": 1450.0,  # scene B, flooded: maximum of the window, RP about 5
    "2021-06-10": 1200.0,  # scene B, held: discharge on the day
}
for day, value in special.items():
    q_main[days == np.datetime64(day)] = value
cell_index = np.arange(nrows * ncols).reshape(nrows, ncols)
discharge = np.empty((len(days), nrows, ncols), np.float32)
discharge[:] = (5.0 + cell_index % 7)[None, :, :]
discharge[:, main_rc[0], main_rc[1]] = q_main
discharge[:, decoy_rc[0], decoy_rc[1]] = 0.9 * q_main


def write_netcdf(path, variable, values, time_values=None):
    ds = gdal.GetDriverByName("netCDF").CreateMultiDimensional(str(path))
    group = ds.GetRootGroup()
    float64 = gdal.ExtendedDataType.Create(gdal.GDT_Float64)
    text = gdal.ExtendedDataType.CreateString()
    dimensions = []
    axes = [
        ("lat", lat_centers, "degrees_north", "latitude", gdal.DIM_TYPE_HORIZONTAL_Y)
    ]
    axes.append(
        ("lon", lon_centers, "degrees_east", "longitude", gdal.DIM_TYPE_HORIZONTAL_X)
    )
    if time_values is not None:
        axes.insert(
            0,
            (
                "time",
                time_values,
                "days since 1981-01-01",
                "time",
                gdal.DIM_TYPE_TEMPORAL,
            ),
        )
    for name, coordinates, units, standard_name, kind in axes:
        dimension = group.CreateDimension(name, kind, None, len(coordinates))
        coordinate = group.CreateMDArray(name, [dimension], float64)
        coordinate.Write(np.asarray(coordinates, np.float64))
        for key, value in [("units", units), ("standard_name", standard_name)]:
            coordinate.CreateAttribute(key, [], text).Write(value)
        dimensions.append(dimension)
    array = group.CreateMDArray(
        variable, dimensions, gdal.ExtendedDataType.Create(gdal.GDT_Float32)
    )
    array.Write(np.ascontiguousarray(values, np.float32))
    ds = None


# Daily values valid at the end of their day (time = day + 1 in "days since"), as glofas_valid_time_ends_day says.
elapsed = (days - np.datetime64("1981-01-01")).astype(int) + 1.0
split = int(np.searchsorted(days, np.datetime64("2001-01-01")))
write_netcdf(data / "glofas_1981_2000.nc", "dis24", discharge[:split], elapsed[:split])
write_netcdf(data / "glofas_2001_2021.nc", "dis24", discharge[split:], elapsed[split:])
write_netcdf(data / "uparea.nc", "uparea", uparea)
wgs84 = osr.SpatialReference()
wgs84.ImportFromEPSG(4326)
for year in (2019, 2020, 2021):
    selected = np.flatnonzero(years == year)
    mem = gdal.GetDriverByName("MEM").Create(
        "", ncols, nrows, len(selected), gdal.GDT_Float32
    )
    mem.SetGeoTransform([grid_west, 0.05, 0, grid_north, 0, -0.05])
    mem.SetProjection(wgs84.ExportToWkt())
    options = [
        "DISCIPLINE=1",
        "PDS_PDTN=0",
        "PDS_TEMPLATE_NUMBERS=0 7 2 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 255 0 0 0 0 0",
        "INPUT_UNIT=K",  # no unit conversion on writing
        "DATA_ENCODING=IEEE_FLOATING_POINT",
    ]
    for b, k in enumerate(selected):
        mem.GetRasterBand(b + 1).WriteArray(discharge[k])
        valid_time = pd.Timestamp(days[k]) + pd.Timedelta(days=1)
        options.append(
            f"BAND_{b + 1}_IDS=CENTER=98 SIGNF_REF_TIME=1 "
            f"REF_TIME={valid_time:%Y-%m-%dT%H:%M:%SZ} PROD_STATUS=0 TYPE=1"
        )
    gdal.GetDriverByName("GRIB").CreateCopy(
        str(data / f"glofas_{year}.grib"), mem, options=options
    )
    mem = None

# --- Reference protected areas for validate.py: ring interiors with reference RP, one area without compartment ---
reference = gpd.GeoDataFrame(
    {"rp_ref": [300.0, 20.0, 100.0, 50.0]},
    geometry=[world(box(*interiors[k])) for k in ["R1", "R2", "R3"]]
    + [world(box(1600, 1300, 1700, 1400))],
    crs=EPSG,
)
reference.to_file(data / "reference.gpkg", driver="GPKG")

# --- Configuration: the Stage B test configuration plus the Stage C keys ---
base_config = (root / "test.toml").read_text()
country_keys = f"""reference_protected_areas = "{data / 'reference.gpkg'}"
reference_rp_field = "rp_ref"
validation_extent = []
"""
source_keys = f"""gfm_stac_url = "http://127.0.0.1:{port}/stac"
gfm_collection = "GFM_TEST"
gfm_assets = {{ flood = "ensemble_flood_extent", exclusion = "exclusion_mask", reference_water = "reference_water_mask", likelihood = "ensemble_likelihood" }}
gfm_date_from = "2015-01-01"
gfm_date_to = "2022-12-31"
gfm_page_limit = 2
gfm_token_env = "DEFENSES_TEST_GFM_TOKEN"
gfm_screen_resolution_m = 100.0
gfm_download_threads = 2
gfm_flood_values = [1]
gfm_exclusion_values = [1]
gfm_reference_water_values = [1]
gfm_nodata_values = [255]
glofas_enabled = true
glofas_dataset_id = "cems-glofas-historical"
glofas_request = {{}}
glofas_time_keys = {{ year = "hyear", month = "hmonth", day = "hday" }}
glofas_area_key = "area"
glofas_first_year = 1981
glofas_last_year = 2021
glofas_file_extension = "nc"
glofas_files = ["{data / 'glofas_1981_2000.nc'}", "{data / 'glofas_2001_2021.nc'}"]
glofas_discharge_variable = "dis24"
glofas_upstream_area_file = "{data / 'uparea.nc'}"
glofas_upstream_area_url = ""
glofas_upstream_area_variable = "uparea"
glofas_upstream_area_to_km2 = 1e-6
glofas_min_upstream_area_km2 = 5.0
glofas_min_days_per_year = 330
glofas_valid_time_ends_day = true
"""
evidence = """
[evidence]
comparison_distance_m = 200.0
min_valid_pixels = 500
f_wet_comparison = 0.5
f_dry_inside = 0.1
f_wet_inside = 0.5
flooded_window_days = 3
glofas_search_radius_km = 5.0
upstream_ratio_min = 0.7
upstream_ratio_max = 1.3
river_level_percentile = 90.0
hint_margin_m = 0.5
hint_dry_share = 0.8
hint_min_scenes = 2
hint_exclusion_m = 50.0
"""
config = base_config.replace("[countries.T4]\n", "[countries.T4]\n" + country_keys, 1)
config = config.replace("[sources]\n", "[sources]\n" + source_keys, 1) + evidence
(root / "test_c.toml").write_text(config)
grib_files = ", ".join(f'"{data / f"glofas_{y}.grib"}"' for y in (2019, 2020, 2021))
grib_config = (
    config.replace(
        f'output_dir = "{root / "output"}"', f'output_dir = "{root / "output_grib"}"'
    )
    .replace(
        f"""glofas_files = ["{data / 'glofas_1981_2000.nc'}", "{data / 'glofas_2001_2021.nc'}"]""",
        f"glofas_files = [{grib_files}]",
    )
    .replace('glofas_discharge_variable = "dis24"', 'glofas_discharge_variable = ""')
    .replace('glofas_file_extension = "nc"', 'glofas_file_extension = "grib"')
)
(root / "test_c_grib.toml").write_text(grib_config)
(root / "output_grib" / "T4").mkdir(parents=True, exist_ok=True)
shutil.copy(
    out / "dtm_footprint.gpkg", root / "output_grib" / "T4" / "dtm_footprint.gpkg"
)

timings = {}
no_proxy = ",".join(
    v for v in [os.environ.get("NO_PROXY", ""), "127.0.0.1", "localhost"] if v
)


def run(module, label, config="test_c.toml", **variables):
    env = dict(
        os.environ,
        DEFENSES_CONFIG_PATH=str(root / config),
        DEFENSES_COUNTRY="T4",
        DEFENSES_TEST_GFM_TOKEN=TOKEN,
        NO_PROXY=no_proxy,
        no_proxy=no_proxy,
        **{f"DEFENSES_{k}": str(v) for k, v in variables.items()},
    )
    start = time.perf_counter()
    proc = subprocess.run(
        [sys.executable, "-m", module],
        env=env,
        capture_output=True,
        text=True,
        cwd=repo,
    )
    timings[label] = time.perf_counter() - start
    (root / f"{label}.log").write_text(proc.stdout + proc.stderr)
    if proc.returncode != 0:
        print(proc.stderr[-3000:])
    return proc


def outcome(table, scene, label):
    rows = table[
        (table.item_id == f"GFM_TEST_{scene}") & (table.compartment_id == ids[label])
    ]
    return set(rows.outcome)


results = []

# --- download_gfm.py ---
proc = run("download.download_gfm", "download_gfm", FORCE="1")
index = (
    pd.read_csv(out / "gfm_index.csv", dtype={"unit_id": str}).fillna("")
    if proc.returncode == 0
    else pd.DataFrame(columns=["unit_id", "item_id"])
)
units = pd.read_csv(out / "units.csv", dtype={"unit_id": str}).unit_id
results.append(
    (
        "download_gfm: every unit keeps scenes A, B, D, E; C (a few wet pixels) is screened out",
        proc.returncode == 0
        and all(
            set(index[index.unit_id == u].item_id) == {f"GFM_TEST_{s}" for s in "ABDE"}
            for u in units
        ),
    )
)
results.append(
    (
        f"STAC search with token, once for all units: two GET pages and one POST page ({len(requests)} requests)",
        sorted(requests) == [("GET", 1), ("GET", 2), ("POST", 3)],
    )
)
proc = run("download.download_gfm", "download_gfm_rerun")
results.append(
    (
        "download_gfm rerun: nothing screened twice, same index (checkpoint screening.csv)",
        proc.returncode == 0
        and "; 0 items to screen" in proc.stderr
        and pd.read_csv(out / "gfm_index.csv", dtype={"unit_id": str})
        .fillna("")
        .equals(index),
    )
)
first = index[(index.unit_id == "000_000") & (index.item_id == "GFM_TEST_A")]
copy_ok = False
if len(first):
    clipped = gdal.Open(first.flood_path.iloc[0])
    copy_srs = osr.SpatialReference(wkt=clipped.GetProjection())
    copy_ok = (
        copy_srs.GetAuthorityCode(None) == str(GFM_EPSG)
        and clipped.RasterXSize * clipped.RasterYSize > 0
        and Path(first.likelihood_path.iloc[0]).exists()
        and index[index.item_id == "GFM_TEST_B"].likelihood_path.eq("").all()
    )
results.append(
    (
        "clipped copies keep the GFM CRS; optional likelihood copied where the item has it",
        copy_ok,
    )
)

# --- gfm_evidence.py without hint exclusion, first, so that the regular hints stay in the output folder ---
(root / "test_c_no_exclusion.toml").write_text(
    config.replace("hint_exclusion_m = 50.0", "hint_exclusion_m = 0.0")
)
proc = run(
    "processing.gfm_evidence",
    "gfm_evidence_no_exclusion",
    config="test_c_no_exclusion.toml",
    HINTS="1",
    DATE_FROM="2016-06-01",
)
unexcluded = (
    gpd.read_file(out / "missing_barrier_hints.gpkg")
    if proc.returncode == 0
    else gpd.GeoDataFrame()
)
# --- gfm_evidence.py, run 1: hints, scenes from June 2016 on (A, B, D) ---
proc = run(
    "processing.gfm_evidence", "gfm_evidence_1", HINTS="1", DATE_FROM="2016-06-01"
)
run1_ok = proc.returncode == 0
table = pd.read_csv(out / "gfm_evidence.csv") if run1_ok else pd.DataFrame()
if len(table):
    results.append(
        (
            "scene A: all three rings held (dry inside, wet floodplain around)",
            all(outcome(table, "A", k) == {"held"} for k in probes),
        )
    )
    results.append(
        (
            "scene B: ring 2 flooded, rings 1 and 3 held",
            outcome(table, "B", "R2") == {"flooded"}
            and outcome(table, "B", "R1") == {"held"}
            and outcome(table, "B", "R3") == {"held"},
        )
    )
    r2_b = table[
        (table.item_id == "GFM_TEST_B") & (table.compartment_id == ids["R2"])
    ].iloc[0]
    results.append(
        (
            f"GFM nodata (255 from the config) excluded: ring 2 valid cells in B {r2_b.inside_valid_px} "
            f"of {int(before['R2'].area_m2 / 4)}",
            r2_b.inside_valid_px < 0.6 * before["R2"].area_m2 / 4,
        )
    )
    results.append(
        (
            "scene D (all excluded): no_test everywhere",
            all(outcome(table, "D", k) == {"no_test"} for k in probes),
        )
    )
    results.append(
        (
            "date filter: scene E (February 2016) not tested",
            not (table.item_id == "GFM_TEST_E").any(),
        )
    )
else:
    results.append(("gfm_evidence run 1 produced tests", False))
layer = gpd.read_file(out / "gfm_evidence.gpkg") if run1_ok else gpd.GeoDataFrame()
results.append(
    (
        "gfm_evidence.gpkg: one compartment polygon per test, for QGIS",
        len(layer) == len(table) > 0 and layer.crs.to_epsg() == EPSG,
    )
)
hints = (
    gpd.read_file(out / "missing_barrier_hints.gpkg") if run1_ok else gpd.GeoDataFrame()
)
pocket_area = (POCKET[2] - POCKET[0]) * (POCKET[3] - POCKET[1])
hint_area = float(hints.area.sum()) if len(hints) else 0.0
results.append(
    (
        f"missing-barrier hints: one area at the dry pocket ({len(hints)} found, "
        f"{hint_area:.0f} m2 of {pocket_area} m2)",
        len(hints) == 1
        and bool(hints.covers(world(Point(900, 848))).iloc[0])
        and 0.7 < hint_area / pocket_area < 1.3,
    )
)

results.append(
    (
        f"without hint_exclusion_m the dry bands along the levees become hints too ({len(unexcluded)} areas)",
        len(unexcluded) > len(hints),
    )
)

# --- download_glofas.py: NetCDF files ---
proc = run("download.download_glofas", "download_glofas_netcdf")
glofas = np.load(out / "glofas_discharge.npz") if proc.returncode == 0 else None
if glofas is not None:
    dates = glofas["dates"].astype("datetime64[D]")
    main = int(np.flatnonzero(np.isclose(glofas["upstream_area_km2"], 5000.0))[0])
    a_day = int(np.flatnonzero(dates == np.datetime64("2016-11-15"))[0])
    e_day = int(np.flatnonzero(dates == np.datetime64("2016-02-01"))[0])
    results.append(
        (
            f"GloFAS NetCDF: days 1981-01-01..2021-12-31 attributed to their own day ({len(dates)} days)",
            dates[0] == np.datetime64("1981-01-01")
            and dates[-1] == np.datetime64("2021-12-31")
            and len(dates) == len(days)
            and glofas["discharge"][a_day, main] == np.float32(special["2016-11-15"]),
        )
    )
    results.append(
        (
            "upstream areas below glofas_min_upstream_area_km2 dropped",
            glofas["upstream_area_km2"].min() >= 5.0
            and np.isclose(glofas["upstream_area_km2"], 4600.0).sum() == 1,
        )
    )
    series = pd.Series(glofas["discharge"][:, main].astype(float))
    maxima = series.groupby(dates.astype("datetime64[Y]").astype(int)).max().to_numpy()
    xi, alpha = gis.Stats.gumbel_lmoments(maxima)
    rp_median = float(gis.Stats.gumbel_rp(np.median(maxima), xi, alpha))
    rp_a = float(gis.Stats.gumbel_rp(glofas["discharge"][a_day, main], xi, alpha))
    rp_b_flooded = float(gis.Stats.gumbel_rp(1450.0, xi, alpha))
    rp_e = float(gis.Stats.gumbel_rp(glofas["discharge"][e_day, main], xi, alpha))
    results.append(
        (
            f"Gumbel fit of the GloFAS annual maxima: xi {xi:.0f}, alpha {alpha:.0f}; median annual maximum "
            f"has RP {rp_median:.2f}",
            abs(rp_median - 2.0) < 0.1,
        )
    )
else:
    results.append(("download_glofas (NetCDF) runs", False))
sample = np.random.default_rng(5).gumbel(1000.0, 300.0, 20000)
fit = gis.Stats.gumbel_lmoments(sample)
rp_sample = float(gis.Stats.gumbel_rp(np.median(sample), *fit))
rp_100 = float(gis.Stats.gumbel_rp(1000.0 - 300.0 * np.log(-np.log(0.99)), *fit))
results.append(
    (
        f"Gumbel sanity (20000 values): xi {fit[0]:.1f}, alpha {fit[1]:.1f}, RP(median) {rp_sample:.3f}, "
        f"RP(true 100-year value) {rp_100:.1f}",
        abs(fit[0] - 1000) < 15
        and abs(fit[1] - 300) < 10
        and abs(rp_sample - 2) < 0.05
        and 90 < rp_100 < 110,
    )
)

# --- download_glofas.py: GRIB files of 2019-2021 give the same series ---
proc = run(
    "download.download_glofas", "download_glofas_grib", config="test_c_grib.toml"
)
grib_ok = False
if proc.returncode == 0 and glofas is not None:
    grib = np.load(root / "output_grib" / "T4" / "glofas_discharge.npz")
    grib_dates = grib["dates"].astype("datetime64[D]")
    rows = np.searchsorted(dates, grib_dates)
    grib_ok = (
        grib_dates[0] == np.datetime64("2019-01-01")
        and grib_dates[-1] == np.datetime64("2021-12-31")
        and np.array_equal(grib["lon"], glofas["lon"])
        and np.array_equal(grib["discharge"], glofas["discharge"][rows])
    )
results.append(("GloFAS GRIB (GRIB_VALID_TIME) equals the NetCDF series", grib_ok))

# --- rp_bounds.py, run 1 and an identical rerun ---
proc = run("processing.rp_bounds", "rp_bounds_1")
bounded = (
    gpd.read_file(out / "compartments_rp.gpkg")
    if proc.returncode == 0
    else gpd.GeoDataFrame()
)
if len(bounded) and glofas is not None:
    after = {k: covering(bounded, k).iloc[0] for k in probes}
    results.append(
        (
            "matching: rings 1 and 2 matched to the main GloFAS cell (ratio 1), not the decoy (ratio 0.92)",
            all(
                after[k].glofas_cell == main
                and abs(after[k].glofas_area_ratio - 1.0) < 1e-9
                and after[k].glofas_distance_km < 5.0
                for k in ("R1", "R2")
            ),
        )
    )
    r1 = after["R1"]
    results.append(
        (
            f"ring 1 combined: held at RP {rp_a:.1f} (A) raises rp_low within the model interval "
            f"[{before['R1'].rp_low:.0f}, {before['R1'].rp_high:.0f}]",
            50 < rp_a < 100
            and r1.n_held == 2
            and r1.n_flooded == 0
            and abs(r1.obs_rp_low - rp_a) < 1e-6 * rp_a
            and np.isnan(r1.obs_rp_high)
            and r1.rp == 100
            and abs(r1.rp_low - rp_a) < 1e-6 * rp_a
            and r1.rp_high == 100
            and not r1.model_obs_conflict
            and r1.rp_source == f"{before['R1'].rp_source}+obs",
        )
    )
    r2 = after["R2"]
    results.append(
        (
            f"ring 2: held at RP {rp_a:.1f} but flooded at RP {rp_b_flooded:.1f} (window maximum) = obs_conflict, "
            "model values kept",
            r2.obs_conflict
            and abs(r2.obs_rp_high - rp_b_flooded) < 1e-6 * rp_b_flooded
            and (r2.rp, r2.rp_low, r2.rp_high, r2.rp_source)
            == (
                before["R2"].rp,
                before["R2"].rp_low,
                before["R2"].rp_high,
                before["R2"].rp_source,
            ),
        )
    )
    r3 = after["R3"]
    results.append(
        (
            "ring 3: override kept despite observations",
            r3.n_held == 2
            and r3.rp_source == "override"
            and r3.rp == r3.rp_low == r3.rp_high == 100,
        )
    )
    columns = ["rp", "rp_low", "rp_high", "rp_source", "obs_rp_low", "obs_rp_high"]
    proc = run("processing.rp_bounds", "rp_bounds_2")
    again = (
        gpd.read_file(out / "compartments_rp.gpkg") if proc.returncode == 0 else None
    )
    results.append(
        (
            "rp_bounds rerun gives identical results (model values kept in *_prior columns)",
            again is not None
            and again[columns].equals(bounded[columns])
            and (again.rp_prior.to_numpy() == model.rp.to_numpy()).all(),
        )
    )
else:
    results.append(("rp_bounds run 1", False))

# --- gfm_evidence.py run 2 (all dates) and rp_bounds.py: a held test above the model interval ---
proc = run("processing.gfm_evidence", "gfm_evidence_2")
table2 = (
    pd.read_csv(out / "gfm_evidence.csv") if proc.returncode == 0 else pd.DataFrame()
)
results.append(
    (
        "scene E (only the area around ring 1 wet): ring 1 held, rings 2 and 3 no_test",
        len(table2) > 0
        and outcome(table2, "E", "R1") == {"held"}
        and outcome(table2, "E", "R2") == {"no_test"}
        and outcome(table2, "E", "R3") == {"no_test"},
    )
)
proc = run("processing.rp_bounds", "rp_bounds_3")
final = (
    gpd.read_file(out / "compartments_rp.gpkg")
    if proc.returncode == 0
    else gpd.GeoDataFrame()
)
if len(final) and glofas is not None:
    r1 = covering(final, "R1").iloc[0]
    results.append(
        (
            f"ring 1 held at RP {rp_e:.0f} (E) above its model interval: model_obs_conflict, rp = 200, "
            "interval from the observation up to the next class",
            100 < rp_e < 200
            and r1.model_obs_conflict
            and abs(r1.obs_rp_low - rp_e) < 1e-6 * rp_e
            and r1.rp == 200
            and abs(r1.rp_low - rp_e) < 1e-6 * rp_e
            and r1.rp_high == 500,
        )
    )
else:
    results.append(("rp_bounds run 3", False))

# --- validate.py ---
proc = run("processing.validate", "validate")
report = (
    json.loads((out / "validation.json").read_text()) if proc.returncode == 0 else {}
)
if report:
    union_result = shapely.union_all(final.geometry.values)
    union_reference = shapely.union_all(reference.geometry.values)
    overlap = union_result.intersection(union_reference).area
    coverage = gpd.read_file(out / "validation_coverage.gpkg")
    matrix = pd.read_csv(out / "validation_rp_confusion.csv", index_col=0)
    matrix.columns = matrix.columns.astype(float)
    expected = pd.DataFrame(0, index=matrix.index, columns=matrix.columns)
    for reference_class, rp in [(200.0, 200.0), (20.0, 100.0), (100.0, 100.0)]:
        expected.loc[reference_class, rp] = 1
    results.append(
        (
            f"area recall {report['recall']:.3f} and precision {report['precision']:.3f}; one reference missed",
            abs(report["recall"] - overlap / union_reference.area) < 1e-9
            and abs(report["precision"] - overlap / union_result.area) < 1e-9
            and int((coverage.covered_m2 == 0).sum()) == 1,
        )
    )
    results.append(
        (
            "RP validation: 3 compartments compared, 2/3 within rp_low..rp_high and in the same class, "
            "confusion matrix",
            report["rp_compartments_compared"] == 3
            and abs(report["rp_share_within_interval"] - 2 / 3) < 1e-9
            and abs(report["rp_share_same_class"] - 2 / 3) < 1e-9
            and list(matrix.index) == CLASSES
            and matrix.equals(expected),
        )
    )
else:
    results.append(("validate runs", False))
server.shutdown()

print("Timings (s): " + ", ".join(f"{k} {v:.1f}" for k, v in timings.items()))
for path in sorted((root / "output" / "T4" / "logs").glob("gfm_evidence_*.log")):
    for line in path.read_text().splitlines():
        if "scenes x" in line:
            print("  " + line.split(": ", 1)[-1])
for label, ok in results:
    print(f"{'PASS' if ok else 'FAIL'}  {label}")
sys.exit(0 if all(ok for _, ok in results) else 1)
