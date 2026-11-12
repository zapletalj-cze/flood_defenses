# Author: Jakub Zapletal
# Date: 2026-11-12
# Transport embankments from the DTM (Stage G step 3) on a synthetic country: an embankment with a bridge over the
# river and an underpass, a cutting through a ridge, a road at grade, a road on a bench of a hillside (one-sided),
# and an embankment over a hidden culvert (a valley under a continuous crest) with a mapped stream crossing it.
# Run: pixi run python -m tests.test_transport  (DEFENSES_TEST_DIR sets the working folder, default: temporary)

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import geopandas as gpd
import numpy as np
import shapely
from osgeo import gdal, osr
from shapely import affinity
from shapely.geometry import LineString, Point

gdal.UseExceptions()

EPSG = 32633  # synthetic test data only
X0, Y_TOP = 600000.0, 5602000.0
WIDTH_M, HEIGHT_M, NATIVE = 2000, 1200, 1.0
RIVER_Y = 600.0

root = Path(
    os.environ.get("DEFENSES_TEST_DIR")
    or tempfile.mkdtemp(prefix="defenses_transport_")
)
data = root / "data"
data.mkdir(parents=True, exist_ok=True)
print(f"Test folder: {root}")


def world(geometry):
    return affinity.translate(geometry, X0, Y_TOP - HEIGHT_M)


cols, rows = int(WIDTH_M / NATIVE), int(HEIGHT_M / NATIVE)
xx, yy = np.meshgrid(np.arange(cols) + 0.5, HEIGHT_M - (np.arange(rows) + 0.5))
# valley floor rising 0.3 % away from the river, 2 m deep channel
dtm = (
    100.0
    + 0.003 * np.abs(yy - RIVER_Y)
    + 0.0005 * xx
    - 2.0 * (np.abs(yy - RIVER_Y) <= 6)
)


def embankment(distance, height, top=4.0, slope=2.0):
    # trapezoid fill: full height within top/2 of the axis, side slopes 1:slope
    return np.clip(height - np.maximum(distance - top / 2, 0) / slope, 0, None)


# R1: north-south embankment 2.5 m at x = 400 with a bridge over the river (30 m) and an underpass (8 m) at y = 900
fill = embankment(np.abs(xx - 400), 2.5)
fill[(np.abs(yy - RIVER_Y) <= 15)] = 0.0
fill[(np.abs(yy - 900) <= 4)] = 0.0
dtm += fill
# R5: east-west embankment 2 m at y = 250 over a valley at x = 900 (no break in the crest: a hidden culvert) and a
# mapped stream at x = 700
valley = np.clip(1.2 * (1 - np.abs(xx - 900) / 60), 0, None) * ((yy > 100) & (yy < 400))
base5 = dtm.copy()
dtm -= valley
crest5 = 102.0 + 0.003 * abs(250 - RIVER_Y) + 0.0005 * xx
near5 = (np.abs(yy - 250) <= 12) & (xx > 480) & (xx < 1320)
dtm = np.where(
    near5, np.maximum(dtm, crest5 - np.maximum(np.abs(yy - 250) - 2, 0) / 2), dtm
)
# R2: cutting at x = 1000 through a ridge (4 m) between y = 850 and 1050
ridge = (yy > 850) & (yy < 1050) & (np.abs(xx - 1000) < 120)
dtm += 4.0 * ridge * (np.abs(xx - 1000) > 6)
# R4: road on a bench of a hillside at x = 1700: terrain rises east of it, a 1.5 m fill drops to the west
hill = (yy > 700) & (yy < 1100)
dtm += hill * np.where(xx > 1704, 0.08 * (xx - 1704) + 1.5, 0.0)
dtm += hill * embankment(np.maximum(1700 - xx, 0) + 0.0, 1.5, top=8.0) * (xx <= 1704)
# R3: road at grade at x = 1400 between y = 100 and 450: nothing to add

srs = osr.SpatialReference()
srs.ImportFromEPSG(EPSG)
ds = gdal.GetDriverByName("GTiff").Create(
    str(data / "dtm.tif"), cols, rows, 1, gdal.GDT_Float32, ["COMPRESS=DEFLATE"]
)
ds.SetGeoTransform([X0, NATIVE, 0, Y_TOP, 0, -NATIVE])
ds.SetProjection(srs.ExportToWkt())
ds.GetRasterBand(1).SetNoDataValue(-9999.0)
ds.GetRasterBand(1).WriteArray(dtm.astype(np.float32))
ds = None

roads = {
    "R1": LineString([(400, 100), (400, 1100)]),
    "R2": LineString([(1000, 760), (1000, 1140)]),
    "R3": LineString([(1400, 100), (1400, 450)]),
    "R4": LineString([(1700, 720), (1700, 1080)]),
    "R5": LineString([(500, 250), (1300, 250)]),
}
gpd.GeoDataFrame(
    {"name": list(roads)}, geometry=[world(g) for g in roads.values()], crs=EPSG
).to_file(data / "roads.gpkg")
rivers = {
    "main": LineString([(WIDTH_M, RIVER_Y), (0, RIVER_Y)]),
    "stream": LineString([(700, 50), (700, RIVER_Y)]),
}
gpd.GeoDataFrame(
    {"name": list(rivers)}, geometry=[world(g) for g in rivers.values()], crs=EPSG
).to_file(data / "rivers.gpkg")
gpd.GeoDataFrame({"name": []}, geometry=gpd.GeoSeries([], crs=EPSG)).to_file(
    data / "levees.gpkg"
)

config = f"""[general]
work_dir = {json.dumps(str(root / 'work'))}
output_dir = {json.dumps(str(root / 'output'))}
n_workers = 2
max_ram_gb_per_worker = 4.0
gdal_cachemax_mb = 256
warp_memory_limit_mb = 256
gdal_threads_per_worker = 1
log_level = "INFO"

[benchmark]
n_windows = 2
window_px = 256
sequential_rows = 64
seed = 1

[countries.T7]
dtm = [{json.dumps(str(data / 'dtm.tif'))}]
dsm = []
levees = [{json.dumps(str(data / 'levees.gpkg'))}]
rivers = [{json.dumps(str(data / 'rivers.gpkg'))}]
transport_lines = [{json.dumps(str(data / 'roads.gpkg'))}]
working_resolution_m = 2.0
corridor_buffer_m = 800.0
unit_core_size_m = 1000.0
unit_halo_m = 800.0
convert_dtm = false
dtm_conversion = "A"
dtm_compress = "ZSTD"
footprint_resolution_m = 20.0
bytes_per_cell = 400.0
crest_box_size_m = 100.0
merge_overlap_ratio = 0.5
drop_contradicted_truncated = true
contradiction_ratio = 0.5
id_rounding_m = 10.0
keep_intermediate = false

[countries.T7.protected_areas]
gap_tolerance_m = 20.0
touch_tolerance_m = 0.01
self_close_min_length_factor = 3.0
crest_search_radius_m = 5.0
barrier_removal_width_m = 10.0
min_margin_m = 0.3
max_level_m = 20.0
min_area_m2 = 10000.0
gated_crossings = false
fill_nodata_holes = true
stream_burn_m = 5.0
stream_burn_buffer_m = 20.0
preserve_river_cells = true
cap_river_holes = true
max_fill_distance_m = 100.0

[countries.T7.rp]
levee_adjacent_distance_m = 10.0
"""
(root / "test.toml").write_text(config)


def run(module, **variables):
    env = dict(
        os.environ,
        DEFENSES_CONFIG_PATH=str(root / "test.toml"),
        DEFENSES_COUNTRY="T7",
        DEFENSES_FORCE="1",
        **{f"DEFENSES_{k}": str(v) for k, v in variables.items()},
    )
    proc = subprocess.run(
        [sys.executable, "-m", module],
        env=env,
        capture_output=True,
        text=True,
        cwd=Path(__file__).resolve().parents[1],
    )
    (root / f"{module.split('.')[-1]}.log").write_text(proc.stdout + proc.stderr)
    if proc.returncode != 0:
        print(proc.stderr[-3000:])
    return proc


list_res = []
for module in [
    "processing.prepare_dtm",
    "processing.build_units",
    "processing.protected_areas",
]:
    proc = run(module, RUN_BENCHMARK="0")
    list_res.append((f"{module} runs", proc.returncode == 0))
proc = run("processing.transport_embankments")
out = root / "output" / "T7"
list_res.append(("transport_embankments runs", proc.returncode == 0))
gdf_seg = (
    gpd.read_file(out / "transport_embankments.gpkg")
    if proc.returncode == 0
    else gpd.GeoDataFrame()
)
gdf_op = (
    gpd.read_file(out / "openings.gpkg") if proc.returncode == 0 else gpd.GeoDataFrame()
)


def class_at(name, x, y):
    # class of the segment of a road nearest to a local point
    point = world(Point(x, y))
    own = gdf_seg[gdf_seg.geometry.distance(world(roads[name])) < 1]
    return (
        own.iloc[int(np.argmin(own.geometry.distance(point)))]["class"]
        if len(own)
        else None
    )


def share(name, cls):
    own = gdf_seg[gdf_seg.geometry.distance(world(roads[name])) < 1]
    return own.loc[own["class"] == cls, "length_m"].sum() / max(
        own["length_m"].sum(), 1e-9
    )


def openings_near(x, y, distance=25):
    return gdf_op[gdf_op.geometry.distance(world(Point(x, y))) <= distance]


if len(gdf_seg):
    list_res.append(
        (
            f"R1 embankment ({share('R1', 'embankment'):.0%} of its length), crest height "
            f"{gdf_seg.loc[gdf_seg['class'] == 'embankment', 'height_min_m'].median():.2f} m on the lower side",
            share("R1", "embankment") > 0.85
            and class_at("R1", 400, 300) == "embankment",
        )
    )
    list_res.append(
        (
            f"R2 cutting through the ridge ({class_at('R2', 1000, 950)}), at grade outside it",
            class_at("R2", 1000, 950) == "cutting"
            and class_at("R2", 1000, 790) == "at_grade",
        )
    )
    list_res.append(
        (f"R3 at grade ({share('R3', 'at_grade'):.0%})", share("R3", "at_grade") > 0.95)
    )
    list_res.append(
        (
            f"R4 on a bench: one-sided ({class_at('R4', 1700, 900)})",
            class_at("R4", 1700, 900) == "one_sided",
        )
    )
    list_res.append(
        (
            f"R5 embankment over the valley, the crest continuous ({class_at('R5', 900, 250)})",
            class_at("R5", 900, 250) == "embankment"
            and share("R5", "embankment") > 0.85,
        )
    )
bridge = openings_near(400, RIVER_Y)
under = openings_near(400, 900)
culvert = openings_near(900, 250, 30)
stream = openings_near(700, 250)
list_res.append(
    (
        f"bridge over the river: {bridge[['type', 'width_m', 'confidence']].to_dict('records')}",
        len(bridge) == 1
        and bridge.iloc[0]["type"] == "bridge"
        and 20 <= bridge.iloc[0]["width_m"] <= 50
        and bridge.iloc[0]["confidence"] == "high"
        and bridge.iloc[0]["invert_z"] < bridge.iloc[0]["crest_z"] - 2,
    )
)
list_res.append(
    (
        f"underpass: {under[['type', 'width_m']].to_dict('records')}",
        len(under) == 1
        and under.iloc[0]["type"] == "underpass"
        and under.iloc[0]["width_m"] <= 30,
    )
)
list_res.append(
    (
        f"hidden culvert in the valley under R5: {culvert[['type', 'confidence', 'evidence']].to_dict('records')}",
        len(culvert) == 1
        and culvert.iloc[0]["type"] == "culvert"
        and culvert.iloc[0]["confidence"] in ("medium", "high"),
    )
)
list_res.append(
    (
        f"culvert where the mapped stream crosses R5: {stream[['type', 'confidence']].to_dict('records')}",
        len(stream) == 1
        and stream.iloc[0]["type"] == "culvert"
        and stream.iloc[0]["confidence"] == "high",
    )
)
others = (
    gdf_op[
        ~gdf_op.index.isin(
            bridge.index.union(under.index).union(culvert.index).union(stream.index)
        )
    ]
    if len(gdf_op)
    else gdf_op
)
list_res.append(
    (f"no other openings ({len(others)})", len(gdf_op) > 0 and len(others) == 0)
)
index = out / "plots" / "transport" / "base" / "index.html"
kinds = [d.name for d in index.parent.iterdir() if d.is_dir()] if index.exists() else []
list_res.append(
    (
        f"plots with index.html ({sorted(kinds)})",
        index.exists()
        and {"cross_sections", "crest_profiles", "openings_map", "summary"}
        <= set(kinds),
    )
)

for label, ok in list_res:
    print(f"{'PASS' if ok else 'FAIL'}  {label}")
sys.exit(0 if all(ok for _, ok in list_res) else 1)
