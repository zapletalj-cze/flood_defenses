# Author: Jakub Zapletal
# Date: 2026-11-04
# End-to-end checks of Stage D on the Stage B test country under its own code T4D. Step D1: DSM inspection and
# benchmark, OSM buildings/water/crossings/walls/industrial and candidate wall tag counts, GHSL settlement model and
# population, CORINE legend mapping verified on a known pixel, Urban Atlas preparation, context zones, zones_to_grid.
# Run: pixi run python -m tests.test_stage_d  (DEFENSES_TEST_DIR sets the working folder, default: temporary;
# the Stage B test runs first in the same folder when its outputs are missing)

import os
import re
import subprocess
import sys
import tempfile
import time
import zipfile
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from osgeo import gdal, ogr, osr
from pyproj import Transformer
from shapely import affinity
from shapely.geometry import Point, box

import gis

gdal.UseExceptions()
ogr.UseExceptions()

# The Stage B test country (tests/test_stage_b.py): synthetic test data only.
EPSG = 32633
X0, Y_TOP = 500000.0, 5501600.0
WIDTH_M, HEIGHT_M = 2400, 1600
MOLLWEIDE = "ESRI:54009"
repo = Path(__file__).resolve().parents[1]

root = Path(
    os.environ.get("DEFENSES_TEST_DIR") or tempfile.mkdtemp(prefix="defenses_stage_d_")
)
data = root / "data"
data_d = root / "data_d"
data_d.mkdir(parents=True, exist_ok=True)
print(f"Test folder: {root}")
if not (root / "output" / "T4" / "compartments_rp.gpkg").exists():
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
        sys.exit("The Stage B test failed; Stage D needs its data.")
out = root / "output" / "T4D"


def world(geometry):
    return affinity.translate(geometry, X0, Y_TOP - HEIGHT_M)


def local_xy(transformer, gx, gy):
    x, y = transformer.transform(gx, gy)
    return x - X0, y - (Y_TOP - HEIGHT_M)


def in_box(x, y, bounds):
    return (x >= bounds[0]) & (x < bounds[2]) & (y >= bounds[1]) & (y < bounds[3])


# --- DSM: the DTM plus a building block and a thin wall, tiled (the DTM is in strip layout) ---
dtm_ds = gdal.Open(str(data / "dtm.tif"))
surface = dtm_ds.GetRasterBand(1).ReadAsArray()
gt = dtm_ds.GetGeoTransform()
rows, cols = np.indices(surface.shape)
lx, ly = (cols + 0.5) * gt[1], HEIGHT_M - (rows + 0.5) * abs(gt[5])
surface = surface + np.where(in_box(lx, ly, (1600, 300, 1700, 400)), 12.0, 0.0)
surface = surface + np.where(
    (np.abs(ly - 820) < 0.5) & (lx > 100) & (lx < 900), 1.5, 0.0
)
dsm_path = data_d / "dsm.tif"
dsm_ds = gdal.GetDriverByName("GTiff").Create(
    str(dsm_path),
    surface.shape[1],
    surface.shape[0],
    1,
    gdal.GDT_Float32,
    ["TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256", "COMPRESS=DEFLATE"],
)
dsm_ds.SetGeoTransform(gt)
dsm_ds.SetProjection(dtm_ds.GetProjection())
dsm_ds.GetRasterBand(1).SetNoDataValue(-9999.0)
dsm_ds.GetRasterBand(1).WriteArray(surface.astype(np.float32))
dsm_ds = dtm_ds = None
del surface, rows, cols, lx, ly

# --- OSM extract: buildings, water, crossings, walls, industrial areas ---
to_lonlat = Transformer.from_crs(EPSG, 4326, always_xy=True)
nodes, ways = [], []


def node(x, y):
    lon, lat = to_lonlat.transform(X0 + x, Y_TOP - HEIGHT_M + y)
    nodes.append((len(nodes) + 1, lat, lon))
    return len(nodes)


def closed(bounds):
    ids = [node(bounds[0], bounds[1]), node(bounds[2], bounds[1])]
    ids += [node(bounds[2], bounds[3]), node(bounds[0], bounds[3])]
    return ids + [ids[0]]


def line(points):
    return [node(x, y) for x, y in points]


way_id = 1000
for refs, tags in [
    (closed((300, 1000, 330, 1020)), {"building": "yes"}),
    (closed((340, 1000, 370, 1020)), {"building": "yes"}),
    (closed((1200, 1500, 1230, 1530)), {"building": "industrial"}),
    (closed((1500, 400, 1550, 450)), {"natural": "water"}),
    (line([(400, 700), (400, 1300)]), {"highway": "residential"}),
    (line([(1000, 1500), (2400, 1500)]), {"railway": "rail"}),
    (line([(600, 600), (600, 650)]), {"highway": "footway"}),
    (line([(100, 820), (900, 820)]), {"barrier": "wall", "wall": "flood_wall"}),
    (line([(1000, 830), (1600, 830)]), {"man_made": "dyke"}),
    (line([(50, 1550), (150, 1550)]), {"barrier": "wall"}),
    (closed((1000, 1200, 1150, 1450)), {"landuse": "industrial"}),
    (closed((200, 300, 350, 450)), {"man_made": "works"}),
]:
    ways.append((way_id, refs, tags))
    way_id += 1
xml = [
    "<?xml version='1.0' encoding='UTF-8'?>",
    '<osm version="0.6" generator="synthetic">',
]
xml += [
    f'  <node id="{nid}" lat="{lat:.8f}" lon="{lon:.8f}" version="1"/>'
    for nid, lat, lon in nodes
]
for wid, refs, tags in ways:
    xml.append(
        f'  <way id="{wid}" version="1">'
        + "".join(f'<nd ref="{r}"/>' for r in refs)
        + "".join(f'<tag k="{k}" v="{v}"/>' for k, v in tags.items())
        + "</way>"
    )
xml.append("</osm>")
(data_d / "t4d.osm").write_text("\n".join(xml))

# --- GHSL settlement model and population: zipped Mollweide tiles of the Stage B tiling schema ---
schema = gpd.read_file(f"/vsizip/{data / 'tiling_schema.zip'}/tiling_schema.shp")
from_moll = Transformer.from_crs(MOLLWEIDE, EPSG, always_xy=True)
moll = osr.SpatialReference()
moll.SetFromUserInput(MOLLWEIDE)
for tile, geometry in zip(schema["tile_id"], schema.geometry):
    tx0, ty0, tx1, ty1 = geometry.bounds
    ncols, nrows = int(round((tx1 - tx0) / 100)), int(round((ty1 - ty0) / 100))
    cx, cy = np.meshgrid(
        tx0 + (np.arange(ncols) + 0.5) * 100, ty1 - (np.arange(nrows) + 0.5) * 100
    )
    x, y = local_xy(from_moll, cx, cy)
    smod = np.where(y >= 700, np.where(x < 1600, 30, 22), 11).astype(np.uint8)
    pop = np.where(in_box(x, y, (250, 900, 1550, 1300)), 50.0, 2.0).astype(np.float32)
    for product, values, dtype in [
        ("SMOD", smod, gdal.GDT_Byte),
        ("POP", pop, gdal.GDT_Float32),
    ]:
        name = f"GHS_{product}_E2020_TEST_{tile}"
        tif = data_d / f"{name}.tif"
        ds = gdal.GetDriverByName("GTiff").Create(str(tif), ncols, nrows, 1, dtype)
        ds.SetGeoTransform([tx0, 100.0, 0, ty1, 0, -100.0])
        ds.SetProjection(moll.ExportToWkt())
        ds.GetRasterBand(1).WriteArray(values)
        ds = None
        with zipfile.ZipFile(data_d / f"{name}.zip", "w") as archive:
            archive.write(tif, tif.name)
        tif.unlink()

# --- CORINE: class indices at 100 m in EPSG:3035 with a legend; one pixel with a known code ---
to_laea = Transformer.from_crs(EPSG, 3035, always_xy=True)
from_laea = Transformer.from_crs(3035, EPSG, always_xy=True)
corners = np.array(
    [
        to_laea.transform(X0 + x, Y_TOP - HEIGHT_M + y)
        for x in (-500, WIDTH_M + 500)
        for y in (-500, HEIGHT_M + 500)
    ]
)
lx0, ly0 = np.floor(corners.min(axis=0) / 100) * 100
lx1, ly1 = np.ceil(corners.max(axis=0) / 100) * 100
ncols, nrows = int((lx1 - lx0) / 100), int((ly1 - ly0) / 100)
cx, cy = np.meshgrid(
    lx0 + (np.arange(ncols) + 0.5) * 100, ly1 - (np.arange(nrows) + 0.5) * 100
)
x, y = local_xy(from_laea, cx, cy)
index = np.full(x.shape, 12, np.uint8)  # 211 arable land
index[in_box(x, y, (600, 1200, 1000, 1450))] = 3  # 121 industrial or commercial units
index[in_box(x, y, (1700, 1200, 1900, 1450))] = 1  # 111 continuous urban fabric
# One isolated industrial pixel in the dense urban area: smaller than min_zone_area_m2, so it is dropped.
small = int(np.argmin(np.hypot(x - 300, y - 1100)))
index.ravel()[small] = 3
small_xy = (float(x.ravel()[small]), float(y.ravel()[small]))
known = np.argwhere(index == 3)[0]
check_x, check_y = lx0 + (known[1] + 0.5) * 100, ly1 - (known[0] + 0.5) * 100
laea = osr.SpatialReference()
laea.ImportFromEPSG(3035)
ds = gdal.GetDriverByName("GTiff").Create(
    str(data_d / "corine_index.tif"), ncols, nrows, 1, gdal.GDT_Byte
)
ds.SetGeoTransform([lx0, 100.0, 0, ly1, 0, -100.0])
ds.SetProjection(laea.ExportToWkt())
ds.GetRasterBand(1).SetNoDataValue(48)
ds.GetRasterBand(1).WriteArray(index)
ds = None
pd.DataFrame(
    {
        "GRID_CODE": [1, 2, 3, 12, 48],
        "CLC_CODE": [111, 112, 121, 211, 999],
        "LABEL3": [
            "Continuous urban fabric",
            "Discontinuous urban fabric",
            "Industrial or commercial units",
            "Non-irrigated arable land",
            "NODATA",
        ],
    }
).to_csv(data_d / "corine_legend.csv", index=False)

# --- Urban Atlas: two-layer GeoPackage (land use polygons and a boundary) in EPSG:3035 ---
ua_path = data_d / "urban_atlas_t4d.gpkg"
ua_path.unlink(missing_ok=True)
ua_ds = ogr.GetDriverByName("GPKG").CreateDataSource(str(ua_path))
layer = ua_ds.CreateLayer("T4D001L1_UA2018", laea, ogr.wkbPolygon)
layer.CreateField(ogr.FieldDefn("code_2018", ogr.OFTString))
for code, bounds in [
    ("12100", (2000, 1200, 2200, 1450)),
    ("11100", (300, 1250, 500, 1450)),
]:
    feature = ogr.Feature(layer.GetLayerDefn())
    feature.SetField("code_2018", code)
    projected = shapely.transform(
        world(box(*bounds)),
        lambda xy: np.column_stack(to_laea.transform(xy[:, 0], xy[:, 1])),
    )
    feature.SetGeometry(ogr.CreateGeometryFromWkb(shapely.to_wkb(projected)))
    layer.CreateFeature(feature)
boundary = ua_ds.CreateLayer("T4D001L1_UA2018_Boundary", laea, ogr.wkbPolygon)
feature = ogr.Feature(boundary.GetLayerDefn())
feature.SetGeometry(ogr.CreateGeometryFromWkt("POLYGON ((0 0, 1 0, 1 1, 0 0))"))
boundary.CreateFeature(feature)
feature = layer = boundary = ua_ds = None

# --- Configuration: the Stage B configuration plus the country T4D (copy of T4 with the Stage D inputs) ---
base = (root / "test.toml").read_text()
start, end = base.index("[countries.T4]\n"), base.index("[sources]\n")
t4d = base[start:end].replace("[countries.T4", "[countries.T4D")
t4d = t4d.replace(
    "dsm = []\n",
    f"""dsm = "{dsm_path}"
corine = "{data_d / 'corine_index.tif'}"
corine_legend = "{data_d / 'corine_legend.csv'}"
corine_legend_index_field = "GRID_CODE"
corine_legend_code_field = "CLC_CODE"
corine_check = [[{check_x}, {check_y}, 121]]
""",
    1,
)
t4d = t4d.replace('osm_pbf = ""', f'osm_pbf = "{data_d / "t4d.osm"}"', 1)
t4d = t4d.replace(
    "[countries.T4D.lulc]\nenabled = false\n",
    f"""[countries.T4D.lulc]
enabled = false
urban_atlas = "{ua_path}"
urban_atlas_layer = "T4D001L1_UA2018"
urban_atlas_class_field = "code_2018"
""",
    1,
)
t4d = re.sub(
    r"\[countries\.T4D\.context\]\n(?:[^\[\n][^\n]*\n|\n)*",
    """[countries.T4D.context]
order = ["industrial", "dense_urban", "urban", "rural"]
corine_industrial_codes = [121, 123]
corine_dense_urban_codes = [111]
urban_atlas_industrial_codes = [12100]
urban_atlas_dense_urban_codes = [11100]
ghsl_smod_dense_urban = [30]
ghsl_smod_urban = [21, 22, 23]
osm_industrial = { landuse = ["industrial"], man_made = ["works"] }
min_zone_area_m2 = 20000.0

""",
    t4d,
    count=1,
)
config = base.replace("[sources]\n", t4d + "[sources]\n", 1)
config = re.sub(
    r"(ghsl_url_templates = \{[^\n]*?)\s*\}\n",
    lambda m: m.group(1)
    + f', smod = "file://{data_d}/GHS_SMOD_E{{epoch}}_TEST_{{tile}}.zip", '
    + f'pop = "file://{data_d}/GHS_POP_E{{epoch}}_TEST_{{tile}}.zip" }}\n',
    config,
    count=1,
)
config = config.replace(
    "osm_tags_walls = {}",
    'osm_tags_walls = { wall = ["flood_wall"], man_made = ["dyke"] }',
    1,
)
(root / "test_d.toml").write_text(config)
(root / "test_d_bad_check.toml").write_text(
    config.replace(
        f"corine_check = [[{check_x}, {check_y}, 121]]",
        f"corine_check = [[{check_x}, {check_y}, 112]]",
        1,
    )
)

timings = {}


def run(module, label, config_name="test_d.toml", **variables):
    env = dict(
        os.environ,
        DEFENSES_CONFIG_PATH=str(root / config_name),
        DEFENSES_COUNTRY="T4D",
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
    (root / f"d_{label}.log").write_text(proc.stdout + proc.stderr)
    if proc.returncode != 0 and "bad" not in label:
        print(proc.stderr[-3000:])
    return proc


results = []

# --- prepare_dtm.py with a DSM, build_units.py ---
proc = run("processing.prepare_dtm", "prepare_dtm")
inspection = (
    pd.read_csv(out / "dsm_inspection.csv")
    if (out / "dsm_inspection.csv").exists()
    else pd.DataFrame()
)
bench = (
    pd.read_csv(out / "dsm_benchmark.csv")
    if (out / "dsm_benchmark.csv").exists()
    else pd.DataFrame()
)
results.append(
    (
        f"prepare_dtm: DSM inspected (tiled {inspection.get('block_x', pd.Series([0])).iloc[0]} px blocks) and "
        f"benchmarked ({len(bench)} windows) like the DTM",
        proc.returncode == 0
        and len(inspection) == 1
        and not bool(inspection["strip_layout"].iloc[0])
        and len(bench) == 2
        and (out / "dtm_benchmark.csv").exists()
        and "Benchmark DSM" in proc.stderr,
    )
)
proc = run("processing.build_units", "build_units")
results.append(("build_units for T4D", proc.returncode == 0))

# --- download_osm.py: Stage D outputs and candidate counts ---
proc = run("download.download_osm", "download_osm", FORCE="1")


def count(name):
    path = out / f"{name}.gpkg"
    return len(gpd.read_file(path)) if path.exists() else -1


osm_counts = {
    n: count(n)
    for n in [
        "osm_buildings",
        "osm_water",
        "osm_crossings",
        "osm_walls",
        "osm_industrial",
    ]
}
results.append(
    (
        f"download_osm: buildings, water, crossings (no footway), walls, industrial areas {osm_counts}",
        proc.returncode == 0
        and osm_counts
        == {
            "osm_buildings": 3,
            "osm_water": 1,
            "osm_crossings": 2,
            "osm_walls": 2,
            "osm_industrial": 2,
        },
    )
)
results.append(
    (
        "candidate wall tags counted in the extract (barrier=wall 2, wall=flood_wall 1, man_made=dyke 1)",
        all(
            t in proc.stderr
            for t in ["barrier=wall: 2", "wall=flood_wall: 1", "man_made=dyke: 1"]
        ),
    )
)

# --- download_ghsl.py: settlement model and population of the same epoch ---
proc = run("download.download_ghsl", "download_ghsl")
smod_ok = False
if (
    proc.returncode == 0
    and (out / "ghsl_smod.vrt").exists()
    and (out / "ghsl_pop.vrt").exists()
):
    values = gdal.Open(str(out / "ghsl_smod.vrt")).ReadAsArray()
    smod_ok = set(np.unique(values)) >= {11, 22, 30}
results.append(("download_ghsl: SMOD and GHS-POP mosaics of the probed epoch", smod_ok))

# --- download_clms.py: a wrong known pixel stops it; the right one passes ---
proc_bad = run(
    "download.download_clms", "download_clms_bad", config_name="test_d_bad_check.toml"
)
proc = run("download.download_clms", "download_clms")
codes_ok = False
if proc.returncode == 0 and (out / "corine_codes.tif").exists():
    codes_ds = gdal.Open(str(out / "corine_codes.tif"))
    codes_gt = codes_ds.GetGeoTransform()
    col = int((check_x - codes_gt[0]) / codes_gt[1])
    row = int((check_y - codes_gt[3]) / codes_gt[5])
    present = set(np.unique(codes_ds.ReadAsArray()))
    codes_ok = (
        int(codes_ds.ReadAsArray(col, row, 1, 1)[0, 0]) == 121
        and {111, 121, 211} <= present
    )
results.append(
    (
        "CORINE legend mapping verified on a known pixel: wrong code stops the script, codes 111/121/211 written",
        proc_bad.returncode != 0
        and "legend check failed" in proc_bad.stderr
        and codes_ok,
    )
)
ua = (
    gpd.read_file(out / "urban_atlas.gpkg")
    if (out / "urban_atlas.gpkg").exists()
    else gpd.GeoDataFrame()
)
results.append(
    (
        "Urban Atlas: layer of the two-layer file merged in the DTM CRS with numeric codes",
        len(ua) == 2
        and set(ua["ua_code"]) == {12100, 11100}
        and ua.crs.to_epsg() == EPSG,
    )
)

# --- context_zones.py ---
proc = run("processing.context_zones", "context_zones")
zones = (
    gpd.read_file(out / "zones.gpkg") if proc.returncode == 0 else gpd.GeoDataFrame()
)
probes = {
    "industrial (CORINE 121)": ((800, 1300), "industrial"),
    "industrial (OSM landuse)": ((1075, 1300), "industrial"),
    "industrial (Urban Atlas over SMOD urban)": ((2100, 1300), "industrial"),
    "industrial (OSM works in the rural part)": ((275, 375), "industrial"),
    "small industrial pixel dropped -> dense_urban": (small_xy, "dense_urban"),
    "dense_urban (SMOD 30)": ((400, 1000), "dense_urban"),
    "dense_urban (CORINE 111 over SMOD urban)": ((1800, 1300), "dense_urban"),
    "urban (SMOD 22)": ((2300, 1000), "urban"),
    "rural": ((1000, 400), "rural"),
}
if len(zones):
    found = {}
    for label, (xy, _) in probes.items():
        hit = zones[zones.covers(world(Point(*xy)))]
        found[label] = hit["context"].iloc[0] if len(hit) == 1 else f"{len(hit)} zones"
    results.append(
        (
            "context zones: order wins where sources overlap; small zones fall through; rural is the rest",
            all(found[k] == v for k, (_, v) in probes.items()),
        )
    )
    for label, (_, expected) in probes.items():
        if found[label] != expected:
            print(f"  {label}: {found[label]} (expected {expected})")
    corridor = shapely.union_all(
        gpd.read_file(out / "corridor.gpkg").to_crs(EPSG).geometry.values
    )
    overlap = sum(
        shapely.intersection(a, b).area
        for i, a in enumerate(zones.geometry)
        for b in zones.geometry.iloc[i + 1 :]
    )
    results.append(
        (
            f"zones do not overlap ({overlap:.1f} m2) and tile the corridor ({zones.area.sum() / corridor.area:.6f})",
            overlap < 1.0 and abs(zones.area.sum() / corridor.area - 1) < 1e-6,
        )
    )
    # zones_to_grid on a unit grid: earlier contexts carry lower values.
    units = pd.read_csv(out / "units.csv", dtype={"unit_id": str})
    unit = units[units["unit_id"] == "000_000"].iloc[0]
    res = 2.0
    grid_gt = (unit.win_xmin, res, 0.0, unit.win_ymax, 0.0, -res)
    shape = (
        int(round((unit.win_ymax - unit.win_ymin) / res)),
        int(round((unit.win_xmax - unit.win_xmin) / res)),
    )
    grid = gis.Vector.zones_to_grid(
        zones,
        ["industrial", "dense_urban", "urban", "rural"],
        shape,
        grid_gt,
        gdal.Open(str(root / "work" / "T4D" / "dtm" / "dtm.vrt")).GetProjection(),
    )
    checks = []
    for (px, py), expected in [((800, 1300), 1), ((400, 1000), 2), ((1000, 400), 4)]:
        c = int((X0 + px - grid_gt[0]) / res)
        r = int((grid_gt[3] - (Y_TOP - HEIGHT_M + py)) / res)
        checks.append(int(grid[r, c]) == expected)
    results.append(("zones_to_grid: context codes on a unit grid", all(checks)))
else:
    results.append(("context_zones runs", False))

print("Timings (s): " + ", ".join(f"{k} {v:.1f}" for k, v in timings.items()))
for label, ok in results:
    print(f"{'PASS' if ok else 'FAIL'}  {label}")
sys.exit(0 if all(ok for _, ok in results) else 1)
