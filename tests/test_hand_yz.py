# Author: Jakub Zapletal
# Date: 2026-11-04
# End-to-end check of HAND_YZ on the Stage B test country: WorldCover tiles found through the tile grid and the URL
# template, Urban Atlas from a two-layer file, land cover value per protected area against an independent computation,
# checkpoint hash, and the sum of building and land cover value in rp_index.py.
# Run: pixi run python -m tests.test_hand_yz  (DEFENSES_TEST_DIR sets the working folder, default: temporary;
# the Stage B test runs first in the same folder when its outputs are missing)

import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from osgeo import gdal, ogr, osr
from pyproj import Transformer
from shapely import affinity
from shapely.geometry import Point, box

gdal.UseExceptions()
ogr.UseExceptions()

# The Stage B test country (tests/test_stage_b.py): synthetic test data only.
EPSG = 32633
UA_EPSG = 3035  # the synthetic Urban Atlas file uses another CRS than the DTM
X0, Y_TOP = 500000.0, 5501600.0
WIDTH_M, HEIGHT_M = 2400, 1600
probes = {"R1": (500, 1100), "R2": (1200, 1050), "R3": (2000, 1100)}
interiors = {
    "R1": (302, 952, 698, 1248),
    "R2": (1102, 952, 1498, 1248),
    "R3": (1802, 952, 2198, 1248),
}
WC_RES = 1.0 / 12000.0  # WorldCover pixel size in degrees
UA_SPLIT_X = 1300.0  # Urban Atlas covers ring 2 west of this line with a class that carries no value
VALUES = {
    "agriculture": 2.0,
    "built_up": 30.0,
    "transport": 50.0,
}  # per m2 of land, synthetic
DEPTHS = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0]
CURVES = {
    "agriculture": [0.0, 0.3, 0.55, 0.8, 1.0, 1.0],
    "built_up": [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
    "transport": [0.0, 0.1, 0.25, 0.5, 0.75, 1.0],
}
repo = Path(__file__).resolve().parents[1]

root = Path(
    os.environ.get("DEFENSES_TEST_DIR") or tempfile.mkdtemp(prefix="defenses_hand_yz_")
)
out = root / "output" / "T4"
data = root / "data_hand_yz"
(data / "wc").mkdir(parents=True, exist_ok=True)
print(f"Test folder: {root}")
if not (out / "protected_areas_rp.gpkg").exists():
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
        sys.exit("The Stage B test failed; HAND_YZ needs its outputs.")
# rp_index.py rewrites protected_areas_rp.gpkg: the Stage B result is kept and restored at the end.
stage_b = root / "protected_areas_rp_stage_b.gpkg"
if not stage_b.exists():
    shutil.copy(out / "protected_areas_rp.gpkg", stage_b)
gdf_model = gpd.read_file(stage_b)
gdf_pa = gpd.read_file(out / "protected_areas.gpkg")


def world(geometry):
    return affinity.translate(geometry, X0, Y_TOP - HEIGHT_M)


def grown(bounds, margin):
    return box(
        bounds[0] - margin, bounds[1] - margin, bounds[2] + margin, bounds[3] + margin
    )


member = {
    k: gdf_pa[gdf_pa.covers(world(Point(*xy)))].iloc[0] for k, xy in probes.items()
}

# --- Synthetic WorldCover: two tiles of the official grid (the corridor crosses 15 degrees east) ---
to_local = Transformer.from_crs(4326, EPSG, always_xy=True)
srs = osr.SpatialReference()
srs.ImportFromEPSG(4326)
list_wc = [
    (grown(interiors["R1"], 20), 40),  # cropland -> agriculture
    (grown(interiors["R2"], 20), 50),  # built-up -> built_up
    (grown(interiors["R3"], 20), 10),  # tree cover: no value, replaced by Urban Atlas
]
tiles = {"N48E012": (14.99, 15.0), "N48E015": (15.0, 15.05)}
lat_south, lat_north = 49.64, 49.68
for tile, (west, east) in tiles.items():
    ncols = int(round((east - west) / WC_RES))
    nrows = int(round((lat_north - lat_south) / WC_RES))
    lon, lat = np.meshgrid(
        west + (np.arange(ncols) + 0.5) * WC_RES,
        lat_north - (np.arange(nrows) + 0.5) * WC_RES,
    )
    x, y = to_local.transform(lon, lat)
    x, y = x - X0, y - (Y_TOP - HEIGHT_M)
    classes = np.full(lon.shape, 30, np.uint8)  # grassland: no value
    for geometry, code in list_wc:
        classes[shapely.contains_xy(geometry, x, y)] = code
    path = data / "wc" / f"ESA_WorldCover_10m_2021_v200_{tile}_Map.tif"
    ds = gdal.GetDriverByName("GTiff").Create(
        str(path), ncols, nrows, 1, gdal.GDT_Byte, ["COMPRESS=DEFLATE"]
    )
    ds.SetGeoTransform([west, WC_RES, 0, lat_north, 0, -WC_RES])
    ds.SetProjection(srs.ExportToWkt())
    ds.SetMetadata({"license": "CC-BY 4.0 (synthetic test tile)"})
    ds.GetRasterBand(1).SetNoDataValue(0)
    ds.GetRasterBand(1).WriteArray(classes)
    ds = None
gdf_grid = gpd.GeoDataFrame(
    {"ll_tile": ["N48E012", "N48E015", "N45E015", "N51E015"]},
    geometry=[
        box(12, 48, 15, 51),
        box(15, 48, 18, 51),
        box(15, 45, 18, 48),
        box(15, 51, 18, 54),
    ],
    crs=4326,
)
gdf_grid.to_file(data / "wc" / "worldcover_grid.geojson", driver="GeoJSON")

# --- Synthetic Urban Atlas: a two-layer GeoPackage (land use polygons and a boundary) in EPSG:3035 ---
to_ua = Transformer.from_crs(EPSG, UA_EPSG, always_xy=True)
ua_srs = osr.SpatialReference()
ua_srs.ImportFromEPSG(UA_EPSG)
list_ua = [
    (
        "14100",
        box(1080, 930, UA_SPLIT_X, 1270),
    ),  # no category: ring 2 west of the line has no value
    ("12210", grown(interiors["R3"], 20)),  # transport over ring 3
    ("12210", box(100, 100, 200, 200)),  # touches no protected area: filtered out
]
ua_path = data / "urban_atlas_test.gpkg"
if ua_path.exists():
    ua_path.unlink()
ua_ds = ogr.GetDriverByName("GPKG").CreateDataSource(str(ua_path))
layer = ua_ds.CreateLayer("TEST001L1_UA2018", ua_srs, ogr.wkbPolygon)
layer.CreateField(ogr.FieldDefn("code_2018", ogr.OFTString))
for code, geometry in list_ua:
    feature = ogr.Feature(layer.GetLayerDefn())
    feature.SetField("code_2018", code)
    projected = shapely.transform(
        world(geometry), lambda xy: np.column_stack(to_ua.transform(xy[:, 0], xy[:, 1]))
    )
    feature.SetGeometry(ogr.CreateGeometryFromWkb(shapely.to_wkb(projected)))
    layer.CreateFeature(feature)
boundary = ua_ds.CreateLayer("TEST001L1_UA2018_Boundary", ua_srs, ogr.wkbPolygon)
feature = ogr.Feature(boundary.GetLayerDefn())
feature.SetGeometry(ogr.CreateGeometryFromWkt("POLYGON ((0 0, 1 0, 1 1, 0 0))"))
boundary.CreateFeature(feature)
feature = layer = boundary = ua_ds = None

# --- Damage tables per land cover category ---
pd.DataFrame(
    [{"country": "T4", **VALUES}, {"country": "XX", **{k: 1.0 for k in VALUES}}]
).to_csv(data / "lulc_max_damage.csv", index=False)
pd.DataFrame({"depth_m": DEPTHS, **CURVES}).to_csv(
    data / "lulc_depth_damage.csv", index=False
)

# --- Configuration: the Stage B test configuration with the land cover table and WorldCover sources ---
lulc_table = f"""[countries.T4.lulc]
enabled = true
worldcover = []
urban_atlas = "{ua_path}"
urban_atlas_layer = "TEST001L1_UA2018"
urban_atlas_class_field = "code_2018"
worldcover_categories = {{ "40" = "agriculture", "50" = "built_up" }}
urban_atlas_categories = {{ "12210" = "transport" }}
max_damage_csv = "{data / 'lulc_max_damage.csv'}"
depth_damage_csv = "{data / 'lulc_depth_damage.csv'}"
"""
source_keys = f"""worldcover_url_template = "file://{data / 'wc'}/ESA_WorldCover_10m_{{year}}_{{version}}_{{tile}}_Map.tif"
worldcover_grid_url = "file://{data / 'wc' / 'worldcover_grid.geojson'}"
worldcover_tile_field = "ll_tile"
worldcover_version = "v200"
worldcover_year = 2021
"""
base = (root / "test.toml").read_text()
assert (
    "[countries.T4.lulc]\nenabled = false\n" in base
), "Stage B configuration without the lulc table"
config = base.replace("[countries.T4.lulc]\nenabled = false\n", lulc_table, 1)
config = config.replace("[sources]\n", "[sources]\n" + source_keys, 1)
(root / "test_hand_yz.toml").write_text(config)


def run(module, label, config_name="test_hand_yz.toml", **variables):
    env = dict(
        os.environ,
        DEFENSES_CONFIG_PATH=str(root / config_name),
        DEFENSES_COUNTRY="T4",
        **{f"DEFENSES_{k}": str(v) for k, v in variables.items()},
    )
    proc = subprocess.run(
        [sys.executable, "-m", module],
        env=env,
        capture_output=True,
        text=True,
        cwd=repo,
    )
    (root / f"{label}.log").write_text(proc.stdout + proc.stderr)
    if proc.returncode != 0:
        print(proc.stderr[-3000:])
    return proc


def expected_value(label, factor=1.0, urban_atlas=True):
    # Independent computation from the unit rasters: category of each protected area cell from the synthetic layout.
    row = member[label]
    unit_dir = root / "work" / "T4" / "units" / row.unit_id
    labels_ds = gdal.Open(str(unit_dir / "labels.tif"))
    gt = labels_ds.GetGeoTransform()
    cells = labels_ds.ReadAsArray() == row.unit_label
    hand = gdal.Open(str(unit_dir / "hand_undef.tif")).ReadAsArray()
    rows, cols = np.nonzero(cells)
    x = gt[0] + (cols + 0.5) * gt[1] - X0
    depth = np.clip(row.protection_level_m - hand[rows, cols], 0.0, None)
    if label == "R1":
        category = np.full(len(x), "agriculture")
    elif label == "R2":
        category = np.where(urban_atlas & (x < UA_SPLIT_X), "", "built_up")
    else:
        category = np.full(len(x), "transport" if urban_atlas else "")
    total = 0.0
    for name in VALUES:
        chosen = category == name
        total += float(
            np.sum(
                VALUES[name]
                * factor
                * gt[1] ** 2
                * np.interp(depth[chosen], DEPTHS, CURVES[name])
            )
        )
    return total


list_res = []

# --- download_worldcover.py ---
proc = run("download.download_worldcover", "download_worldcover")
record = (
    json.loads((out / "worldcover_sources.json").read_text())
    if proc.returncode == 0
    else {}
)
list_res.append(
    (
        f"download_worldcover: tiles of the grid that touch the corridor ({record.get('tiles')}), VRT, licence",
        proc.returncode == 0
        and record.get("tiles") == ["N48E012", "N48E015"]
        and (out / "worldcover.vrt").exists()
        and "CC-BY" in record.get("licence", ""),
    )
)

# --- HAND_YZ.py ---
proc = run("processing.HAND_YZ", "hand_yz")
df_exp = (
    pd.read_csv(out / "exposure_lulc.csv").set_index("protected_area_id")
    if proc.returncode == 0
    else pd.DataFrame()
)
if len(df_exp):
    got = {
        k: float(df_exp.loc[member[k].protected_area_id, "value_lulc"]) for k in probes
    }
    want = {k: expected_value(k) for k in probes}
    list_res.append(
        (
            "land cover value per protected area equals the independent computation "
            + ", ".join(f"{k} {got[k]:.0f}/{want[k]:.0f}" for k in probes),
            all(abs(got[k] - want[k]) <= 1e-6 * max(want[k], 1.0) for k in probes)
            and all(want[k] > 0 for k in probes),
        )
    )
    list_res.append(
        (
            "WorldCover classes: cropland -> agriculture (ring 1), built-up -> built_up (ring 2)",
            df_exp.loc[member["R1"].protected_area_id, "lulc_agriculture"] == got["R1"]
            and df_exp.loc[member["R2"].protected_area_id, "lulc_built_up"]
            == got["R2"],
        )
    )
    list_res.append(
        (
            "Urban Atlas takes precedence: ring 3 (tree cover) valued as transport, west of ring 2 without value",
            df_exp.loc[member["R3"].protected_area_id, "lulc_transport"] == got["R3"]
            and df_exp.loc[member["R3"].protected_area_id, "lulc_built_up"] == 0
            and 0 < got["R2"] < 0.8 * expected_value("R2", urban_atlas=False),
        )
    )
else:
    list_res.append(("HAND_YZ runs and writes exposure_lulc.csv", False))
list_copy = list((root / "work" / "T4" / "lulc").glob("*_layer.gpkg"))
gdf_kept = (
    gpd.read_file(root / "work" / "T4" / "lulc" / "urban_atlas_categories.gpkg")
    if (root / "work" / "T4" / "lulc" / "urban_atlas_categories.gpkg").exists()
    else gpd.GeoDataFrame()
)
list_res.append(
    (
        f"Urban Atlas layer copied from the two-layer file; {len(gdf_kept)} polygons touch protected areas",
        len(list_copy) == 1 and len(gdf_kept) == 2,
    )
)

# --- Checkpoints: an unchanged rerun skips all units, new category values recompute them ---
proc = run("processing.HAND_YZ", "hand_yz_rerun")
rerun_skipped = proc.returncode == 0 and "0 to process" in proc.stderr
pd.DataFrame(
    [
        {"country": "T4", **{k: 2 * v for k, v in VALUES.items()}},
        {"country": "XX", **{k: 1.0 for k in VALUES}},
    ]
).to_csv(data / "lulc_max_damage.csv", index=False)
proc = run("processing.HAND_YZ", "hand_yz_new_values")
df_dbl = (
    pd.read_csv(out / "exposure_lulc.csv").set_index("protected_area_id")
    if proc.returncode == 0
    else pd.DataFrame()
)
list_res.append(
    (
        "unchanged rerun skips every unit; changed values are detected and recomputed (value doubles)",
        rerun_skipped
        and "0 to process" not in proc.stderr
        and len(df_dbl) == len(df_exp)
        and np.allclose(df_dbl["value_lulc"], 2 * df_exp["value_lulc"]),
    )
)

# --- rp_index.py: building value plus land cover value ---
proc = run("processing.rp_index", "rp_index_lulc")
gdf_comb = (
    gpd.read_file(out / "protected_areas_rp.gpkg")
    if proc.returncode == 0
    else gpd.GeoDataFrame()
)
if len(gdf_comb):
    gdf_mrg = gdf_comb.set_index("protected_area_id")
    gdf_before = gdf_model.set_index("protected_area_id")
    list_res.append(
        (
            "rp_index: value_exposed = value_ghsl + value_lulc, GHSL part unchanged from Stage B",
            np.allclose(gdf_mrg.value_exposed, gdf_mrg.value_ghsl + gdf_mrg.value_lulc)
            and np.allclose(
                gdf_mrg.loc[gdf_before.index, "value_ghsl"], gdf_before.value_exposed
            )
            and np.allclose(
                gdf_mrg.loc[df_dbl.index, "value_lulc"], df_dbl["value_lulc"]
            )
            and "lulc_transport" in gdf_mrg.columns,
        )
    )
else:
    list_res.append(("rp_index runs with land cover exposure", False))
shutil.copy(stage_b, out / "protected_areas_rp.gpkg")

for label, ok in list_res:
    print(f"{'PASS' if ok else 'FAIL'}  {label}")
sys.exit(0 if all(ok for _, ok in list_res) else 1)
