# Author: Jakub Zapletal
# Date: 2026-11-04
# End-to-end check of Stage B on a synthetic country: OSM download and extraction, water sources (gate at a levee
# crossing, pumping station at a downstream end), pumping stations, GHSL download and mosaic, value warp totals,
# exposure, score, quantile mapping, cap, interval and overrides.
# Run: pixi run python -m tests.test_stage_b  (DEFENSES_TEST_DIR sets the working folder, default: temporary)

import json
import os
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from osgeo import gdal, osr
from pyproj import Transformer
from shapely import affinity
from shapely.geometry import LineString, Point, box

import gis

gdal.UseExceptions()

EPSG = 32633  # synthetic test data only
MOLLWEIDE = "ESRI:54009"
X0, Y_TOP = 500000.0, 5501600.0
WIDTH_M, HEIGHT_M, NATIVE = 2400, 1600, 1.0
RIVER_Y = 800.0

root = Path(
    os.environ.get("DEFENSES_TEST_DIR") or tempfile.mkdtemp(prefix="defenses_stage_b_")
)
data = root / "data"
data.mkdir(parents=True, exist_ok=True)
print(f"Test folder: {root}")


def world(geometry):
    return affinity.translate(geometry, X0, Y_TOP - HEIGHT_M)


def ring(xmin, ymin, xmax, ymax):
    return LineString(
        [(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax), (xmin, ymin)]
    )


levees = {
    "R1_gate": ring(300, 950, 700, 1250),
    "R2_pump": ring(1100, 950, 1500, 1250),
    "R3_plain": ring(1800, 950, 2200, 1250),
}
rivers = {
    "main": (LineString([(WIDTH_M, RIVER_Y), (0, RIVER_Y)]), 5000.0),
    "trib_gate": (LineString([(500, 1200), (500, RIVER_Y + 5)]), 15.0),
    # digitised against the flow: the downstream end must come from the DTM
    "trib_pump": (LineString([(1300, 960), (1300, 1200)]), 3.0),
    "trib_free": (LineString([(2300, 1500), (2300, RIVER_Y + 5)]), 40.0),
}
probes = {"R1_gate": (500, 1100), "R2_pump": (1200, 1050), "R3_plain": (2000, 1100)}

# --- Terrain: valley, channels for the tributaries, embankments on the rings ---
cols, rows = int(WIDTH_M / NATIVE), int(HEIGHT_M / NATIVE)
xx, yy = np.meshgrid(
    (np.arange(cols) + 0.5) * NATIVE, HEIGHT_M - (np.arange(rows) + 0.5) * NATIVE
)
dtm = (
    100.0
    + 0.001 * xx
    + 0.005 * np.abs(yy - RIVER_Y)
    - 2.0 * (np.abs(yy - RIVER_Y) <= 5)
).astype(np.float32)
for name, (line, _) in rivers.items():
    if name != "main":
        xmin, ymin, xmax, ymax = line.bounds
        dtm[(np.abs(xx - xmin) <= 2) & (yy >= ymin) & (yy <= ymax)] -= 1.0
for line in levees.values():
    xmin, ymin, xmax, ymax = line.bounds
    c0, c1 = int(max(xmin - 3, 0)), int(min(xmax + 3, WIDTH_M))
    r0, r1 = int(max(HEIGHT_M - ymax - 3, 0)), int(min(HEIGHT_M - ymin + 3, HEIGHT_M))
    sub = shapely.points(xx[r0:r1, c0:c1].ravel(), yy[r0:r1, c0:c1].ravel())
    dtm[r0:r1, c0:c1][shapely.dwithin(line, sub, 1.5).reshape(r1 - r0, c1 - c0)] += 2.5
srs = osr.SpatialReference()
srs.ImportFromEPSG(EPSG)
ds = gdal.GetDriverByName("GTiff").Create(
    str(data / "dtm.tif"), cols, rows, 1, gdal.GDT_Float32, ["COMPRESS=DEFLATE"]
)
ds.SetGeoTransform([X0, NATIVE, 0, Y_TOP, 0, -NATIVE])
ds.SetProjection(srs.ExportToWkt())
ds.GetRasterBand(1).SetNoDataValue(-9999.0)
ds.GetRasterBand(1).WriteArray(dtm)
ds = None

gpd.GeoDataFrame(
    {"name": list(levees)}, geometry=[world(g) for g in levees.values()], crs=EPSG
).to_file(data / "levees.gpkg", driver="GPKG")
gpd.GeoDataFrame(
    {"name": list(rivers), "upstream_km2": [a for _, a in rivers.values()]},
    geometry=[world(g) for g, _ in rivers.values()],
    crs=EPSG,
).to_file(data / "rivers.gpkg", driver="GPKG")
gpd.GeoDataFrame(
    {"rp": [100], "note": ["local knowledge"]},
    geometry=[world(Point(2000, 1100))],
    crs=EPSG,
).to_file(data / "overrides.gpkg", driver="GPKG")

# --- OSM extract (synthetic .osm served through a file:// template) ---
to_lonlat = Transformer.from_crs(EPSG, 4326, always_xy=True)
list_node, list_way = [], []


def node(x, y, tags=None):
    lon, lat = to_lonlat.transform(X0 + x, Y_TOP - HEIGHT_M + y)
    list_node.append((len(list_node) + 1, lat, lon, tags or {}))
    return len(list_node)


node(500, 953, {"waterway": "sluice_gate"})
node(1300, 955, {"man_made": "pumping_station", "name": "PS ring 2"})
node(2300, 1000, {"waterway": "weir"})
list_far = [node(100, 1490), node(110, 1490), node(110, 1500), node(100, 1500)]
list_way.append(
    (1000, list_far + [list_far[0]], {"man_made": "pumping_station", "building": "yes"})
)
list_way.append((1001, [node(900, 1400), node(1000, 1400)], {"waterway": "ditch"}))
list_xml = [
    "<?xml version='1.0' encoding='UTF-8'?>",
    '<osm version="0.6" generator="synthetic">',
]
for nid, lat, lon, tags in list_node:
    list_xml.append(
        f'  <node id="{nid}" lat="{lat:.8f}" lon="{lon:.8f}" version="1">'
        + "".join(f'<tag k="{k}" v="{v}"/>' for k, v in tags.items())
        + "</node>"
    )
for wid, list_ref, tags in list_way:
    list_xml.append(
        f'  <way id="{wid}" version="1">'
        + "".join(f'<nd ref="{r}"/>' for r in list_ref)
        + "".join(f'<tag k="{k}" v="{v}"/>' for k, v in tags.items())
        + "</way>"
    )
list_xml.append("</osm>")
(data / "test-region.osm").write_text("\n".join(list_xml))

# --- GHSL: two zipped Mollweide tiles for epoch 2020, zipped tiling schema ---
to_moll = Transformer.from_crs(EPSG, MOLLWEIDE, always_xy=True)
from_moll = Transformer.from_crs(MOLLWEIDE, EPSG, always_xy=True)
corners = np.array(
    [
        to_moll.transform(X0 + x, Y_TOP - HEIGHT_M + y)
        for x in (0, WIDTH_M)
        for y in (0, HEIGHT_M)
    ]
)
mx0, my0 = np.floor(corners.min(axis=0) / 100) * 100 - 500
mx1, my1 = np.ceil(corners.max(axis=0) / 100) * 100 + 500
split = np.round((mx0 + mx1) / 200) * 100
town = box(250, 900, 1550, 1300)  # local coordinates: buildings around rings 1 and 2
moll_srs = osr.SpatialReference()
moll_srs.SetFromUserInput(MOLLWEIDE)
tiles = {"R1_C1": (mx0, split), "R1_C2": (split, mx1)}
for tile, (tx0, tx1) in tiles.items():
    ncols, nrows = int(round((tx1 - tx0) / 100)), int(round((my1 - my0) / 100))
    cx, cy = np.meshgrid(
        tx0 + (np.arange(ncols) + 0.5) * 100, my1 - (np.arange(nrows) + 0.5) * 100
    )
    ux, uy = from_moll.transform(cx, cy)
    inside = shapely.contains_xy(world(town), ux, uy)
    for product, volume in [("BUILT_V", 30000.0), ("BUILT_V_NRES", 8000.0)]:
        name = f"GHS_{product}_E2020_TEST_{tile}"
        tif = data / f"{name}.tif"
        ds = gdal.GetDriverByName("GTiff").Create(
            str(tif), ncols, nrows, 1, gdal.GDT_Float32
        )
        ds.SetGeoTransform([tx0, 100.0, 0, my1, 0, -100.0])
        ds.SetProjection(moll_srs.ExportToWkt())
        ds.GetRasterBand(1).WriteArray(np.where(inside, volume, 0.0).astype(np.float32))
        ds = None
        with zipfile.ZipFile(data / f"{name}.zip", "w") as archive:
            archive.write(tif, tif.name)
        tif.unlink()
gdf_sch = gpd.GeoDataFrame(
    {"tile_id": list(tiles)},
    geometry=[box(t0, my0, t1, my1) for t0, t1 in tiles.values()],
    crs=MOLLWEIDE,
)
gdf_sch.to_file(data / "tiling_schema.shp", driver="ESRI Shapefile")
with zipfile.ZipFile(data / "tiling_schema.zip", "w") as archive:
    for part in data.glob("tiling_schema.*"):
        if part.suffix != ".zip":
            archive.write(part, part.name)
pd.DataFrame(
    {
        "country": ["T4", "XX"],
        "residential": [800.0, 1.0],
        "non_residential": [600.0, 1.0],
    }
).to_csv(data / "max_damage.csv", index=False)
pd.DataFrame(
    {
        "depth_m": [0, 0.5, 1, 2, 3, 6],
        "res": [0, 0.25, 0.4, 0.6, 0.75, 1.0],
        "nonres": [0, 0.2, 0.35, 0.55, 0.7, 1.0],
    }
).to_csv(data / "depth_damage.csv", index=False)

config = f"""placeholders = ["country.rp.max_rp"]
[general]
work_dir = "{root / 'work'}"
output_dir = "{root / 'output'}"
n_workers = 2
max_ram_gb_per_worker = 4.0
gdal_cachemax_mb = 256
warp_memory_limit_mb = 256
gdal_threads_per_worker = 1
log_level = "INFO"

[benchmark]
n_windows = 2
window_px = 256
sequential_rows = 32
seed = 1

[countries.T4]
dtm = "{data / 'dtm.tif'}"
dsm = []
building_footprints = []
water_polygons = []
levees = "{data / 'levees.gpkg'}"
rivers = "{data / 'rivers.gpkg'}"
river_upstream_area_field = "upstream_km2"
upstream_area_raster = ""
overrides = "{data / 'overrides.gpkg'}"
working_resolution_m = 2.0
corridor_buffer_m = 600.0
unit_core_size_m = 1200.0
unit_halo_m = 600.0
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
geofabrik_region = "test-region"
osm_pbf = ""
endpoint_elevation_radius_m = 3.0

[countries.T4.compartments]
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

[countries.T4.rp]
classes = [10, 20, 50, 100, 200, 500]
distribution_classes = [20, 50, 100, 200]
distribution_shares = [0.2, 0.3, 0.4, 0.1]
max_rp = 200
cap_levels_m = [0.5, 1.0, 5.6]
cap_rp = [10, 50, 100]
score_weight_value = 0.5
score_weight_upstream = 0.5
levee_adjacent_distance_m = 10.0
storey_height_m = 3.0
max_damage_csv = "{data / 'max_damage.csv'}"
depth_damage_csv = "{data / 'depth_damage.csv'}"
value_total_tolerance = 0.01

[countries.T4.lulc]
enabled = false

[countries.T4.context]
osm_industrial = {{ landuse = ["industrial"], man_made = ["works"] }}

[sources]
request_pause_s = 0.1
max_retries = 1
geofabrik_url_template = "file://{data}/{{region}}.osm"
ghsl_url_templates = {{ total = "file://{data}/GHS_BUILT_V_E{{epoch}}_TEST_{{tile}}.zip", non_residential = "file://{data}/GHS_BUILT_V_NRES_E{{epoch}}_TEST_{{tile}}.zip" }}
ghsl_tile_schema = "file://{data}/tiling_schema.zip"
ghsl_tile_field = "tile_id"
ghsl_epoch = 0
ghsl_nodata_as_zero = true
osm_tags_hydraulic = {{ waterway = ["sluice_gate", "lock_gate", "weir", "dam"], man_made = ["pumping_station"] }}
osm_tags_drains = {{ waterway = ["drain", "ditch", "canal"] }}
osm_tags_closing = {{ waterway = ["sluice_gate", "lock_gate"] }}
osm_tags_pumping = {{ man_made = ["pumping_station"] }}
osm_tags_water = {{ natural = ["water"], waterway = ["riverbank"], landuse = ["reservoir", "basin"] }}

[barriers]
osm_tags_walls = {{}}
osm_tags_crossings = {{ highway = ["residential", "service"], railway = ["rail"] }}
osm_wall_candidates = {{ barrier = ["wall"], wall = ["flood_wall"], man_made = ["dyke", "embankment"] }}
"""
(root / "test.toml").write_text(config)


def run(module, **variables):
    env = dict(
        os.environ,
        DEFENSES_CONFIG_PATH=str(root / "test.toml"),
        DEFENSES_COUNTRY="T4",
        **{f"DEFENSES_{k}": v for k, v in variables.items()},
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


def hits(frame, label):
    return int(frame.covers(world(Point(*probes[label]))).sum())


out = root / "output" / "T4"
list_res = []
for module in ["processing.prepare_dtm", "processing.build_units"]:
    proc = run(module, RUN_BENCHMARK="0")
    list_res.append((f"{module} runs", proc.returncode == 0))

proc = run("processing.compartments")
gdf_before = (
    gpd.read_file(out / "compartments.gpkg")
    if proc.returncode == 0
    else gpd.GeoDataFrame()
)
list_res.append(
    (
        "without OSM: tributaries inside rings 1 and 2 are water sources (no protection), ring 3 protected",
        proc.returncode == 0
        and hits(gdf_before, "R1_gate") == 0
        and hits(gdf_before, "R2_pump") == 0
        and hits(gdf_before, "R3_plain") == 1,
    )
)

proc = run("download.download_osm")
gdf_pts = (
    gpd.read_file(out / "osm_hydraulic_points.gpkg")
    if proc.returncode == 0
    else pd.DataFrame()
)
gdf_drain = (
    gpd.read_file(out / "osm_drains.gpkg") if proc.returncode == 0 else pd.DataFrame()
)
list_res.append(
    (
        f"download_osm: 4 hydraulic points (nodes and a way as point), 1 drain; found {len(gdf_pts)}, {len(gdf_drain)}",
        proc.returncode == 0
        and len(gdf_pts) == 4
        and len(gdf_drain) == 1
        and set(gdf_pts.get("osm_element", [])) == {"node", "area"},
    )
)
list_res.append(
    ("download_osm logs counts per tag value", "waterway=sluice_gate: 1" in proc.stderr)
)

proc = run("processing.osm_hydraulics")
gdf_src = (
    gpd.read_file(out / "rivers_sources.gpkg").set_index("name")
    if proc.returncode == 0
    else None
)
expected_sources = {
    "main": True,
    "trib_gate": False,
    "trib_pump": False,
    "trib_free": True,
}
list_res.append(
    (
        f"osm_hydraulics: water sources {expected_sources}",
        gdf_src is not None
        and all(
            bool(gdf_src.loc[k, "is_source"]) == v for k, v in expected_sources.items()
        ),
    )
)
if gdf_src is not None:
    list_res.append(
        (
            "closed by a gate at the levee crossing",
            "closing structure" in gdf_src.loc["trib_gate", "reason"],
        )
    )
    pump_end = shapely.Point(
        gdf_src.loc["trib_pump", "downstream_x"],
        gdf_src.loc["trib_pump", "downstream_y"],
    )
    list_res.append(
        (
            "downstream end from the DTM (line digitised against the flow) next to the pumping station",
            "pumping station" in gdf_src.loc["trib_pump", "reason"]
            and gdf_src.loc["trib_pump", "downstream_method"] == "dtm elevation"
            and pump_end.distance(world(Point(1300, 960))) < 1.0,
        )
    )
gdf_pump = (
    gpd.read_file(out / "osm_pumping.gpkg") if proc.returncode == 0 else pd.DataFrame()
)
list_res.append(
    ("osm_pumping: only the pumping station next to a levee", len(gdf_pump) == 1)
)

proc = run("processing.compartments", FORCE="1")
gdf_after = (
    gpd.read_file(out / "compartments.gpkg")
    if proc.returncode == 0
    else gpd.GeoDataFrame()
)
list_res.append(
    (
        "with closed tributaries all three rings are protected; ring 2 has pumping",
        proc.returncode == 0
        and all(hits(gdf_after, k) == 1 for k in probes)
        and bool(
            gdf_after[
                gdf_after.covers(world(Point(*probes["R2_pump"])))
            ].has_pumping.iloc[0]
        )
        and not bool(
            gdf_after[
                gdf_after.covers(world(Point(*probes["R1_gate"])))
            ].has_pumping.iloc[0]
        ),
    )
)

proc = run("download.download_ghsl")
record = (
    json.loads((out / "ghsl_sources.json").read_text()) if proc.returncode == 0 else {}
)
list_res.append(
    (
        f"download_ghsl: latest available epoch found ({record.get('epoch')}), two tiles per product mosaicked",
        proc.returncode == 0
        and record.get("epoch") == 2020
        and all(len(v["sources"]) == 2 for v in record.get("products", {}).values()),
    )
)

proc = run("processing.rp_index")
list_ratio = [
    r
    for d in (root / "work" / "T4" / "units").iterdir()
    if (d / "rp_index.json").exists()
    for r in json.loads((d / "rp_index.json").read_text())["total_ratio"].values()
]
list_res.append(
    (
        f"value warp preserves totals within 1 % (ratios {', '.join(f'{r:.4f}' for r in list_ratio)})",
        proc.returncode == 0
        and list_ratio
        and all(abs(r - 1) <= 0.01 for r in list_ratio),
    )
)
gdf_rp = (
    gpd.read_file(out / "compartments_rp.gpkg")
    if proc.returncode == 0
    else gpd.GeoDataFrame()
)
if len(gdf_rp) and all(
    int(gdf_rp.covers(world(Point(*xy))).sum()) == 1 for xy in probes.values()
):
    list_cls = [10, 20, 50, 100, 200, 500]
    by_probe = {
        k: gdf_rp[gdf_rp.covers(world(Point(*xy)))].iloc[0] for k, xy in probes.items()
    }
    list_col = [
        "v_rel",
        "levee_length_m",
        "upstream_area_km2",
        "score",
        "rp_model",
        "rp_cap",
        "rp",
        "rp_low",
        "rp_high",
        "rp_source",
        "obs_rp_low",
        "obs_rp_high",
        "n_held",
        "n_flooded",
        "obs_conflict",
        "model_obs_conflict",
    ]
    list_res.append(
        ("compartments_rp has all columns", all(c in gdf_rp.columns for c in list_col))
    )
    list_res.append(
        (
            "rings 1 and 2 lie in the town (exposed value > 0), ring 3 has none",
            by_probe["R1_gate"].value_exposed > 0
            and by_probe["R2_pump"].value_exposed > 0
            and by_probe["R3_plain"].value_exposed == 0,
        )
    )
    gdf_model = gdf_rp[gdf_rp.rp_source != "override"]
    list_cap = [
        np.nan if i >= 3 else [10, 50, 100][i]
        for i in np.searchsorted(
            [0.5, 1.0, 5.6], gdf_model.protection_level_m, side="right"
        )
    ]
    list_res.append(
        (
            "cap from protection level; rp = min(model, cap, max_rp) and a class value",
            np.allclose(gdf_model.rp_cap.to_numpy(), list_cap, equal_nan=True)
            and all(r in list_cls for r in gdf_model.rp)
            and (
                gdf_model.rp
                <= np.fmin(
                    np.fmin(
                        gdf_model.rp_model.fillna(1e9), gdf_model.rp_cap.fillna(1e9)
                    ),
                    200,
                )
            ).all(),
        )
    )
    list_res.append(
        (
            "interval rp_low <= rp <= rp_high, bounded by cap and max_rp",
            (gdf_model.rp_low <= gdf_model.rp).all()
            and (gdf_model.rp <= gdf_model.rp_high).all()
            and (gdf_model.rp_high <= 200).all(),
        )
    )
    list_res.append(
        (
            "override sets rp, rp_low and rp_high of ring 3",
            by_probe["R3_plain"].rp_source == "override"
            and by_probe["R3_plain"].rp
            == by_probe["R3_plain"].rp_low
            == by_probe["R3_plain"].rp_high
            == 100,
        )
    )
    list_res.append(
        (
            "v_rel = value / levee length",
            np.allclose(
                gdf_model.v_rel, gdf_model.value_exposed / gdf_model.levee_length_m
            ),
        )
    )

# Quantile mapping and percentile ranks on many scores.
scores = np.random.default_rng(3).random(997)
mapped = gis.Stats.quantile_classes(scores, [20, 50, 100, 200], [0.2, 0.3, 0.4, 0.1])
list_share = [float(np.mean(mapped == c)) for c in [20, 50, 100, 200]]
list_res.append(
    (
        f"quantile mapping reproduces the shares within rounding ({list_share})",
        all(abs(s - t) <= 1 / 997 for s, t in zip(list_share, [0.2, 0.3, 0.4, 0.1])),
    )
)
ties = gis.Stats.quantile_classes([1.0, 1.0, 2.0, 3.0], [20, 50], [0.5, 0.5])
list_res.append(("ties get the same class", ties[0] == ties[1]))

if not (
    len(gdf_rp)
    and all(int(gdf_rp.covers(world(Point(*xy))).sum()) == 1 for xy in probes.values())
):
    list_res.append(("compartments_rp contains the three rings", False))
for label, ok in list_res:
    print(f"{'PASS' if ok else 'FAIL'}  {label}")
sys.exit(0 if all(ok for _, ok in list_res) else 1)
