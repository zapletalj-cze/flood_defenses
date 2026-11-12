# Author: Jakub Zapletal
# Date: 2026-11-04
# End-to-end check of Stage B on a synthetic country: tributaries inside levee rings open (gated_crossings off) and
# closed (on), GHSL download and mosaic, value warp totals, exposure, score, quantile mapping, cap, interval and
# overrides.
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
    # drainage ditch inside ring 1, closed by the levee (it ends at the levee)
    "trib_gate": (LineString([(500, 1200), (500, 958)]), 15.0),
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

data_url = data.as_uri()
config = f"""placeholders = ["country.rp.max_rp"]
[general]
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
sequential_rows = 32
seed = 1

[countries.T4]
dtm = {json.dumps(str(data / 'dtm.tif'))}
dsm = []
building_footprints = []
water_polygons = []
levees = {json.dumps(str(data / 'levees.gpkg'))}
rivers = {json.dumps(str(data / 'rivers.gpkg'))}
river_upstream_area_field = "upstream_km2"
upstream_area_raster = ""
overrides = {json.dumps(str(data / 'overrides.gpkg'))}
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

[countries.T4.protected_areas]
gap_tolerance_m = 20.0
touch_tolerance_m = 0.01
self_close_min_length_factor = 3.0
crest_search_radius_m = 5.0
barrier_removal_width_m = 10.0
min_margin_m = 0.3
max_level_m = 20.0
min_area_m2 = 10000.0
gated_crossings = true
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
max_damage_csv = {json.dumps(str(data / 'max_damage.csv'))}
depth_damage_csv = {json.dumps(str(data / 'depth_damage.csv'))}
value_total_tolerance = 0.01

[countries.T4.lulc]
enabled = false

[sources]
request_pause_s = 0.1
max_retries = 1
ghsl_url_templates = {{ total = {json.dumps(data_url + "/GHS_BUILT_V_E{epoch}_TEST_{tile}.zip")}, non_residential = {json.dumps(data_url + "/GHS_BUILT_V_NRES_E{epoch}_TEST_{tile}.zip")} }}
ghsl_tile_schema = {json.dumps(data_url + "/tiling_schema.zip")}
ghsl_tile_field = "tile_id"
ghsl_epoch = 0
ghsl_nodata_as_zero = true

"""
(root / "test.toml").write_text(config)
(root / "test_open.toml").write_text(
    config.replace("gated_crossings = true", "gated_crossings = false")
)


def run(module, config="test.toml", **variables):
    env = dict(
        os.environ,
        DEFENSES_CONFIG_PATH=str(root / config),
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

proc = run("processing.protected_areas", config="test_open.toml")
gdf_before = (
    gpd.read_file(out / "protected_areas.gpkg")
    if proc.returncode == 0
    else gpd.GeoDataFrame()
)
list_res.append(
    (
        "gated_crossings off: tributaries inside rings 1 and 2 are water sources (no protection), ring 3 protected",
        proc.returncode == 0
        and hits(gdf_before, "R1_gate") == 0
        and hits(gdf_before, "R2_pump") == 0
        and hits(gdf_before, "R3_plain") == 1,
    )
)
gdf_leak = (
    gpd.read_file(out / "leaks.gpkg")
    if (out / "leaks.gpkg").exists()
    else gpd.GeoDataFrame(columns=["cause", "geometry"], geometry="geometry")
)


def leaks_in(name):
    inside = gdf_leak.within(world(shapely.Polygon(levees[name].coords)))
    return sorted(gdf_leak.loc[inside, "cause"])


list_res.append(
    (
        "leaks.gpkg, crossings open: rings 1 and 2 leak through the stream inside them (water_source_behind), "
        f"ring 3 does not ({leaks_in('R1_gate')}, {leaks_in('R2_pump')}, {leaks_in('R3_plain')})",
        leaks_in("R1_gate") == ["water_source_behind"]
        and leaks_in("R2_pump") == ["water_source_behind"]
        and leaks_in("R3_plain") == [],
    )
)
gdf_byg = (
    gpd.read_file(out / "protected_areas_by_group.gpkg")
    if (out / "protected_areas_by_group.gpkg").exists()
    else gpd.GeoDataFrame(columns=["group_id", "geometry"], geometry="geometry")
)
gdf_rg = (
    gpd.read_file(out / "river_groups.gpkg")
    if (out / "river_groups.gpkg").exists()
    else None
)


def hits_group(frame, label, group):
    return int(
        frame[frame["group_id"] == group].covers(world(Point(*probes[label]))).sum()
    )


list_res.append(
    (
        "per river group (protected_areas_by_group.gpkg), crossings open: the main river (group 1) is held back by "
        "all three rings, the small streams together (group 0) only by ring 3",
        len(gdf_byg) > 0
        and gdf_rg is not None
        and sorted(gdf_rg["group_id"]) == [0] * (len(rivers) - 1) + [1]
        and all(hits_group(gdf_byg, k, 1) == 1 for k in probes)
        and hits_group(gdf_byg, "R3_plain", 0) == 1
        and hits_group(gdf_byg, "R1_gate", 0) == 0,
    )
)
gdf_channel = (
    gpd.read_file(out / "rivers_channel.gpkg")
    if (out / "rivers_channel.gpkg").exists()
    else None
)
list_res.append(
    (
        "rivers snapped once per country (rivers_channel.gpkg; without guide lines the lines stay where they are)",
        gdf_channel is not None
        and len(gdf_channel) == len(rivers)
        and all(
            g.equals(world(rivers[n][0]))
            for n, g in zip(gdf_channel["name"], gdf_channel.geometry)
        ),
    )
)

proc = run("processing.protected_areas", FORCE="1")
gdf_after = (
    gpd.read_file(out / "protected_areas.gpkg")
    if proc.returncode == 0
    else gpd.GeoDataFrame()
)
list_res.append(
    (
        "gated_crossings on: tributaries enclosed by the rings are closed, all three rings are protected",
        proc.returncode == 0 and all(hits(gdf_after, k) == 1 for k in probes),
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
    gpd.read_file(out / "protected_areas_rp.gpkg")
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
        (
            "protected_areas_rp has all columns",
            all(c in gdf_rp.columns for c in list_col),
        )
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
    list_res.append(("protected_areas_rp contains the three rings", False))
# File B downstream: rp_index with LAYER = B values the protected areas of every river group with the labels and HAND
# of that group (gated crossings: all three rings in group 1)
proc = run("processing.rp_index", LAYER="B")
gdf_b = gpd.read_file(out / "protected_areas_by_group.gpkg")
gdf_rpb = (
    gpd.read_file(out / "protected_areas_by_group_rp.gpkg")
    if (out / "protected_areas_by_group_rp.gpkg").exists()
    else gpd.GeoDataFrame(
        columns=["group_id", "value_exposed", "rp", "geometry"], geometry="geometry"
    )
)
rasters_ok = all(
    Path(gis.Run.unit_layer_path(root / "work" / "T4" / "units" / u, name, g)).exists()
    for u, g in gdf_b[["unit_id", "group_id"]].drop_duplicates().itertuples(index=False)
    for name in ("labels.tif", "hand_undef.tif")
)


def value_b(label, group):
    hit = gdf_rpb[
        (gdf_rpb["group_id"] == group) & gdf_rpb.covers(world(Point(*probes[label])))
    ]
    return float(hit["value_exposed"].iloc[0]) if len(hit) else float("nan")


list_res.append(
    (
        f"rp_index LAYER=B: {len(gdf_rpb)} of {len(gdf_b)} protected areas per river group valued with the rasters "
        f"of their group (group rasters {'present' if rasters_ok else 'missing'}); rings 1 and 2 in the town have "
        f"value ({value_b('R1_gate', 1):.3g}, {value_b('R2_pump', 1):.3g}), ring 3 none ({value_b('R3_plain', 1):.3g})",
        proc.returncode == 0
        and rasters_ok
        and len(gdf_rpb) == len(gdf_b) > 0
        and gdf_rpb["value_exposed"].notna().all()
        and gdf_rpb["rp"].notna().all()
        and value_b("R1_gate", 1) > 0
        and value_b("R2_pump", 1) > 0
        and value_b("R3_plain", 1) == 0,
    )
)

for label, ok in list_res:
    print(f"{'PASS' if ok else 'FAIL'}  {label}")
sys.exit(0 if all(ok for _, ok in list_res) else 1)
