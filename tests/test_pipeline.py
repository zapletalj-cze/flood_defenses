# Author: Jakub Zapletal
# Date: 2026-11-04
# End-to-end check of prepare_dtm.py, build_units.py and compartments.py on a synthetic country: closed ring across
# a core boundary (seam), bridged and wide gaps, diagonal levee, nested rings, compartment wider than a window,
# nodata hole, empty levee input, DTM conversion options A and B, unit selection and checkpointing.
# Run: pixi run python -m tests.test_pipeline  (DEFENSES_TEST_DIR sets the working folder, default: temporary)

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

import gis

gdal.UseExceptions()

EPSG = 32633  # synthetic test data only
X0, Y_TOP = 500000.0, 5502000.0
WIDTH_M, HEIGHT_M, NATIVE = 3000, 2000, 1.0
RIVER_Y = 1000.0

root = Path(
    os.environ.get("DEFENSES_TEST_DIR") or tempfile.mkdtemp(prefix="defenses_pipeline_")
)
data = root / "data"
data.mkdir(parents=True, exist_ok=True)
print(f"Test folder: {root}")


def world(geometry):
    # Local coordinates (origin lower-left, metres) to map coordinates.
    return affinity.translate(geometry, X0, Y_TOP - HEIGHT_M)


def ring(xmin, ymin, xmax, ymax):
    return LineString(
        [(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax), (xmin, ymin)]
    )


# Levees: (geometry in local coordinates, crest height above ground). The river runs along y = 1000; every
# feature lies within the 500 m corridor. Cores are 1000 m, the halo 500 m (= corridor buffer).
levees = {
    "A_seam": (ring(870, 1150, 1170, 1350), 2.5),
    "B_bridged_1": (
        LineString([(1700, 1245), (1700, 1150), (1900, 1150), (1900, 1245)]),
        2.5,
    ),
    "B_bridged_2": (
        LineString([(1900, 1255), (1900, 1350), (1700, 1350), (1700, 1255)]),
        2.5,
    ),
    "C_wide_1": (
        LineString([(2300, 1230), (2300, 1150), (2500, 1150), (2500, 1230)]),
        2.5,
    ),
    "C_wide_2": (
        LineString([(2500, 1270), (2500, 1350), (2300, 1350), (2300, 1270)]),
        2.5,
    ),
    "D_diamond": (
        LineString([(281, 1251), (401, 1371), (521, 1251), (401, 1131), (281, 1251)]),
        2.5,
    ),
    "E_outer": (ring(1250, 600, 1600, 900), 2.0),
    "E_inner": (ring(1350, 680, 1500, 830), 3.0),
    "F_outer": (ring(1900, 600, 2250, 900), 3.0),
    "F_inner": (ring(2000, 680, 2150, 830), 1.5),
    "G_long": (ring(100, 1400, 2900, 1480), 2.5),
    # B2: levee on the river bank; its removal zone reaches the river cells.
    "H_bank": (LineString([(100, 1008), (500, 1008)]), 2.5),
    # B3/B4: U-shaped levee open towards a strip without DTM coverage (x >= 2960).
    "U_open": (LineString([(2960, 600), (2700, 600), (2700, 800), (2960, 800)]), 2.5),
}
probes = {  # local point -> expected number of compartments covering it
    "A_seam": ((1020, 1250), 1),
    "B_bridged": ((1800, 1200), 1),
    "C_wide": ((2400, 1250), 0),
    "D_diamond": ((401, 1251), 1),
    "E_inner": ((1425, 755), 1),
    "E_outer": ((1290, 880), 1),
    "F_inner": ((2075, 755), 1),
    "F_outer": ((1950, 880), 1),
    "G_long": ((1500, 1440), 1),
    "U_open": ((2830, 700), 1),
    "L1_hole_artifact": ((2700, 1550), 0),
}

# Terrain: valley rising 0.5 % away from the river and 0.1 % downstream to upstream, 2 m deep channel,
# embankments 3 m wide on every levee line; a nodata hole over the channel (lidar water surface).
cols, rows = int(WIDTH_M / NATIVE), int(HEIGHT_M / NATIVE)
x = (np.arange(cols) + 0.5) * NATIVE
y = HEIGHT_M - (np.arange(rows) + 0.5) * NATIVE
xx, yy = np.meshgrid(x, y)
dtm = (
    100.0
    + 0.001 * xx
    + 0.005 * np.abs(yy - RIVER_Y)
    - 2.0 * (np.abs(yy - RIVER_Y) <= 5)
).astype(np.float32)
for name, (line, height) in levees.items():
    xmin, ymin, xmax, ymax = line.bounds
    c0, c1 = int(max(xmin - 3, 0)), int(min(xmax + 3, WIDTH_M))
    r0, r1 = int(max(HEIGHT_M - ymax - 3, 0)), int(min(HEIGHT_M - ymin + 3, HEIGHT_M))
    sub = shapely.points(xx[r0:r1, c0:c1].ravel(), yy[r0:r1, c0:c1].ravel())
    mask = shapely.dwithin(line, sub, 1.5).reshape(r1 - r0, c1 - c0)
    dtm[r0:r1, c0:c1][mask] += height
# L1: water-surface hole over the river; interpolation from the banks would raise the river cells.
dtm[(np.abs(yy - RIVER_Y) <= 4) & (xx > 2600) & (xx < 2800)] = -9999.0
# B3/B4: no DTM coverage east of x = 2960.
dtm[xx >= 2960] = -9999.0

srs = osr.SpatialReference()
srs.ImportFromEPSG(EPSG)
dtm_files, dtm_files_dry = [], []
# T3: the whole river channel without DTM data (lidar water surface along the entire river).
dtm_dry = dtm.copy()
dtm_dry[np.abs(yy - RIVER_Y) <= 4] = -9999.0
for name, array, files in [
    ("dtm_block", dtm, dtm_files),
    ("dtm_water", dtm_dry, dtm_files_dry),
]:
    for k, (c0, c1) in enumerate([(0, cols // 2), (cols // 2, cols)]):
        path = data / f"{name}_{k}.tif"
        ds = gdal.GetDriverByName("GTiff").Create(
            str(path),
            c1 - c0,
            rows,
            1,
            gdal.GDT_Float32,
            ["COMPRESS=DEFLATE", "PREDICTOR=3"],
        )  # strip layout, no overviews
        ds.SetGeoTransform([X0 + c0 * NATIVE, NATIVE, 0, Y_TOP, 0, -NATIVE])
        ds.SetProjection(srs.ExportToWkt())
        ds.GetRasterBand(1).SetNoDataValue(-9999.0)
        ds.GetRasterBand(1).WriteArray(array[:, c0:c1])
        ds = None
        files.append(str(path))

levee_names = list(levees)
half = len(levee_names) // 2
for file_name, names in [
    ("levees_ml.gpkg", levee_names[:half]),
    ("levees_osm.gpkg", levee_names[half:]),
]:
    gpd.GeoDataFrame(
        {"name": names}, geometry=[world(levees[n][0]) for n in names], crs=EPSG
    ).to_file(data / file_name, driver="GPKG")
gpd.GeoDataFrame({"name": []}, geometry=gpd.GeoSeries([], crs=EPSG)).to_file(
    data / "levees_empty.gpkg", driver="GPKG"
)
river = LineString([(WIDTH_M, RIVER_Y), (0, RIVER_Y)])
gpd.GeoDataFrame({"name": ["main"]}, geometry=[world(river)], crs=EPSG).to_file(
    data / "rivers.gpkg", driver="GPKG"
)
# B4: tributary entering from the area without DTM coverage into the U-shaped levee.
gpd.GeoDataFrame(
    {"name": ["trib_east"]},
    geometry=[world(LineString([(2995, 700), (2750, 700)]))],
    crs=EPSG,
).to_file(data / "rivers_trib.gpkg", driver="GPKG")


def country_table(
    code, levee_files, convert, option, river_files=None, gated=False, dtm_list=None
):
    river_files = river_files or [str(data / "rivers.gpkg")]
    return f"""
[countries.{code}]
dtm = {json.dumps(dtm_list or dtm_files)}
dsm = []
levees = {json.dumps(levee_files)}
rivers = {json.dumps(river_files)}
working_resolution_m = 2.0
corridor_buffer_m = 500.0
unit_core_size_m = 1000.0
unit_halo_m = 500.0
convert_dtm = {str(convert).lower()}
dtm_conversion = "{option}"
dtm_compress = "ZSTD"
footprint_resolution_m = 20.0
bytes_per_cell = 400.0
crest_box_size_m = 100.0
merge_overlap_ratio = 0.5
drop_contradicted_truncated = true
contradiction_ratio = 0.5
id_rounding_m = 10.0
keep_intermediate = true

[countries.{code}.compartments]
gap_tolerance_m = 20.0
touch_tolerance_m = 0.01
self_close_min_length_factor = 3.0
crest_search_radius_m = 5.0
barrier_removal_width_m = 10.0
min_margin_m = 0.3
max_level_m = 20.0
min_area_m2 = 10000.0
gated_crossings = {str(gated).lower()}
fill_nodata_holes = true
stream_burn_m = 5.0
stream_burn_buffer_m = 20.0
preserve_river_cells = true
cap_river_holes = true
max_fill_distance_m = 100.0

[countries.{code}.rp]
levee_adjacent_distance_m = 10.0
"""


def write_config(convert):
    text = f"""placeholders = ["country.bytes_per_cell"]
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
n_windows = 3
window_px = 512
sequential_rows = 64
seed = 1
"""
    text += country_table(
        "T1",
        [str(data / "levees_ml.gpkg"), str(data / "levees_osm.gpkg")],
        convert,
        "B",
    )
    text += country_table("T0", str(data / "levees_empty.gpkg"), convert, "A")
    text += country_table(
        "T2",
        [str(data / "levees_ml.gpkg"), str(data / "levees_osm.gpkg")],
        False,
        "A",
        river_files=[str(data / "rivers.gpkg"), str(data / "rivers_trib.gpkg")],
        gated=True,
    )
    text += country_table(
        "T3",
        [str(data / "levees_ml.gpkg"), str(data / "levees_osm.gpkg")],
        False,
        "A",
        dtm_list=dtm_files_dry,
    )
    (root / "test.toml").write_text(text)


def run(script, country, **variables):
    env = dict(
        os.environ,
        DEFENSES_CONFIG_PATH=str(root / "test.toml"),
        DEFENSES_COUNTRY=country,
        **{f"DEFENSES_{k}": v for k, v in variables.items()},
    )
    proc = subprocess.run(
        [sys.executable, "-m", f"processing.{script}"],
        env=env,
        capture_output=True,
        text=True,
        cwd=Path(__file__).resolve().parents[1],
    )
    (root / f"{script}_{country}.log").write_text(proc.stdout + proc.stderr)
    if proc.returncode != 0:
        print(proc.stderr[-3000:])
    return proc


results = []
write_config(convert=False)
for country in ["T1", "T0", "T2", "T3"]:
    proc = run("prepare_dtm", country)
    results.append(
        (
            f"{country}: prepare_dtm (inspection, footprint, benchmark)",
            proc.returncode == 0,
        )
    )
    proc = run("build_units", country)
    results.append((f"{country}: build_units", proc.returncode == 0))
write_config(convert=True)
for country, option in [("T1", "B"), ("T0", "A")]:
    proc = run("prepare_dtm", country, RUN_BENCHMARK="0")
    converted = list((root / "work" / country / "dtm").glob(f"*_{option}.tif"))
    info = gdal.Info(str(converted[0]), format="json") if converted else {}
    band = info.get("bands", [{}])[0]
    results.append(
        (
            f"{country}: conversion option {option} gives tiled files with 5 overviews",
            proc.returncode == 0
            and len(converted) == 2
            and band.get("block") == [512, 512]
            and len(band.get("overviews", [])) == 5,
        )
    )
print(
    open(root / "prepare_dtm_T1.log").read().split("Overview use")[1][:200]
    if "Overview use" in open(root / "prepare_dtm_T1.log").read()
    else "no overview message"
)

# Compartments: first two units only, then the rest (checkpoints), then a rerun that must skip everything.
units = sorted(gpd.read_file(root / "output" / "T1" / "units_core.gpkg").unit_id)
proc = run("compartments", "T1", UNITS=",".join(units[:2]))
results.append(("T1: compartments on two selected units", proc.returncode == 0))
proc = run("compartments", "T1")
results.append(
    (
        "T1: remaining units processed, finished units skipped",
        proc.returncode == 0
        and f"{len(units) - 2} to process, 2 already finished" in proc.stderr,
    )
)
proc = run("compartments", "T1")
results.append(
    (
        "T1: rerun skips every unit",
        proc.returncode == 0 and "0 to process" in proc.stderr,
    )
)

final = gpd.read_file(root / "output" / "T1" / "compartments.gpkg")
for label, ((px, py), expected) in probes.items():
    hits = int(final.covers(world(Point(px, py))).sum())
    results.append(
        (
            f"{label}: {expected} compartment(s) at the probe, found {hits}",
            hits == expected,
        )
    )
results.append((f"T1: total compartments 9, found {len(final)}", len(final) == 9))
u_open = final[final.covers(world(Point(2830, 700)))]
results.append(
    (
        "B3: compartment closed by missing DTM coverage kept and flagged touches_dtm_edge",
        len(u_open) == 1
        and bool(u_open.touches_dtm_edge.iloc[0])
        and not bool(u_open.truncated.iloc[0]),
    )
)
others_flagged = int(final.touches_dtm_edge.astype(bool).sum()) - len(u_open)
results.append(
    ("B3: no other compartment flagged touches_dtm_edge", others_flagged == 0)
)

# B2: river cells under the removal zone of the bank levee keep their elevation.
unit = root / "work" / "T1" / "units" / "000_000"
ds = gdal.Open(str(unit / "dtm.tif"))
unit_gt, unit_dtm = ds.GetGeoTransform(), ds.ReadAsArray()
unit_undef = gdal.Open(str(unit / "dtm_undefended.tif")).ReadAsArray()
unit_stream = (
    gis.Raster.rasterize_geometries(
        [world(river)],
        [1],
        unit_dtm.shape,
        unit_gt,
        ds.GetProjection(),
        False,
        np.uint8,
    )
    > 0
)
cell_x = unit_gt[0] + (np.arange(unit_dtm.shape[1]) + 0.5) * unit_gt[1] - X0
cell_y = (
    unit_gt[3] + (np.arange(unit_dtm.shape[0]) + 0.5) * unit_gt[5] - (Y_TOP - HEIGHT_M)
)
along_h = (cell_x > 150) & (cell_x < 450)
river_mask = unit_stream & along_h[None, :]
between = (np.abs(cell_y - 1004) < 1.5)[:, None] & along_h[None, :]
results.append(
    (
        "B2: river cells next to a levee keep their elevation; cells beside them are interpolated",
        river_mask.any()
        and float(np.nanmax(np.abs(unit_undef[river_mask] - unit_dtm[river_mask])))
        < 1e-6
        and float(np.nanmax(np.abs(unit_undef[between] - unit_dtm[between]))) > 0.01,
    )
)
along = sum(
    json.loads((d / "done.json").read_text()).get("river_cells_interpolated_along", 0)
    for d in (root / "work" / "T1" / "units").iterdir()
)
hole_unit = root / "work" / "T1" / "units" / "002_000"
ds = gdal.Open(str(hole_unit / "dtm.tif"))
hole_gt, hole_dtm = ds.GetGeoTransform(), ds.ReadAsArray()
hole_stream = (
    gis.Raster.rasterize_geometries(
        [world(river)],
        [1],
        hole_dtm.shape,
        hole_gt,
        ds.GetProjection(),
        False,
        np.uint8,
    )
    > 0
)
hx = hole_gt[0] + (np.arange(hole_dtm.shape[1]) + 0.5) * hole_gt[1] - X0
hy = hole_gt[3] + (np.arange(hole_dtm.shape[0]) + 0.5) * hole_gt[5] - (Y_TOP - HEIGHT_M)
hxx, hyy = np.meshgrid(hx, hy)
in_hole = hole_stream & (hxx > 2610) & (hxx < 2790)
bed = 100.0 + 0.001 * hxx - 2.0
deviation = (
    float(np.abs(hole_dtm[in_hole] - bed[in_hole]).max()) if in_hole.any() else 99.0
)
results.append(
    (
        f"L1: river cells in the gap follow the river profile ({along} cells, max deviation "
        f"from the bed {deviation:.2f} m)",
        along > 0 and deviation < 0.1,
    )
)
long = final[final.covers(world(Point(1500, 1440)))]
results.append(
    (
        "G_long: truncated pieces merged into one compartment spanning the valley",
        len(long) == 1
        and bool(long.truncated.iloc[0])
        and long.geometry.iloc[0].length > 5000,
    )
)
results.append(("compartment ids are unique", final.compartment_id.is_unique))
river_world = world(river)
source_ok = (
    shapely.distance(shapely.points(final.source_x, final.source_y), river_world) <= 4.0
)
entry_ok = (
    shapely.distance(
        shapely.points(final.entry_x, final.entry_y), final.geometry.values
    )
    <= 4.0
)
results.append(
    (
        "source cells lie on the river, entry cells on their compartment",
        source_ok.all() and entry_ok.all(),
    )
)

barriers = gpd.read_file(root / "output" / "T1" / "barriers.gpkg")
crest_hand = {
    name: barriers[
        barriers.intersects(world(levees[name][0]).buffer(0.5))
        & ~barriers.bridged.astype(bool)
    ].crest_hand_min_m.min()
    for name in ["E_outer", "E_inner", "F_outer", "F_inner"]
}
for pair in ["E", "F"]:
    inner = final[final.covers(world(Point(*probes[f"{pair}_inner"][0])))]
    expected = max(crest_hand[f"{pair}_outer"], crest_hand[f"{pair}_inner"])
    level = float(inner.protection_level_m.iloc[0]) if len(inner) else float("nan")
    results.append(
        (
            f"nested rings {pair}: inner protection {level:.2f} m = higher ring crest {expected:.2f} m",
            abs(level - expected) <= 0.05,
        )
    )
results.append(
    (
        "bridged gaps recorded in barriers.gpkg",
        int(barriers.bridged.astype(bool).sum()) >= 2,
    )
)

unit_dirs = sorted((root / "work" / "T1" / "units").iterdir())
diagonal = 0
for unit_dir in unit_dirs:
    barrier = gdal.Open(str(unit_dir / "barrier.tif")).ReadAsArray() == 1
    diagonal += int(
        (
            (
                barrier[:-1, :-1]
                & barrier[1:, 1:]
                & ~barrier[:-1, 1:]
                & ~barrier[1:, :-1]
            )
            | (
                barrier[:-1, 1:]
                & barrier[1:, :-1]
                & ~barrier[:-1, :-1]
                & ~barrier[1:, 1:]
            )
        ).sum()
    )
results.append(
    (
        f"barrier rasters are 4-connected (diagonal-only contacts: {diagonal})",
        diagonal == 0,
    )
)
filled = sum(
    json.loads((d / "done.json").read_text()).get("filled_gap_cells", 0)
    for d in unit_dirs
)
results.append(
    (f"nodata hole over the river interpolated ({filled} cells)", filled > 0)
)

# B4: in gated mode a tributary entering from missing DTM coverage is a water source.
proc = run("compartments", "T2")
gated_final = gpd.read_file(root / "output" / "T2" / "compartments.gpkg")
u_gated = int(gated_final.covers(world(Point(2830, 700))).sum())
others_same = all(
    int(gated_final.covers(world(Point(*xy))).sum()) == expected
    for label, (xy, expected) in probes.items()
    if label != "U_open"
)
results.append(
    (
        f"B4: tributary from missing coverage floods the U in gated mode ({u_gated} compartments "
        f"at the U, other probes unchanged: {others_same})",
        proc.returncode == 0 and u_gated == 0 and others_same,
    )
)

# T3: the whole river is nodata; gaps are filled and the river cells capped at the waterline.
proc = run("compartments", "T3")
water_final = gpd.read_file(root / "output" / "T3" / "compartments.gpkg")
water_hits = {
    label: int(water_final.covers(world(Point(*xy))).sum())
    for label, (xy, expected) in probes.items()
}
water_ok = all(
    water_hits[label] == expected for label, (xy, expected) in probes.items()
)
results.append(
    (
        f"T3 river entirely nodata: same compartments as T1 ({len(water_final)}; probes "
        f"{'ok' if water_ok else water_hits})",
        proc.returncode == 0 and water_ok and len(water_final) == 9,
    )
)

proc = run("compartments", "T0")
margins = [
    gdal.Open(str(d / "margin.tif")).ReadAsArray()
    for d in (root / "work" / "T0" / "units").iterdir()
]
max_margin = max(
    float(np.nanmax(np.abs(m))) if np.isfinite(m).any() else 0.0 for m in margins
)
empty_final = gpd.read_file(root / "output" / "T0" / "compartments.gpkg")
results.append(
    (
        f"T0 empty levees: margin zero everywhere (max |margin| {max_margin}), no compartments",
        proc.returncode == 0 and max_margin == 0.0 and len(empty_final) == 0,
    )
)

for label, ok in results:
    print(f"{'PASS' if ok else 'FAIL'}  {label}")
sys.exit(0 if all(ok for _, ok in results) else 1)
