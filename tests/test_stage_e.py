# Author: Jakub Zapletal
# Date: 2026-11-04
# Checks of Stage E (E.7), step by step. E1: least-cost path, Chaikin and hillshade kernels; snapping on a synthetic
# country T5 (valley minimum 15 m off the line, deeper ditch beyond max_offset_m, reversed line, wide river with a
# nodata water surface, quality-check fallback, pilot units, chainage and tangents), plots and index.html.
# Run: pixi run python -m tests.test_stage_e  (DEFENSES_TEST_DIR sets the working folder, default: temporary)

import json
import os
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from osgeo import gdal, osr
from shapely import affinity

import gis

gdal.UseExceptions()

EPSG = 32633  # synthetic test data only
X0, Y0 = 600000.0, 5600000.0  # lower-left corner of the test world
WIDTH_M, HEIGHT_M, NATIVE = 8000, 3200, 1.0
RES = 2.0  # implied.snap.resolution_m of the test configuration
WORK_RES = 5.0  # working_resolution_m of stage A
# Levee rings (local x0, y0, x1, y1) on the north bank of the main river: A partly covers the town of the E3/E4
# checks; B, 10 m high, protects nearly all of a second built reach.
RINGS = [(4600.0, 1600.0, 5100.0, 1880.0), (6950.0, 1180.0, 7650.0, 2150.0)]
RING_HEIGHTS = [3.0, 10.0]  # B protects its whole envelope (HAND up to envelope_hand_m)
repo = Path(__file__).resolve().parents[1]
root = Path(
    os.environ.get("DEFENSES_TEST_DIR") or tempfile.mkdtemp(prefix="defenses_stage_e_")
)
data = root / "data_e"
out = root / "output" / "T5"
for folder in (data, out, root / "work" / "T5" / "dtm"):
    folder.mkdir(parents=True, exist_ok=True)
print(f"Test folder: {root}")
srs = osr.SpatialReference()
srs.ImportFromEPSG(EPSG)
wkt = srs.ExportToWkt()
list_res, timings = [], {}


def thalweg(x):
    # Main valley: meandering thalweg, local y in m.
    return 1200.0 + 150.0 * np.sin(2 * np.pi * x / 4000.0)


def world(geometry):
    return affinity.translate(geometry, X0, Y0)


# --- E1 kernels ---
cost = np.full((40, 60), 10.0)
cost[5, 0:50] = 0.1
cost[5:35, 49] = 0.1
cost[34, 10:50] = 0.1
cost[20, :] = np.nan
cost[20, 49] = 0.1
path, total = gis.Hydro.least_cost_path(cost, (5, 0), (34, 10))
blocked = gis.Hydro.least_cost_path(
    np.where(np.arange(60)[None, :] == 30, np.nan, 1.0) * np.ones((40, 1)),
    (0, 0),
    (0, 59),
)
list_res.append(
    (
        "least_cost_path: follows the cheap channel through the only gap; a closed wall gives no path",
        len(path) > 0
        and bool(np.all(cost.ravel()[path] == 0.1))
        and blocked[1] == float("inf")
        and len(blocked[0]) == 0,
    )
)
x = np.linspace(0, 100, 51)
y = np.where(np.arange(51) % 2, 1.5, 0.0)
y[25] = 20.0
zigzag = shapely.LineString(np.c_[x, y])
smooth = gis.Vector.chaikin(zigzag, 2, 2.0)
list_res.append(
    (
        "chaikin: end points kept, no point moves more than max_shift, a sharper corner is kept",
        smooth.coords[0] == zigzag.coords[0]
        and smooth.coords[-1] == zigzag.coords[-1]
        and zigzag.hausdorff_distance(smooth) <= 2.0
        and (50.0, 20.0) in list(smooth.coords),
    )
)
z_test = np.add.outer(np.linspace(0, 50, 100), np.linspace(0, 20, 120)).astype(
    np.float32
)
z_test[:5, :5] = np.nan
shade = gis.Raster.hillshade(z_test, (0.0, 2.0, 0.0, 200.0, 0.0, -2.0), wkt)
list_res.append(
    (
        "hillshade: uint8 with 0 at nodata",
        shade.dtype == np.uint8 and shade[0, 0] == 0 and shade[50, 50] > 0,
    )
)

g, used, br = gis.Stats.merge_levels(
    np.array([0, 1, 2, 3, 4, 7, 8]),
    np.array([True, True, False, True, True, True, True]),
    np.array([2.0, 2.1, np.nan, 2.2, 2.0, 2.1, 3.5]),
    0.5,
    1,
)
list_res.append(
    (
        f"merging: three similar subsegments around one park gap form one group with the gap bridged "
        f"(groups {g.tolist()}, bridged {np.flatnonzero(br).tolist()}); a gap of two and a level step start new groups",
        g[:5].tolist() == [0, 0, 0, 0, 0]
        and br.tolist() == [False, False, True, False, False, False, False]
        and abs(used[2] - 2.15) < 1e-9
        and g[5] == 1
        and g[6] == 2,
    )
)
edges_t = np.arange(0, 10.01, 0.25)
env_t = np.full((2, 40), 4000.0)
built_t = np.zeros((2, 40))
built_t[0, 16:] = 0.5 * env_t[0, 16:]
built_t[1, 8:] = np.linspace(0.02, 0.5, 32) * env_t[1, 8:]
ob, sm = gis.Stats.onset(built_t / env_t, 3, 0.5, 4)
sh = gis.Stats.sharpness(sm, ob, 4)
list_res.append(
    (
        f"sharpness: a sharp step {sh[0]:.2f}, a gradual rise {sh[1]:.2f} (sharp_high 0.5)",
        sh[0] >= 0.5 and sh[1] < 0.5,
    )
)

# --- Synthetic country T5 (1 m DTM written in row blocks) ---
# Main river: valley along x, thalweg y = thalweg(x), falling 0.5 m/km to the east; its line runs 15 m north of the
# thalweg. A ditch 3 m deep (deeper than the 2 m channel) runs 40 m south of the line, outside max_offset_m = 30 m.
# Tributary: valley along x = 2000 rising to the north; its line is drawn from the confluence upstream (reversed).
# Wide river: water surface 70 m wide (nodata in the DTM, water polygon) along y = 2600; line 12 m off its centre.
# Zigzag river: channel zigzagging 28 m around a straight line; the snapped path is about twice as long (fallback).
dtm_path = data / "dtm_t5.tif"
if not dtm_path.exists():
    cols, rows = int(WIDTH_M / NATIVE), int(HEIGHT_M / NATIVE)
    ds = gdal.GetDriverByName("GTiff").Create(
        str(dtm_path),
        cols,
        rows,
        1,
        gdal.GDT_Float32,
        ["TILED=YES", "COMPRESS=DEFLATE", "PREDICTOR=3", "BIGTIFF=IF_SAFER"],
    )
    ds.SetGeoTransform((X0, NATIVE, 0.0, Y0 + HEIGHT_M, 0.0, -NATIVE))
    ds.SetProjection(wkt)
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(-9999.0)
    rng = np.random.default_rng(5)
    gx = (np.arange(cols) + 0.5) * NATIVE
    for r0 in range(0, rows, 400):
        gy = HEIGHT_M - (np.arange(r0, min(r0 + 400, rows)) + 0.5) * NATIVE
        xx, yy = np.meshgrid(gx, gy)
        d = np.abs(yy - thalweg(xx))
        z = 100.0 + 0.0005 * (8000.0 - xx) + 0.01 * d
        z -= 2.0 * np.maximum(0.0, 1.0 - d / 6.0)
        z -= 3.0 * np.maximum(0.0, 1.0 - np.abs(yy - (thalweg(xx) - 25.0)) / 2.0)
        dt = np.abs(xx - 2000.0)
        z_trib = 103.0 + 0.002 * (yy - 1200.0) + 0.01 * dt
        z_trib -= 2.0 * np.maximum(0.0, 1.0 - dt / 6.0)
        z = np.where(yy > 1150.0, np.minimum(z, z_trib), z)
        dz = np.abs(yy - (400.0 + 28.0 * np.sin(2 * np.pi * xx / 60.0)))
        z_zig = 100.6 + 0.0005 * (8000.0 - xx) + 0.01 * np.abs(yy - 400.0)
        z_zig -= 2.0 * np.maximum(0.0, 1.0 - dz / 6.0)
        z = np.where((xx > 4300) & (xx < 6700) & (yy < 700), np.minimum(z, z_zig), z)
        for (x0, y0, x1, y1), height in zip(RINGS, RING_HEIGHTS):
            # embankments 6 m wide along the ring rectangles
            on_x = (
                ((np.abs(xx - x0) <= 3) | (np.abs(xx - x1) <= 3))
                & (yy >= y0 - 3)
                & (yy <= y1 + 3)
            )
            on_y = (
                ((np.abs(yy - y0) <= 3) | (np.abs(yy - y1) <= 3))
                & (xx >= x0 - 3)
                & (xx <= x1 + 3)
            )
            z = np.where(on_x | on_y, z + height, z)
        z += rng.normal(0.0, 0.03, z.shape)
        wet = (np.abs(yy - 2600.0) <= 35.0) & (xx > 3900.0) & (xx < 7700.0)
        z[wet] = -9999.0
        band.WriteArray(z.astype(np.float32), 0, r0)
    ds = None
gdal.BuildVRT(str(root / "work" / "T5" / "dtm" / "dtm.vrt"), [str(dtm_path)])

xs = np.arange(100.0, 7900.0 + 1e-9, 10.0)
lines = {
    "main": shapely.LineString(np.c_[xs, thalweg(xs) + 15.0]),
    "trib": shapely.LineString([(2000.0, 1260.0), (2000.0, 3150.0)]),
    "wide": shapely.LineString([(4000.0, 2612.0), (7600.0, 2612.0)]),
    "zigzag": shapely.LineString([(4500.0, 400.0), (6500.0, 400.0)]),
    "closed": shapely.LineString([(500.0, 2900.0), (1500.0, 2900.0)]),
}
gdf_rn = gpd.GeoDataFrame(
    {
        "river_id": np.arange(len(lines)),
        "name": list(lines),
        "is_source": [True, True, True, True, False],
    },
    geometry=[world(g) for g in lines.values()],
    crs=f"EPSG:{EPSG}",
)
# the rivers input holds the water sources only (the closed line is drawn in the DTM but is not a river)
gis.Vector.write_gpkg(gdf_rn[gdf_rn["is_source"]], data / "rivers.gpkg")
pd.DataFrame(
    {"country": ["T5"], "residential": [1000.0], "non_residential": [800.0]}
).to_csv(data / "max_damage.csv", index=False)
gis.Vector.write_gpkg(
    gpd.GeoDataFrame(
        {"name": [f"ring_{k}" for k in range(len(RINGS))]},
        geometry=[
            world(shapely.LineString(shapely.box(*r).exterior.coords)) for r in RINGS
        ],
        crs=f"EPSG:{EPSG}",
    ),
    data / "levees.gpkg",
)
gis.Vector.write_gpkg(
    gpd.GeoDataFrame(
        {"name": ["wide"]},
        geometry=[world(shapely.box(3900.0, 2565.0, 7700.0, 2635.0))],
        crs=f"EPSG:{EPSG}",
    ),
    data / "water.gpkg",
)
unit_box = world(shapely.box(3000.0, 0.0, 6000.0, HEIGHT_M))
for name in ("units_window.gpkg", "units_core.gpkg"):
    gis.Vector.write_gpkg(
        gpd.GeoDataFrame({"unit_id": ["U1"]}, geometry=[unit_box], crs=f"EPSG:{EPSG}"),
        out / name,
    )


def write_config(name, **kw):
    v = {
        "units": "[]",
        "plot_all": "[]",
        "max_per_type": 200,
        "ratio": 1.5,
        "length": 500.0,
        "recompute": "false",
        "footprints": "[]",
        "built_source": "footprints",
        "include_implied": "false",
        "reference": "[]",
    }
    v.update(kw)
    (root / name).write_text(f"""[general]
work_dir = "{(root / 'work').as_posix()}"
output_dir = "{(root / 'output').as_posix()}"
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

[countries.T5]
dtm = "{dtm_path.as_posix()}"
dsm = []
levees = "{(data / 'levees.gpkg').as_posix()}"
rivers = "{(data / 'rivers.gpkg').as_posix()}"
water_polygons = "{(data / 'water.gpkg').as_posix()}"
building_footprints = {v["footprints"]}
river_upstream_area_field = ""
reference_protected_areas = {v["reference"]}
working_resolution_m = {WORK_RES}
corridor_buffer_m = 1000.0
unit_core_size_m = 2000.0
unit_halo_m = 1000.0
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
keep_intermediate = true

[countries.T5.protected_areas]
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

[countries.T5.rp]
levee_adjacent_distance_m = 10.0
classes = [10, 20, 50, 100, 200, 500]
storey_height_m = 3.0
max_damage_csv = "{(data / 'max_damage.csv').as_posix()}"

[countries.T5.context]
order = ["industrial", "dense_urban", "urban", "rural"]

[countries.T5.context.params.industrial]
score_weight_value = 0.5
score_weight_population = 0.3
score_weight_upstream = 0.2

[countries.T5.context.params.dense_urban]
score_weight_value = 0.5
score_weight_population = 0.3
score_weight_upstream = 0.2

[countries.T5.context.params.urban]
score_weight_value = 0.5
score_weight_population = 0.3
score_weight_upstream = 0.2

[countries.T5.context.params.rural]
score_weight_value = 0.5
score_weight_population = 0.3
score_weight_upstream = 0.2

[implied]
units = {v["units"]}
recompute_hand = {v["recompute"]}
envelope_hand_m = 10.0
subsegment_length_m = {v["length"]}
subsegment_length_grid = [500.0, 1000.0]
built_source = "{v["built_source"]}"
exclude_building_tags = {{ building = ["boathouse", "houseboat"], man_made = ["wastewater_plant", "water_works"] }}
built_share_min = 0.3
built_share_min_grid = [0.2, 0.3]
built_area_min_m2 = 20000.0
hand_bin_m = 0.25
min_bin_area_m2 = 2000.0
smooth_bins = 3
onset_fraction = 0.5
onset_persist_bins = 4
quantiles = [0.02, 0.05, 0.10]
primary_level_method = "onset"
sharp_window_m = 1.0
sharp_high = 0.5
component_built_min_m2 = 5000.0
merge_max_dlevel_m = 0.5
merge_max_gap = 1
corroboration_distance_m = 50.0
corroboration_min = 0.3
explicit_overlap_max = 0.8
regional_sop_rp = 10
implied_corroborated_max_rp = 100
historic_low_core_uncorroborated_cap = true
ghsl_epochs = [1990, 2020]
epoch_cutoff = 1990
growth_low_share_min = 0.1
min_cells_temporal = 20
rp_distribution_classes = [20, 50, 100]
rp_distribution_shares = [0.4, 0.4, 0.2]
include_implied = {v["include_implied"]}
include_in_review = true

[evidence]
comparison_distance_m = 100.0
min_valid_pixels = 10
f_wet_comparison = 0.5
f_dry_inside = 0.1
f_wet_inside = 0.5

[sources]
gfm_flood_values = [1]
gfm_exclusion_values = [1]
gfm_reference_water_values = [1]
gfm_nodata_values = [255]

[implied.snap]
step_m = 5.0
resolution_m = {RES}
max_offset_m = 30.0
min_window_m = 15.0
cost_epsilon = 0.01
chunk_length_m = 5000.0
chunk_overlap_m = 200.0
max_length_ratio = {v["ratio"]}
wide_river_m = 40.0
direction_check_m = 200.0
max_nodata_share = 0.5
chainage_origin = "downstream"

[implied.plots]
enabled = true
dpi = 80
max_per_type = {v["max_per_type"]}
sample_seed = 42
show_max = 2
plot_all_in_units = {v["plot_all"]}
""")


def run(label, config, module="processing.snap_rivers", **variables):
    env = dict(
        os.environ,
        DEFENSES_CONFIG_PATH=str(root / config),
        DEFENSES_COUNTRY="T5",
        **{f"DEFENSES_{k}": str(v) for k, v in variables.items()},
    )
    env.pop("DISPLAY", None)
    env.pop("WAYLAND_DISPLAY", None)
    begin = time.perf_counter()
    proc = subprocess.run(
        [sys.executable, "-m", module],
        env=env,
        capture_output=True,
        text=True,
        cwd=repo,
    )
    timings[label] = time.perf_counter() - begin
    (root / f"e_{label}.log").write_text(proc.stdout + proc.stderr)
    if proc.returncode != 0:
        print(proc.stderr[-3000:])
    return proc


# --- Pilot run: unit U1 only, plot_all_in_units forces every chunk map despite max_per_type = 1 ---
write_config("test_e_pilot.toml", units='["U1"]', plot_all='["U1"]', max_per_type=1)
proc = run("pilot", "test_e_pilot.toml", RUN_ID="pilot", SHOW_PLOTS=1)
if proc.returncode == 0:
    gdf_snap = gpd.read_file(out / "rivers_snapped.gpkg")
    maps = sorted((out / "plots" / "implied" / "pilot" / "snap_map").glob("*.png"))
    n_chunks = int(gdf_snap["n_chunks"].sum())
    inside = bool(len(gdf_snap)) and bool(
        gdf_snap.geometry.within(unit_box.buffer(RES)).all()
    )
    list_res.append(
        (
            f"pilot units: rivers clipped to the U1 window ({len(gdf_snap)} lines, "
            f"{gdf_snap['length_original_m'].sum() / 1000:.2f} km, all inside: {inside})",
            inside and len(gdf_snap) >= 3,
        )
    )
    list_res.append(
        (
            f"plot_all_in_units: all {n_chunks} chunk maps of U1 plotted with max_per_type = 1 ({len(maps)})",
            len(maps) == n_chunks,
        )
    )
    list_res.append(
        (
            "SHOW_PLOTS without a display: warning logged, figures saved",
            "no display available" in proc.stderr and len(maps) > 0,
        )
    )
else:
    list_res.append(("snap_rivers runs on pilot units", False))

# --- Full run ---
write_config("test_e.toml")
proc = run("full", "test_e.toml")
if proc.returncode == 0:
    gdf_snap = gpd.read_file(out / "rivers_snapped.gpkg").set_index("name")
    gdf_pts = gpd.read_file(out / "river_points.gpkg")
    perf = json.loads((out / "snap_rivers_performance.json").read_text())
    list_res.append(
        (
            f"all rivers of the input are snapped: {sorted(gdf_snap.index)}",
            sorted(gdf_snap.index) == ["main", "trib", "wide", "zigzag"],
        )
    )
    main = gdf_snap.loc["main"]
    xs_true = np.arange(300.0, 7700.0, 2.0)
    true_line = world(shapely.LineString(np.c_[xs_true, thalweg(xs_true)]))
    # Points every 2 m along the snapped line (a straight stretch keeps few vertices after simplification).
    xy = shapely.get_coordinates(
        shapely.line_interpolate_point(
            main.geometry, np.arange(0.0, main.geometry.length, RES)
        )
    )
    middle = (xy[:, 0] > X0 + 300.0) & (xy[:, 0] < X0 + 7700.0)
    dist_true = shapely.distance(true_line, shapely.points(xy[middle]))
    list_res.append(
        (
            f"valley minimum 15 m off the line found within one cell: max {dist_true.max():.2f} m, "
            f"mean {dist_true.mean():.2f} m from the thalweg (cell {RES} m); mean offset from the line "
            f"{main.mean_offset_m:.1f} m",
            dist_true.max() <= RES and 13.0 <= main.mean_offset_m <= 16.0,
        )
    )
    ditch = world(shapely.LineString(np.c_[xs_true, thalweg(xs_true) - 25.0]))
    gap = float(shapely.distance(ditch, main.geometry))
    list_res.append(
        (
            f"deeper ditch beyond max_offset_m not chosen: {gap:.1f} m between the snapped line and the ditch, "
            f"max offset {main.max_offset_m:.1f} m",
            gap > 20.0 and main.max_offset_m <= 30.0 + RES,
        )
    )
    list_res.append(
        (
            f"main river: {int(main.n_chunks)} chunks stitched into one line, not flipped",
            int(main.n_chunks) == 2
            and main.geometry.geom_type == "LineString"
            and not bool(main.flipped),
        )
    )
    trib = gdf_snap.loc["trib"]
    trib_xy = shapely.get_coordinates(trib.geometry)
    list_res.append(
        (
            f"reversed tributary flipped (flipped = {bool(trib.flipped)}; now from y {trib_xy[0, 1] - Y0:.0f} "
            f"to {trib_xy[-1, 1] - Y0:.0f})",
            bool(trib.flipped) and trib_xy[0, 1] > trib_xy[-1, 1],
        )
    )
    wide = gdf_snap.loc["wide"]
    wide_xy = shapely.get_coordinates(
        shapely.line_interpolate_point(
            wide.geometry, np.arange(0.0, wide.geometry.length, RES)
        )
    )
    mid_w = (wide_xy[:, 0] > X0 + 4300.0) & (wide_xy[:, 0] < X0 + 7300.0)
    off_c = np.abs(wide_xy[mid_w, 1] - (Y0 + 2600.0))
    list_res.append(
        (
            f"wide river (nodata water surface) follows the polygon centre: {int(wide.wide_river_chunks)} wide "
            f"chunk, max {off_c.max():.2f} m from the centre line",
            int(wide.wide_river_chunks) == 1 and off_c.max() <= RES,
        )
    )
    zig = gdf_snap.loc["zigzag"]
    list_res.append(
        (
            f"quality check: a path {zig.n_kept_original}/{zig.n_chunks} chunks longer than max_length_ratio keeps "
            f"the original line (Hausdorff {zig.geometry.hausdorff_distance(world(lines['zigzag'])):.2f} m)",
            int(zig.n_kept_original) == int(zig.n_chunks) == 1
            and zig.geometry.hausdorff_distance(world(lines["zigzag"])) <= RES,
        )
    )
    p_main = gdf_pts[gdf_pts["river_id"] == gdf_snap.loc["main", "river_id"]]
    p_main = p_main.sort_values("chainage_m")
    first = shapely.get_coordinates(p_main.geometry.iloc[[0]])[0]
    end = shapely.get_coordinates(main.geometry)[-1]
    spacing = np.diff(p_main["chainage_m"].to_numpy())
    p_trib = gdf_pts[gdf_pts["river_id"] == gdf_snap.loc["trib", "river_id"]]
    list_res.append(
        (
            f"river_points: chainage 0 at the downstream end ({np.hypot(*(first - end)):.2f} m), spacing "
            f"{spacing.min():.1f}-{spacing.max():.1f} m, tangents downstream (main dx {p_main['tangent_dx'].mean():.2f}, "
            f"tributary dy {p_trib['tangent_dy'].mean():.2f})",
            np.hypot(*(first - end)) < 0.01
            and np.allclose(spacing, 5.0)
            and p_main["tangent_dx"].mean() > 0.9
            and p_trib["tangent_dy"].mean() < -0.9
            and np.allclose(
                np.hypot(gdf_pts["tangent_dx"], gdf_pts["tangent_dy"]), 1.0
            ),
        )
    )
    run_dir = out / "plots" / "implied" / "base"
    profiles = sorted((run_dir / "snap_profile").glob("*.png"))
    maps = sorted((run_dir / "snap_map").glob("*.png"))
    index = (
        (run_dir / "index.html").read_text()
        if (run_dir / "index.html").exists()
        else ""
    )
    listed = all(p.name in index for p in profiles + maps)
    list_res.append(
        (
            f"plots: {len(profiles)} snap_profile and {len(maps)} snap_map figures; index.html lists all with metrics",
            len(profiles) == 4
            and len(maps) == int(gdf_snap["n_chunks"].sum())
            and listed
            and "max_offset_m" in index
            and "kept_original" in index,
        )
    )
    list_res.append(
        (
            "kept-original chunk reported in the log and flagged first in the plots",
            re.search(r"chunk 0 kept original: length ratio", proc.stderr) is not None,
        )
    )
    print(
        "snap_rivers timings: "
        + ", ".join(
            f"{k} {v}"
            for k, v in perf.items()
            if k.endswith("_s") or "rss" in k or k in ("river_km", "chunks")
        )
    )
else:
    list_res.append(("snap_rivers runs", False))

# --- Stage A on T5 (inputs of E2-E4): DTM preparation, units, protected areas behind the levee rings ---
write_config("test_e.toml")
stage_a_ok = True
for module in ("prepare_dtm", "build_units", "protected_areas"):
    proc = run(module, "test_e.toml", f"processing.{module}")
    if proc.returncode != 0:
        list_res.append((f"stage A {module} runs for T5", False))
        stage_a_ok = False
        break
if stage_a_ok:
    gdf_pa = gpd.read_file(out / "protected_areas.gpkg")
    list_res.append(
        (
            f"stage A on T5: {len(gdf_pa)} explicit protected areas behind the levee rings",
            len(gdf_pa) >= 2,
        )
    )
    # --- E2: subsegments of 1000 m directly, then with HAND recomputed, then the base of 500 m ---
    write_config("test_e_1000.toml", length=1000.0)
    proc_1000 = run(
        "segments_1000", "test_e_1000.toml", "processing.implied_segments", FORCE=1
    )
    df_1000 = (
        gpd.read_file(out / "subsegments.gpkg").drop(columns="geometry")
        if proc_1000.returncode == 0
        else None
    )
    write_config("test_e_rh.toml", recompute="true")
    proc_rh = run(
        "segments_rehand", "test_e_rh.toml", "processing.implied_segments", FORCE=1
    )
    proc = run("segments", "test_e.toml", "processing.implied_segments", FORCE=1)
    proc_again = run("segments_again", "test_e.toml", "processing.implied_segments")
    if proc.returncode == 0 and df_1000 is not None:
        gdf_seg = gpd.read_file(out / "subsegments.gpkg")
        gdf_snap = gpd.read_file(out / "rivers_snapped.gpkg").set_index("name")
        main_id = int(gdf_snap.loc["main", "river_id"])
        seg_main = gdf_seg[gdf_seg["river_id"] == main_id]
        rep = shapely.get_coordinates(seg_main.geometry.representative_point())
        north = rep[:, 1] - Y0 > thalweg(rep[:, 0] - X0)
        left = (seg_main["side"] == "left").to_numpy()
        list_res.append(
            (
                f"bank side on the meander: {int((north == left).sum())} of {len(seg_main)} subsegments of the "
                "main river on the correct bank (left = north, looking downstream to the east)",
                len(seg_main) > 10 and bool(np.all(north == left)),
            )
        )
        union = shapely.union_all(gdf_seg.geometry.values).area
        list_res.append(
            (
                f"core cells counted once: union {union / 1e6:.3f} km2 = sum of envelope areas "
                f"{gdf_seg['envelope_area_m2'].sum() / 1e6:.3f} km2 = sum of polygon areas",
                abs(union - gdf_seg["envelope_area_m2"].sum()) < 1.0
                and abs(gdf_seg.area.sum() - union) < 1.0,
            )
        )
        df_agg = (
            gdf_seg.assign(k=gdf_seg["k"] // 2)
            .groupby(["river_id", "side", "k"])["envelope_area_m2"]
            .sum()
        )
        df_direct = df_1000.set_index(["river_id", "side", "k"])["envelope_area_m2"]
        both = df_agg.to_frame("agg").join(df_direct.to_frame("direct"), how="outer")
        list_res.append(
            (
                f"aggregation: {len(both)} subsegments of 1000 m from 500 m equal the direct computation "
                f"(max difference {float((both['agg'] - both['direct']).abs().max()):.3f} m2)",
                bool(both.notna().all().all())
                and float((both["agg"] - both["direct"]).abs().max()) < 1e-6,
            )
        )
        perf_again = json.loads((out / "implied_segments_performance.json").read_text())
        list_res.append(
            (
                f"checkpoints: a second run processes {perf_again['processed']} of {perf_again['units']} units",
                perf_again["processed"] == 0 and proc_again.returncode == 0,
            )
        )
        unit_dirs = sorted(
            (root / "work" / "T5" / "units").glob("*/hand_undef_snapped.tif")
        )
        zero = []
        for path in unit_dirs:
            ds = gdal.Open(str(path))
            hand = ds.GetRasterBand(1).ReadAsArray()
            cells = (
                gis.Raster.rasterize_geometries(
                    gdf_snap.geometry.values,
                    np.ones(len(gdf_snap)),
                    hand.shape,
                    ds.GetGeoTransform(),
                    ds.GetProjection(),
                    True,
                    np.uint8,
                )
                > 0
            )
            vals = hand[cells & np.isfinite(hand)]
            zero.append(float(np.abs(vals).max()) if len(vals) else 0.0)
        list_res.append(
            (
                f"recompute_hand: hand_undef_snapped.tif in {len(unit_dirs)} units through the stage A code path, "
                f"HAND 0 on the snapped rivers (max |HAND| {max(zero or [np.nan]):.3f} m)",
                proc_rh.returncode == 0 and len(unit_dirs) > 0 and max(zero) == 0.0,
            )
        )
        perf = json.loads(
            (root / "output" / "T5" / "implied_segments_performance.json").read_text()
        )
        log_text = (root / "e_segments.log").read_text()
        print(
            "implied_segments: "
            + (re.search(r"Timings: (.*)", log_text) or [None, ""])[1]
        )
    else:
        list_res.append(("implied_segments runs", False))

# --- E3: buildings laid out from the subsegments of E2, GHSL epochs, zones; implied_protection parts A-B ---
if stage_a_ok and (out / "subsegments.gpkg").exists():
    gdf_seg = gpd.read_file(out / "subsegments.gpkg")
    gdf_snap = gpd.read_file(out / "rivers_snapped.gpkg").set_index("name")
    main_id = int(gdf_snap.loc["main", "river_id"])
    left = gdf_seg[
        (gdf_seg["river_id"] == main_id) & (gdf_seg["side"] == "left")
    ].copy()
    cx = left.geometry.centroid.x.to_numpy() - X0
    town = left[(cx > 2600) & (cx < 5400)].sort_values("k")
    park_k = int(town["k"].iloc[len(town) // 2])
    mills_k = int(town["k"].iloc[1])
    reach_b = left[(cx > 7000) & (cx < 7600)]
    gx, gy = np.meshgrid(
        np.arange(0.0, WIDTH_M, 20.0) + 10.0, np.arange(0.0, HEIGHT_M, 20.0) + 10.0
    )
    gx, gy = gx.ravel(), gy.ravel()
    d = gy - thalweg(gx)
    pts = shapely.points(gx + X0, gy + Y0)
    list_b, list_tag = [], []

    def in_segments(segs):
        tree = shapely.STRtree(segs.geometry.values)
        hit = np.zeros(len(pts), bool)
        hit[np.unique(tree.query(pts, predicate="within")[0])] = True
        return hit

    town_cells = in_segments(town[town["k"] != park_k]) & (d >= 200) & (d <= 750)
    mills_cells = (
        in_segments(town[town["k"] == mills_k])
        & (d >= 10)
        & (d < 190)
        & (np.arange(len(gx)) % 5 == 0)
    )
    reach_cells = in_segments(reach_b) & (d >= 150) & (d <= 750)
    for mask, half in ((town_cells, 8.0), (mills_cells, 8.0), (reach_cells, 8.0)):
        list_b.append(
            shapely.box(
                gx[mask] + X0 - half,
                gy[mask] + Y0 - half,
                gx[mask] + X0 + half,
                gy[mask] + Y0 + half,
            )
        )
    boat_x = float(town["geometry"].iloc[2].centroid.x)
    boat = shapely.box(
        boat_x - 6,
        Y0 + thalweg(boat_x - X0) + 12,
        boat_x + 6,
        Y0 + thalweg(boat_x - X0) + 20,
    )
    geoms_b = np.concatenate(list_b + [[boat]])
    gdf_b = gpd.GeoDataFrame(
        {
            "building": ["yes"] * (len(geoms_b) - 1) + ["boathouse"],
            "man_made": [None] * len(geoms_b),
        },
        geometry=geoms_b,
        crs=f"EPSG:{EPSG}",
    )
    gis.Vector.write_gpkg(gdf_b, data / "buildings.gpkg")
    # GHSL built-up surface per epoch (m2 per 100 m cell): 1990 without the mills and the second reach.
    cells_x = np.arange(0.0, WIDTH_M, 100.0)
    cells_y = np.arange(0.0, HEIGHT_M, 100.0)
    for epoch, parts in ((1990, [list_b[0]]), (2020, list_b)):
        union = shapely.union_all(np.concatenate(parts))
        boxes = shapely.box(
            cells_x[None, :] + X0,
            (HEIGHT_M - cells_y[:, None] - 100.0) + Y0,
            cells_x[None, :] + X0 + 100.0,
            (HEIGHT_M - cells_y[:, None]) + Y0,
        )
        area = shapely.area(shapely.intersection(boxes, union)).astype(np.float32)
        tif = data / f"ghsl_built_s_{epoch}.tif"
        gis.Raster.write_array(
            tif, area, (X0, 100.0, 0.0, Y0 + HEIGHT_M, 0.0, -100.0), wkt, nodata=-1.0
        )
        gdal.BuildVRT(str(out / f"ghsl_built_s_{epoch}.vrt"), [str(tif)])
    zone_town = world(shapely.box(2500.0, 1250.0, 5500.0, 2100.0))
    corridor = shapely.union_all(gpd.read_file(out / "corridor.gpkg").geometry.values)
    gis.Vector.write_gpkg(
        gpd.GeoDataFrame(
            {
                "zone_id": [0, 1],
                "context": ["dense_urban", "rural"],
                "sources": ["test", "test"],
            },
            geometry=[zone_town.intersection(corridor), corridor.difference(zone_town)],
            crs=f"EPSG:{EPSG}",
        ),
        out / "zones.gpkg",
    )
    write_config("test_e.toml", footprints=f'"{(data / "buildings.gpkg").as_posix()}"')
    proc = run("protection", "test_e.toml", "processing.implied_protection")
    if proc.returncode == 0:
        gdf_imp = gpd.read_file(
            out / "implied" / "base" / "implied_subsegments.gpkg"
        ).set_index("subsegment_id")
        # Expected onset: stage A HAND at the lower edge of the town (200 m from the thalweg), sampled from the
        # unit rasters. Stage A measures HAND from the river line as drawn, 15 m off the channel.
        list_h = []
        for path in sorted((root / "work" / "T5" / "units").glob("*/hand_undef.tif")):
            ds = gdal.Open(str(path))
            hgt, harr = ds.GetGeoTransform(), ds.GetRasterBand(1).ReadAsArray()
            ex = np.arange(2700.0, 5300.0, 50.0)
            col = ((ex + X0 - hgt[0]) / hgt[1]).astype(int)
            row = ((hgt[3] - (thalweg(ex) + 200.0 + Y0)) / -hgt[5]).astype(int)
            ok = (row >= 0) & (row < harr.shape[0]) & (col >= 0) & (col < harr.shape[1])
            list_h.append(harr[row[ok], col[ok]])
        edge = float(np.nanmedian(np.concatenate(list_h)))
        t_ids = [f"{main_id}_left_{k}" for k in town["k"] if k != park_k]
        t = gdf_imp.loc[t_ids]
        list_res.append(
            (
                f"onset at the town edge: {len(t)} town subsegments eligible {int(t['eligible'].sum())}, onset "
                f"{t['level_onset_m'].min():.2f}-{t['level_onset_m'].max():.2f} m (edge at HAND {edge:.2f} m)",
                bool(t["eligible"].all())
                and bool((t["level_onset_m"] - edge).abs().max() <= 0.5),
            )
        )
        mill = gdf_imp.loc[f"{main_id}_left_{mills_k}"]
        list_res.append(
            (
                f"isolated low buildings: onset {mill['level_onset_m']:.2f} m ignores them, q0.05 "
                f"{mill['level_q0_05_m']:.2f} m and q0.02 {mill['level_q0_02_m']:.2f} m follow them",
                abs(mill["level_onset_m"] - edge) <= 0.5
                and mill["level_q0_05_m"] < edge - 0.25,
            )
        )
        list_res.append(
            (
                f"sharp town edge: sharpness {t['sharpness'].min():.2f}-{t['sharpness'].max():.2f} (sharp_high 0.5)",
                bool((t["sharpness"] >= 0.5).all()),
            )
        )
        park = gdf_imp.loc[f"{main_id}_left_{park_k}"]
        right = gdf_imp[(gdf_imp["river_id"] == main_id) & (gdf_imp["side"] == "right")]
        list_res.append(
            (
                f"eligibility: the park subsegment not eligible ({park['reason']}); {int(right['eligible'].sum())} "
                "eligible subsegments on the unbuilt right bank",
                not bool(park["eligible"]) and int(right["eligible"].sum()) == 0,
            )
        )
        list_res.append(
            (
                f"context: town subsegments labelled {sorted(set(t['context']))}",
                set(t["context"]) == {"dense_urban"},
            )
        )
        log_text = proc.stderr
        curves = np.load(out / "implied" / "base" / "curves.npz")
        list_res.append(
            (
                "exclusions and curves: the boathouse excluded and logged; curves.npz with "
                f"{curves['envelope_area'].shape} envelope bins and {curves['built_area_epochs'].shape[0]} GHSL epochs",
                "1 with building in ['boathouse', 'houseboat'] excluded" in log_text
                and curves["envelope_area"].shape == (len(gdf_imp), 40)
                and curves["built_area_epochs"].shape[0] == 2,
            )
        )
        run_dir = out / "plots" / "implied" / "base"
        n_curve = len(list((run_dir / "curve").glob("*.png")))
        index = (run_dir / "index.html").read_text()
        list_res.append(
            (
                f"plots: {n_curve} curve figures and the summary; index.html lists curve, summary and the snap types",
                n_curve == min(200, len(gdf_imp))
                and (run_dir / "summary" / "summary.png").exists()
                and all(
                    k in index for k in ("curve", "summary", "snap_map", "level_onset")
                ),
            )
        )
        print(
            "implied_protection: "
            + (re.search(r"Timings: (.*)", log_text) or [None, ""])[1]
        )
    else:
        list_res.append(("implied_protection (A-B) runs", False))

# --- E4: implied protected areas, precedence, temporal signal, corroboration, RP ---
imp_path = out / "implied" / "base" / "implied_protected_areas.gpkg"
if stage_a_ok and imp_path.exists():
    gdf_ic = gpd.read_file(imp_path)
    gdf_pa = gpd.read_file(out / "protected_areas.gpkg")
    ring_a = gdf_pa[gdf_pa.intersects(world(shapely.box(*RINGS[0])))]
    town_c = gdf_ic[
        (gdf_ic["river_id"] == main_id)
        & (gdf_ic["side"] == "left")
        & (gdf_ic["k_from"] <= park_k)
        & (gdf_ic["k_to"] >= park_k)
    ]
    list_res.append(
        (
            f"merging: the town forms {len(town_c)} implied protected area over subsegments "
            f"{town_c['k_from'].tolist()}-{town_c['k_to'].tolist()}, {town_c['n_bridged'].tolist()} bridged (park)",
            len(town_c) == 1 and int(town_c["n_bridged"].iloc[0]) >= 1,
        )
    )
    inside_a = (
        float(
            shapely.intersection(
                town_c.geometry.iloc[0], shapely.union_all(ring_a.geometry.values)
            ).area
        )
        if len(town_c) and len(ring_a)
        else -1
    )
    list_res.append(
        (
            (
                f"precedence: the town clipped by ring A (overlap share {town_c['explicit_overlap_share'].iloc[0]:.3f}, "
                f"{inside_a:.0f} m2 left inside it)"
                if len(town_c)
                else "precedence: town protected area missing"
            ),
            len(town_c) == 1
            and 0 < town_c["explicit_overlap_share"].iloc[0] < 0.8
            and 0 <= inside_a < 1.0,
        )
    )
    dropped = pd.read_csv(out / "implied" / "base" / "implied_dropped.csv")
    reach_b_c = gdf_ic[
        (gdf_ic["river_id"] == main_id)
        & (gdf_ic["side"] == "left")
        & (gdf_ic["k_to"] < town["k"].min())
    ]
    list_res.append(
        (
            f"precedence: the reach inside ring B dropped ({'; '.join(dropped['reason'])})",
            len(reach_b_c) == 0
            and dropped["reason"].str.contains("explicit_overlap_share").any(),
        )
    )
    if len(town_c):
        r = town_c.iloc[0]
        list_res.append(
            (
                f"temporal: onset per epoch {r['onset_epochs']}, growth below the onset {r['growth_low_share']:.3f}, "
                f"historic low core {bool(r['historic_low_core'])}",
                bool(r["historic_low_core"])
                and not bool(r["growth_into_floodplain"])
                and abs(r["delta_onset_m"]) <= 0.5,
            )
        )
        list_res.append(
            (
                f"RP: crest support {r['crest_support']:.2f}, corroborated {bool(r['corroborated'])}, rp {r['rp']:.0f} "
                f"({r['rp_low']:.0f}-{r['rp_high']:.0f}), cap {r['cap_reason']}, confidence {r['confidence']}",
                r["rp_source"] == "implied"
                and r["confidence"] == "low"
                and r["rp"] <= 20
                and r["rp_low"] >= 10
                and r["rp_low"] <= r["rp"] <= r["rp_high"],
            )
        )
    run_dir = out / "plots" / "implied" / "base"
    counts = {
        k: len(list((run_dir / k).glob("*.png")))
        for k in ("longitudinal", "temporal", "protected_area_map")
    }
    index = (run_dir / "index.html").read_text()
    list_res.append(
        (
            f"plots: {counts} figures, listed in index.html with their metrics",
            counts["longitudinal"] == 4
            and counts["temporal"] == len(gdf_ic)
            and counts["protected_area_map"] == len(gdf_ic)
            and "crest_support" in index,
        )
    )
else:
    list_res.append(
        ("implied_protection (C-F) writes implied_protected_areas.gpkg", False)
    )

# --- E5: GFM scenes and event RPs, reference areas, integration, sensitivity grid ---
if stage_a_ok and imp_path.exists() and len(town_c):
    town_id = town_c["implied_id"].iloc[0]
    scene_gt = (X0 + 2500.0, 10.0, 0.0, Y0 + 2200.0, 0.0, -10.0)
    sx = (np.arange(310) + 0.5) * 10.0 + 2500.0
    sy = 2200.0 - (np.arange(110) + 0.5) * 10.0
    wet = (
        (sx[None, :] > 3200)
        & (sx[None, :] < 3600)
        & (sy[:, None] > 1400)
        & (sy[:, None] < 1700)
    ).astype(np.uint8)
    list_idx = []
    for item, when in (
        ("SCENE_A", "2020-05-01T10:00:00Z"),
        ("SCENE_B", "2021-06-01T10:00:00Z"),
        ("SCENE_C", "2022-07-01T10:00:00Z"),
    ):
        paths = {}
        for role, arr in (
            ("flood", wet),
            ("exclusion", np.zeros_like(wet)),
            ("reference_water", np.zeros_like(wet)),
            ("likelihood", wet * 90),
        ):
            paths[role] = data / f"{item}_{role}.tif"
            gis.Raster.write_array(paths[role], arr, scene_gt, wkt, nodata=255)
        for u in pd.read_csv(out / "units.csv", dtype={"unit_id": str}).itertuples():
            list_idx.append(
                {
                    "unit_id": u.unit_id,
                    "item_id": item,
                    "datetime": when,
                    **{f"{k}_path": str(v) for k, v in paths.items()},
                }
            )
    pd.DataFrame(list_idx).to_csv(out / "gfm_index.csv", index=False)
    # event RPs per test as rp_bounds.py combines them over the sources (Stage F); EFAS higher than GloFAS on
    # SCENE_A, so that conservative_both (flooded: the higher RP) is what the Stage E scripts read
    pd.DataFrame(
        {
            "protected_area_id": [town_id, town_id],
            "item_id": ["SCENE_A", "SCENE_B"],
            "outcome": ["flooded", "flooded"],
            "rp_combined": [10.0, 50.0],
            "used": ["both", "glofas"],
            "source_disagreement": [False, False],
            "rp_glofas": [6.0, 50.0],
            "rp_efas": [10.0, np.nan],
        }
    ).to_csv(out / "event_rp_combined.csv", index=False)
    gis.Vector.write_gpkg(
        gpd.GeoDataFrame(
            {"name": ["town"]},
            geometry=[world(shapely.box(2700.0, 1250.0, 5300.0, 2000.0))],
            crs=f"EPSG:{EPSG}",
        ),
        data / "reference.gpkg",
    )
    gdf_clip = gpd.read_file(out / "protected_areas.gpkg")[
        ["protected_area_id", "geometry"]
    ].assign(rp=50.0, rp_source="model")
    gis.Vector.write_gpkg(gdf_clip, out / "defenses_clip.gpkg")
    write_config(
        "test_e.toml",
        footprints=f'"{(data / "buildings.gpkg").as_posix()}"',
        include_implied="true",
        reference=f'"{(data / "reference.gpkg").as_posix()}"',
    )
    proc = run("protection_e5", "test_e.toml", "processing.implied_protection")
    clip_ok = False
    if proc.returncode == 0 and (out / "defenses_clip_with_implied.gpkg").exists():
        gdf_cwi = gpd.read_file(out / "defenses_clip_with_implied.gpkg")
        kept = gdf_cwi[gdf_cwi["protected_area_id"].notna()]
        added = (
            gdf_cwi[gdf_cwi["implied_id"].notna()]
            if "implied_id" in gdf_cwi
            else gdf_cwi.iloc[:0]
        )
        clip_ok = (
            len(kept) == len(gdf_clip)
            and all(
                a.equals(b)
                for a, b in zip(
                    kept.sort_values("protected_area_id").geometry,
                    gdf_clip.sort_values("protected_area_id").geometry,
                )
            )
            and len(added) > 0
            and float(
                shapely.intersection(
                    shapely.union_all(added.geometry.values),
                    shapely.union_all(gdf_clip.geometry.values),
                ).area
            )
            < 1.0
        )
    list_res.append(
        (
            "include_implied: defenses_clip.gpkg features unchanged, implied parts added outside them only",
            clip_ok,
        )
    )
    rev = out / "implied" / "base" / "review_implied.gpkg"
    list_res.append(
        (
            f"include_in_review: review_implied.gpkg with reasons {sorted(set(gpd.read_file(rev)['reasons'])) if rev.exists() else None}",
            rev.exists()
            and "uncorroborated" in ",".join(gpd.read_file(rev)["reasons"]),
        )
    )
    proc = run("evaluate", "test_e.toml", "processing.implied_evaluate")
    if proc.returncode == 0:
        sens = pd.read_csv(out / "sensitivity.csv")
        tests_e = pd.read_csv(out / "sensitivity_tests.csv")
        base_r = sens[
            (sens["method"] == "onset")
            & (sens["built_share_min"] == 0.3)
            & (sens["length_m"] == 500.0)
        ].iloc[0]
        t_town = tests_e[
            (tests_e["run_id"] == base_r["run_id"]) & (tests_e["implied_id"] == town_id)
        ].set_index("item_id")
        list_res.append(
            (
                f"sensitivity grid: {len(sens)} runs (4 methods x 2 shares x 2 lengths), "
                f"{int((sens['length_m'] == 1000.0).sum())} at 1000 m with {int(sens.loc[sens['length_m'] == 1000.0, 'implied_protected_areas'].sum())} protected areas",
                len(sens) == 16
                and int(
                    sens.loc[
                        sens["length_m"] == 1000.0, "implied_protected_areas"
                    ].sum()
                )
                > 0,
            )
        )
        list_res.append(
            (
                f"positive metric: recall_built {base_r['recall_built']:.3f} with implied, {base_r['recall_built_explicit']:.3f} "
                f"explicit only, implied_gain {base_r['implied_gain']:.3f}",
                base_r["implied_gain"] > 0
                and base_r["recall_built"] >= base_r["recall_built_explicit"],
            )
        )
        statuses = t_town["status"].to_dict() if len(t_town) else {}
        list_res.append(
            (
                f"evaluation: wet built cells inside the town (rp {t_town['rp'].iloc[0] if len(t_town) else 'n/a'}): {statuses}; "
                f"wrongly_protected_share {base_r['wrongly_protected_share']:.3f}, potential {base_r['potential_share']:.3f}",
                statuses
                == {"SCENE_A": "error", "SCENE_B": "ok", "SCENE_C": "potential"}
                and abs(base_r["wrongly_protected_share"] - 1 / 3) < 0.01,
            )
        )
        run_dir = out / "plots" / "implied" / "sensitivity"
        list_res.append(
            (
                "sensitivity plot listed in its index.html",
                (run_dir / "sensitivity" / "sensitivity.png").exists()
                and "sensitivity" in (run_dir / "index.html").read_text(),
            )
        )
        print(
            "implied_evaluate: "
            + (re.search(r"Timings: (.*)", proc.stderr) or [None, ""])[1]
        )
    else:
        list_res.append(("implied_evaluate runs", False))
    proc = run(
        "protected_areas_force", "test_e.toml", "processing.protected_areas", FORCE=1
    )
    list_res.append(
        (
            "protected_areas.py snaps the rivers again with FORCE (rivers_channel.gpkg) and recomputes all units",
            proc.returncode == 0
            and "Rivers snapped into their channel (rivers_channel.gpkg)" in proc.stderr
            and "Units: 8 to process" in proc.stderr,
        )
    )

print("Timings (s): " + ", ".join(f"{k} {v:.1f}" for k, v in timings.items()))
for label, ok in list_res:
    print(f"{'PASS' if ok else 'FAIL'}  {label}")
sys.exit(0 if all(ok for _, ok in list_res) else 1)
