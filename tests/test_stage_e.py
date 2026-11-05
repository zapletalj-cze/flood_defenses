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
import shapely
from osgeo import gdal, osr
from shapely import affinity

import gis

gdal.UseExceptions()

EPSG = 32633  # synthetic test data only
X0, Y0 = 600000.0, 5600000.0  # lower-left corner of the test world
WIDTH_M, HEIGHT_M, NATIVE = 8000, 3200, 1.0
RES = 2.0  # implied.snap.resolution_m of the test configuration
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
gis.Vector.write_gpkg(gdf_rn, out / "rivers_sources.gpkg")
gis.Vector.write_gpkg(gdf_rn, data / "rivers.gpkg")
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


def write_config(name, units="[]", plot_all="[]", max_per_type=200, ratio=1.5):
    (root / name).write_text(f"""[general]
work_dir = "{(root / 'work').as_posix()}"
output_dir = "{(root / 'output').as_posix()}"
n_workers = 2
gdal_cachemax_mb = 256
warp_memory_limit_mb = 256
gdal_threads_per_worker = 1
log_level = "INFO"

[countries.T5]
rivers = "{(data / 'rivers.gpkg').as_posix()}"
water_polygons = "{(data / 'water.gpkg').as_posix()}"

[implied]
units = {units}

[implied.snap]
step_m = 5.0
resolution_m = {RES}
max_offset_m = 30.0
min_window_m = 15.0
cost_epsilon = 0.01
chunk_length_m = 5000.0
chunk_overlap_m = 200.0
max_length_ratio = {ratio}
wide_river_m = 40.0
direction_check_m = 200.0
max_nodata_share = 0.5
chainage_origin = "downstream"

[implied.plots]
enabled = true
dpi = 80
max_per_type = {max_per_type}
sample_seed = 42
show_max = 2
plot_all_in_units = {plot_all}
""")


def run(label, config, **variables):
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
        [sys.executable, "-m", "processing.snap_rivers"],
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
            f"only water sources are snapped (is_source): {sorted(gdf_snap.index)}",
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

print("Timings (s): " + ", ".join(f"{k} {v:.1f}" for k, v in timings.items()))
for label, ok in list_res:
    print(f"{'PASS' if ok else 'FAIL'}  {label}")
sys.exit(0 if all(ok for _, ok in list_res) else 1)
