# Author: Jakub Zapletal
# Date: 2026-11-04
# Checks of the numba kernels in gis.py on small synthetic grids: closed ring, gap, diagonal leak,
# nested rings, seed propagation and breaching; river snapping into a channel of a small DTM file; levee lines onto
# the crest, gap profiles, flood paths, river groups.
# Run: pixi run python -m tests.test_kernels

import sys
import tempfile
from pathlib import Path

import numpy as np
import shapely
from osgeo import osr
from scipy import ndimage

import gis

list_res = []

# Floodplain rising away from a river in row 0; HAND equals the row number times 0.1 m.
rows, cols = 60, 60
hand = np.repeat((np.arange(rows) * 0.1)[:, None], cols, axis=1).astype(np.float32)
seeds = np.zeros((rows, cols), bool)
seeds[0] = True


def ring(r0, r1, c0, c1):
    mask = np.zeros((rows, cols), bool)
    mask[r0, c0 : c1 + 1] = mask[r1, c0 : c1 + 1] = True
    mask[r0 : r1 + 1, c0] = mask[r0 : r1 + 1, c1] = True
    return mask


def defended(barriers):
    surface = hand.copy()
    for mask, crest in barriers:
        surface[mask] = np.maximum(surface[mask], crest)
    return surface


spill_undef, _, _ = gis.Hydro.priority_flood(hand, seeds)
list_res.append(
    (
        "no barriers: spill equals HAND on a monotonic floodplain",
        np.allclose(spill_undef, hand),
    )
)

closed = ring(20, 40, 20, 40)
spill_def, seed_def, order_def = gis.Hydro.priority_flood(
    defended([(closed, 5.0)]), seeds
)
inside = np.zeros((rows, cols), bool)
inside[21:40, 21:40] = True
margin = spill_def - spill_undef
labels, n = ndimage.label((margin >= 0.3) & ~closed)
list_res.append(
    ("closed ring: one protected area", n == 1 and labels[inside].min() > 0)
)
list_res.append(
    (
        "closed ring: protection level equals the lowest crest",
        np.isclose(spill_def[inside].min(), 5.0),
    )
)

gap = closed.copy()
gap[20, 29:32] = False
spill_gap, _, _ = gis.Hydro.priority_flood(defended([(gap, 5.0)]), seeds)
list_res.append(
    (
        "gap in the ring removes protection",
        np.allclose(spill_gap[inside], spill_undef[inside]),
    )
)

# Diamond whose boundary cells touch only diagonally; a 4-connected flood must not pass between them.
rr, cc = np.mgrid[0:rows, 0:cols]
diamond_distance = np.abs(rr - 30) + np.abs(cc - 30)
stair = diamond_distance == 15
stair_inside = diamond_distance < 15
spill_stair, _, _ = gis.Hydro.priority_flood(defended([(stair, 6.0)]), seeds)
list_res.append(
    (
        "no diagonal leak through a staircase barrier",
        np.all(spill_stair[stair_inside] >= 6.0 - 1e-6),
    )
)

# Nested rings: the inner protected area is protected up to the higher of the two crests.
outer = ring(10, 50, 10, 50)
inner = ring(25, 35, 25, 35)
core = np.zeros((rows, cols), bool)
core[26:35, 26:35] = True
for outer_crest, inner_crest in [(5.5, 7.0), (7.0, 5.5)]:
    spill_nested, _, _ = gis.Hydro.priority_flood(
        defended([(outer, outer_crest), (inner, inner_crest)]), seeds
    )
    list_res.append(
        (
            f"nested rings ({outer_crest}, {inner_crest}): inner level = higher crest",
            np.isclose(spill_nested[core].min(), max(outer_crest, inner_crest)),
        )
    )

# Seed propagation and entry: two rivers (left column low, right column), ring entered from the side.
seeds2 = np.zeros((rows, cols), bool)
seeds2[:, 0] = True
seeds2[:, -1] = True
side = (
    np.abs(np.arange(cols) - 0.0)[None, :].repeat(rows, axis=0).astype(np.float32) * 0.1
)
surface2 = side.copy()
spill2, seed2, order2 = gis.Hydro.priority_flood(surface2, seeds2)
seed_cols = seed2 % cols
list_res.append(
    ("seed index points to a seed cell", np.all(seeds2.ravel()[seed2.ravel()]))
)
list_res.append(
    (
        "cells are reached in non-decreasing spill order",
        np.all(np.diff(spill2.ravel()[np.argsort(order2.ravel())]) >= -1e-6),
    )
)
list_res.append(
    (
        "every cell keeps the seed of its minimal path (left river)",
        np.all(seed_cols[:, :30] == 0),
    )
)
entry = np.unravel_index(
    np.argmin(np.where(inside, order_def, np.iinfo(order_def.dtype).max)), (rows, cols)
)
list_res.append(
    (
        "entry cell of the ring is adjacent to its lowest barrier cell",
        spill_def[entry] == spill_def[inside].min() and seeds.ravel()[seed_def[entry]],
    )
)

# Breaching: a pit in a slope drains to the stream through the carved path; elevations are unchanged.
dem = (np.arange(rows)[:, None] * 0.1 + np.zeros((1, cols))).astype(np.float32)
dem[30:33, 30:33] = -5.0
dem_copy = dem.copy()
stream = np.zeros((rows, cols), bool)
stream[0, :] = True
reached = gis.Hydro.breach_to_streams(dem, stream)
list_res.append(
    (
        "breaching: every cell reaches a stream cell",
        np.all(reached >= 0) and np.all(stream.ravel()[reached.ravel()]),
    )
)
list_res.append(
    ("breaching does not modify the input DTM", np.array_equal(dem, dem_copy))
)
dem_nan = dem.copy()
dem_nan[:, 45] = np.nan
reached_nan = gis.Hydro.breach_to_streams(dem_nan, stream)
list_res.append(("nodata cells are not routed", np.all(reached_nan[:, 45] == -1)))

# Bug B1: a flat-bottomed channel sloping west with the river rasterized on its centre line. Without burning,
# paths run along channel cells parallel to the river and reach it only at the lowest (west) end.
yy, xx = np.mgrid[0:41, 0:200].astype(float)
valley = (
    100 + 0.002 * xx + 0.01 * np.abs(yy - 20) - 2.0 * (np.abs(yy - 20) <= 3)
).astype(np.float32)
valley[np.abs(yy - 20) <= 3] = (98 + 0.002 * xx)[np.abs(yy - 20) <= 3]
river = np.zeros(valley.shape, bool)
river[20] = True
far_col = lambda r: np.abs(r % 200 - np.arange(200)[None, :])[5:15, 50:150].max()
list_res.append(
    (
        "B1 reproduced: without burning paths reach the river far downstream",
        far_col(gis.Hydro.breach_to_streams(valley, river, 0.0)) > 20,
    )
)
list_res.append(
    (
        "B1 fixed: with burning every cell reaches the river within 2 columns",
        far_col(gis.Hydro.breach_to_streams(valley, river, 5.0, 5)) <= 2,
    )
)

# River snapping: a 2 m DTM with a channel along a sine; the river line lies 60 m beside it.
tmp = Path(tempfile.mkdtemp(prefix="defenses_kernels_"))
nx, ny, res = 1000, 500, 2.0
gx, gy = np.meshgrid((np.arange(nx) + 0.5) * res, 1000 - (np.arange(ny) + 0.5) * res)
centre = 500 + 40 * np.sin(gx / 150)
z = 100 + 0.002 * np.abs(gy - centre) - 2.0 * (np.abs(gy - centre) < 6) + 0.0001 * gx
srs = osr.SpatialReference()
srs.ImportFromEPSG(32633)
gis.Raster.write_array(
    tmp / "dtm.tif",
    z.astype(np.float32),
    (0, res, 0, 1000, 0, -res),
    srs.ExportToWkt(),
    nodata=-9999,
)
xs = np.linspace(5, 1995, 60)
river = shapely.linestrings(np.c_[xs, 500 + 40 * np.sin(xs / 150) + 60])
guide = shapely.linestrings(np.c_[xs, 500 + 40 * np.sin(xs / 150) + 3])
snap_p = {
    k.split(".")[1]: v
    for k, v in gis.Config.DEFAULTS.items()
    if k.startswith("river_snap.")
}


def channel_offset(line):
    xy = shapely.get_coordinates(
        shapely.line_interpolate_point(line, np.arange(0, line.length, 2.0))
    )
    return np.abs(xy[:, 1] - (500 + 40 * np.sin(xy[:, 0] / 150)))


out_g, st_g = gis.Hydro.snap_to_channel([river], tmp / "dtm.tif", [guide], snap_p)
out_n, _ = gis.Hydro.snap_to_channel([river], tmp / "dtm.tif", None, snap_p)
out_d, _ = gis.Hydro.snap_to_channel(
    [river], tmp / "dtm.tif", None, {**snap_p, "max_shift_without_guide_m": 80.0}
)
off_g, off_d = channel_offset(out_g[0]), channel_offset(out_d[0])
list_res.append(
    (
        f"snap_to_channel: with a guide line the river lies in the channel (median {np.median(off_g):.1f} m, max "
        f"{off_g.max():.1f} m from 60 m); without one it stays put; DTM only (80 m) median {np.median(off_d):.1f} m",
        np.median(off_g) < 3.0
        and off_g.max() < 7.0
        and out_n[0].equals(river)
        and np.median(off_d) < 10.0
        and st_g["guided_share"] == 1.0,
    )
)

# Levee crest: 1 m DTM with a 1.8 m high embankment (4 m top) along a sine and a 40 m gap without it; a levee
# line 9 m beside the crest, one on the crest axis and one on flat ground.
nx, ny = 600, 400
gx, gy = np.meshgrid(np.arange(nx) + 0.5, 400 - (np.arange(ny) + 0.5))
axis_y = 200 + 30 * np.sin(gx / 80)
across = np.abs(gy - axis_y) * np.cos(np.arctan(30 / 80 * np.cos(gx / 80)))
z = 100 + 0.002 * gx + np.clip(1.8 - np.maximum(across - 2, 0) / 2, 0, None)
z[:, 280:320] = 100 + 0.002 * gx[:, 280:320]
gis.Raster.write_array(
    tmp / "levee.tif",
    z.astype(np.float32),
    (0, 1.0, 0, 400, 0, -1.0),
    srs.ExportToWkt(),
    nodata=-9999,
)
xs = np.arange(10, 590, 5.0)
crest_xy = np.c_[xs, 200 + 30 * np.sin(xs / 80)]
crest_p = {
    k.split(".")[1]: v
    for k, v in gis.Config.DEFAULTS.items()
    if k.startswith("levee_snap.")
}
moved, df_ridge, _ = gis.Hydro.snap_to_crest(
    [
        shapely.linestrings(crest_xy - [0, 9]),
        shapely.linestrings(crest_xy),
        shapely.linestrings([[50, 50], [550, 60]]),
    ],
    str(tmp / "levee.tif"),
    crest_p,
)
pts = shapely.get_coordinates(shapely.segmentize(moved[0], 2.0))
dev = np.abs(pts[:, 1] - (200 + 30 * np.sin(pts[:, 0] / 80)))
on_axis = shapely.get_coordinates(shapely.segmentize(moved[1], 2.0))
dev_axis = np.abs(on_axis[:, 1] - (200 + 30 * np.sin(on_axis[:, 0] / 80)))
list_res.append(
    (
        f"snap_to_crest: a levee line 9 m beside the crest moves onto it (median {np.median(dev):.1f} m, max "
        f"{dev.max():.1f} m from the axis, top 2 m either side), one on the axis stays within "
        f"{dev_axis.max():.1f} m, one on flat ground stays (ridge share {df_ridge['ridge_share'].iloc[2]:.0%})",
        np.median(dev) <= 2.0
        and dev.max() <= 3.0
        and dev_axis.max() <= 1.0
        and moved[2].equals(shapely.linestrings([[50, 50], [550, 60]]))
        and df_ridge["ridge_share"].iloc[2] == 0,
    )
)
end_low, along_min = gis.Hydro.bridge_profile(
    [shapely.linestrings(crest_xy[[52, 64]]), shapely.linestrings(crest_xy[[18, 30]])],
    str(tmp / "levee.tif"),
    15.0,
)
list_res.append(
    (
        f"bridge_profile: across the gap without the embankment the terrain drops {end_low[0] - along_min[0]:.1f} m "
        f"below the crest, along the embankment {end_low[1] - along_min[1]:.2f} m",
        end_low[0] - along_min[0] > 1.5 and end_low[1] - along_min[1] < 0.1,
    )
)

# Flood paths: water from column 0 reaches the cells behind a wall in column 10 through its gap in row 15.
surface = np.zeros((20, 20), np.float32)
surface[:, 10] = 5.0
surface[15, 10] = 0.0
left = np.zeros((20, 20), bool)
left[:, 0] = True
_, _, order = gis.Hydro.priority_flood(surface, left)
event = np.zeros((20, 20), np.int8)
event[:, 10] = 2
event[15, 10] = 3
code, cell = gis.Hydro.first_event_on_path(
    order, left, event, [5 * 20 + 18, 5 * 20 + 5]
)
list_res.append(
    (
        f"first_event_on_path: the path behind the wall leads back through the gap ({code[0]}, cell "
        f"{divmod(int(cell[0]), 20)}), the one in front to its seed ({code[1]})",
        code[0] == 3
        and cell[0] == 15 * 20 + 10
        and code[1] == -1
        and cell[1] % 20 == 0,
    )
)

# River groups: a main stem of five 10 km reaches, two tributaries and a minor stream.
reaches = [
    shapely.linestrings(c)
    for c in [
        [[0, 0], [10e3, 0]],
        [[10e3, 0], [20e3, 0]],
        [[20e3, 0], [30e3, 0]],
        [[30e3, 0], [40e3, 0]],
        [[40e3, 0], [50e3, 0]],
        [[20e3, 5e3], [20e3, 0]],
        [[30e3, 8e3], [30e3, 0]],
        [[40e3, 3e3], [40e3, 0]],
    ]
]
groups = gis.Vector.river_groups(
    reaches, [5000, 3000, 2900, 1200, 1100, 900, 150, 50], 2.0, 30000.0, 100.0, 0.01
)
list_res.append(
    (
        f"river_groups: the main stem joins up to 30 km, a jump in upstream area above 2 or a tributary starts a "
        f"group, a minor stream is group 0 ({groups.tolist()})",
        groups.tolist() == [1, 1, 1, 2, 2, 3, 4, 0],
    )
)

for label, ok in list_res:
    print(f"{'PASS' if ok else 'FAIL'}  {label}")
sys.exit(0 if all(ok for _, ok in list_res) else 1)
