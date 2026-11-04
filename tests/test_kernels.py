# Author: Jakub Zapletal
# Date: 2026-11-04
# Checks of the numba kernels in gis.py on small synthetic grids: closed ring, gap, diagonal leak,
# nested rings, seed propagation and breaching. Run: pixi run python -m tests.test_kernels

import sys

import numpy as np
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
    ("closed ring: one protected compartment", n == 1 and labels[inside].min() > 0)
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

# Nested rings: the inner compartment is protected up to the higher of the two crests.
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

for label, ok in list_res:
    print(f"{'PASS' if ok else 'FAIL'}  {label}")
sys.exit(0 if all(ok for _, ok in list_res) else 1)
