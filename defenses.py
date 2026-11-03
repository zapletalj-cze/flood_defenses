# Author: Jakub Zapletal
# Date: 2026-11-03
# Derives areas protected by flood levees ("Defenses") from terrain, levee lines and the river network.
# Protection is the difference between priority-flood spill levels above the river with and without levees.

import argparse
import heapq
import logging
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from numba import njit
from osgeo import gdal, ogr, osr
from scipy import ndimage
from shapely.geometry import LineString, Point, box
from shapely.ops import nearest_points

gdal.UseExceptions()
ogr.UseExceptions()
osr.UseExceptions()

NODATA = -9999.0
# Endpoints closer than this to another line are considered already connected.
TOUCH_TOLERANCE_M = 0.01

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("defenses")


def read_vectors(paths, crs):
    frames = []
    for path in paths:
        ds = ogr.Open(str(path))
        if ds is None or ds.GetLayerCount() != 1:
            raise SystemExit(f"{path}: expected exactly one layer.")
        ds = None
        gdf = gpd.read_file(path)
        if gdf.crs is None:
            log.warning("%s has no CRS, assuming the DTM CRS.", path)
            gdf = gdf.set_crs(crs)
        gdf = gdf.to_crs(crs)
        gdf["source"] = Path(path).stem
        frames.append(gdf[["source", "geometry"]])
    gdf = gpd.GeoDataFrame(pd.concat(frames, ignore_index=True), geometry="geometry", crs=crs)
    return gdf[gdf.geometry.notna() & ~gdf.geometry.is_empty].reset_index(drop=True)


def array_to_ds(array, geotransform, wkt, path="", driver="MEM", nodata=NODATA):
    gdal_type = gdal.GDT_Int32 if np.issubdtype(array.dtype, np.integer) else gdal.GDT_Float32
    options = ["COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=IF_SAFER"] if driver == "GTiff" else []
    if driver == "GTiff" and gdal_type == gdal.GDT_Float32:
        options.append("PREDICTOR=3")
    ds = gdal.GetDriverByName(driver).Create(str(path), array.shape[1], array.shape[0], 1, gdal_type, options)
    ds.SetGeoTransform(geotransform)
    ds.SetProjection(wkt)
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(nodata)
    band.WriteArray(np.where(np.isnan(array), nodata, array) if gdal_type == gdal.GDT_Float32 else array)
    return ds


def rasterize(geometries, values, shape, geotransform, wkt, all_touched):
    vector_ds = gdal.GetDriverByName("MEM").Create("", 0, 0, 0, gdal.GDT_Unknown)
    srs = osr.SpatialReference(wkt=wkt)
    layer = vector_ds.CreateLayer("features", srs, ogr.wkbUnknown)
    layer.CreateField(ogr.FieldDefn("value", ogr.OFTInteger))
    for geometry, value in zip(geometries, values):
        feature = ogr.Feature(layer.GetLayerDefn())
        feature.SetGeometry(ogr.CreateGeometryFromWkb(shapely.to_wkb(geometry)))
        feature.SetField("value", int(value))
        layer.CreateFeature(feature)
    raster_ds = array_to_ds(np.zeros(shape, np.int32), geotransform, wkt, nodata=-1)
    gdal.RasterizeLayer(raster_ds, [1], layer, options=["ATTRIBUTE=value", f"ALL_TOUCHED={str(all_touched).upper()}"])
    return raster_ds.ReadAsArray()


@njit(cache=True)
def downstream_stream_cell(dem, stream):
    # Complete breaching: priority flood from stream cells, depressions are carved along the flood path
    # (flats get an epsilon gradient), then D8 steepest descent on the conditioned copy.
    rows, cols = dem.shape
    z = dem.ravel().copy()
    is_stream = stream.ravel()
    n = rows * cols
    parent = np.full(n, -1, np.int64)
    seen = np.zeros(n, np.bool_)
    eps = 1e-5
    di = np.array([-1, -1, -1, 0, 0, 1, 1, 1])
    dj = np.array([-1, 0, 1, -1, 1, -1, 0, 1])
    dist = np.array([1.4142135623730951, 1.0, 1.4142135623730951, 1.0, 1.0, 1.4142135623730951, 1.0,
                     1.4142135623730951])
    heap = [(0.0, np.int64(0))]
    heap.pop()
    for c in range(n):
        if is_stream[c] and not np.isnan(z[c]):
            seen[c] = True
            heap.append((z[c], np.int64(c)))
    heapq.heapify(heap)
    while len(heap) > 0:
        _, c = heapq.heappop(heap)
        ci = c // cols
        cj = c % cols
        for k in range(8):
            ni = ci + di[k]
            nj = cj + dj[k]
            if ni < 0 or nj < 0 or ni >= rows or nj >= cols:
                continue
            nb = ni * cols + nj
            if seen[nb] or np.isnan(z[nb]):
                continue
            seen[nb] = True
            parent[nb] = c
            if z[nb] <= z[c]:
                if z[c] - z[nb] < eps:
                    z[nb] = z[c] + eps
                else:
                    target = z[nb] - eps
                    p = c
                    while p >= 0 and z[p] > target:
                        z[p] = target
                        target -= eps
                        p = parent[p]
            heapq.heappush(heap, (z[nb], nb))

    receiver = np.full(n, -1, np.int64)
    for c in range(n):
        if not seen[c] or is_stream[c]:
            continue
        ci = c // cols
        cj = c % cols
        best = 0.0
        for k in range(8):
            ni = ci + di[k]
            nj = cj + dj[k]
            if ni < 0 or nj < 0 or ni >= rows or nj >= cols:
                continue
            nb = ni * cols + nj
            if not seen[nb]:
                continue
            slope = (z[c] - z[nb]) / dist[k]
            if slope > best:
                best = slope
                receiver[c] = nb
        if receiver[c] < 0:
            receiver[c] = parent[c]

    reached = np.full(n, -1, np.int64)
    for c in np.argsort(z):
        if not seen[c]:
            continue
        reached[c] = c if is_stream[c] else reached[receiver[c]]
    return reached.reshape(rows, cols)


@njit(cache=True)
def priority_flood(surface, seeds):
    # Lowest level at which water from the seeds reaches each cell (4-connected).
    rows, cols = surface.shape
    spill = np.full((rows, cols), np.inf)
    di = np.array([-1, 1, 0, 0])
    dj = np.array([0, 0, -1, 1])
    heap = [(0.0, 0, 0)]
    heap.pop()
    for i in range(rows):
        for j in range(cols):
            if seeds[i, j]:
                spill[i, j] = 0.0
                heap.append((0.0, i, j))
    heapq.heapify(heap)
    while len(heap) > 0:
        level, i, j = heapq.heappop(heap)
        if level > spill[i, j]:
            continue
        for k in range(4):
            ni = i + di[k]
            nj = j + dj[k]
            if ni < 0 or nj < 0 or ni >= rows or nj >= cols or np.isnan(surface[ni, nj]):
                continue
            new_level = max(level, surface[ni, nj])
            if new_level < spill[ni, nj]:
                spill[ni, nj] = new_level
                heapq.heappush(heap, (new_level, ni, nj))
    return spill


parser = argparse.ArgumentParser(description="Derive areas protected by flood levees (Defenses layer).")
parser.add_argument("--dtm", nargs="+", required=True, help="DTM GeoTIFF tile(s) in one CRS.")
parser.add_argument("--levees", nargs="+", required=True, help="Levee line file(s), one layer each.")
parser.add_argument("--rivers", nargs="+", required=True, help="Modelled river line file(s), one layer each.")
parser.add_argument("--extent", nargs="+", help="Optional pilot extent polygon file(s).")
parser.add_argument("--output-dir", required=True, type=Path)
parser.add_argument("--gap-tolerance-m", type=float, default=20.0)
parser.add_argument("--crest-search-radius-m", type=float, default=5.0)
parser.add_argument("--barrier-removal-width-m", type=float, default=10.0)
parser.add_argument("--min-margin-m", type=float, default=0.3)
parser.add_argument("--max-level-m", type=float, default=20.0)
parser.add_argument("--min-area-m2", type=float, default=10000.0)
parser.add_argument("--resolution", type=float, help="Working resolution in DTM units (default: native).")
parser.add_argument("--keep-intermediate", action="store_true", help="Write intermediate rasters as GeoTIFF.")
args = parser.parse_args()
args.output_dir.mkdir(parents=True, exist_ok=True)

# --- DTM ---
reference_srs = None
for tile in args.dtm:
    srs = osr.SpatialReference(wkt=gdal.Open(str(tile)).GetProjection())
    if reference_srs is None:
        reference_srs = srs
    elif not srs.IsSame(reference_srs):
        raise SystemExit(f"DTM tile {tile} has a different CRS than {args.dtm[0]}. Reproject the tiles first.")
vrt = gdal.BuildVRT("/vsimem/dtm.vrt", [str(t) for t in args.dtm])
wkt = vrt.GetProjection()
crs = reference_srs.ExportToWkt()
gt_vrt = vrt.GetGeoTransform()
native_res = abs(gt_vrt[1])

if args.extent:
    xmin, ymin, xmax, ymax = read_vectors(args.extent, crs).total_bounds
    col0 = max(int(np.floor((xmin - gt_vrt[0]) / gt_vrt[1])), 0)
    col1 = min(int(np.ceil((xmax - gt_vrt[0]) / gt_vrt[1])), vrt.RasterXSize)
    row0 = max(int(np.floor((ymax - gt_vrt[3]) / gt_vrt[5])), 0)
    row1 = min(int(np.ceil((ymin - gt_vrt[3]) / gt_vrt[5])), vrt.RasterYSize)
    if col1 <= col0 or row1 <= row0:
        raise SystemExit("The extent does not overlap the DTM.")
    native_ds = gdal.Translate("", vrt, format="MEM", srcWin=[col0, row0, col1 - col0, row1 - row0])
else:
    native_ds = gdal.Translate("", vrt, format="MEM")
gt_native = native_ds.GetGeoTransform()
band = native_ds.GetRasterBand(1)
dtm_native = band.ReadAsArray().astype(np.float64)
if band.GetNoDataValue() is not None:
    dtm_native[dtm_native == band.GetNoDataValue()] = np.nan
log.info("Native DTM grid: %d x %d cells, %.3f m.", dtm_native.shape[1], dtm_native.shape[0], native_res)

coarse = args.resolution is not None and args.resolution > native_res * (1 + 1e-9)
if args.resolution is not None and args.resolution < native_res * (1 - 1e-9):
    raise SystemExit(f"Working resolution {args.resolution} is finer than the native {native_res}.")
if coarse:
    work_ds = gdal.Warp("", array_to_ds(dtm_native, gt_native, wkt), format="MEM", xRes=args.resolution,
                        yRes=args.resolution, resampleAlg="average", srcNodata=NODATA, dstNodata=NODATA)
    dtm = work_ds.ReadAsArray().astype(np.float64)
    dtm[dtm == NODATA] = np.nan
    gt = work_ds.GetGeoTransform()
    log.info("Working grid: %d x %d cells, %.3f m (average resampling).", dtm.shape[1], dtm.shape[0], gt[1])
else:
    dtm = dtm_native
    gt = gt_native
res = abs(gt[1])
cell_area = abs(gt[1] * gt[5])
shape = dtm.shape
valid = ~np.isnan(dtm)
grid_box = box(gt[0], gt[3] + gt[5] * shape[0], gt[0] + gt[1] * shape[1], gt[3])

# --- Barriers: merge, bridge gaps ---
levees = read_vectors(args.levees, crs)
levees = levees[levees.intersects(grid_box)].explode(index_parts=False).reset_index(drop=True)
levees = levees[levees.geom_type == "LineString"].reset_index(drop=True)
levees["bridged"] = False
log.info("Levee lines in the processing extent: %d.", len(levees))

bridges = {}
for line_idx, line in enumerate(levees.geometry):
    for endpoint in (Point(line.coords[0]), Point(line.coords[-1])):
        candidates = [c for c in levees.sindex.query(endpoint, predicate="dwithin", distance=args.gap_tolerance_m)
                      if c != line_idx]
        if not candidates:
            continue
        distances = levees.geometry.iloc[candidates].distance(endpoint).to_numpy()
        if distances.min() <= TOUCH_TOLERANCE_M:
            continue
        other_idx = candidates[int(np.argmin(distances))]
        other = levees.geometry.iloc[other_idx]
        other_ends = [Point(other.coords[0]), Point(other.coords[-1])]
        end_distances = [endpoint.distance(e) for e in other_ends]
        # Prefer endpoint-to-endpoint, otherwise connect to the nearest point on the other line.
        target = other_ends[int(np.argmin(end_distances))] if min(end_distances) <= args.gap_tolerance_m \
            else nearest_points(other, endpoint)[0]
        key = tuple(sorted([(round(endpoint.x, 3), round(endpoint.y, 3)), (round(target.x, 3), round(target.y, 3))]))
        sources = sorted({levees.at[line_idx, "source"], levees.at[other_idx, "source"]})
        bridges[key] = ("|".join(sources), LineString([endpoint, target]))
bridge_lines = gpd.GeoDataFrame({"source": [b[0] for b in bridges.values()], "bridged": True},
                                geometry=[b[1] for b in bridges.values()], crs=crs)
log.info("Bridged gaps: %d.", len(bridge_lines))

# --- Barriers: rasterize, crest elevation ---
footprint_radius = int(np.ceil(args.crest_search_radius_m / native_res))
yy, xx = np.mgrid[-footprint_radius:footprint_radius + 1, -footprint_radius:footprint_radius + 1]
footprint = np.hypot(yy, xx) * native_res <= args.crest_search_radius_m
crest_native = ndimage.maximum_filter(np.where(np.isnan(dtm_native), -np.inf, dtm_native), footprint=footprint,
                                      mode="constant", cval=-np.inf)
crest_native[~np.isfinite(crest_native)] = np.nan

levee_cells = rasterize(levees.geometry, np.ones(len(levees)), shape, gt, wkt, True) > 0
crest = np.full(shape, np.nan)
if coarse:
    levee_native = rasterize(levees.geometry, np.ones(len(levees)), dtm_native.shape, gt_native, wkt, True) > 0
    crest_ds = gdal.Warp("", array_to_ds(np.where(levee_native, crest_native, np.nan), gt_native, wkt),
                         format="MEM", outputBounds=grid_box.bounds, width=shape[1], height=shape[0],
                         resampleAlg="max", srcNodata=NODATA, dstNodata=NODATA)
    crest_coarse = crest_ds.ReadAsArray().astype(np.float64)
    crest[levee_cells] = np.where(crest_coarse == NODATA, np.nan, crest_coarse)[levee_cells]
else:
    crest[levee_cells] = crest_native[levee_cells]

bridge_cells = np.zeros(shape, bool)
if len(bridge_lines):
    bridge_id = rasterize(bridge_lines.geometry, np.arange(1, len(bridge_lines) + 1), shape, gt, wkt, True)
    bridge_cells = (bridge_id > 0) & ~levee_cells
    ends = np.array([[line.coords[0], line.coords[-1]] for line in bridge_lines.geometry])
    end_cols = np.clip(((ends[..., 0] - gt_native[0]) / gt_native[1]).astype(int), 0, dtm_native.shape[1] - 1)
    end_rows = np.clip(((ends[..., 1] - gt_native[3]) / gt_native[5]).astype(int), 0, dtm_native.shape[0] - 1)
    end_crest = crest_native[end_rows, end_cols]
    # A missing crest at one end takes the value of the other end.
    end_crest = np.where(np.isnan(end_crest), end_crest[:, ::-1], end_crest)
    rows, cols = np.nonzero(bridge_cells)
    seg = bridge_id[rows, cols] - 1
    px = gt[0] + (cols + 0.5) * gt[1]
    py = gt[3] + (rows + 0.5) * gt[5]
    a = ends[seg, 0]
    ab = ends[seg, 1] - a
    t = np.clip(((px - a[:, 0]) * ab[:, 0] + (py - a[:, 1]) * ab[:, 1]) / (ab ** 2).sum(axis=1), 0, 1)
    crest[rows, cols] = end_crest[seg, 0] + t * (end_crest[seg, 1] - end_crest[seg, 0])

barrier = levee_cells | bridge_cells
# Close diagonal-only contacts so that water cannot pass between two barrier cells touching at a corner.
added = 0
while True:
    diag = barrier[:-1, :-1] & barrier[1:, 1:] & ~barrier[:-1, 1:] & ~barrier[1:, :-1]
    anti = barrier[:-1, 1:] & barrier[1:, :-1] & ~barrier[:-1, :-1] & ~barrier[1:, 1:]
    if not diag.any() and not anti.any():
        break
    height = np.where(valid, dtm, -np.inf)
    i, j = np.nonzero(diag)
    first = height[i, j + 1] >= height[i + 1, j]
    barrier[np.where(first, i, i + 1), np.where(first, j + 1, j)] = True
    i, j = np.nonzero(anti)
    first = height[i, j] >= height[i + 1, j + 1]
    barrier[np.where(first, i, i + 1), np.where(first, j, j + 1)] = True
    added += diag.sum() + anti.sum()
log.info("Barrier cells: %d (%d added for 4-connectivity).", barrier.sum(), added)

missing = barrier & np.isnan(crest)
if missing.any() and (~np.isnan(crest)).any():
    nearest = ndimage.distance_transform_edt(np.isnan(crest), return_distances=False, return_indices=True)
    crest[missing] = crest[nearest[0][missing], nearest[1][missing]]
crest[barrier & np.isnan(crest)] = dtm[barrier & np.isnan(crest)]

# --- Rivers ---
rivers = read_vectors(args.rivers, crs)
rivers = rivers[rivers.intersects(grid_box)].reset_index(drop=True)
stream = (rasterize(rivers.geometry, np.ones(len(rivers)), shape, gt, wkt, False) > 0) & valid
if not stream.any():
    raise SystemExit("No river cells inside the processing extent.")
crossings = barrier & stream
barrier_full = barrier.copy()
barrier &= ~stream
log.info("Stream cells: %d, river-barrier crossing cells: %d.", stream.sum(), crossings.sum())

# --- Undefended DTM ---
dtm_undef = dtm.copy()
if barrier_full.any():
    removal = (ndimage.distance_transform_edt(~barrier_full, sampling=res) <= args.barrier_removal_width_m) & valid
    fill_ds = array_to_ds(np.where(removal, np.nan, dtm), gt, wkt)
    mask_ds = array_to_ds((~removal & valid).astype(np.int32), gt, wkt, nodata=-1)
    max_search = int(np.ceil(4 * args.barrier_removal_width_m / res)) + 10
    gdal.FillNodata(fill_ds.GetRasterBand(1), mask_ds.GetRasterBand(1), maxSearchDist=max_search, smoothingIterations=0)
    filled = fill_ds.ReadAsArray().astype(np.float64)
    unfilled = removal & (filled == NODATA)
    if unfilled.any():
        log.warning("%d removal cells could not be filled, keeping the original DTM there.", unfilled.sum())
    dtm_undef[removal & ~unfilled] = filled[removal & ~unfilled]

# --- HAND on the undefended DTM ---
reached = downstream_stream_cell(dtm_undef, stream)
unreached = valid & (reached < 0)
if unreached.any():
    log.warning("%d cells are not connected to any stream cell, using the nearest stream cell.", unreached.sum())
    nearest = ndimage.distance_transform_edt(~stream, return_distances=False, return_indices=True)
    reached[unreached] = (nearest[0] * shape[1] + nearest[1])[unreached]
river_ref = np.where(valid, dtm_undef.ravel()[np.maximum(reached, 0)].reshape(shape), np.nan)
hand_undef = dtm_undef - river_ref

# --- Spill levels ---
surface_def = hand_undef.copy()
surface_def[barrier] = np.maximum(crest[barrier] - river_ref[barrier], hand_undef[barrier])
spill_undef = priority_flood(hand_undef, stream)
spill_def = priority_flood(surface_def, stream)
spill_undef[~np.isfinite(spill_undef)] = np.nan
spill_def[~np.isfinite(spill_def)] = np.nan
margin = spill_def - spill_undef

# Stream cells in areas fully enclosed by barriers make those areas unprotected.
enclosed, _ = ndimage.label(~barrier_full & valid)
edge_labels = np.unique(np.concatenate([enclosed[0], enclosed[-1], enclosed[:, 0], enclosed[:, -1]]))
for label in np.setdiff1d(np.unique(enclosed[stream]), np.append(edge_labels, 0)):
    region = enclosed == label
    log.warning("Stream cells inside a levee-enclosed area: %d stream cells, area %.0f m2, centre %.1f %.1f; "
                "this area is treated as unprotected.", (region & stream).sum(), region.sum() * cell_area,
                gt[0] + (np.nonzero(region)[1].mean() + 0.5) * gt[1],
                gt[3] + (np.nonzero(region)[0].mean() + 0.5) * gt[5])

# --- Protected compartments ---
protected = (margin >= args.min_margin_m) & (spill_undef <= args.max_level_m) & ~barrier
labels, n_labels = ndimage.label(protected)
counts = np.bincount(labels.ravel())
keep = np.nonzero(counts * cell_area >= args.min_area_m2)[0]
keep = keep[keep > 0]
relabel = np.zeros(n_labels + 1, np.int32)
relabel[keep] = np.arange(1, len(keep) + 1)
labels = relabel[labels]
ids = np.arange(1, len(keep) + 1)
log.info("Protected compartments: %d of %d candidates kept (min area %.0f m2).", len(ids), n_labels, args.min_area_m2)

polygons = gdal.GetDriverByName("MEM").Create("", 0, 0, 0, gdal.GDT_Unknown)
layer = polygons.CreateLayer("compartments", osr.SpatialReference(wkt=wkt), ogr.wkbPolygon)
layer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))
label_ds = array_to_ds(labels, gt, wkt, nodata=0)
mask_ds = array_to_ds((labels > 0).astype(np.int32), gt, wkt, nodata=-1)
gdal.Polygonize(label_ds.GetRasterBand(1), mask_ds.GetRasterBand(1), layer, 0)
features = [(f.GetField("id"), shapely.from_wkb(bytes(f.GetGeometryRef().ExportToWkb()))) for f in layer]
compartments = gpd.GeoDataFrame({"id": [f[0] for f in features]}, geometry=[f[1] for f in features], crs=crs)
compartments = compartments.dissolve(by="id").reset_index()
stats = pd.DataFrame({
    "id": ids,
    "area_m2": counts[keep] * cell_area,
    "protection_level_m": ndimage.minimum(spill_def, labels, ids) if len(ids) else [],
    "margin_median_m": ndimage.median(margin, labels, ids) if len(ids) else [],
    "margin_max_m": ndimage.maximum(margin, labels, ids) if len(ids) else [],
    "river_ref_median_m": ndimage.median(river_ref, labels, ids) if len(ids) else [],
})
compartments = compartments.merge(stats, on="id", how="right")
compartments["rp"] = pd.Series([None] * len(compartments), dtype="Float64")
compartments = compartments[["id", "area_m2", "protection_level_m", "margin_median_m", "margin_max_m",
                             "river_ref_median_m", "rp", "geometry"]]
compartments.to_file(args.output_dir / "protected_areas.gpkg", driver="GPKG")

# --- Barrier lines with crest attributes ---
barrier_lines = gpd.GeoDataFrame(pd.concat([levees[["source", "bridged", "geometry"]], bridge_lines],
                                           ignore_index=True), geometry="geometry", crs=crs)
crest_attr = np.where(barrier, crest, np.nan)
crest_hand = crest_attr - river_ref
crest_min = []
crest_hand_min = []
for line in barrier_lines.geometry:
    points = shapely.line_interpolate_point(line, np.append(np.arange(0, line.length, res / 2), line.length))
    xy = shapely.get_coordinates(points)
    cols = ((xy[:, 0] - gt[0]) / gt[1]).astype(int)
    rows = ((xy[:, 1] - gt[3]) / gt[5]).astype(int)
    inside = (cols >= 0) & (rows >= 0) & (cols < shape[1]) & (rows < shape[0])
    values = crest_attr[rows[inside], cols[inside]]
    hand_values = crest_hand[rows[inside], cols[inside]]
    crest_min.append(np.nanmin(values) if np.isfinite(values).any() else np.nan)
    crest_hand_min.append(np.nanmin(hand_values) if np.isfinite(hand_values).any() else np.nan)
barrier_lines["crest_min_m"] = crest_min
barrier_lines["crest_hand_min_m"] = crest_hand_min
barrier_lines.to_file(args.output_dir / "barriers.gpkg", driver="GPKG")

if args.keep_intermediate:
    for name, array in [("dtm_undefended", dtm_undef), ("hand_undef", hand_undef), ("river_ref", river_ref),
                        ("spill_undef", spill_undef), ("spill_def", spill_def), ("margin", margin)]:
        array_to_ds(array, gt, wkt, path=args.output_dir / f"{name}.tif", driver="GTiff")
log.info("Done. Outputs written to %s.", args.output_dir)
