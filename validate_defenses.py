# Author: Jakub Zapletal
# Date: 2026-11-03
# Compares derived protected areas with reference protected-area polygons.
# Reports area-based recall and precision and writes the covered share of each reference polygon.

import argparse
import logging
from pathlib import Path

import geopandas as gpd
import pandas as pd
from osgeo import ogr

ogr.UseExceptions()

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("validate_defenses")


def read_vectors(paths, crs):
    frames = []
    for path in paths:
        ds = ogr.Open(str(path))
        if ds is None or ds.GetLayerCount() != 1:
            raise SystemExit(f"{path}: expected exactly one layer.")
        ds = None
        gdf = gpd.read_file(path)
        if gdf.crs is None:
            log.warning("%s has no CRS, assuming %s.", path, crs)
            gdf = gdf.set_crs(crs)
        gdf = gdf.to_crs(crs)
        gdf["source"] = Path(path).stem
        frames.append(gdf)
    gdf = gpd.GeoDataFrame(pd.concat(frames, ignore_index=True), geometry="geometry", crs=crs)
    return gdf[gdf.geometry.notna() & ~gdf.geometry.is_empty].reset_index(drop=True)


parser = argparse.ArgumentParser(description="Validate protected areas against reference polygons.")
parser.add_argument("--result", nargs="+", required=True, help="protected_areas.gpkg file(s).")
parser.add_argument("--reference", nargs="+", required=True, help="Reference protected-area polygon file(s).")
parser.add_argument("--extent", nargs="+", help="Optional extent polygon file(s); both layers are clipped to it.")
parser.add_argument("--output", required=True, type=Path, help="Output GeoPackage with per-reference coverage.")
args = parser.parse_args()

crs = gpd.read_file(args.result[0], rows=0).crs
if crs is None or not crs.is_projected:
    raise SystemExit("The result must be in a projected CRS to compute areas.")
result = read_vectors(args.result, crs)
reference = read_vectors(args.reference, crs)
if args.extent:
    extent = read_vectors(args.extent, crs).union_all()
    result = gpd.clip(result, extent)
    reference = gpd.clip(reference, extent)
reference = reference[reference.area > 0].reset_index(drop=True)

result_union = result.union_all()
reference_union = reference.union_all()
intersection = result_union.intersection(reference_union).area
recall = intersection / reference_union.area if reference_union.area > 0 else float("nan")
precision = intersection / result_union.area if result_union.area > 0 else float("nan")
log.info("Reference area: %.0f m2, result area: %.0f m2, intersection: %.0f m2.",
         reference_union.area, result_union.area, intersection)
log.info("Recall: %.3f, precision: %.3f.", recall, precision)

reference["area_m2"] = reference.area
reference["covered_m2"] = reference.geometry.intersection(result_union).area
reference["covered_share"] = reference["covered_m2"] / reference["area_m2"]
log.info("Reference polygons: %d, fully missed: %d, covered over 50 %%: %d.", len(reference),
         (reference["covered_m2"] == 0).sum(), (reference["covered_share"] > 0.5).sum())
args.output.parent.mkdir(parents=True, exist_ok=True)
reference.to_file(args.output, driver="GPKG")
