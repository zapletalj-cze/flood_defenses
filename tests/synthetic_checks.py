# Author: Jakub Zapletal
# Date: 2026-11-03
# Runs defenses.py and validate_defenses.py on synthetic terrain and checks the expected behaviour
# (empty levees, closed ring, bridged and wide gaps, diagonal leakage, coarse grid, river inside a ring, CRS check).

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import geopandas as gpd
import numpy as np
from osgeo import gdal, osr
from shapely.geometry import LineString, box

gdal.UseExceptions()

REPO = Path(__file__).resolve().parents[1]
RES = 2.0
COLS, ROWS = 300, 200
X0, Y0 = 500000.0, 5500000.0
EPSG = 32633  # synthetic test data only

parser = argparse.ArgumentParser(description="Synthetic checks of the Defenses scripts.")
parser.add_argument("--workdir", type=Path, help="Keep test data and outputs here (default: temporary folder).")
args = parser.parse_args()
workdir = args.workdir or Path(tempfile.mkdtemp(prefix="defenses_checks_"))
workdir.mkdir(parents=True, exist_ok=True)
print(f"Working directory: {workdir}")

# Local coordinates: x along the valley, y across it, origin at the lower-left corner.
x = (np.arange(COLS) + 0.5) * RES
y = ROWS * RES - (np.arange(ROWS) + 0.5) * RES
xx, yy = np.meshgrid(x, y)
river_y = 100.0
base = 100.0 + 0.002 * xx + 0.005 * np.abs(yy - river_y) - 2.0 * (np.abs(yy - river_y) <= 4.0)

# Square ring levee (local coords), the embankment crest is 2 m off the digitised line.
ring = [(220, 180), (380, 180), (380, 340), (220, 340), (220, 180)]
ring_line = LineString(ring)
# Diamond levee with vertices on cell centres; the embankment occupies only the exactly diagonal cells.
diamond = LineString([(221, 261), (301, 341), (381, 261), (301, 181), (221, 261)])


def inward(line, offset=2.0):
    # Embankment crest 2 m inside the digitised levee line.
    shift = {220: offset, 380: -offset, 180: offset, 340: -offset}
    return LineString([(px + shift.get(px, 0), py + shift.get(py, 0)) for px, py in line.coords])


def to_world(geometry):
    return gpd.GeoSeries([geometry]).affine_transform([1, 0, 0, 1, X0, Y0 - ROWS * RES]).iloc[0]


def embankment(lines, height=2.5, half_width=3.0):
    dtm = base.copy()
    for line in lines:
        distance = gpd.GeoSeries(gpd.points_from_xy(xx.ravel(), yy.ravel())).distance(line).to_numpy()
        dtm += height * (distance.reshape(dtm.shape) <= half_width)
    return dtm


def write_dtm(dtm, name, epsg=EPSG):
    paths = []
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    half = COLS // 2
    for k, (c0, c1) in enumerate([(0, half), (half, COLS)]):
        path = workdir / f"{name}_tile{k}.tif"
        ds = gdal.GetDriverByName("GTiff").Create(str(path), c1 - c0, ROWS, 1, gdal.GDT_Float32)
        ds.SetGeoTransform([X0 + c0 * RES, RES, 0, Y0, 0, -RES])
        ds.SetProjection(srs.ExportToWkt())
        ds.GetRasterBand(1).SetNoDataValue(-9999)
        ds.GetRasterBand(1).WriteArray(dtm[:, c0:c1])
        ds = None
        paths.append(str(path))
    return paths


def write_lines(lines, name):
    path = workdir / f"{name}.gpkg"
    gpd.GeoDataFrame(geometry=gpd.GeoSeries([to_world(g) for g in lines], crs=EPSG)).to_file(path, driver="GPKG")
    return str(path)


def run(name, dtm, levees, extra=(), rivers=None):
    out = workdir / name
    command = [sys.executable, str(REPO / "defenses.py"), "--dtm", *dtm, "--levees", levees,
               "--rivers", rivers or river_file, "--output-dir", str(out), "--keep-intermediate", *extra]
    proc = subprocess.run(command, capture_output=True, text=True)
    (workdir / f"{name}.log").write_text(proc.stdout + proc.stderr)
    if proc.returncode != 0:
        return proc, None
    return proc, gpd.read_file(out / "protected_areas.gpkg")


river_file = write_lines([LineString([(600, river_y), (0, river_y)])], "river")
dtm_ring = write_dtm(embankment([inward(ring_line)]), "ring")
results = []

proc, areas = run("empty", dtm_ring, write_lines([], "levees_empty"))
margin = gdal.Open(str(workdir / "empty" / "margin.tif")).ReadAsArray()
margin = margin[margin != -9999]
results.append(("empty levee input gives zero margin", proc.returncode == 0 and np.abs(margin).max() == 0
                and len(areas) == 0))

proc, areas = run("ring", dtm_ring, write_lines([ring_line], "levees_ring"))
results.append(("closed ring gives one compartment", areas is not None and len(areas) == 1))
if areas is not None and len(areas) == 1:
    print(areas.drop(columns="geometry").to_string(index=False))

# Ring split into two lines with two gaps in the middle of the left and right sides.
def split_ring(gap):
    half = gap / 2
    return [LineString([(220, 260 - half), (220, 180), (380, 180), (380, 260 - half)]),
            LineString([(380, 260 + half), (380, 340), (220, 340), (220, 260 + half)])]


for gap, expected, label in [(10, 1, "gap below tolerance is bridged"), (30, 0, "gap above tolerance removes protection")]:
    dtm = write_dtm(embankment([inward(g) for g in split_ring(gap)]), f"gap{gap}")
    proc, areas = run(f"gap{gap}", dtm, write_lines(split_ring(gap), f"levees_gap{gap}"))
    results.append((f"{label} ({gap} m)", areas is not None and len(areas) == expected))

dtm_diamond = embankment([])
on_diamond = np.abs(np.abs(xx - 301) + np.abs(yy - 261) - 80) < 1e-6
dtm_diamond[on_diamond] += 2.5
proc, areas = run("diagonal", write_dtm(dtm_diamond, "diamond"), write_lines([diamond], "levees_diamond"))
added = [line for line in (proc.stdout + proc.stderr).splitlines() if "added for 4-connectivity" in line]
print(added[0].split("INFO ")[-1] if added else "no barrier log line")
results.append(("no diagonal leakage through a staircase levee", areas is not None and len(areas) == 1))

proc, areas = run("ring_coarse", dtm_ring, str(workdir / "levees_ring.gpkg"), ["--resolution", "4"])
results.append(("closed ring at 4 m working resolution", areas is not None and len(areas) == 1))

extent_file = workdir / "extent.gpkg"
gpd.GeoDataFrame(geometry=[to_world(box(100, 0, 500, 400))], crs=EPSG).to_file(extent_file, driver="GPKG")
proc, areas = run("ring_extent", dtm_ring, str(workdir / "levees_ring.gpkg"), ["--extent", str(extent_file)])
results.append(("closed ring with pilot extent", areas is not None and len(areas) == 1))

rivers_inside = write_lines([LineString([(600, river_y), (0, river_y)]), LineString([(260, 300), (340, 220)])],
                            "river_inside")
proc, areas = run("river_inside", dtm_ring, str(workdir / "levees_ring.gpkg"), rivers=rivers_inside)
warned = "inside a levee-enclosed area" in proc.stderr
results.append(("stream inside the ring: unprotected and warned", areas is not None and len(areas) == 0 and warned))

wrong_crs = write_dtm(base, "wrong_crs", epsg=32634)
proc, _ = run("crs_mismatch", [dtm_ring[0], wrong_crs[1]], str(workdir / "levees_ring.gpkg"))
results.append(("different tile CRS stops with a message", proc.returncode != 0 and "different CRS" in proc.stderr))

reference = workdir / "reference.gpkg"
gpd.GeoDataFrame(geometry=[to_world(box(224, 184, 376, 336))], crs=EPSG).to_file(reference, driver="GPKG")
proc = subprocess.run([sys.executable, str(REPO / "validate_defenses.py"), "--result",
                       str(workdir / "ring" / "protected_areas.gpkg"), "--reference", str(reference),
                       "--output", str(workdir / "validation.gpkg")], capture_output=True, text=True)
recall_line = [line for line in proc.stderr.splitlines() if "Recall" in line]
print(recall_line[0].split("INFO ")[-1] if recall_line else proc.stderr)
coverage = gpd.read_file(workdir / "validation.gpkg")["covered_share"] if proc.returncode == 0 else [0]
results.append(("validation recall of the ring above 0.95", proc.returncode == 0 and min(coverage) > 0.95))

for label, ok in results:
    print(f"{'PASS' if ok else 'FAIL'}  {label}")
sys.exit(0 if all(ok for _, ok in results) else 1)
