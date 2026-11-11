# Author: Jakub Zapletal
# Date: 2026-11-03
# Downloads layers of the USACE National Levee Database (NLD) from its public OGC WFS service with GDAL.
# Every layer is written to its own single-layer GeoPackage, optionally restricted to an extent.

import argparse
import logging
import re
import time
from pathlib import Path

import geopandas as gpd
import pandas as pd
from osgeo import gdal, ogr

gdal.UseExceptions()
ogr.UseExceptions()

NLD_WFS_URL = "https://levees.sec.usace.army.mil/mapserver/public/ows"

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("download_nld")

parser = argparse.ArgumentParser(
    description="Download National Levee Database layers from the NLD WFS service."
)
parser.add_argument("--output-dir", required=True, type=Path)
parser.add_argument(
    "--url", default=NLD_WFS_URL, help="WFS endpoint (default: public NLD OGC service)."
)
parser.add_argument(
    "--layers",
    nargs="+",
    help="Layer names or case-insensitive substrings to download (default: all layers).",
)
parser.add_argument(
    "--extent",
    nargs="+",
    help="Optional polygon file(s); only features in their bounding box.",
)
parser.add_argument(
    "--page-size", type=int, default=1000, help="Features per WFS request."
)
parser.add_argument("--retries", type=int, default=4)
parser.add_argument(
    "--list", action="store_true", help="List the available layers and exit."
)
args = parser.parse_args()

gdal.SetConfigOption("OGR_WFS_PAGING_ALLOWED", "ON")
gdal.SetConfigOption("OGR_WFS_PAGE_SIZE", str(args.page_size))
gdal.SetConfigOption("GDAL_HTTP_TIMEOUT", "120")

source = None
for attempt in range(args.retries + 1):
    try:
        source = gdal.OpenEx(f"WFS:{args.url}", gdal.OF_VECTOR)
        break
    except RuntimeError as error:
        if attempt == args.retries:
            raise SystemExit(f"Cannot open the WFS service {args.url}: {error}")
        log.warning("Opening the WFS service failed (%s), retrying.", error)
        time.sleep(2 ** (attempt + 1))

names = [source.GetLayer(i).GetName() for i in range(source.GetLayerCount())]
log.info("Layers offered by the service: %d.", len(names))
if args.list:
    for name in names:
        print(name)
    raise SystemExit(0)

if args.layers:
    selected = [
        n for n in names if any(p == n or p.lower() in n.lower() for p in args.layers)
    ]
    if not selected:
        raise SystemExit(
            f"No layer matches {args.layers}. Use --list to see the available layers."
        )
else:
    selected = names

translate_options = {}
if args.extent:
    frames = [gpd.read_file(p) for p in args.extent]
    extent = gpd.GeoDataFrame(
        pd.concat([f.to_crs(frames[0].crs) for f in frames]), crs=frames[0].crs
    )
    translate_options = {
        "spatFilter": list(extent.total_bounds),
        "spatSRS": extent.crs.to_wkt(),
    }

args.output_dir.mkdir(parents=True, exist_ok=True)
for name in selected:
    stem = re.sub(r"[^A-Za-z0-9_-]+", "_", name.split(":")[-1]).strip("_")
    path = args.output_dir / f"{stem}.gpkg"
    for attempt in range(args.retries + 1):
        try:
            result = gdal.VectorTranslate(
                str(path),
                source,
                format="GPKG",
                layers=[name],
                layerName=stem,
                accessMode="overwrite",
                **translate_options,
            )
            count = result.GetLayer(0).GetFeatureCount()
            result = None
            break
        except RuntimeError as error:
            if attempt == args.retries:
                log.error("Layer %s failed: %s", name, error)
                count = None
                break
            log.warning("Layer %s failed (%s), retrying.", name, error)
            time.sleep(2 ** (attempt + 1))
    if count == 0:
        path.unlink(missing_ok=True)
        log.info("%s: no features, skipped.", name)
    elif count:
        log.info("%s: %d features -> %s", name, count, path)
log.info("Done.")
