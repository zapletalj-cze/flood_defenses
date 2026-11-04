# Author: Jakub Zapletal
# Date: 2026-11-04
# Prepares Copernicus Land Monitoring Service data for the context zones from local files (the CLMS portal needs a
# login and no download source could be verified): CORINE class index -> CLC code through the legend, verified on known
# pixels and clipped to the corridor; Urban Atlas merged into one single-layer GeoPackage in the DTM CRS.

import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from osgeo import gdal

import gis

# --- Run variables (no command-line parameters; DEFENSES_<NAME> environment variables override them) ---
CONFIG_PATH = "config/defenses.toml"
COUNTRY = "XX"

KEYS = [
    "general.work_dir",
    "general.output_dir",
    "general.log_level",
    "country.corine",
    "country.corine_legend",
    "country.corine_legend_index_field",
    "country.corine_legend_code_field",
    "country.corine_check",
    "country.lulc.urban_atlas",
    "country.lulc.urban_atlas_layer",
    "country.lulc.urban_atlas_class_field",
]

if __name__ == "__main__":
    run = gis.Run.env_overrides(CONFIG_PATH=CONFIG_PATH, COUNTRY=COUNTRY)
    cfg = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], KEYS)
    general, country = cfg["general"], cfg["country"]
    work_dir = Path(general["work_dir"]) / run["COUNTRY"]
    out_dir = Path(general["output_dir"]) / run["COUNTRY"]
    log = gis.Run.logger("download_clms", out_dir / "logs", general["log_level"])
    gis.Config.warn_placeholders(cfg)
    lulc = country["lulc"]
    if not country["corine"] and not lulc["urban_atlas"]:
        raise SystemExit(
            "Neither countries.<code>.corine nor countries.<code>.lulc.urban_atlas is set; nothing to prepare. "
            "Download the data from the CLMS portal (login needed) and list the local files."
        )
    corridor_path, vrt_path = out_dir / "corridor.gpkg", work_dir / "dtm" / "dtm.vrt"
    if not corridor_path.exists() or not vrt_path.exists():
        raise SystemExit(
            "corridor.gpkg or dtm.vrt is missing; run prepare_dtm.py and build_units.py first."
        )
    wkt = gdal.Open(str(vrt_path)).GetProjection()
    gdf_corr = gpd.read_file(corridor_path).to_crs(wkt)
    record = {"corine": {}, "urban_atlas": {}}

    # --- CORINE: legend mapping, verified on known pixels, clipped to the corridor ---
    if country["corine"]:
        index_field = country["corine_legend_index_field"]
        code_field = country["corine_legend_code_field"]
        if not country["corine_legend"] or not index_field or not code_field:
            raise SystemExit(
                "Set countries.<code>.corine_legend and its index and code fields: CORINE rasters store class "
                "indices, not three-digit codes."
            )
        df_leg = pd.read_csv(country["corine_legend"], sep=None, engine="python")
        list_abs = [f for f in [index_field, code_field] if f not in df_leg.columns]
        if list_abs:
            raise SystemExit(
                f"Legend columns {list_abs} not found; columns: {list(df_leg.columns)}"
            )
        df_leg = df_leg[[index_field, code_field]].dropna()
        mapping = dict(
            zip(df_leg[index_field].astype(int), df_leg[code_field].astype(int))
        )
        vrt = gis.Raster.build_vrt(
            country["corine"], work_dir / "clms" / "corine_index.vrt"
        )
        ds = gdal.Open(vrt)
        band, gt = ds.GetRasterBand(1), ds.GetGeoTransform()
        if not country["corine_check"]:
            raise SystemExit(
                "countries.<code>.corine_check is empty: give at least one pixel with a known CLC code "
                "[[x, y, code]] in the raster CRS, so that the legend mapping can be verified."
            )
        list_chk = []
        for x, y, code in country["corine_check"]:
            col = int(np.floor((x - gt[0]) / gt[1]))
            row = int(np.floor((y - gt[3]) / gt[5]))
            if not (0 <= col < ds.RasterXSize and 0 <= row < ds.RasterYSize):
                raise SystemExit(
                    f"Check pixel ({x}, {y}) lies outside the CORINE raster."
                )
            value = int(band.ReadAsArray(col, row, 1, 1)[0, 0])
            mapped = mapping.get(value)
            list_chk.append(
                {"x": x, "y": y, "value": value, "mapped": mapped, "expected": code}
            )
            if mapped != int(code):
                raise SystemExit(
                    f"CORINE legend check failed at ({x}, {y}): raster value {value} maps to {mapped}, "
                    f"expected code {code}. Check corine_legend and its fields."
                )
        log.info(
            "CORINE legend verified on %d known pixel(s): %s.",
            len(list_chk),
            ", ".join(f"{c['value']} -> {c['mapped']}" for c in list_chk),
        )
        # Corridor bounds in the raster CRS, one pixel of margin.
        area = gpd.GeoSeries(
            [shapely.segmentize(gdf_corr.union_all(), 1000.0)], crs=wkt
        ).to_crs(ds.GetProjection())
        xmin, ymin, xmax, ymax = area.total_bounds
        bounds = (xmin - gt[1], ymin + gt[5], xmax + gt[1], ymax - gt[5])
        values, window_gt = gis.Raster.read_native(ds, bounds)
        codes = gis.Raster.map_codes(values, mapping)
        unmapped = np.isfinite(values) & (codes == 0)
        if unmapped.any():
            log.warning(
                "%d CORINE pixels have values without a legend entry: %s.",
                int(unmapped.sum()),
                sorted(set(values[unmapped].astype(int).tolist()))[:20],
            )
        gis.Raster.write_array(
            out_dir / "corine_codes.tif",
            codes.astype(np.uint16),
            window_gt,
            ds.GetProjection(),
            nodata=0,
        )
        present, counts = np.unique(codes[codes > 0], return_counts=True)
        log.info(
            "corine_codes.tif: %d x %d pixels; most frequent codes: %s.",
            codes.shape[1],
            codes.shape[0],
            ", ".join(
                f"{c} ({n})"
                for c, n in sorted(zip(present, counts), key=lambda t: -t[1])[:10]
            ),
        )
        record["corine"] = {
            "files": gis.Config.as_list(country["corine"]),
            "legend": country["corine_legend"],
            "checks": list_chk,
            "output": str(out_dir / "corine_codes.tif"),
        }

    # --- Urban Atlas: single-layer copies, merged, reprojected, polygons touching the corridor ---
    if lulc["urban_atlas"]:
        field = lulc["urban_atlas_class_field"]
        if not field:
            raise SystemExit("Set countries.<code>.lulc.urban_atlas_class_field.")
        list_path = []
        for path in gis.Config.as_list(lulc["urban_atlas"]):
            list_name = gis.Vector.layer_names(path)
            if len(list_name) == 1:
                list_path.append(path)
            elif lulc["urban_atlas_layer"]:
                list_path.append(
                    gis.Vector.single_layer_copy(
                        path,
                        lulc["urban_atlas_layer"],
                        work_dir / "clms" / f"{Path(path).stem}_layer.gpkg",
                    )
                )
            else:
                raise SystemExit(
                    f"{path} has layers {list_name}; set countries.<code>.lulc.urban_atlas_layer."
                )
        gdf_ua = gis.Vector.read_merged(
            list_path, wkt, bbox=tuple(gdf_corr.total_bounds), keep_columns=[field]
        )
        if field not in gdf_ua.columns:
            raise SystemExit(f"Field {field!r} not found in the Urban Atlas files.")
        touching, _ = gdf_corr.sindex.query(
            gdf_ua.geometry.values, predicate="intersects"
        )
        gdf_ua = gdf_ua.iloc[np.unique(touching)].reset_index(drop=True)
        codes = pd.to_numeric(gdf_ua[field].astype(str).str.strip(), errors="coerce")
        if codes.isna().any():
            log.warning(
                "%d Urban Atlas polygons have a non-numeric class code; they are kept with code 0.",
                int(codes.isna().sum()),
            )
        gdf_ua["ua_code"] = codes.fillna(0).astype(int)
        gis.Vector.write_gpkg(
            gdf_ua[["ua_code", "source", "geometry"]], out_dir / "urban_atlas.gpkg"
        )
        log.info(
            "urban_atlas.gpkg: %d polygons touching the corridor; classes: %s.",
            len(gdf_ua),
            ", ".join(
                f"{c} ({n})"
                for c, n in gdf_ua["ua_code"].value_counts().head(12).items()
            ),
        )
        record["urban_atlas"] = {
            "files": gis.Config.as_list(lulc["urban_atlas"]),
            "field": field,
            "polygons": int(len(gdf_ua)),
            "output": str(out_dir / "urban_atlas.gpkg"),
        }
    (out_dir / "clms_sources.json").write_text(
        json.dumps(record, indent=2, default=str)
    )
