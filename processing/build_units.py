# Author: Jakub Zapletal
# Date: 2026-11-04
# Builds the river corridor (rivers buffered and clipped to the DTM footprint) and the processing units:
# square cores aligned to the working grid, windows with halo, per-unit domain and memory estimate.

import math
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
    "general.n_workers",
    "general.max_ram_gb_per_worker",
    "general.log_level",
    "country.rivers",
    "country.working_resolution_m",
    "country.corridor_buffer_m",
    "country.unit_core_size_m",
    "country.unit_halo_m",
    "country.bytes_per_cell",
]

if __name__ == "__main__":
    run = gis.Run.env_overrides(CONFIG_PATH=CONFIG_PATH, COUNTRY=COUNTRY)
    cfg = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], KEYS)
    general, country = cfg["general"], cfg["country"]
    work_dir = Path(general["work_dir"]) / run["COUNTRY"]
    out_dir = Path(general["output_dir"]) / run["COUNTRY"]
    log = gis.Run.logger("build_units", out_dir / "logs", general["log_level"])
    gis.Config.warn_placeholders(cfg)

    vrt_path = work_dir / "dtm" / "dtm.vrt"
    footprint_path = out_dir / "dtm_footprint.gpkg"
    if not vrt_path.exists() or not footprint_path.exists():
        raise SystemExit(
            "dtm.vrt or dtm_footprint.gpkg is missing; run prepare_dtm.py first."
        )
    ds = gdal.Open(str(vrt_path))
    wkt = ds.GetProjection()
    gt = ds.GetGeoTransform()
    native_res = abs(gt[1])
    resolution = country["working_resolution_m"]
    if abs(resolution / native_res - round(resolution / native_res)) > 1e-6:
        log.warning(
            "Working resolution %.3f is not a multiple of the DTM pixel %.3f; working cells do not align "
            "with native pixels (averaging then weights partial pixels).",
            resolution,
            native_res,
        )

    # --- Corridor ---
    footprint = gpd.read_file(footprint_path).to_crs(wkt).union_all()
    rivers = gis.Vector.read_merged(country["rivers"], wkt)
    corridor = shapely.make_valid(
        shapely.union_all(
            shapely.buffer(rivers.geometry.values, country["corridor_buffer_m"])
        )
    ).intersection(footprint)
    raster_km2 = ds.RasterXSize * ds.RasterYSize * native_res * abs(gt[5]) / 1e6
    log.info(
        "Corridor: %.1f km2 = %.1f %% of the valid DTM area (%.1f km2) and %.1f %% of the raster extent.",
        corridor.area / 1e6,
        100 * corridor.area / footprint.area,
        footprint.area / 1e6,
        100 * corridor.area / 1e6 / raster_km2,
    )
    gis.Vector.write_gpkg(
        gpd.GeoDataFrame(
            {"area_km2": [corridor.area / 1e6]}, geometry=[corridor], crs=wkt
        ),
        out_dir / "corridor.gpkg",
    )

    # --- Units on the working grid anchored at the DTM origin ---
    n_cols = math.ceil(ds.RasterXSize * native_res / resolution)
    n_rows = math.ceil(ds.RasterYSize * abs(gt[5]) / resolution)
    core_cells = max(2, int(round(country["unit_core_size_m"] / resolution)))
    halo_cells = int(math.ceil(country["unit_halo_m"] / resolution))
    if abs(core_cells * resolution - country["unit_core_size_m"]) > 1e-6:
        log.info(
            "Core size rounded to %.1f m (whole working cells).",
            core_cells * resolution,
        )
    records = []
    for ci in range(math.ceil(n_cols / core_cells)):
        for ri in range(math.ceil(n_rows / core_cells)):
            c0, c1 = ci * core_cells, min((ci + 1) * core_cells, n_cols)
            r0, r1 = ri * core_cells, min((ri + 1) * core_cells, n_rows)
            wc0, wc1 = max(c0 - halo_cells, 0), min(c1 + halo_cells, n_cols)
            wr0, wr1 = max(r0 - halo_cells, 0), min(r1 + halo_cells, n_rows)
            records.append(
                {
                    "unit_id": f"{ci:03d}_{ri:03d}",
                    "core_xmin": gt[0] + c0 * resolution,
                    "core_xmax": gt[0] + c1 * resolution,
                    "core_ymin": gt[3] - r1 * resolution,
                    "core_ymax": gt[3] - r0 * resolution,
                    "win_xmin": gt[0] + wc0 * resolution,
                    "win_xmax": gt[0] + wc1 * resolution,
                    "win_ymin": gt[3] - wr1 * resolution,
                    "win_ymax": gt[3] - wr0 * resolution,
                    "win_col0": wc0,
                    "win_row0": wr0,
                    "win_ncols": wc1 - wc0,
                    "win_nrows": wr1 - wr0,
                    "native_col0": int(round(wc0 * resolution / native_res)),
                    "native_row0": int(round(wr0 * resolution / native_res)),
                    "native_ncols": int(round((wc1 - wc0) * resolution / native_res)),
                    "native_nrows": int(round((wr1 - wr0) * resolution / native_res)),
                }
            )
    units = pd.DataFrame(records)
    core_geoms = shapely.box(
        units.core_xmin, units.core_ymin, units.core_xmax, units.core_ymax
    )
    shapely.prepare(corridor)
    keep = shapely.intersects(corridor, core_geoms)
    units, core_geoms = units[keep].reset_index(drop=True), core_geoms[keep]
    window_geoms = shapely.box(
        units.win_xmin, units.win_ymin, units.win_xmax, units.win_ymax
    )
    log.info(
        "Units: %d of %d cores intersect the corridor (core %.0f m, halo %.0f m, working cell %.1f m).",
        len(units),
        len(records),
        core_cells * resolution,
        halo_cells * resolution,
        resolution,
    )

    # Processing domain: corridor widened by the halo, clipped to each window; cells outside stay nodata.
    domain = shapely.buffer(corridor, country["unit_halo_m"])
    domain_geoms = shapely.intersection(window_geoms, domain)
    units["est_window_cells"] = units.win_ncols * units.win_nrows
    units["est_domain_cells"] = np.round(
        shapely.area(domain_geoms) / resolution**2
    ).astype(np.int64)
    units["est_ram_gb"] = units.est_window_cells * country["bytes_per_cell"] / 1e9
    units = units.sort_values("est_domain_cells", ascending=False).reset_index(
        drop=True
    )
    no_river = ~shapely.intersects(
        shapely.union_all(rivers.geometry.values),
        shapely.box(units.win_xmin, units.win_ymin, units.win_xmax, units.win_ymax),
    )
    if no_river.any():
        log.warning(
            "%d unit windows contain no river; they cannot produce compartments: %s",
            int(no_river.sum()),
            ", ".join(units.unit_id[no_river]),
        )

    max_cells = general["max_ram_gb_per_worker"] * 1e9 / country["bytes_per_cell"]
    if units.est_window_cells.max() > max_cells:
        side = math.floor(math.sqrt(max_cells)) - 2 * halo_cells
        raise SystemExit(
            f"The largest unit window needs {units.est_ram_gb.max():.1f} GB at {country['bytes_per_cell']:.0f} "
            f"bytes per cell, more than max_ram_gb_per_worker = {general['max_ram_gb_per_worker']} GB. "
            f"Use unit_core_size_m <= {max(side, 2) * resolution:.0f} m."
        )
    log.info(
        "Largest window %d x %d cells needs about %.1f GB (bytes_per_cell %.0f); %d workers need about %.1f GB "
        "of RAM in total. Domain cells to process: %.1f million.",
        units.win_ncols.max(),
        units.win_nrows.max(),
        units.est_ram_gb.max(),
        country["bytes_per_cell"],
        general["n_workers"],
        general["n_workers"] * units.est_ram_gb.max(),
        units.est_domain_cells.sum() / 1e6,
    )

    attributes = units[
        ["unit_id", "est_window_cells", "est_domain_cells", "est_ram_gb"]
    ]
    sorted_cores = shapely.box(
        units.core_xmin, units.core_ymin, units.core_xmax, units.core_ymax
    )
    sorted_windows = shapely.box(
        units.win_xmin, units.win_ymin, units.win_xmax, units.win_ymax
    )
    gis.Vector.write_gpkg(
        gpd.GeoDataFrame(attributes, geometry=sorted_cores, crs=wkt),
        out_dir / "units_core.gpkg",
    )
    gis.Vector.write_gpkg(
        gpd.GeoDataFrame(attributes, geometry=sorted_windows, crs=wkt),
        out_dir / "units_window.gpkg",
    )
    gis.Vector.write_gpkg(
        gpd.GeoDataFrame(
            units[["unit_id"]],
            geometry=shapely.intersection(sorted_windows, domain),
            crs=wkt,
        ),
        work_dir / "units_domain.gpkg",
    )
    units.to_csv(out_dir / "units.csv", index=False)
    log.info(
        "Written corridor.gpkg, units_core.gpkg, units_window.gpkg and units.csv to %s.",
        out_dir,
    )
