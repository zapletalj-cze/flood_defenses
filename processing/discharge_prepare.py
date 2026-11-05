# Author: Jakub Zapletal
# Date: 2026-11-05
# Stage F: daily discharge per river cell of each source (time convention per source, EFAS 6-hourly aggregated to
# days), annual maxima of complete hydrological years, Gumbel fits by L-moments; compact outputs and fit plots.

import json
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import shapely
from osgeo import gdal

import gis
from plots import COLOURS, LINE_WIDTH, SOURCE_NAMES, Plots

gdal.UseExceptions()

# --- Run variables (no command-line parameters; DEFENSES_<NAME> environment variables override them) ---
CONFIG_PATH = "config/defenses.toml"
COUNTRY = "XX"
RUN_ID = "prepare"  # plots in <output_dir>/<country>/plots/discharge/<RUN_ID>
SHOW_PLOTS = False  # show up to implied.plots.show_max figures per type while running (needs a display)

KEYS = [
    "general.work_dir",
    "general.output_dir",
    "general.log_level",
    "discharge.sources",
    "discharge.raw_dir",
    "discharge.aggregation",
    "discharge.hydrological_year_start_month",
    "discharge.min_years",
    "discharge.bbox_buffer_deg",
    *[
        f"discharge.{s}.{k}"
        for s in ("glofas", "efas")
        for k in (
            "step_hours",
            "valid_time_marks",
            "netcdf_variable",
            "upstream_area_variable",
            "upstream_area_to_km2",
            "min_upstream_area_km2",
        )
    ],
    *[
        f"implied.plots.{k}"
        for k in (
            "enabled",
            "dpi",
            "max_per_type",
            "sample_seed",
            "show_max",
            "plot_all_in_units",
        )
    ],
]

if __name__ == "__main__":
    run = gis.Run.env_overrides(
        CONFIG_PATH=CONFIG_PATH, COUNTRY=COUNTRY, RUN_ID=RUN_ID, SHOW_PLOTS=SHOW_PLOTS
    )
    cfg = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], KEYS)
    general, dis = cfg["general"], cfg["discharge"]
    work_dir = Path(general["work_dir"]) / run["COUNTRY"]
    out_dir = Path(general["output_dir"]) / run["COUNTRY"]
    log = gis.Run.logger("discharge_prepare", out_dir / "logs", general["log_level"])
    gis.Config.warn_placeholders(cfg)
    raw_dir = Path(dis["raw_dir"]) if dis["raw_dir"] else work_dir / "discharge" / "raw"
    fp_path = out_dir / "dtm_footprint.gpkg"
    if not fp_path.exists():
        raise SystemExit("dtm_footprint.gpkg missing; run prepare_dtm.py first.")
    west, south, east, north = gpd.read_file(fp_path).to_crs(4326).total_bounds
    b = dis["bbox_buffer_deg"]
    bbox = (west - b, south - b, east + b, north + b)
    perf, fits_for_plots = {}, {}

    for src in dis["sources"]:
        t0 = time.perf_counter()
        c = dis[src]
        src_dir = raw_dir / src
        list_file = sorted(
            p
            for p in src_dir.glob("*")
            if p.suffix in (".grib", ".nc", ".grb", ".grib2")
            and not p.stem.startswith("upstream_area")
        )
        up_files = sorted(src_dir.glob("upstream_area.*"))
        if not list_file or not up_files:
            log.warning(
                "%s: no discharge files or no upstream area file in %s; source skipped.",
                src,
                src_dir,
            )
            continue
        if c["upstream_area_to_km2"] <= 0:
            raise SystemExit(
                f"discharge.{src}.upstream_area_to_km2 is not set; check the units of the upstream area file."
            )

        # --- River cells: upstream area at least min_upstream_area_km2 inside the area ---
        up_name = (
            f'NETCDF:"{up_files[0]}":{c["upstream_area_variable"]}'
            if c["upstream_area_variable"]
            else str(up_files[0])
        )
        ds_up = gdal.Open(up_name)
        ugt = ds_up.GetGeoTransform()
        area = (
            ds_up.GetRasterBand(1).ReadAsArray().astype(np.float64)
            * c["upstream_area_to_km2"]
        )
        nodata = ds_up.GetRasterBand(1).GetNoDataValue()
        if nodata is not None:
            area[ds_up.GetRasterBand(1).ReadAsArray() == nodata] = np.nan
        ds_up = None
        lon_c = ugt[0] + (np.arange(area.shape[1]) + 0.5) * ugt[1]
        lat_c = ugt[3] + (np.arange(area.shape[0]) + 0.5) * ugt[5]
        inside = (
            (lon_c[None, :] >= bbox[0])
            & (lon_c[None, :] <= bbox[2])
            & (lat_c[:, None] >= bbox[1])
            & (lat_c[:, None] <= bbox[3])
        )
        rows, cols = np.nonzero(inside & (area >= c["min_upstream_area_km2"]))
        if not len(rows):
            log.warning(
                "%s: no cell reaches min_upstream_area_km2 in the area; source skipped.",
                src,
            )
            continue
        lon, lat = lon_c[cols], lat_c[rows]
        cell_id = rows.astype(np.int64) * area.shape[1] + cols

        # --- Discharge of these cells, aggregated to days file by file (days by the time convention of the source);
        # the steps of the last day of a file wait for the next file, since a day can span two files ---
        steps = 24 // c["step_hours"]
        how = dis["aggregation"]
        if steps == 1 and how == "daily_max":
            log.info(
                "%s is daily: daily_max is not available, daily_mean is used.", src
            )
            how = "daily_mean"
        list_day, list_q, list_stamp = [], [], []
        carry_day = np.empty(0, "datetime64[D]")
        carry_q = np.empty((0, len(lon)), np.float32)
        for path in list_file + [None]:
            if path is None:
                # end of the files: the carried day is aggregated as it is
                step_day, step_q = carry_day, carry_q
                done = np.ones(len(step_day), bool)
            else:
                name = (
                    f'NETCDF:"{path}":{c["netcdf_variable"]}'
                    if c["netcdf_variable"] and path.suffix == ".nc"
                    else str(path)
                )
                ds = gdal.Open(name)
                gt = ds.GetGeoTransform()
                fc = np.floor((lon - gt[0]) / gt[1] + 1e-6).astype(int)
                fr = np.floor((lat - gt[3]) / gt[5] + 1e-6).astype(int)
                if (
                    (fc < 0).any()
                    or (fr < 0).any()
                    or (fc >= ds.RasterXSize).any()
                    or (fr >= ds.RasterYSize).any()
                ):
                    raise SystemExit(
                        f"{path.name} does not cover all river cells of the area; download a larger area."
                    )
                x0, y0 = int(fc.min()), int(fr.min())
                w, h = int(fc.max() - x0 + 1), int(fr.max() - y0 + 1)
                stamps = gis.Discharge.band_times(ds)
                list_stamp += [stamps.min(), stamps.max()]
                values = np.empty((ds.RasterCount, len(lon)), np.float32)
                for k in range(ds.RasterCount):
                    band = ds.GetRasterBand(k + 1)
                    block = band.ReadAsArray(x0, y0, w, h).astype(np.float32)
                    if band.GetNoDataValue() is not None:
                        block[block == np.float32(band.GetNoDataValue())] = np.nan
                    values[k] = block[fr - y0, fc - x0]
                ds = None
                step_day = np.concatenate(
                    [
                        carry_day,
                        gis.Discharge.period_days(stamps, c["valid_time_marks"]),
                    ]
                )
                step_q = np.concatenate([carry_q, values])
                done = step_day < step_day.max()
            if done.any():
                d_, q_ = gis.Discharge.aggregate_daily(
                    step_day[done], step_q[done], how, steps
                )
                list_day.append(d_)
                list_q.append(q_)
            carry_day, carry_q = step_day[~done], step_q[~done]
        days, q = np.concatenate(list_day), np.concatenate(list_q)
        # the first and last valid times with their days, to check the time convention on real files
        t_first, t_last = min(list_stamp), max(list_stamp)
        log.info(
            "%s: valid times %s to %s UTC -> days %s to %s (valid_time_marks = %s, %d steps per day).",
            src,
            np.datetime64(int(t_first), "s"),
            np.datetime64(int(t_last), "s"),
            *gis.Discharge.period_days([t_first, t_last], c["valid_time_marks"]),
            c["valid_time_marks"],
            steps,
        )
        if len(np.unique(days)) < len(days):
            log.warning(
                "%s: %d days occur in more than one group of files (overlapping files?); one value per day is kept.",
                src,
                len(days) - len(np.unique(days)),
            )
        # complete calendar of days; missing days stay NaN
        full = np.arange(days[0], days[-1] + np.timedelta64(1, "D"))
        series = np.full((len(full), len(lon)), np.nan, np.float32)
        series[(days - full[0]).astype(int)] = q
        n_missing = int(np.isnan(series).all(axis=1).sum())
        if n_missing:
            log.warning(
                "%s: %d days without data between %s and %s.",
                src,
                n_missing,
                full[0],
                full[-1],
            )

        # --- Annual maxima of complete hydrological years and Gumbel fits ---
        years, maxima = gis.Discharge.annual_maxima(
            full, series, dis["hydrological_year_start_month"]
        )
        n_years, xi, alpha, rmse = gis.Discharge.gumbel_fits(maxima, dis["min_years"])
        keep = n_years >= dis["min_years"]
        log.info(
            "%s: %d river cells, %d with at least %d complete years (of %d hydrological years); aggregation %s.",
            src,
            len(lon),
            int(keep.sum()),
            dis["min_years"],
            len(years),
            how,
        )
        np.savez_compressed(
            out_dir / f"discharge_{src}.npz",
            cell_id=cell_id[keep],
            lon=lon[keep],
            lat=lat[keep],
            upstream_area_km2=area[rows, cols][keep],
            dates=full.astype(np.int64),
            values=series[:, keep],
            years=years,
            maxima=maxima[:, keep],
        )
        gdf_cells = gpd.GeoDataFrame(
            {
                "cell_id": cell_id[keep],
                "upstream_area_km2": area[rows, cols][keep],
                "n_years": n_years[keep],
                "gumbel_xi": xi[keep],
                "gumbel_alpha": alpha[keep],
                "fit_rmse": rmse[keep],
                "lon": lon[keep],
                "lat": lat[keep],
            },
            geometry=shapely.points(lon[keep], lat[keep]),
            crs="EPSG:4326",
        )
        gis.Vector.write_gpkg(gdf_cells, out_dir / f"discharge_cells_{src}.gpkg")
        fits_for_plots[src] = (gdf_cells, years, maxima[:, keep])
        perf[src] = {
            "files": len(list_file),
            "cells": int(keep.sum()),
            "days": len(full),
            "seconds": round(time.perf_counter() - t0, 2),
            "peak_rss_mb": round(gis.Run.peak_rss_mb()),
        }
        log.info(
            "%s: discharge_%s.npz and discharge_cells_%s.gpkg in %.1f s.",
            src,
            src,
            src,
            perf[src]["seconds"],
        )

    if not fits_for_plots:
        raise SystemExit("No source prepared.")

    # --- Fit plots per cell: annual maxima at Gringorten positions and the Gumbel curve ---
    t = time.perf_counter()
    settings = cfg["implied"]["plots"]
    if settings["enabled"]:
        plots = Plots(
            out_dir / "plots" / "discharge" / run["RUN_ID"],
            settings,
            run["SHOW_PLOTS"],
            log,
        )
        plt = plots.plt
        for src, (gdf_cells, years, maxima) in fits_for_plots.items():
            ids = [f"{src}_{int(k)}" for k in gdf_cells["cell_id"]]
            # flagged: the worst fits relative to the scale, and short records
            score = gdf_cells["fit_rmse"] / gdf_cells["gumbel_alpha"]
            flagged = [
                i
                for i, s_, n in zip(ids, score, gdf_cells["n_years"])
                if s_ > score.quantile(0.9) or n < 1.5 * dis["min_years"]
            ]
            for name in plots.select(f"fit_{src}", ids, flagged):
                j = ids.index(name)
                r = gdf_cells.iloc[j]
                x = np.sort(maxima[:, j][np.isfinite(maxima[:, j])])
                n = len(x)
                rp_obs = 1.0 / (1.0 - (np.arange(1, n + 1) - 0.44) / (n + 0.12))
                rp_max = max(1000.0, rp_obs.max() * 2)
                rp_line = np.logspace(np.log10(1.01), np.log10(rp_max), 200)
                fig, ax = plt.subplots(figsize=(7.5, 4.8), layout="constrained")
                ax.plot(
                    rp_line,
                    gis.Stats.gumbel_q(rp_line, r["gumbel_xi"], r["gumbel_alpha"]),
                    color=COLOURS["source"][src],
                    lw=LINE_WIDTH,
                    label="Gumbel fit (L-moments)",
                )
                ax.plot(
                    rp_obs,
                    x,
                    "o",
                    color=COLOURS["source"][src],
                    mfc="none",
                    ms=5,
                    label="annual maxima (Gringorten positions)",
                )
                Plots.rp_axis(ax, "x", rp_max)
                ax.set_xlabel("Return period (years)")
                ax.set_ylabel("Discharge (m³/s)")
                ax.legend(loc="upper left")
                ax.set_title(
                    f"fit · {SOURCE_NAMES[src]} cell {int(r['cell_id'])} · run {run['RUN_ID']}\n"
                    f"upstream area {r['upstream_area_km2']:.0f} km², n_years {int(r['n_years'])}, "
                    f"fit_rmse {r['fit_rmse']:.1f} m³/s"
                )
                plots.save(
                    fig,
                    f"fit_{src}",
                    name,
                    {
                        k: r[k]
                        for k in (
                            "upstream_area_km2",
                            "n_years",
                            "gumbel_xi",
                            "gumbel_alpha",
                            "fit_rmse",
                        )
                    },
                )
        log.info("Plots: %s", plots.finish())
    perf["plots_s"] = round(time.perf_counter() - t, 2)
    (out_dir / "discharge_prepare_performance.json").write_text(
        json.dumps(perf, indent=1)
    )
    log.info("Timings: %s", json.dumps(perf))
