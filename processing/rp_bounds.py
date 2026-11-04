# Author: Jakub Zapletal
# Date: 2026-11-04
# Bounds the return period of each compartment with GFM tests and GloFAS discharge: a held test is a lower bound, a
# flooded test an upper bound (Gumbel fit of annual maxima, L-moments); combines them with the model interval.

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

import gis

# --- Run variables (no command-line parameters; DEFENSES_<NAME> environment variables override them) ---
CONFIG_PATH = "config/defenses.toml"
COUNTRY = "XX"

KEYS = [
    "general.output_dir",
    "general.log_level",
    "country.rp.classes",
    "sources.glofas_enabled",
    "sources.glofas_min_days_per_year",
    "evidence.glofas_search_radius_km",
    "evidence.upstream_ratio_min",
    "evidence.upstream_ratio_max",
    "evidence.flooded_window_days",
]
MODEL_COLUMNS = ["rp", "rp_low", "rp_high", "rp_source"]

if __name__ == "__main__":
    run = gis.Run.env_overrides(CONFIG_PATH=CONFIG_PATH, COUNTRY=COUNTRY)
    cfg = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], KEYS)
    general, evidence = cfg["general"], cfg["evidence"]
    out_dir = Path(general["output_dir"]) / run["COUNTRY"]
    log = gis.Run.logger("rp_bounds", out_dir / "logs", general["log_level"])
    gis.Config.warn_placeholders(cfg)
    if not cfg["sources"]["glofas_enabled"]:
        log.info("glofas_enabled is false; no observed bounds.")
        raise SystemExit(0)
    needed = [
        out_dir / "compartments_rp.gpkg",
        out_dir / "gfm_evidence.csv",
        out_dir / "glofas_discharge.npz",
    ]
    if not all(path.exists() for path in needed):
        raise SystemExit(
            f"Missing inputs {[str(p) for p in needed if not p.exists()]}; run rp_index.py, "
            "gfm_evidence.py and download_glofas.py first."
        )
    result = gpd.read_file(needed[0])
    # The model values of rp_index.py are kept as *_prior, so that reruns start from them again.
    for column in MODEL_COLUMNS:
        if f"{column}_prior" not in result.columns:
            result[f"{column}_prior"] = result[column]
        result[column] = result[f"{column}_prior"]
    classes = np.asarray(cfg["country"]["rp"]["classes"], float)

    # --- Matching compartments to GloFAS cells: within the radius, upstream area ratio in range, closest ratio ---
    glofas = np.load(needed[2])
    cells = gpd.GeoSeries(
        gpd.points_from_xy(glofas["lon"], glofas["lat"]), crs=4326
    ).to_crs(result.crs)
    tree = cKDTree(np.column_stack([cells.x, cells.y]))
    cell_area = glofas["upstream_area_km2"]
    match_cell = np.full(len(result), -1, np.int64)
    match_distance = np.full(len(result), np.nan)
    match_ratio = np.full(len(result), np.nan)
    for i, (x, y, area) in enumerate(
        zip(result["source_x"], result["source_y"], result["upstream_area_km2"])
    ):
        if not (np.isfinite(x) and np.isfinite(area) and area > 0):
            continue
        candidates = np.array(
            tree.query_ball_point(
                [x, y], r=evidence["glofas_search_radius_km"] * 1000.0
            ),
            int,
        )
        if not len(candidates):
            continue
        ratio = cell_area[candidates] / area
        ok = (ratio >= evidence["upstream_ratio_min"]) & (
            ratio <= evidence["upstream_ratio_max"]
        )
        if not ok.any():
            continue
        best = candidates[ok][np.argmin(np.abs(np.log(ratio[ok])))]
        match_cell[i] = best
        match_distance[i] = np.hypot(cells.x[best] - x, cells.y[best] - y) / 1000.0
        match_ratio[i] = cell_area[best] / area
    result["glofas_cell"] = match_cell
    result["glofas_distance_km"] = match_distance
    result["glofas_area_ratio"] = match_ratio
    matched = result["glofas_cell"] >= 0
    log.info(
        "Matched to a GloFAS cell: %d of %d compartments.",
        int(matched.sum()),
        len(result),
    )

    # --- Gumbel fit (L-moments) of the annual maxima of each matched cell ---
    dates = glofas["dates"].astype("datetime64[D]")
    discharge = glofas["discharge"]
    years = dates.astype("datetime64[Y]").astype(int) + 1970
    fits = {}
    for cell in np.unique(result.loc[matched, "glofas_cell"]).astype(int):
        series = pd.Series(discharge[:, cell]).groupby(years)
        complete = series.count() >= cfg["sources"]["glofas_min_days_per_year"]
        maxima = series.max()[complete]
        fits[cell] = gis.Stats.gumbel_lmoments(maxima.to_numpy())
        if len(maxima) < 10:
            log.warning(
                "GloFAS cell %d: only %d complete years for the Gumbel fit.",
                cell,
                len(maxima),
            )

    # --- Event return periods of the GFM tests ---
    tests = pd.read_csv(needed[1])
    tests = tests[tests["outcome"].isin(["held", "flooded"])].copy()
    tests = tests.merge(
        result[["compartment_id", "glofas_cell"]], on="compartment_id", how="inner"
    )
    tests = tests[tests["glofas_cell"] >= 0]
    day_index = {d: i for i, d in enumerate(dates)}
    window = evidence["flooded_window_days"]
    event_rp = []
    acquisition_days = pd.to_datetime(
        tests["datetime"], utc=True, format="ISO8601"
    ).dt.date
    for outcome, day, cell in zip(
        tests["outcome"], acquisition_days, tests["glofas_cell"].astype(int)
    ):
        k = day_index.get(np.datetime64(day, "D"))
        xi, alpha = fits.get(cell, (np.nan, np.nan))
        if k is None or not np.isfinite(alpha):
            event_rp.append(np.nan)
            continue
        if outcome == "held":
            q = discharge[
                k, cell
            ]  # discharge on the acquisition day: a valid lower bound
        else:
            q = np.nanmax(
                discharge[max(k - window + 1, 0) : k + 1, cell]
            )  # peak before and on the day: upper bound
        event_rp.append(float(gis.Stats.gumbel_rp(q, xi, alpha)))
    tests["event_rp"] = event_rp
    tests = tests[np.isfinite(tests["event_rp"])]
    log.info(
        "GFM tests with an event return period: %d (held %d, flooded %d).",
        len(tests),
        int((tests["outcome"] == "held").sum()),
        int((tests["outcome"] == "flooded").sum()),
    )

    held = tests[tests["outcome"] == "held"].groupby("compartment_id")["event_rp"]
    flooded = tests[tests["outcome"] == "flooded"].groupby("compartment_id")["event_rp"]
    result["n_held"] = result["compartment_id"].map(held.count()).fillna(0).astype(int)
    result["n_flooded"] = (
        result["compartment_id"].map(flooded.count()).fillna(0).astype(int)
    )
    result["obs_rp_low"] = result["compartment_id"].map(held.max())
    result["obs_rp_high"] = result["compartment_id"].map(flooded.min())
    result["obs_conflict"] = (result["obs_rp_low"] > result["obs_rp_high"]).fillna(
        False
    )
    result["model_obs_conflict"] = False

    # --- Combination with the model interval ---
    cases = {
        "combined": 0,
        "model_obs_conflict": 0,
        "obs_conflict": 0,
        "override_kept": 0,
    }
    log_classes = np.log(classes)
    for i in result.index[(result["n_held"] + result["n_flooded"]) > 0]:
        row = result.loc[i]
        if row["rp_source_prior"] == "override":
            cases["override_kept"] += 1
            continue
        if row["obs_conflict"]:
            cases["obs_conflict"] += 1
            continue
        obs_low = row["obs_rp_low"] if np.isfinite(row["obs_rp_low"]) else -np.inf
        obs_high = row["obs_rp_high"] if np.isfinite(row["obs_rp_high"]) else np.inf
        low = max(row["rp_low_prior"], obs_low)
        high = min(row["rp_high_prior"], obs_high)
        if low > high:
            # The observations contradict the model interval: only the observed interval is used, a side
            # without observations stays open.
            result.loc[i, "model_obs_conflict"] = True
            cases["model_obs_conflict"] += 1
            low, high = obs_low, obs_high
        else:
            cases["combined"] += 1
        value = float(np.clip(row["rp_prior"], low, high))
        # Snap to the class inside [low, high] closest to the value (log scale), else to the closest class.
        inside = (classes >= low) & (classes <= high)
        pool = np.flatnonzero(inside) if inside.any() else np.arange(len(classes))
        k = pool[np.argmin(np.abs(log_classes[pool] - np.log(value)))]
        snapped = classes[k]
        # An open side of the interval becomes the neighbouring class of the result.
        if not np.isfinite(low):
            low = classes[max(k - 1, 0)]
        if not np.isfinite(high):
            high = classes[min(k + 1, len(classes) - 1)]
        result.loc[i, "rp"] = snapped
        result.loc[i, "rp_low"] = min(low, snapped)
        result.loc[i, "rp_high"] = max(high, snapped)
        result.loc[i, "rp_source"] = f"{row['rp_source_prior']}+obs"
    for column in ["rp", "rp_low", "rp_high"]:
        result[column] = result[column].astype(float)
    gis.Vector.write_gpkg(result, out_dir / "compartments_rp.gpkg")
    log.info(
        "Observed bounds: %s. compartments_rp.gpkg updated (%d compartments).",
        ", ".join(f"{k} {v}" for k, v in cases.items()),
        len(result),
    )
