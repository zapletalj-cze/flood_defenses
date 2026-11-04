# Author: Jakub Zapletal
# Date: 2026-11-03
# Compares the compartments with reference protected-area polygons: area-based recall and precision, covered share of
# each reference polygon and, with reference return periods, a confusion matrix of RP classes.

import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely

import gis

# --- Run variables (no command-line parameters; DEFENSES_<NAME> environment variables override them) ---
CONFIG_PATH = "config/defenses.toml"
COUNTRY = "XX"

KEYS = [
    "general.output_dir",
    "general.log_level",
    "country.reference_protected_areas",
    "country.reference_rp_field",
    "country.validation_extent",
    "country.rp.classes",
]

if __name__ == "__main__":
    run = gis.Run.env_overrides(CONFIG_PATH=CONFIG_PATH, COUNTRY=COUNTRY)
    cfg = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], KEYS)
    general, country = cfg["general"], cfg["country"]
    out_dir = Path(general["output_dir"]) / run["COUNTRY"]
    log = gis.Run.logger("validate", out_dir / "logs", general["log_level"])
    gis.Config.warn_placeholders(cfg)
    if not country["reference_protected_areas"]:
        raise SystemExit(
            "countries.<code>.reference_protected_areas is empty; nothing to validate."
        )
    result_path = out_dir / "compartments_rp.gpkg"
    if not result_path.exists():
        result_path = out_dir / "compartments.gpkg"
    if not result_path.exists():
        raise SystemExit("No compartments found; run compartments.py first.")
    result = gpd.read_file(result_path)
    crs = result.crs
    if crs is None or not crs.is_projected:
        raise SystemExit(
            "The compartments must be in a projected CRS to compute areas."
        )
    field = country["reference_rp_field"]
    reference = gis.Vector.read_merged(
        country["reference_protected_areas"],
        crs.to_wkt(),
        keep_columns=[field] if field else None,
    )
    if country["validation_extent"]:
        extent = gis.Vector.read_merged(
            country["validation_extent"], crs.to_wkt()
        ).union_all()
        result = gpd.clip(result, extent)
        reference = gpd.clip(reference, extent)
    reference = reference[reference.area > 0].reset_index(drop=True)
    result = result[result.area > 0].reset_index(drop=True)
    log.info(
        "Validating %s against %d reference polygons.", result_path.name, len(reference)
    )

    # --- Area-based recall and precision ---
    result_union = result.union_all()
    reference_union = reference.union_all()
    intersection = result_union.intersection(reference_union).area
    recall = (
        intersection / reference_union.area
        if reference_union.area > 0
        else float("nan")
    )
    precision = (
        intersection / result_union.area if result_union.area > 0 else float("nan")
    )
    log.info(
        "Reference area: %.0f m2, result area: %.0f m2, intersection: %.0f m2.",
        reference_union.area,
        result_union.area,
        intersection,
    )
    log.info("Recall: %.3f, precision: %.3f.", recall, precision)
    reference["area_m2"] = reference.area
    reference["covered_m2"] = reference.geometry.intersection(result_union).area
    reference["covered_share"] = reference["covered_m2"] / reference["area_m2"]
    log.info(
        "Reference polygons: %d, fully missed: %d, covered over 50 %%: %d.",
        len(reference),
        int((reference["covered_m2"] == 0).sum()),
        int((reference["covered_share"] > 0.5).sum()),
    )
    gis.Vector.write_gpkg(reference, out_dir / "validation_coverage.gpkg")
    report = {
        "result": result_path.name,
        "reference_polygons": len(reference),
        "recall": recall,
        "precision": precision,
        "reference_area_m2": reference_union.area,
        "result_area_m2": result_union.area,
    }

    # --- Return periods: each compartment against the reference polygon it overlaps most ---
    if field and "rp" in result.columns:
        if field not in reference.columns:
            raise SystemExit(
                f"Reference field {field!r} not found; columns: {list(reference.columns)}"
            )
        classes = np.asarray(country["rp"]["classes"], float)
        left, right = reference.sindex.query(
            result.geometry.values, predicate="intersects"
        )
        pairs = pd.DataFrame({"compartment": left, "reference": right})
        pairs["overlap"] = shapely.area(
            shapely.intersection(
                result.geometry.values[left], reference.geometry.values[right]
            )
        )
        best = pairs.sort_values("overlap", ascending=False).drop_duplicates(
            "compartment"
        )
        best = best[best["overlap"] > 0]
        reference_rp = pd.to_numeric(reference[field], errors="coerce").to_numpy()[
            best["reference"]
        ]
        known = np.isfinite(reference_rp) & (reference_rp > 0)
        best, reference_rp = best[known], reference_rp[known]
        compared = result.iloc[best["compartment"].to_numpy()]
        # Reference values snapped to the nearest class (log scale) for the confusion matrix.
        snapped = classes[
            np.argmin(
                np.abs(np.log(reference_rp)[:, None] - np.log(classes)[None, :]), axis=1
            )
        ]
        matrix = pd.crosstab(
            pd.Series(snapped, name="reference_rp"),
            pd.Series(compared["rp"].to_numpy(), name="rp"),
            dropna=False,
        ).reindex(index=classes, columns=classes, fill_value=0)
        matrix.to_csv(out_dir / "validation_rp_confusion.csv")
        within = (compared["rp_low"].to_numpy() <= reference_rp) & (
            reference_rp <= compared["rp_high"].to_numpy()
        )
        exact = snapped == compared["rp"].to_numpy()
        report.update(
            {
                "rp_compartments_compared": int(len(best)),
                "rp_share_within_interval": (
                    float(within.mean()) if len(best) else float("nan")
                ),
                "rp_share_same_class": (
                    float(exact.mean()) if len(best) else float("nan")
                ),
            }
        )
        log.info(
            "RP: %d compartments with a reference RP; reference RP within rp_low..rp_high for %.1f %%, same "
            "class for %.1f %%. Confusion matrix in validation_rp_confusion.csv.",
            len(best),
            100 * report["rp_share_within_interval"],
            100 * report["rp_share_same_class"],
        )
    elif field:
        log.warning(
            "The compartments have no rp column (run rp_index.py); RP validation skipped."
        )
    (out_dir / "validation.json").write_text(json.dumps(report, indent=2))
