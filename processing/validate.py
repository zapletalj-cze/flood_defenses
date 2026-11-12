# Author: Jakub Zapletal
# Date: 2026-11-03
# Compares the protected areas with reference protected-area polygons: area-based recall and precision, covered share of
# each reference polygon and, with reference return periods, a confusion matrix of RP classes. LAYER "A" or "B" (per
# river group; overlapping areas count once, outputs with the suffix _b).

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
LAYER = "A"  # "A": protected areas of all water sources together, "B": per river group

KEYS = [
    "general.output_dir",
    "general.log_level",
    "country.reference_protected_areas",
    "country.reference_rp_field",
    "country.validation_extent",
    "country.rp.classes",
]

if __name__ == "__main__":
    run = gis.Run.env_overrides(CONFIG_PATH=CONFIG_PATH, COUNTRY=COUNTRY, LAYER=LAYER)
    files = gis.Run.layer_files(run["LAYER"])
    sfx = files["suffix"]
    cfg = gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], KEYS)
    general, country = cfg["general"], cfg["country"]
    out_dir = Path(general["output_dir"]) / run["COUNTRY"]
    log = gis.Run.logger("validate", out_dir / "logs", general["log_level"])
    gis.Config.warn_placeholders(cfg)
    if not country["reference_protected_areas"]:
        raise SystemExit(
            "countries.<code>.reference_protected_areas is empty; nothing to validate."
        )
    result_path = out_dir / files["rp"]
    if not result_path.exists():
        result_path = out_dir / files["areas"]
    if not result_path.exists():
        raise SystemExit("No protected areas found; run protected_areas.py first.")
    gdf_pa = gpd.read_file(result_path)
    crs = gdf_pa.crs
    if crs is None or not crs.is_projected:
        raise SystemExit(
            "The protected areas must be in a projected CRS to compute areas."
        )
    field = country["reference_rp_field"]
    gdf_ref = gis.Vector.read_merged(
        country["reference_protected_areas"],
        crs.to_wkt(),
        keep_columns=[field] if field else None,
    )
    if country["validation_extent"]:
        extent = gis.Vector.read_merged(
            country["validation_extent"], crs.to_wkt()
        ).union_all()
        gdf_pa = gpd.clip(gdf_pa, extent)
        gdf_ref = gpd.clip(gdf_ref, extent)
    gdf_ref = gdf_ref[gdf_ref.area > 0].reset_index(drop=True)
    gdf_pa = gdf_pa[gdf_pa.area > 0].reset_index(drop=True)
    log.info(
        "Validating %s against %d reference polygons.", result_path.name, len(gdf_ref)
    )

    # --- Area-based recall and precision ---
    result_union = gdf_pa.union_all()
    ref_union = gdf_ref.union_all()
    intersection = result_union.intersection(ref_union).area
    recall = intersection / ref_union.area if ref_union.area > 0 else float("nan")
    precision = (
        intersection / result_union.area if result_union.area > 0 else float("nan")
    )
    log.info(
        "Reference area: %.0f m2, result area: %.0f m2, intersection: %.0f m2.",
        ref_union.area,
        result_union.area,
        intersection,
    )
    log.info("Recall: %.3f, precision: %.3f.", recall, precision)
    gdf_ref["area_m2"] = gdf_ref.area
    gdf_ref["covered_m2"] = gdf_ref.geometry.intersection(result_union).area
    gdf_ref["covered_share"] = gdf_ref["covered_m2"] / gdf_ref["area_m2"]
    log.info(
        "Reference polygons: %d, fully missed: %d, covered over 50 %%: %d.",
        len(gdf_ref),
        int((gdf_ref["covered_m2"] == 0).sum()),
        int((gdf_ref["covered_share"] > 0.5).sum()),
    )
    gis.Vector.write_gpkg(gdf_ref, out_dir / f"validation_coverage{sfx}.gpkg")
    report = {
        "result": result_path.name,
        "reference_polygons": len(gdf_ref),
        "recall": recall,
        "precision": precision,
        "reference_area_m2": ref_union.area,
        "result_area_m2": result_union.area,
    }

    # --- Return periods: each protected area against the reference polygon it overlaps most ---
    if field and "rp" in gdf_pa.columns:
        if field not in gdf_ref.columns:
            raise SystemExit(
                f"Reference field {field!r} not found; columns: {list(gdf_ref.columns)}"
            )
        classes = np.asarray(country["rp"]["classes"], float)
        left, right = gdf_ref.sindex.query(
            gdf_pa.geometry.values, predicate="intersects"
        )
        df_pairs = pd.DataFrame({"protected area": left, "reference": right})
        df_pairs["overlap"] = shapely.area(
            shapely.intersection(
                gdf_pa.geometry.values[left], gdf_ref.geometry.values[right]
            )
        )
        df_best = df_pairs.sort_values("overlap", ascending=False).drop_duplicates(
            "protected area"
        )
        df_best = df_best[df_best["overlap"] > 0]
        reference_rp = pd.to_numeric(gdf_ref[field], errors="coerce").to_numpy()[
            df_best["reference"]
        ]
        known = np.isfinite(reference_rp) & (reference_rp > 0)
        df_best, reference_rp = df_best[known], reference_rp[known]
        gdf_cmp = gdf_pa.iloc[df_best["protected area"].to_numpy()]
        # Reference values snapped to the nearest class (log scale) for the confusion matrix.
        snapped = classes[
            np.argmin(
                np.abs(np.log(reference_rp)[:, None] - np.log(classes)[None, :]), axis=1
            )
        ]
        df_cm = pd.crosstab(
            pd.Series(snapped, name="reference_rp"),
            pd.Series(gdf_cmp["rp"].to_numpy(), name="rp"),
            dropna=False,
        ).reindex(index=classes, columns=classes, fill_value=0)
        df_cm.to_csv(out_dir / f"validation_rp_confusion{sfx}.csv")
        within = (gdf_cmp["rp_low"].to_numpy() <= reference_rp) & (
            reference_rp <= gdf_cmp["rp_high"].to_numpy()
        )
        exact = snapped == gdf_cmp["rp"].to_numpy()
        report.update(
            {
                "rp_protected_areas_compared": int(len(df_best)),
                "rp_share_within_interval": (
                    float(within.mean()) if len(df_best) else float("nan")
                ),
                "rp_share_same_class": (
                    float(exact.mean()) if len(df_best) else float("nan")
                ),
            }
        )
        log.info(
            "RP: %d protected areas with a reference RP; reference RP within rp_low..rp_high for %.1f %%, same "
            "class for %.1f %%. Confusion matrix in validation_rp_confusion.csv.",
            len(df_best),
            100 * report["rp_share_within_interval"],
            100 * report["rp_share_same_class"],
        )
    elif field:
        log.warning(
            "The protected areas have no rp column (run rp_index.py); RP validation skipped."
        )
    (out_dir / f"validation{sfx}.json").write_text(json.dumps(report, indent=2))
