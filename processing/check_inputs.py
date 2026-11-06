# Author: Jakub Zapletal
# Date: 2026-11-05
# Checks the configuration and the inputs before a run: the keys every script needs (core and optional steps),
# the input files (open, one layer, CRS, the upstream area column), the damage tables and the output folders.
# Prints one overview and writes nothing except a test file in the output and work folders.

import importlib
import os
import sys
from pathlib import Path

import pyogrio
from osgeo import gdal, osr

import gis

gdal.UseExceptions()

# --- Run variables (no command-line parameters; DEFENSES_<NAME> environment variables override them) ---
CONFIG_PATH = "config/defenses.toml"
COUNTRY = "XX"  # "XX": the country of [general] country

CORE = [
    "processing.prepare_dtm",
    "processing.build_units",
    "processing.compartments",
    "download.download_ghsl",
    "processing.rp_index",
]
OPTIONAL = [
    "download.download_osm",
    "processing.osm_hydraulics",
    "download.download_gfm",
    "processing.gfm_evidence",
    "download.download_discharge",
    "processing.discharge_prepare",
    "processing.rp_bounds",
    "processing.validate",
]

if __name__ == "__main__":
    run = gis.Run.env_overrides(CONFIG_PATH=CONFIG_PATH, COUNTRY=COUNTRY)
    gis.log.setLevel("ERROR")
    list_row = []  # (group, item, ok, message)

    def add(group, item, ok, message=""):
        list_row.append((group, item, ok, message))

    # --- Configuration: the keys of every script ---
    for group, modules in (("core", CORE), ("optional", OPTIONAL)):
        for name in modules:
            keys = importlib.import_module(name).KEYS
            try:
                gis.Config.load(run["CONFIG_PATH"], run["COUNTRY"], keys)
                add(f"configuration ({group})", name.split(".")[1], True)
            except SystemExit as exc:
                lines = str(exc).splitlines()
                add(
                    f"configuration ({group})",
                    name.split(".")[1],
                    False,
                    "; ".join(l.strip(" -") for l in lines[1:]) or lines[0],
                )
    raw = gis.Config.read_merged(run["CONFIG_PATH"])
    if run["COUNTRY"] not in raw.get("countries", {}):
        sys.exit(
            f"Country {run['COUNTRY']!r} is not in {run['CONFIG_PATH']} (set [general] country)."
        )
    country = gis.Config.merge(
        raw.get("country_defaults", {}), raw["countries"][run["COUNTRY"]], "country"
    )

    # --- Inputs ---
    def crs_name(wkt):
        srs = osr.SpatialReference(wkt=wkt) if wkt else None
        return (
            (srs.GetAuthorityName(None) or "")
            + ":"
            + (srs.GetAuthorityCode(None) or "?")
            if srs
            else "none"
        )

    dtm_crs = None
    for path in gis.Config.as_list(country.get("dtm")) or [""]:
        if not path:
            add("inputs", "dtm", False, "not set")
            continue
        try:
            ds = gdal.Open(path)
            gt = ds.GetGeoTransform()
            dtm_crs = ds.GetProjection()
            add(
                "inputs",
                f"dtm {Path(path).name}",
                True,
                f"{ds.RasterXSize} x {ds.RasterYSize}, {abs(gt[1]):g} m, {crs_name(dtm_crs)}",
            )
        except Exception as exc:  # noqa: BLE001 - reported
            add("inputs", f"dtm {Path(path).name}", False, str(exc)[:200])
    field = country.get("river_upstream_area_field", "")
    for kind in ("levees", "rivers"):
        for path in gis.Config.as_list(country.get(kind)) or [""]:
            if not path:
                add("inputs", kind, False, "not set")
                continue
            try:
                layers = pyogrio.list_layers(path)
                if len(layers) != 1:
                    add(
                        "inputs",
                        f"{kind} {Path(path).name}",
                        False,
                        f"{len(layers)} layers; save one layer per file",
                    )
                    continue
                info = pyogrio.read_info(path)
                problems = []
                if kind == "rivers" and field and field not in list(info["fields"]):
                    problems.append(
                        f"no column {field!r} for river_upstream_area_field; columns: {list(info['fields'])}"
                    )
                if not info["crs"]:
                    problems.append("no CRS")
                elif dtm_crs:
                    srs = osr.SpatialReference()
                    srs.SetFromUserInput(info["crs"])
                    if not srs.IsSame(osr.SpatialReference(wkt=dtm_crs)):
                        # not an error: the scripts reproject vectors to the DTM CRS
                        info["crs"] = f"{info['crs']} (reprojected to the DTM CRS)"
                add(
                    "inputs",
                    f"{kind} {Path(path).name}",
                    not problems,
                    "; ".join(problems)
                    or f"{info['features']} features, {info['crs']}",
                )
            except Exception as exc:  # noqa: BLE001 - reported
                add("inputs", f"{kind} {Path(path).name}", False, str(exc)[:200])
    if not field and not country.get("upstream_area_raster"):
        add(
            "inputs",
            "upstream area",
            False,
            "set river_upstream_area_field (column name) or upstream_area_raster",
        )
    pbf, region = country.get("osm_pbf", ""), country.get("geofabrik_region", "")
    if pbf:
        add(
            "inputs (optional)",
            "osm_pbf",
            Path(pbf).exists(),
            (
                f"{Path(pbf).stat().st_size / 1e9:.1f} GB"
                if Path(pbf).exists()
                else "file not found"
            ),
        )
    else:
        add(
            "inputs (optional)",
            "osm_pbf / geofabrik_region",
            bool(region),
            region or "not set: no OSM water sources and pumping stations (optional)",
        )
    for key in ("max_damage_csv", "depth_damage_csv"):
        path = country.get("rp", {}).get(key, "")
        add(
            "inputs",
            f"rp.{key}",
            bool(path) and Path(path).exists(),
            path or "not set (needed by rp_index)",
        )

    # --- Output folders ---
    for key in ("work_dir", "output_dir"):
        folder = raw.get("general", {}).get(key, "")
        try:
            target = Path(folder) / run["COUNTRY"]
            target.mkdir(parents=True, exist_ok=True)
            (target / ".write_test").write_text("ok")
            (target / ".write_test").unlink()
            add("folders", key, bool(folder), str(target))
        except Exception as exc:  # noqa: BLE001 - reported
            add("folders", key, False, f"{folder!r}: {exc}")

    # --- Overview ---
    print(f"\nCheck of {run['CONFIG_PATH']} for country {run['COUNTRY']}\n")
    group = None
    for g, item, ok, message in list_row:
        if g != group:
            print(f"{g}:")
            group = g
        print(
            f"  {'OK     ' if ok else 'PROBLEM'}  {item}{'  ' + message if message else ''}"
        )
    core_bad = [r for r in list_row if not r[2] and not r[0].endswith("(optional)")]
    print(
        "\nCore steps can run."
        if not core_bad
        else f"\n{len(core_bad)} problem(s) block the core steps (see PROBLEM above)."
    )
    sys.exit(0 if not core_bad else 1)
