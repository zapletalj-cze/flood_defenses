# Author: Jakub Zapletal
# Date: 2026-11-06
# Compact report of a country's inputs, outputs and run logs, written into the repository (reports/<country>/) so
# that it can be committed and reviewed elsewhere: report.md (a few tens of KB at most) and overview.png (one small
# map). Reads samples and summaries only: DTM windows, vector statistics within the DTM footprint, output tables,
# the newest log of every script. Never copies data.

import json
import re
import sys
import time
from collections import Counter
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pyogrio
import shapely
from osgeo import gdal, osr

import gis
from plots import CATEGORICAL, INK

gdal.UseExceptions()

# --- Run variables (no command-line parameters; DEFENSES_<NAME> environment variables override them) ---
CONFIG_PATH = "config/defenses.toml"
COUNTRY = "XX"  # "XX": the country of [general] country
REPORT_DIR = (
    "reports"  # folder in the repository; the report goes to <REPORT_DIR>/<country>/
)
PLOT = True  # overview.png with the DTM footprint, levees and protected areas

MAX_REPORT_KB = 60  # report.md is cut at this size
DTM_WINDOWS = 24  # random 256 x 256 windows read for the DTM statistics
TOP_ROWS = 10  # largest features listed per output layer
TOP_VALUES = 8  # values listed per text column
LOG_LINES = 12  # distinct warning / error messages listed per log
PLOT_PX = 1400  # long side of overview.png in pixels
PERCENTILES = [0, 5, 25, 50, 75, 95, 100]

# Output layers summarised when present, with the columns shown in the list of the largest features
OUTPUTS = {
    "protected_areas.gpkg": [
        "protected_area_id",
        "area_m2",
        "protection_level_m",
        "margin_median_m",
        "truncated",
        "unit_id",
    ],
    "protected_areas_rp.gpkg": [
        "protected_area_id",
        "area_m2",
        "rp",
        "rp_low",
        "rp_high",
        "rp_source",
    ],
    "protected_areas_by_group.gpkg": [
        "protected_area_id",
        "group_id",
        "area_m2",
        "protection_level_m",
        "unit_id",
    ],
    "protected_areas_by_group_rp.gpkg": [
        "protected_area_id",
        "group_id",
        "area_m2",
        "rp",
        "rp_low",
        "rp_high",
    ],
    "leaks.gpkg": [
        "leak_id",
        "cause",
        "x",
        "y",
        "levee_length_m",
        "lost_margin_median_m",
        "line_ids",
    ],
    "transport_embankments.gpkg": [
        "segment_id",
        "class",
        "length_m",
        "height_min_m",
        "openings",
    ],
    "openings.gpkg": ["opening_id", "type", "x", "y", "width_m", "confidence"],
    "events.gpkg": [
        "event_id",
        "group_id",
        "source",
        "peak_date",
        "peak_rp",
        "n_peaks",
    ],
    "barriers.gpkg": [],
    "levees_crest.gpkg": [],
    "enclosed_streams.gpkg": [],
    "rivers_channel.gpkg": [],
    "river_groups.gpkg": [],
    "units_core.gpkg": [],
    "gfm_evidence.gpkg": [],
    "implied/base/implied_protected_areas.gpkg": [],
    "defenses_clip.gpkg": [],
}
TABLES = [
    "units.csv",
    "gfm_evidence.csv",
    "event_rp.csv",
    "validation_rp_confusion.csv",
]
JSONS = [
    "performance.json",
    "dtm_estimates.json",
    "validation.json",
    "rp_bounds_performance.json",
]


def fmt(v):
    if isinstance(v, (float, np.floating)):
        if not np.isfinite(v):
            return str(v)
        return f"{v:.4g}" if abs(v) < 1e5 else f"{v:,.0f}"
    return str(v)


def pct_row(values):
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if not len(v):
        return "no finite values"
    return ", ".join(
        f"p{p} {fmt(x)}" for p, x in zip(PERCENTILES, np.percentile(v, PERCENTILES))
    )


def profile(gdf, list_line, show_cols=()):
    """Columns, geometry types, sizes, numeric percentiles, flags, text values and the largest features"""
    geom = gdf.geometry
    types = Counter(geom.geom_type.fillna("None"))
    list_line.append(
        f"- features: {len(gdf):,}; geometry types: {dict(types)}; empty or invalid: "
        f"{int(geom.is_empty.sum())} / {int((~geom.is_valid).sum())}"
    )
    if len(gdf):
        x0, y0, x1, y1 = gdf.total_bounds
        list_line.append(f"- bounds: {fmt(x0)}, {fmt(y0)}, {fmt(x1)}, {fmt(y1)}")
        dims = shapely.get_dimensions(geom.values)
        if (dims == 1).any():
            length = shapely.length(geom.values[dims == 1]) / 1000
            list_line.append(
                f"- line length km: total {fmt(length.sum())}; per feature {pct_row(length)}"
            )
        if (dims == 2).any():
            area = shapely.area(geom.values[dims == 2]) / 1e6
            list_line.append(
                f"- polygon area km2: total {fmt(area.sum())}; per feature {pct_row(area)}"
            )
    for col in [c for c in gdf.columns if c != gdf.geometry.name][:40]:
        s = gdf[col]
        if pd.api.types.is_bool_dtype(s):
            list_line.append(f"- `{col}` (bool): true {int(s.sum()):,}")
        elif pd.api.types.is_numeric_dtype(s):
            list_line.append(
                f"- `{col}` ({s.dtype}): null {int(s.isna().sum()):,}; {pct_row(s)}"
            )
        else:
            vc = s.astype(str).value_counts()
            top = ", ".join(f"{k[:30]!r} {v}" for k, v in vc.head(TOP_VALUES).items())
            list_line.append(f"- `{col}` (text): {len(vc):,} distinct; {top}")
    cols = [c for c in show_cols if c in gdf.columns]
    if cols and len(gdf):
        key = "area_m2" if "area_m2" in gdf.columns else None
        order = gdf.sort_values(key, ascending=False) if key else gdf
        list_line.append(f"- largest {min(TOP_ROWS, len(gdf))}:")
        list_line.append("")
        list_line.append("| " + " | ".join(cols) + " |")
        list_line.append("|" + "---|" * len(cols))
        for _, r in order.head(TOP_ROWS).iterrows():
            list_line.append("| " + " | ".join(fmt(r[c]) for c in cols) + " |")
        list_line.append("")


if __name__ == "__main__":
    run = gis.Run.env_overrides(
        CONFIG_PATH=CONFIG_PATH, COUNTRY=COUNTRY, REPORT_DIR=REPORT_DIR, PLOT=PLOT
    )
    raw = gis.Config.read_merged(run["CONFIG_PATH"])
    code = run["COUNTRY"]
    if code not in raw.get("countries", {}):
        sys.exit(
            f"Country {code!r} is not in {run['CONFIG_PATH']} (set [general] country)."
        )
    country = gis.Config.merge(
        raw.get("country_defaults", {}), raw["countries"][code], "country"
    )
    general = raw.get("general", {})
    out_dir = Path(general.get("output_dir", "")) / code
    rep_dir = Path(run["REPORT_DIR"]) / code
    rep_dir.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    rng = np.random.default_rng(1)
    L = [
        f"# Data report {code}",
        "",
        f"Created {time.strftime('%Y-%m-%d %H:%M')}; configuration {run['CONFIG_PATH']}; "
        f"output folder {out_dir}",
        "",
    ]

    # --- Configuration in use (the user's keys and the country table) ---
    L += ["## Configuration", ""]
    for key in ("n_workers", "max_ram_gb_per_worker"):
        L.append(f"- general.{key} = {general.get(key)!r}")
    for key in (
        "dtm",
        "levees",
        "rivers",
        "river_upstream_area_field",
        "upstream_area_raster",
        "river_guides",
        "working_resolution_m",
        "corridor_buffer_m",
        "unit_core_size_m",
        "unit_halo_m",
        "bytes_per_cell",
    ):
        L.append(f"- country.{key} = {country.get(key)!r}")
    pa_cfg = country.get("protected_areas", {})
    L.append(
        "- country.protected_areas: "
        + ", ".join(f"{k} = {v!r}" for k, v in pa_cfg.items())
    )
    L.append("")

    # --- DTM: storage and statistics of random windows ---
    L += ["## Inputs", "", "### DTM", ""]
    dtm_wkt, footprint = None, None
    for path in gis.Config.as_list(country.get("dtm")):
        try:
            ds = gdal.Open(path)
            band = ds.GetRasterBand(1)
            gt, dtm_wkt = ds.GetGeoTransform(), ds.GetProjection()
            srs = osr.SpatialReference(wkt=dtm_wkt)
            bx, by = band.GetBlockSize()
            L.append(
                f"- {Path(path).name}: {ds.RasterXSize:,} x {ds.RasterYSize:,} px, {abs(gt[1]):g} m, "
                f"{srs.GetAuthorityName(None)}:{srs.GetAuthorityCode(None)}, {gdal.GetDataTypeName(band.DataType)}, "
                f"nodata {band.GetNoDataValue()}, block {bx} x {by}, overviews {band.GetOverviewCount()}, "
                f"compression {ds.GetMetadataItem('COMPRESSION', 'IMAGE_STRUCTURE')}"
            )
            nodata, list_v, n_nod = band.GetNoDataValue(), [], 0
            for _ in range(DTM_WINDOWS):
                w, h = min(256, ds.RasterXSize), min(256, ds.RasterYSize)
                x = int(rng.integers(0, ds.RasterXSize - w + 1))
                y = int(rng.integers(0, ds.RasterYSize - h + 1))
                a = band.ReadAsArray(x, y, w, h).astype(float)
                bad = ~np.isfinite(a) | ((a == nodata) if nodata is not None else False)
                n_nod += int(bad.sum())
                list_v.append(a[~bad])
            v = np.concatenate(list_v)
            L.append(
                f"  - {DTM_WINDOWS} random windows: nodata share {n_nod / (n_nod + len(v)):.1%}; "
                f"elevation {pct_row(v)}"
            )
        except Exception as exc:  # noqa: BLE001 - reported
            L.append(f"- {path}: cannot read ({str(exc)[:150]})")
    fp_path = out_dir / "dtm_footprint.gpkg"
    area_name = "DTM footprint"
    if fp_path.exists():
        footprint = gpd.read_file(fp_path)
        L.append(f"- dtm_footprint.gpkg: area {fmt(footprint.area.sum() / 1e6)} km2")
    elif (out_dir / "corridor.gpkg").exists():
        # without the footprint the corridor limits the inputs and the map
        footprint = gpd.read_file(out_dir / "corridor.gpkg")
        area_name = "corridor"
        L.append(
            f"- no dtm_footprint.gpkg: inputs and map limited to corridor.gpkg ({fmt(footprint.area.sum() / 1e6)} km2)"
        )
    L.append("")

    # --- Levees and rivers within the DTM footprint ---
    gdf_lv = None
    for kind in ("levees", "rivers"):
        L += [f"### {kind}", ""]
        for path in gis.Config.as_list(country.get(kind)):
            try:
                info = pyogrio.read_info(path)
                L.append(
                    f"- {Path(path).name}: {info['features']:,} features in the file, {info['geometry_type']}, "
                    f"{info['crs']}, columns {list(info['fields'])[:30]}"
                )
                gdf = (
                    gpd.read_file(path, bbox=footprint)
                    if footprint is not None
                    else gpd.read_file(path)
                )
                if dtm_wkt and len(gdf):
                    gdf = gdf.to_crs(dtm_wkt)
                if footprint is not None and len(gdf):
                    gdf = gdf[gdf.intersects(footprint.to_crs(gdf.crs).union_all())]
                L.append(
                    f"- inside the {area_name}{'' if footprint is not None else ' (no footprint yet: whole file)'}:"
                )
                profile(gdf, L)
                if kind == "levees":
                    gdf_lv = gdf if gdf_lv is None else pd.concat([gdf_lv, gdf])
            except Exception as exc:  # noqa: BLE001 - reported
                L.append(f"- {path}: cannot read ({str(exc)[:150]})")
        L.append("")

    # --- Outputs ---
    L += ["## Outputs", ""]
    if out_dir.exists():
        files = sorted(
            p
            for p in out_dir.rglob("*")
            if p.is_file() and "logs" not in p.parts and "plots" not in p.parts
        )
        L.append(
            f"- {len(files)} files, {sum(p.stat().st_size for p in files) / 1e9:.2f} GB (without logs and plots); "
            "top level: "
            + ", ".join(
                f"{p.name} {p.stat().st_size / 1e6:.1f} MB"
                for p in files
                if p.parent == out_dir
            )[:1500]
        )
    L.append("")
    layers = {}
    for name, cols in OUTPUTS.items():
        path = out_dir / name
        if not path.exists():
            continue
        L += [f"### {name}", ""]
        try:
            gdf = gpd.read_file(path)
            layers[name] = gdf
            profile(gdf, L, cols)
            if name == "protected_areas_by_group.gpkg" and len(gdf):
                # areas of different river groups overlap: the union counts every place once
                union_b = shapely.union_all(gdf.geometry.values)
                gdf_a = layers.get("protected_areas.gpkg")
                union_a = (
                    shapely.union_all(gdf_a.geometry.values)
                    if gdf_a is not None and len(gdf_a)
                    else shapely.Polygon()
                )
                L.append(
                    f"- union of all groups (overlaps counted once): {fmt(union_b.area / 1e6)} km2; file A "
                    f"{fmt(union_a.area / 1e6)} km2; in B but not in A {fmt(union_b.difference(union_a).area / 1e6)} "
                    f"km2; in A but not in B {fmt(union_a.difference(union_b).area / 1e6)} km2"
                )
            if name == "leaks.gpkg" and len(gdf):
                by = gdf.groupby("cause")["levee_length_m"].agg(["size", "sum"])
                L.append(
                    "- levee length per cause: "
                    + ", ".join(
                        f"{c} {int(r['size']):,} places {fmt(r['sum'] / 1000)} km"
                        for c, r in by.sort_values("sum", ascending=False).iterrows()
                    )
                )
        except Exception as exc:  # noqa: BLE001 - reported
            L.append(f"- cannot read ({str(exc)[:150]})")
        L.append("")
    for name in TABLES:
        path = out_dir / name
        if path.exists():
            df = pd.read_csv(path)
            L += [f"### {name}", "", f"- rows {len(df):,}"]
            for col in df.columns[:25]:
                s = df[col]
                if pd.api.types.is_numeric_dtype(s):
                    L.append(f"- `{col}`: {pct_row(s)}")
                else:
                    L.append(
                        f"- `{col}`: "
                        + ", ".join(
                            f"{k[:25]!r} {v}"
                            for k, v in s.astype(str)
                            .value_counts()
                            .head(TOP_VALUES)
                            .items()
                        )
                    )
            L.append("")
    for name in JSONS:
        path = out_dir / name
        if path.exists():
            text = json.dumps(json.loads(path.read_text()), separators=(", ", ": "))
            L += [f"### {name}", "", "```", text[:3000], "```", ""]

    # --- Logs: the newest log of every script, its end and its distinct warnings and errors ---
    L += ["## Logs (newest per script)", ""]
    newest = {}
    for path in (
        sorted((out_dir / "logs").glob("*.log")) if (out_dir / "logs").exists() else []
    ):
        newest[re.sub(r"_\d{8}_\d{6}$", "", path.stem)] = path
    for script, path in sorted(newest.items()):
        lines = path.read_text(errors="replace").splitlines()
        msgs = Counter()
        for ln in lines:
            m = re.search(r" (WARNING|ERROR|CRITICAL) [\w.]+: (.*)", ln)
            if m:
                msgs[
                    f"{m.group(1)} " + re.sub(r"\d+(\.\d+)?", "#", m.group(2))[:220]
                ] += 1
        L.append(f"### {script} ({path.name}, {len(lines):,} lines)")
        L.append("")
        for msg, n in msgs.most_common(LOG_LINES):
            L.append(f"- {n} x {msg}")
        L.append("- last lines:")
        L += ["```"] + [ln[:300] for ln in lines[-6:]] + ["```", ""]

    # --- Overview map ---
    if run["PLOT"] and (footprint is not None or "protected_areas.gpkg" in layers):
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(
            figsize=(PLOT_PX / 100, PLOT_PX / 100 * 0.8), dpi=100, layout="constrained"
        )
        if footprint is not None:
            footprint.boundary.plot(ax=ax, color=INK["muted"], lw=0.6)
        gdf_b_map = layers.get("protected_areas_by_group.gpkg")
        if gdf_b_map is not None and len(gdf_b_map):
            gdf_b_map.plot(ax=ax, color=CATEGORICAL[2], alpha=0.35, lw=0)
        gdf_pa_map = layers.get(
            "protected_areas_rp.gpkg", layers.get("protected_areas.gpkg")
        )
        if gdf_pa_map is not None and len(gdf_pa_map):
            # fill plus a point per protected area, so that small ones stay visible at country scale
            gdf_pa_map.plot(
                ax=ax, color=CATEGORICAL[0], alpha=0.6, edgecolor=CATEGORICAL[0], lw=0.6
            )
            pts = gdf_pa_map.representative_point()
            ax.scatter(pts.x, pts.y, s=6, color=CATEGORICAL[0], lw=0)
        if gdf_lv is not None and len(gdf_lv):
            gdf_lv.simplify(50).plot(ax=ax, color=CATEGORICAL[1], lw=0.8)
        ax.set_title(
            f"{code}: {area_name} (grey), levees in it (orange), protected areas A (blue), per river group B (green)",
            loc="left",
            fontsize=9,
        )
        ax.set_aspect("equal")
        if footprint is not None:
            # the map shows the processed area, not the whole extent of the input files
            x0, y0, x1, y1 = footprint.total_bounds
            pad = 0.03 * max(x1 - x0, y1 - y0)
            ax.set_xlim(x0 - pad, x1 + pad)
            ax.set_ylim(y0 - pad, y1 + pad)
        ax.ticklabel_format(style="plain", useOffset=False)
        ax.tick_params(labelsize=7)
        fig.savefig(rep_dir / "overview.png", dpi=100)
        plt.close(fig)
        L.append(
            f"Map: overview.png ({(rep_dir / 'overview.png').stat().st_size / 1e3:.0f} KB)"
        )

    # --- Write, cut to the size limit ---
    text = "\n".join(L) + "\n"
    limit = MAX_REPORT_KB * 1024
    if len(text.encode()) > limit:
        text = (
            text.encode()[:limit].decode(errors="ignore")
            + f"\n\n(cut at {MAX_REPORT_KB} KB)\n"
        )
    (rep_dir / "report.md").write_text(text, encoding="utf-8")
    print(
        f"Report: {rep_dir / 'report.md'} ({len(text.encode()) / 1e3:.0f} KB) in {time.perf_counter() - start:.0f} s. "
        f"Commit and push the folder {rep_dir}."
    )
