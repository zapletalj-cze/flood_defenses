# Author: Jakub Zapletal
# Date: 2026-11-12
# Explanatory figures of the README (docs/img/): the protection principle on a cross-section, river snapping into the
# channel, the rasters of one synthetic unit and open versus gated river crossings. Synthetic data only: the Stage B
# test country is built (or reused) in TEST_DIR and protected_areas.py runs twice on it, with intermediate rasters.
# Run: pixi run make-figures  (DEFENSES_TEST_DIR keeps the test folder; default: temporary)

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import geopandas as gpd
import matplotlib
import numpy as np
import shapely
from osgeo import gdal, osr

import gis
from matplotlib.colors import LinearSegmentedColormap

from plots import CATEGORICAL, HILLSHADE_GREYS, INK, LINE_WIDTH

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

gdal.UseExceptions()

# --- Run variables (no command-line parameters; DEFENSES_<NAME> environment variables override them) ---
TEST_DIR = ""  # folder of the Stage B test country; empty = temporary
OUT_DIR = "docs/img"
DPI = 110
UNIT = "000_000"  # unit whose rasters are shown
EPSG = 32633  # CRS of the synthetic test data

plt.rcParams.update(
    {
        "figure.facecolor": INK["surface"],
        "axes.facecolor": INK["surface"],
        "savefig.facecolor": INK["surface"],
        "axes.edgecolor": INK["axis"],
        "axes.labelcolor": INK["secondary"],
        "axes.titlecolor": INK["primary"],
        "axes.titlesize": 10,
        "axes.titlelocation": "left",
        "axes.grid": True,
        "axes.axisbelow": True,
        "grid.color": INK["grid"],
        "grid.linewidth": 0.6,
        "xtick.color": INK["axis"],
        "ytick.color": INK["axis"],
        "xtick.labelcolor": INK["secondary"],
        "ytick.labelcolor": INK["secondary"],
        "legend.frameon": False,
        "legend.labelcolor": INK["secondary"],
        "font.size": 9,
    }
)
BLUE, ORANGE, AQUA, YELLOW = (
    CATEGORICAL[0],
    CATEGORICAL[1],
    CATEGORICAL[2],
    CATEGORICAL[3],
)
# elevation in a light grey ramp (low = dark), so that coloured lines and areas stay readable
GREYS = LinearSegmentedColormap.from_list("elevation", HILLSHADE_GREYS)


def read(path):
    ds = gdal.Open(str(path))
    arr = ds.ReadAsArray().astype(float)
    nodata = ds.GetRasterBand(1).GetNoDataValue()
    if nodata is not None:
        arr[arr == nodata] = np.nan
    gt = ds.GetGeoTransform()
    extent = (
        gt[0],
        gt[0] + ds.RasterXSize * gt[1],
        gt[3] + ds.RasterYSize * gt[5],
        gt[3],
    )
    return arr, gt, ds.GetProjection(), extent


def lines(ax, gdf, origin, **kw):
    for g in shapely.get_parts(gdf.geometry.values):
        if shapely.get_type_id(g) in (1, 2):
            xy = shapely.get_coordinates(g) - origin
            ax.plot(xy[:, 0], xy[:, 1], **kw)


def areas(ax, gdf, origin, **kw):
    for g in shapely.get_parts(gdf.geometry.values):
        xy = shapely.get_coordinates(g.exterior) - origin
        ax.fill(xy[:, 0], xy[:, 1], **kw)


if __name__ == "__main__":
    run = gis.Run.env_overrides(TEST_DIR=TEST_DIR, OUT_DIR=OUT_DIR, DPI=DPI, UNIT=UNIT)
    out = Path(run["OUT_DIR"])
    out.mkdir(parents=True, exist_ok=True)
    repo = Path(__file__).resolve().parents[1]

    # --- 1. Principle on a cross-section: spill levels with and without the levee ---
    res = 2.0
    x = np.arange(0.0, 800.0, res)
    terrain = np.where(x < 20, 0.0, 1.0 + 0.004 * x)
    terrain -= np.where(
        (x > 250) & (x < 400), 0.8 * np.sin(np.pi * (x - 250) / 150), 0.0
    )
    levee = (x >= 120) & (x < 136)
    defended = np.where(levee, 3.5, terrain)
    seeds = (x < 20)[None, :]
    spill_undef = gis.Hydro.priority_flood(terrain[None, :].astype(np.float32), seeds)[
        0
    ][0]
    spill_def = gis.Hydro.priority_flood(defended[None, :].astype(np.float32), seeds)[
        0
    ][0]
    margin = spill_def - spill_undef
    protected = (margin >= 0.3) & ~levee
    fig, (ax, ax2) = plt.subplots(
        2,
        1,
        figsize=(9, 5.6),
        sharex=True,
        height_ratios=[3, 1.3],
        layout="constrained",
    )
    ax.fill_between(
        x, -0.5, defended, color="#d9d8d2", lw=0, label="terrain (levee included)"
    )
    ax.plot(
        x,
        spill_undef,
        color=INK["muted"],
        ls=(0, (4, 2)),
        lw=LINE_WIDTH,
        label="spill level without the levee",
    )
    ax.plot(
        x,
        spill_def,
        color=BLUE,
        lw=LINE_WIDTH + 0.4,
        label="spill level with the levee",
    )
    ax.fill_between(
        x,
        spill_undef,
        spill_def,
        where=protected,
        color=BLUE,
        alpha=0.15,
        lw=0,
        label="margin of protected cells",
    )
    ax.annotate("river", (8, 0.1), fontsize=8, color=INK["secondary"])
    ax.annotate(
        "levee crest",
        (128, 3.5),
        (150, 4.35),
        fontsize=8,
        color=INK["secondary"],
        arrowprops={"arrowstyle": "-", "color": INK["muted"], "lw": 0.6},
    )
    ax.annotate(
        "water must rise to the crest\nto reach these cells",
        (320, 3.55),
        fontsize=8,
        color=INK["secondary"],
    )
    ax.annotate(
        "upslope: reached at its own\nheight anyway, not protected",
        (640, 2.0),
        fontsize=8,
        color=INK["secondary"],
    )
    ax.set_ylabel("Height above the river (m)")
    ax.set_ylim(-0.5, 6.2)
    ax.legend(loc="upper left", fontsize=8, ncols=2)
    ax.set_title(
        "Principle: a cell is protected when the levee raises the level at which river water reaches it"
    )
    ax2.plot(x, margin, color=BLUE, lw=LINE_WIDTH)
    ax2.axhline(0.3, color=INK["muted"], lw=0.8, ls=(0, (4, 2)))
    ax2.annotate("min_margin_m", (700, 0.4), fontsize=8, color=INK["secondary"])
    ax2.fill_between(
        x,
        0,
        1,
        where=protected,
        transform=ax2.get_xaxis_transform(),
        color=BLUE,
        alpha=0.08,
        lw=0,
    )
    ax2.set_ylabel("Margin (m)")
    ax2.set_xlabel("Distance from the river (m)")
    fig.savefig(out / "principle_profile.png", dpi=run["DPI"])
    plt.close(fig)

    # --- 2. River snapping: a channel along a sine, the river line 60 m beside it ---
    tmp = Path(tempfile.mkdtemp(prefix="defenses_figures_"))
    nx, ny, res = 1000, 500, 2.0
    gx, gy = np.meshgrid(
        (np.arange(nx) + 0.5) * res, 1000 - (np.arange(ny) + 0.5) * res
    )
    centre = 500 + 40 * np.sin(gx / 150)
    z = (
        100
        + 0.002 * np.abs(gy - centre)
        - 2.0 * (np.abs(gy - centre) < 6)
        + 0.0001 * gx
    )
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(EPSG)
    gt_s = (0, res, 0, 1000, 0, -res)
    gis.Raster.write_array(
        tmp / "dtm.tif", z.astype(np.float32), gt_s, srs.ExportToWkt(), nodata=-9999
    )
    xs = np.linspace(5, 1995, 60)
    river = shapely.linestrings(np.c_[xs, 500 + 40 * np.sin(xs / 150) + 60])
    guide = shapely.linestrings(np.c_[xs, 500 + 40 * np.sin(xs / 150) + 3])
    snap_p = {
        k.split(".")[1]: v
        for k, v in gis.Config.DEFAULTS.items()
        if k.startswith("river_snap.")
    }
    out_g, _ = gis.Hydro.snap_to_channel([river], tmp / "dtm.tif", [guide], snap_p)
    out_d, _ = gis.Hydro.snap_to_channel(
        [river], tmp / "dtm.tif", None, {**snap_p, "max_shift_without_guide_m": 80.0}
    )
    fig, (ax, ax2) = plt.subplots(
        1, 2, figsize=(11, 4.4), width_ratios=[1.6, 1], layout="constrained"
    )
    ax.imshow(z, extent=(0, 2000, 0, 1000), cmap=GREYS, interpolation="nearest")
    for geom, kw in (
        (
            river,
            {
                "color": ORANGE,
                "ls": (0, (4, 2)),
                "lw": LINE_WIDTH + 0.4,
                "label": "river line (coarse network)",
            },
        ),
        (
            guide,
            {"color": AQUA, "lw": LINE_WIDTH, "label": "guide line (river_guides)"},
        ),
        (
            out_d[0],
            {
                "color": YELLOW,
                "lw": LINE_WIDTH + 0.4,
                "label": "snapped, DTM only (80 m cap)",
            },
        ),
        (
            out_g[0],
            {
                "color": BLUE,
                "lw": LINE_WIDTH + 0.4,
                "label": "snapped with the guide line",
            },
        ),
    ):
        xy = shapely.get_coordinates(geom)
        ax.plot(xy[:, 0], xy[:, 1], **kw)
    ax.set_xlim(0, 1000)
    ax.set_ylim(350, 650)
    ax.set_aspect("equal")
    ax.grid(False)
    ax.legend(loc="lower left", fontsize=8)
    ax.set_title(
        "Snapping into the channel (grey: elevation, the channel is the dark band)"
    )
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    for geom, colour, label in (
        (river, ORANGE, "river line"),
        (out_d[0], YELLOW, "DTM only"),
        (out_g[0], BLUE, "with guide"),
    ):
        pts = shapely.get_coordinates(
            shapely.line_interpolate_point(geom, np.arange(0, geom.length, 2.0))
        )
        off = np.abs(pts[:, 1] - (500 + 40 * np.sin(pts[:, 0] / 150)))
        ax2.plot(
            pts[:, 0],
            off,
            color=colour,
            lw=LINE_WIDTH,
            label=f"{label}: median {np.median(off):.1f} m",
        )
    ax2.set_xlabel("x (m)")
    ax2.set_ylabel("Distance from the channel axis (m)")
    ax2.legend(loc="upper right", fontsize=8)
    ax2.set_title("Distance from the channel along the line")
    fig.savefig(out / "river_snapping.png", dpi=run["DPI"])
    plt.close(fig)

    # --- 3. Stage B test country: protected_areas.py with open and gated crossings, intermediate rasters kept ---
    root = Path(run["TEST_DIR"] or tempfile.mkdtemp(prefix="defenses_figures_b_"))
    env = dict(os.environ, DEFENSES_TEST_DIR=str(root))
    if not (root / "test.toml").exists():
        subprocess.run(
            [sys.executable, "-m", "tests.test_stage_b"],
            env=env,
            cwd=repo,
            check=False,
            capture_output=True,
        )
    variants = {}
    for name, config in (("open", "test_open.toml"), ("gated", "test.toml")):
        text = (
            (root / config)
            .read_text()
            .replace("keep_intermediate = false", "keep_intermediate = true")
        )
        (root / f"figure_{name}.toml").write_text(text)
        proc = subprocess.run(
            [sys.executable, "-m", "processing.protected_areas"],
            env=dict(
                env,
                DEFENSES_CONFIG_PATH=str(root / f"figure_{name}.toml"),
                DEFENSES_COUNTRY="T4",
                DEFENSES_FORCE="1",
                DEFENSES_UNITS="",
            ),
            cwd=repo,
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            sys.exit(
                f"protected_areas.py failed for the {name} variant:\n{proc.stderr[-2000:]}"
            )
        keep = root / f"figure_{name}"
        shutil.rmtree(keep, ignore_errors=True)
        shutil.copytree(root / "work" / "T4" / "units" / run["UNIT"], keep)
        shutil.copy(
            root / "output" / "T4" / "protected_areas.gpkg", keep / "merged.gpkg"
        )
        variants[name] = keep
    gdf_lv = gpd.read_file(root / "data" / "levees.gpkg")
    gdf_rn = gpd.read_file(root / "output" / "T4" / "rivers_channel.gpkg")
    dtm, gt_u, wkt, ext = read(variants["gated"] / "dtm.tif")
    full, _, _, ext_f = read(root / "work" / "T4" / "dtm" / "dtm.vrt")

    def base(ax, title, elevation, extent):
        origin = np.array([extent[0], extent[2]])
        ext0 = (0, extent[1] - extent[0], 0, extent[3] - extent[2])
        im = ax.imshow(elevation, extent=ext0, cmap=GREYS, interpolation="nearest")
        lines(ax, gdf_rn, origin, color=BLUE, lw=LINE_WIDTH + 0.4)
        lines(ax, gdf_lv, origin, color=INK["primary"], lw=LINE_WIDTH + 0.4)
        ax.set_xlim(ext0[0], ext0[1])
        ax.set_ylim(ext0[2], ext0[3])
        ax.set_aspect("equal")
        ax.grid(False)
        ax.set_title(title)
        return im, origin

    origin = np.array([ext[0], ext[2]])
    ext0 = (0, ext[1] - ext[0], 0, ext[3] - ext[2])

    # 3a. The rasters of one unit (gated variant)
    hand, *_ = read(variants["gated"] / "hand_undef.tif")
    margin_u, *_ = read(variants["gated"] / "margin.tif")
    labels, *_ = read(variants["gated"] / "labels.tif")
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.2), layout="constrained")
    im, _ = base(
        axes[0, 0],
        "1. Terrain, levee rings (black), rivers in the channel (blue)",
        dtm,
        ext,
    )
    fig.colorbar(im, ax=axes[0, 0], shrink=0.8, label="Elevation (m)")
    im = axes[0, 1].imshow(
        np.clip(hand, 0, 8), extent=ext0, cmap="Blues_r", vmin=0, vmax=8
    )
    lines(axes[0, 1], gdf_lv, origin, color=INK["primary"], lw=LINE_WIDTH)
    axes[0, 1].set_title("2. HAND on the undefended DTM (levees removed)")
    fig.colorbar(im, ax=axes[0, 1], shrink=0.8, label="Height above the river (m)")
    im = axes[1, 0].imshow(
        np.where(margin_u > 0, margin_u, np.nan),
        extent=ext0,
        cmap="Oranges",
        vmin=0,
        vmax=3,
    )
    lines(axes[1, 0], gdf_lv, origin, color=INK["primary"], lw=LINE_WIDTH)
    axes[1, 0].set_title("3. Margin = spill level with minus without the levees")
    fig.colorbar(im, ax=axes[1, 0], shrink=0.8, label="Margin (m)")
    base(
        axes[1, 1],
        "4. Protected areas (blue): margin ≥ min_margin_m, area ≥ min_area_m2",
        dtm,
        ext,
    )
    areas(
        axes[1, 1],
        gpd.read_file(variants["gated"] / "merged.gpkg"),
        origin,
        color=BLUE,
        alpha=0.45,
        lw=0,
    )
    for ax in axes.ravel():
        ax.set_xlim(ext0[0], ext0[1])
        ax.set_ylim(ext0[2], ext0[3])
        ax.set_aspect("equal")
        ax.grid(False)
        ax.tick_params(labelsize=7)
    fig.suptitle(
        "From terrain to protected areas on one unit of the synthetic test country",
        x=0.01,
        ha="left",
    )
    fig.savefig(out / "unit_steps.png", dpi=run["DPI"])
    plt.close(fig)

    # 3b. Open versus gated crossings
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), layout="constrained")
    for ax, (name, title) in zip(
        axes,
        (
            (
                "open",
                "gated_crossings = false: the streams inside rings 1 and 2 are water\nsources and flood them; only ring 3 is protected (strips north of the\nrings: valley floors that the rings shield from those streams)",
            ),
            (
                "gated",
                "gated_crossings = true: barriers stay closed at crossings and enclosed\nstreams are not water sources; all three rings are protected",
            ),
        ),
    ):
        _, org = base(ax, title, full, ext_f)
        areas(
            ax,
            gpd.read_file(variants[name] / "merged.gpkg"),
            org,
            color=BLUE,
            alpha=0.45,
            lw=0,
        )
        ax.tick_params(labelsize=7)
    fig.savefig(out / "gated_crossings.png", dpi=run["DPI"])
    plt.close(fig)
    sizes = {p.name: round(p.stat().st_size / 1e3) for p in sorted(out.glob("*.png"))}
    print(f"Figures in {out} (KB): {sizes}")
