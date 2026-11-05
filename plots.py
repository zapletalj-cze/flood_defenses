# Author: Jakub Zapletal
# Date: 2026-11-04
# Shared plotting helpers of the Stage E scripts (matplotlib): backend, figure budget per type, saving with an
# optional non-blocking display, fixed colours per entity, and the static index.html of a plot run folder.

import html
import json
import logging
import os
import shutil
import sys
from pathlib import Path

import matplotlib
import numpy as np

import gis

# Validated categorical palette (light surface), fixed per entity and never assigned by rank.
CATEGORICAL = [
    "#2a78d6",
    "#eb6834",
    "#1baf7a",
    "#eda100",
    "#e87ba4",
    "#008300",
    "#6250d6",
    "#e34948",
]
DASHES = ["-", (0, (4, 2)), (0, (1, 1.5)), (0, (6, 2, 1, 2)), (0, (2, 2))]
COLOURS = {
    "water": CATEGORICAL[0],
    "snapped": CATEGORICAL[1],
    "original": CATEGORICAL[2],
    # contexts of zones.gpkg, slots in the fixed order; no zone is neutral
    "context": {
        "none": "#c3c2b7",
        "industrial": CATEGORICAL[0],
        "dense_urban": CATEGORICAL[1],
        "urban": CATEGORICAL[2],
        "rural": CATEGORICAL[3],
    },
}
# Ink and chrome: text never takes a data colour.
INK = {
    "primary": "#0b0b0b",
    "secondary": "#52514e",
    "muted": "#898781",
    "grid": "#e1e0d9",
    "axis": "#c3c2b7",
    "surface": "#fcfcfb",
}
LINE_WIDTH = 1.2  # points; 2 px at 120 dpi
HILLSHADE_GREYS = ("#a3a29c", "#f6f6f3")  # light grey ramp under coloured lines
HILLSHADE_STRETCH = (
    2,
    98,
)  # percentiles mapped to the ends of the ramp (flat floodplains stay readable)
MAP_MIN_ASPECT = 1 / 3  # short side of a map window at least a third of the long side
MAP_MAX_CELLS = (
    1500  # cells along the long side of a map background; coarser cells beyond
)
DASHED = (
    0,
    (4, 2),
)  # secondary encoding of the original line (its colour is below 3:1 on the surface)
PAUSE_S = 0.2  # display time of a shown figure before processing continues


class Plots:
    """
    Figures of one run folder <output_dir>/<country>/plots/implied/<run_id>: every figure
    type has its own folder with the PNG files and a manifest.json of the plotted items
    and their metrics; index.html of the run folder lists all types.
    """

    def __init__(self, run_dir, settings, show=False, logger=None):
        """
        :param run_dir: plot run folder
        :param settings: the [implied.plots] table (dpi, max_per_type, sample_seed,
                         show_max, plot_all_in_units)
        :param show: display figures while processing (SHOW_PLOTS); ignored without a display
        :param logger: logging.Logger
        """
        self.run_dir = Path(run_dir)
        self.settings = settings
        self.log = logger or logging.getLogger("plots")
        self.show = bool(show) and Plots.display_available()
        if show and not self.show:
            self.log.warning(
                "SHOW_PLOTS: no display available; figures are only saved."
            )
        if not self.show:
            matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt

        self.plt = plt
        plt.rcParams.update(
            {
                "figure.facecolor": INK["surface"],
                "axes.facecolor": INK["surface"],
                "savefig.facecolor": INK["surface"],
                "axes.edgecolor": INK["axis"],
                "axes.linewidth": 0.6,
                "axes.labelcolor": INK["secondary"],
                "axes.titlecolor": INK["primary"],
                "axes.titlesize": 10,
                "axes.titlelocation": "left",
                "axes.spines.top": False,
                "axes.spines.right": False,
                "axes.grid": True,
                "grid.color": INK["grid"],
                "grid.linewidth": 0.6,
                "grid.linestyle": "-",
                "xtick.color": INK["axis"],
                "ytick.color": INK["axis"],
                "xtick.labelcolor": INK["secondary"],
                "ytick.labelcolor": INK["secondary"],
                "font.size": 9,
                "legend.frameon": False,
                "legend.labelcolor": INK["secondary"],
                "lines.linewidth": LINE_WIDTH,
                "lines.solid_capstyle": "round",
                "lines.solid_joinstyle": "round",
                "axes.formatter.useoffset": False,
                "axes.formatter.limits": [-9, 9],
            }
        )
        self.entries = {}
        self.shown = {}
        self.open = []

    @staticmethod
    def display_available():
        """
        :return: True on Windows and macOS, or with an X11/Wayland display on Linux
        """
        return sys.platform in ("win32", "darwin") or bool(
            os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY")
        )

    @staticmethod
    def method_styles(methods):
        """
        Line style per level method (onset first, then the quantiles): one categorical
        slot and one dash pattern each, in the order given, so a method keeps its style
        in every figure of a run
        :param methods: method names in a fixed order
        :return: dict name -> matplotlib line keyword arguments
        """
        return {
            m: {
                "color": CATEGORICAL[i % len(CATEGORICAL)],
                "ls": DASHES[i % len(DASHES)],
                "lw": LINE_WIDTH + (0.4 if i == 0 else 0.0),
            }
            for i, m in enumerate(methods)
        }

    @staticmethod
    def halo():
        """
        :return: path effects drawing a surface-coloured ring under a line on a busy background
        """
        from matplotlib import patheffects

        return [
            patheffects.Stroke(linewidth=LINE_WIDTH + 2.0, foreground="white"),
            patheffects.Normal(),
        ]

    @staticmethod
    def hillshade_image(ax, shade, extent):
        """
        Draws a hillshade (uint8, 0 = nodata) as a light grey background, stretched
        between the HILLSHADE_STRETCH percentiles; nodata stays transparent
        :param ax: matplotlib axes
        :param shade: hillshade array
        :param extent: (xmin, xmax, ymin, ymax) of the array
        """
        from matplotlib.colors import LinearSegmentedColormap, Normalize

        valid = shade[shade > 0]
        lo, hi = np.percentile(valid, HILLSHADE_STRETCH) if len(valid) else (1.0, 255.0)
        ax.imshow(
            np.ma.masked_equal(shade, 0),
            cmap=LinearSegmentedColormap.from_list("hillshade", HILLSHADE_GREYS),
            norm=Normalize(lo, max(hi, lo + 1.0)),
            extent=extent,
            interpolation="nearest",
        )

    @staticmethod
    def window(bounds, resolution):
        """
        Map window around bounds: the short side padded to MAP_MIN_ASPECT of the long
        side, the cell size coarsened so that the long side has at most MAP_MAX_CELLS
        cells, the bounds aligned to the cell size
        :param bounds: (xmin, ymin, xmax, ymax)
        :param resolution: finest cell size
        :return: (aligned bounds, cell size)
        """
        x0, y0, x1, y1 = bounds
        w, h = x1 - x0, y1 - y0
        if h < w * MAP_MIN_ASPECT:
            pad = (w * MAP_MIN_ASPECT - h) / 2
            y0, y1 = y0 - pad, y1 + pad
        elif w < h * MAP_MIN_ASPECT:
            pad = (h * MAP_MIN_ASPECT - w) / 2
            x0, x1 = x0 - pad, x1 + pad
        cell = max(resolution, max(x1 - x0, y1 - y0) / MAP_MAX_CELLS)
        return (
            np.floor(x0 / cell) * cell,
            np.floor(y0 / cell) * cell,
            np.ceil(x1 / cell) * cell,
            np.ceil(y1 / cell) * cell,
        ), cell

    @staticmethod
    def hillshade_layer(ax, ds, bounds, cell, warp_mb):
        """
        Reads the DTM window (average to cell), draws its hillshade and sets the axes
        to the window
        :param ax: matplotlib axes
        :param ds: DTM dataset
        :param bounds: aligned window bounds
        :param cell: cell size
        :param warp_mb: warp memory in MB
        :return: extent (xmin, xmax, ymin, ymax)
        """
        z, gt = gis.Raster.read_window(ds, bounds, cell, "average", warp_mb)
        extent = (gt[0], gt[0] + z.shape[1] * cell, gt[3] - z.shape[0] * cell, gt[3])
        Plots.hillshade_image(
            ax, gis.Raster.hillshade(z, gt, ds.GetProjection()), extent
        )
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
        ax.set_aspect("equal")
        ax.grid(False)
        return extent

    def select(self, kind, ids, flagged=(), units_of=None):
        """
        Items to plot for one figure type: flagged items first, then a random sample
        (sample_seed), up to max_per_type; plus every item of the units in
        plot_all_in_units. Earlier figures of this type in the run folder are removed.
        :param kind: figure type (folder name)
        :param ids: item names in a stable order
        :param flagged: names of flagged or anomalous items
        :param units_of: dict name -> unit ids the item touches (for plot_all_in_units)
        :return: list of names to plot
        """
        folder = self.run_dir / kind
        if folder.exists():
            shutil.rmtree(folder)
        folder.mkdir(parents=True, exist_ok=True)
        self.entries[kind] = []
        ids = list(ids)
        flagged = set(flagged)
        first = [i for i in ids if i in flagged]
        rest = [i for i in ids if i not in flagged]
        order = np.random.default_rng(self.settings["sample_seed"]).permutation(
            len(rest)
        )
        list_sel = (first + [rest[k] for k in order])[: self.settings["max_per_type"]]
        list_all = []
        wanted = set(self.settings["plot_all_in_units"])
        if wanted and units_of:
            chosen = set(list_sel)
            list_all = [
                i for i in ids if i not in chosen and wanted & set(units_of.get(i, ()))
            ]
        self.log.info(
            "%s: %d of %d items plotted (%d flagged, %d more from plot_all_in_units).",
            kind,
            len(list_sel) + len(list_all),
            len(ids),
            len(first),
            len(list_all),
        )
        return list_sel + list_all

    def save(self, fig, kind, name, metrics):
        """
        Saves a figure as <run_dir>/<kind>/<name>.png and records it with its metrics;
        with SHOW_PLOTS the first show_max figures of a type stay open, shown without
        blocking, all others are closed.
        :return: path of the PNG file
        """
        path = self.run_dir / kind / f"{name}.png"
        fig.savefig(path, dpi=self.settings["dpi"])
        self.entries.setdefault(kind, []).append(
            {"name": name, "file": f"{kind}/{name}.png", "metrics": metrics}
        )
        if self.show and self.shown.get(kind, 0) < self.settings["show_max"]:
            self.shown[kind] = self.shown.get(kind, 0) + 1
            fig.show()
            self.plt.pause(PAUSE_S)
            self.open.append(fig)
        else:
            self.plt.close(fig)
        return path

    def finish(self):
        """
        Writes manifest.json of every plotted type and index.html of the run folder;
        closes the figures kept open for display
        :return: path of index.html
        """
        for kind, list_entry in self.entries.items():
            (self.run_dir / kind / "manifest.json").write_text(
                json.dumps(list_entry, indent=1)
            )
        for fig in self.open:
            self.plt.close(fig)
        self.open = []
        return Plots.write_index(self.run_dir)

    @staticmethod
    def write_index(run_dir):
        """
        Static index.html of a run folder (standard library only): per figure type,
        thumbnails linking to the full images and a table of the metrics per item.
        Types are read from the manifest.json files, so figures of every script that
        wrote into the run folder are listed.
        :param run_dir: plot run folder
        :return: path of index.html
        """
        run_dir = Path(run_dir)
        sections = []
        for manifest in sorted(run_dir.glob("*/manifest.json")):
            kind = manifest.parent.name
            list_entry = json.loads(manifest.read_text())
            columns = sorted({k for e in list_entry for k in e["metrics"]})
            thumbs = "".join(
                f'<figure><a href="{html.escape(e["file"])}"><img src="{html.escape(e["file"])}" '
                f'alt="{html.escape(e["name"])}" loading="lazy"></a>'
                f'<figcaption>{html.escape(e["name"])}</figcaption></figure>'
                for e in list_entry
            )
            head = "".join(f"<th>{html.escape(c)}</th>" for c in ["item"] + columns)
            rows = "".join(
                "<tr>"
                + f'<td><a href="{html.escape(e["file"])}">{html.escape(e["name"])}</a></td>'
                + "".join(
                    f"<td>{html.escape(Plots.format_value(e['metrics'].get(c, '')))}</td>"
                    for c in columns
                )
                + "</tr>"
                for e in list_entry
            )
            sections.append(
                f'<section id="{html.escape(kind)}"><h2>{html.escape(kind)} '
                f"<small>{len(list_entry)} figures</small></h2>"
                f'<div class="thumbs">{thumbs}</div>'
                f'<div class="table"><table><thead><tr>{head}</tr></thead><tbody>{rows}</tbody></table></div>'
                "</section>"
            )
        nav = " · ".join(
            f'<a href="#{html.escape(m.parent.name)}">{html.escape(m.parent.name)}</a>'
            for m in sorted(run_dir.glob("*/manifest.json"))
        )
        page = (
            '<!doctype html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width, initial-scale=1">'
            f"<title>Plots {html.escape(run_dir.name)}</title><style>"
            f"body{{margin:0 auto;max-width:1400px;padding:16px;background:{INK['surface']};"
            f"color:{INK['primary']};font:14px system-ui,-apple-system,'Segoe UI',sans-serif}}"
            f"h1{{font-size:20px}}h2{{font-size:16px;margin-top:32px}}small{{color:{INK['muted']};font-weight:normal}}"
            f"a{{color:{INK['primary']}}}nav{{color:{INK['secondary']}}}"
            ".thumbs{display:flex;flex-wrap:wrap;gap:12px}figure{margin:0;width:240px}"
            f"figure img{{width:240px;height:auto;border:1px solid {INK['grid']}}}"
            f"figcaption{{font-size:12px;color:{INK['secondary']};overflow-wrap:anywhere}}"
            ".table{overflow-x:auto;margin-top:12px}table{border-collapse:collapse;font-size:12px}"
            f"th,td{{padding:4px 8px;border-bottom:1px solid {INK['grid']};text-align:left;"
            "font-variant-numeric:tabular-nums;white-space:nowrap}"
            "</style></head><body>"
            f"<h1>Plots, run {html.escape(run_dir.name)}</h1><nav>{nav}</nav>"
            + "".join(sections)
            + "</body></html>"
        )
        path = run_dir / "index.html"
        run_dir.mkdir(parents=True, exist_ok=True)
        path.write_text(page, encoding="utf-8")
        return path

    @staticmethod
    def format_value(value):
        """
        :return: a metric value as short text for the index table
        """
        if isinstance(value, bool):
            return "yes" if value else "no"
        if isinstance(value, (int, np.integer)):
            return str(int(value))
        if isinstance(value, (float, np.floating)):
            if not np.isfinite(value):
                return ""
            return f"{value:.2f}" if abs(value) < 1000 else f"{value:.0f}"
        return str(value)
