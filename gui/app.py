# Author: Jakub Zapletal
# Date: 2026-11-12
# Desktop GUI of the Defenses package, an alternative to the shell: a local web application (standard library
# server, page in an app window of Edge or Chrome) to edit the configuration, run the steps with their run variables,
# follow the console with progress and look at results and logs. Listens on 127.0.0.1 only.
# Run: pixi run gui

import ast
import glob
import json
import os
import re
import secrets
import shutil
import signal
import subprocess
import sys
import threading
import time
import tomllib
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import pyogrio

import gis

# --- Run variables (DEFENSES_<NAME> environment variables override them) ---
CONFIG_PATH = "config/defenses.toml"
PORT = 0  # 0 = any free port
OPEN_WINDOW = True  # open an app window of Edge or Chrome (else the default browser)

ROOT = Path(__file__).resolve().parents[1]
STATIC = Path(__file__).resolve().parent / "static"
MAX_LINES = 20000  # console lines kept per job
TYPES = {
    ".html": "text/html; charset=utf-8",
    ".css": "text/css; charset=utf-8",
    ".js": "text/javascript; charset=utf-8",
    ".svg": "image/svg+xml",
    ".png": "image/png",
    ".md": "text/plain; charset=utf-8",
}

# Steps in run order: module, group, title, description, run variables shown, outputs (globs in the output folder,
# "work:" = in the work folder, "repo:" = in the repository), log name in output_dir/<country>/logs
GROUPS = [
    (
        "core",
        "Core pipeline",
        "From the DTM to protected areas and their return period",
    ),
    (
        "observations",
        "Observations",
        "Satellite flood maps and river discharge as evidence",
    ),
    ("exposure", "Exposure", "Land cover as a supplement to building volume"),
    (
        "context",
        "Context and implied protection",
        "Stages D and E, experimental",
    ),
    ("tools", "Report and tools", "Summaries to share, figures and self-tests"),
]
STEPS = [
    (
        "processing.check_inputs",
        "core",
        "Check inputs",
        "Keys, files, columns and CRS of every step before a run",
        [],
        [],
        None,
    ),
    (
        "processing.prepare_dtm",
        "core",
        "Prepare DTM",
        "Inspection, footprint, read benchmark, optional conversion",
        ["FORCE", "RUN_BENCHMARK"],
        ["dtm_inspection.csv", "dtm_footprint.gpkg"],
        "prepare_dtm",
    ),
    (
        "processing.build_units",
        "core",
        "Build units",
        "Corridor along the rivers and processing units with a halo",
        [],
        ["units.csv", "units_core.gpkg"],
        "build_units",
    ),
    (
        "processing.protected_areas",
        "core",
        "Protected areas",
        "Rivers and levees onto the DTM, HAND and floods per unit, files A and B, leaks",
        ["UNITS", "FORCE"],
        [
            "protected_areas.gpkg",
            "protected_areas_by_group.gpkg",
            "leaks.gpkg",
            "levees_crest.gpkg",
        ],
        "protected_areas",
    ),
    (
        "processing.transport_embankments",
        "core",
        "Transport embankments",
        "Road and railway embankments from the DTM, their openings (optional: transport_lines)",
        ["RUN_ID"],
        ["transport_embankments.gpkg", "openings.gpkg"],
        "transport_embankments",
    ),
    (
        "download.download_ghsl",
        "core",
        "Building volume (GHSL)",
        "Latest GHSL epoch mosaicked for the country",
        [],
        ["ghsl_total.vrt"],
        "download_ghsl",
    ),
    (
        "processing.rp_index",
        "core",
        "Return period index",
        "Exposed value, score, return period class and interval",
        ["UNITS", "FORCE", "LAYER"],
        ["protected_areas_rp.gpkg"],
        "rp_index",
    ),
    (
        "download.download_gfm",
        "observations",
        "GFM scenes",
        "Copernicus Global Flood Monitoring scenes over the units",
        ["UNITS", "FORCE"],
        ["gfm_index.csv"],
        "download_gfm",
    ),
    (
        "processing.gfm_evidence",
        "observations",
        "GFM evidence",
        "Held, failed or flooded per protected area and scene",
        ["UNITS", "HINTS", "DATE_FROM", "DATE_TO", "LAYER"],
        ["gfm_evidence.gpkg"],
        "gfm_evidence",
    ),
    (
        "download.download_discharge",
        "observations",
        "Discharge download",
        "GloFAS and EFAS reanalysis from the EWDS",
        [],
        ["work:discharge/raw/download_log.csv"],
        "download_discharge",
    ),
    (
        "processing.discharge_prepare",
        "observations",
        "Discharge series",
        "Daily series per river cell",
        ["RUN_ID"],
        ["discharge_*.npz"],
        "discharge_prepare",
    ),
    (
        "processing.rp_bounds",
        "observations",
        "Observed bounds",
        "Event return periods and bounds of the return period",
        ["RUN_ID", "LAYER"],
        ["event_rp_combined.csv"],
        "rp_bounds",
    ),
    (
        "processing.flood_events",
        "observations",
        "Flood events",
        "Events per river group from the discharge peaks (Stage H)",
        ["RUN_ID"],
        ["events.gpkg"],
        "flood_events",
    ),
    (
        "processing.validate",
        "observations",
        "Validation",
        "Coverage and return period confusion matrix",
        ["LAYER"],
        ["validation.json"],
        "validate",
    ),
    (
        "download.download_worldcover",
        "exposure",
        "WorldCover",
        "ESA WorldCover tiles for the corridor",
        [],
        ["worldcover.vrt"],
        "download_worldcover",
    ),
    (
        "processing.HAND_YZ",
        "exposure",
        "Land cover exposure",
        "Exposure from land cover per protected area",
        ["UNITS", "FORCE", "LAYER"],
        ["exposure_lulc.csv"],
        "HAND_YZ",
    ),
    (
        "download.download_clms",
        "context",
        "CORINE and Urban Atlas",
        "Copernicus land monitoring layers",
        [],
        ["corine_codes.tif"],
        "download_clms",
    ),
    (
        "processing.context_zones",
        "context",
        "Context zones",
        "Industrial sites, dense urban areas",
        [],
        ["zones.gpkg"],
        "context_zones",
    ),
    (
        "processing.snap_rivers",
        "context",
        "Thalweg snapping",
        "Rivers onto the thalweg for implied protection",
        ["RUN_ID"],
        ["rivers_snapped.gpkg"],
        "snap_rivers",
    ),
    (
        "processing.implied_segments",
        "context",
        "Bank subsegments",
        "Envelope and bank subsegments along the rivers",
        ["FORCE"],
        ["subsegments.gpkg"],
        "implied_segments",
    ),
    (
        "processing.implied_protection",
        "context",
        "Implied protection",
        "Built-up curves, levels and implied protected areas",
        ["RUN_ID"],
        ["implied/*/implied_protected_areas.gpkg"],
        "implied_protection",
    ),
    (
        "processing.implied_evaluate",
        "context",
        "Implied evaluation",
        "Sensitivity of the implied protection",
        ["FORCE"],
        ["sensitivity.csv"],
        "implied_evaluate",
    ),
    (
        "processing.data_report",
        "tools",
        "Report",
        "Compact report of inputs, outputs and logs to share (reports/<country>)",
        [],
        ["repo:reports/{country}/report.md"],
        None,
    ),
    (
        "docs.make_figures",
        "tools",
        "Explanatory figures",
        "Figures of the README from synthetic data",
        [],
        ["repo:docs/img/unit_steps.png"],
        None,
    ),
    (
        "tests.test_kernels",
        "tools",
        "Self-test",
        "Checks of the numerical kernels on small synthetic grids",
        [],
        [],
        None,
    ),
]
LEAK_TEXT = {
    "river_crossing": "Open river crossing through the levee",
    "low_crest": "Water passes over a low crest",
    "line_end": "Water flows around a free line end",
    "water_source_behind": "Stream behind the levee",
    "terrain_above_crest": "Terrain behind lies above the crest",
}


def script_default(module, name):
    # value of a run variable at the top of a script (NAME = literal), None if absent or not a literal
    path = ROOT / Path(*module.split(".")).with_suffix(".py")
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"))
    except (OSError, SyntaxError):
        return None
    for node in tree.body:
        if (
            isinstance(node, ast.Assign)
            and len(node.targets) == 1
            and isinstance(node.targets[0], ast.Name)
            and node.targets[0].id == name
        ):
            try:
                return ast.literal_eval(node.value)
            except ValueError:
                return None
    return None


class Jobs:
    """Queue of step runs, one at a time, with their console lines and progress"""

    def __init__(self):
        self.lock = threading.Condition()
        self.items = []
        self.next_id = 1
        threading.Thread(target=self.worker, daemon=True).start()

    def add(self, module, title, country, variables):
        with self.lock:
            job = {
                "id": self.next_id,
                "module": module,
                "title": title,
                "country": country,
                "variables": variables,
                "status": "queued",
                "created": time.time(),
                "started": None,
                "ended": None,
                "returncode": None,
                "lines": [],
                "offset": 0,
                "progress": None,
                "eta": "",
                "warnings": 0,
                "errors": 0,
                "proc": None,
            }
            self.next_id += 1
            self.items.append(job)
            self.lock.notify_all()
            return job["id"]

    def get(self, job_id):
        with self.lock:
            return next((j for j in self.items if j["id"] == job_id), None)

    def summary(self, job):
        keys = [
            "id",
            "module",
            "title",
            "country",
            "variables",
            "status",
            "created",
            "started",
            "ended",
            "returncode",
            "progress",
            "eta",
            "warnings",
            "errors",
        ]
        out = {k: job[k] for k in keys}
        out["n_lines"] = job["offset"] + len(job["lines"])
        return out

    def worker(self):
        while True:
            with self.lock:
                while not any(j["status"] == "queued" for j in self.items):
                    self.lock.wait()
                job = next(j for j in self.items if j["status"] == "queued")
                job["status"] = "running"
                job["started"] = time.time()
            self.run(job)

    def run(self, job):
        env = {k: v for k, v in os.environ.items() if not k.startswith("DEFENSES_")}
        env.update(
            {
                "PYTHONUNBUFFERED": "1",
                "PYTHONIOENCODING": "utf-8",
                "DEFENSES_CONFIG_PATH": str(SETTINGS["config_path"]),
                "DEFENSES_COUNTRY": job["country"],
            }
        )
        env.update({f"DEFENSES_{k}": str(v) for k, v in job["variables"].items()})
        kwargs = (
            {"creationflags": subprocess.CREATE_NEW_PROCESS_GROUP}
            if os.name == "nt"
            else {"start_new_session": True}
        )
        try:
            proc = subprocess.Popen(
                [sys.executable, "-m", job["module"]],
                cwd=str(ROOT),
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                encoding="utf-8",
                errors="replace",
                bufsize=1,
                **kwargs,
            )
        except OSError as exc:
            self.append(job, f"Cannot start: {exc}")
            with self.lock:
                job.update(status="failed", ended=time.time(), returncode=-1)
            return
        job["proc"] = proc
        for line in proc.stdout:
            self.append(job, line.rstrip("\n"))
        proc.wait()
        with self.lock:
            job["ended"] = time.time()
            job["returncode"] = proc.returncode
            if job["status"] == "stopping":
                job["status"] = "stopped"
            else:
                job["status"] = "done" if proc.returncode == 0 else "failed"
            job["proc"] = None

    def append(self, job, line):
        with self.lock:
            job["lines"].append(line)
            if len(job["lines"]) > MAX_LINES:
                drop = len(job["lines"]) - MAX_LINES
                del job["lines"][:drop]
                job["offset"] += drop
            m = re.search(r"\[(\d+)/(\d+)\]", line)
            if m and int(m.group(2)) > 0:
                job["progress"] = [int(m.group(1)), int(m.group(2))]
            m = re.search(r"ETA (\d+:\d\d:\d\d)", line)
            if m:
                job["eta"] = m.group(1)
            job["warnings"] += " WARNING " in line
            job["errors"] += (" ERROR " in line) or line.startswith("Traceback")

    def stop(self, job_id):
        with self.lock:
            job = next((j for j in self.items if j["id"] == job_id), None)
            if job is None:
                return False
            if job["status"] == "queued":
                job.update(status="cancelled", ended=time.time())
                return True
            if job["status"] != "running" or job["proc"] is None:
                return False
            job["status"] = "stopping"
            pid = job["proc"].pid
        # the step may run worker processes: the whole process tree is stopped
        if os.name == "nt":
            subprocess.run(
                ["taskkill", "/F", "/T", "/PID", str(pid)], capture_output=True
            )
        else:
            try:
                os.killpg(pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
        return True

    def clear(self):
        with self.lock:
            self.items = [
                j
                for j in self.items
                if j["status"] in ("queued", "running", "stopping")
            ]

    def busy(self):
        with self.lock:
            return any(j["status"] in ("running", "stopping") for j in self.items)


def read_config():
    # the user configuration alone (not merged), and its parse error if any
    path = SETTINGS["config_path"]
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        return "", {}, f"Cannot read {path}: {exc}"
    try:
        return text, tomllib.loads(text), ""
    except tomllib.TOMLDecodeError as exc:
        return text, {}, f"TOML error: {exc}"


def folders(country):
    # output and work folders of a country from the merged configuration
    try:
        general = gis.Config.read_merged(str(SETTINGS["config_path"])).get(
            "general", {}
        )
    except (SystemExit, OSError, tomllib.TOMLDecodeError):
        general = {}
    out = Path(general.get("output_dir", "outputs")) / country
    work = Path(general.get("work_dir", "work")) / country
    return (out if out.is_absolute() else ROOT / out), (
        work if work.is_absolute() else ROOT / work
    )


def find(pattern, out_dir, work_dir, country):
    # newest file matching an output pattern, or None
    if pattern.startswith("work:"):
        base, pattern = work_dir, pattern[5:]
    elif pattern.startswith("repo:"):
        base, pattern = ROOT, pattern[5:].format(country=country)
    else:
        base = out_dir
    hits = [Path(p) for p in glob.glob(str(base / pattern), recursive=True)]
    return max(hits, key=lambda p: p.stat().st_mtime) if hits else None


def step_states(country):
    out_dir, work_dir = folders(country)
    logs = out_dir / "logs"
    with JOBS.lock:
        jobs = [j for j in JOBS.items if j["country"] == country]
    list_state = []
    for module, group, title, about, options, outputs, log_name in STEPS:
        hits = [find(p, out_dir, work_dir, country) for p in outputs]
        found = [h for h in hits if h is not None]
        newest_log = None
        if log_name and logs.exists():
            list_log = sorted(logs.glob(f"{log_name}_*.log"))
            newest_log = list_log[-1] if list_log else None
        times = [p.stat().st_mtime for p in found] + (
            [newest_log.stat().st_mtime] if newest_log else []
        )
        last_job = next((j for j in reversed(jobs) if j["module"] == module), None)
        # the first output marks a finished step; the others are shown when present
        status = "done" if hits and hits[0] is not None else "idle"
        if last_job and last_job["status"] in ("queued", "running", "stopping"):
            status = last_job["status"]
        elif last_job and last_job["status"] in ("failed", "stopped"):
            status = last_job["status"]
        list_state.append(
            {
                "module": module,
                "group": group,
                "title": title,
                "about": about,
                "options": [
                    {"name": n, "default": script_default(module, n)} for n in options
                ],
                "outputs": [
                    {
                        "name": p.split(":")[-1].format(country=country),
                        "found": h is not None,
                    }
                    for p, h in zip(outputs, hits)
                ],
                "status": status,
                "last": max(times) if times else None,
                "log": newest_log.name if newest_log else None,
                "job": last_job["id"] if last_job else None,
            }
        )
    return list_state


def layer_summary(path, columns):
    # attributes of a layer without geometry, cached by modification time
    key = (str(path), path.stat().st_mtime)
    if key not in CACHE:
        CACHE[key] = pyogrio.read_dataframe(
            str(path),
            read_geometry=False,
            columns=[c for c in columns if c in pyogrio.read_info(str(path))["fields"]],
        )
    return CACHE[key]


def results(country):
    out_dir, work_dir = folders(country)
    res = {"output_dir": str(out_dir), "exists": out_dir.exists()}
    if (out_dir / "units.csv").exists():
        n_units = max(0, len((out_dir / "units.csv").read_text().splitlines()) - 1)
        n_done = len(list((work_dir / "units").glob("*/done.json")))
        res["units"] = {"total": n_units, "done": n_done}
    path = out_dir / "protected_areas.gpkg"
    if path.exists():
        df = layer_summary(path, ["area_m2", "protection_level_m", "truncated"])
        res["protected"] = {
            "count": len(df),
            "area_km2": float(df["area_m2"].sum() / 1e6) if len(df) else 0.0,
            "level_median_m": (
                float(df["protection_level_m"].median()) if len(df) else None
            ),
        }
    path = out_dir / "protected_areas_by_group.gpkg"
    if path.exists():
        df = layer_summary(path, ["area_m2", "group_id"])
        res["by_group"] = {
            "count": len(df),
            "groups": int(df["group_id"].nunique()) if len(df) else 0,
            "area_km2": float(df["area_m2"].sum() / 1e6) if len(df) else 0.0,
        }
    path = out_dir / "leaks.gpkg"
    if path.exists():
        df = layer_summary(
            path,
            [
                "leak_id",
                "cause",
                "x",
                "y",
                "levee_length_m",
                "lost_margin_median_m",
                "line_ids",
            ],
        )
        by = (
            df.groupby("cause")["levee_length_m"].agg(["size", "sum"]).reset_index()
            if len(df)
            else None
        )
        res["leaks"] = {
            "count": len(df),
            "by_cause": [
                {
                    "cause": r.cause,
                    "text": LEAK_TEXT.get(r.cause, r.cause),
                    "places": int(r.size),
                    "km": float(r.sum / 1000),
                }
                for r in (
                    by.sort_values("sum", ascending=False).itertuples()
                    if by is not None
                    else []
                )
            ],
            "top": json.loads(df.head(12).to_json(orient="records")),
        }
    path = out_dir / "levees_crest.gpkg"
    if path.exists():
        df = layer_summary(path, ["no_ridge", "ridge_share", "shift_median_m"])
        res["levees"] = {
            "lines": len(df),
            "no_ridge": int(df["no_ridge"].sum()) if "no_ridge" in df else 0,
            "moved": (
                int((df["shift_median_m"] > 0).sum()) if "shift_median_m" in df else 0
            ),
        }
    path = out_dir / "performance.json"
    if path.exists():
        res["performance"] = json.loads(path.read_text())
    overview = ROOT / "reports" / country / "overview.png"
    res["overview"] = f"reports/{country}/overview.png" if overview.exists() else None
    res["files"] = (
        [
            {"name": p.name, "size": p.stat().st_size, "mtime": p.stat().st_mtime}
            for p in sorted(
                out_dir.iterdir(), key=lambda p: p.stat().st_mtime, reverse=True
            )
            if p.is_file()
        ][:60]
        if out_dir.exists()
        else []
    )
    return res


def browse(kind, multiple):
    # native file or folder dialog of the machine (tkinter); (paths, error text)
    try:
        import tkinter
        from tkinter import filedialog
    except ImportError:
        return [], "No file dialog on this computer (tkinter missing); type the path."
    with BROWSE_LOCK:
        try:
            root = tkinter.Tk()
        except tkinter.TclError:
            return [], "No file dialog without a display; type the path."
        root.withdraw()
        root.attributes("-topmost", True)
        if kind == "folder":
            paths = [filedialog.askdirectory(parent=root)]
        elif multiple:
            paths = list(filedialog.askopenfilenames(parent=root))
        else:
            paths = [filedialog.askopenfilename(parent=root)]
        root.destroy()
    return [str(Path(p)).replace("\\", "/") for p in paths if p], ""


def open_folder(path):
    path = Path(path)
    if not path.exists():
        return False
    try:
        if os.name == "nt":
            os.startfile(str(path))  # noqa: S606 - local folder chosen by the GUI
        else:
            subprocess.Popen(["xdg-open", str(path)])
    except OSError:
        return False
    return True


class Handler(BaseHTTPRequestHandler):
    server_version = "Defenses"

    def log_message(self, fmt, *args):
        pass

    def send(self, code, body, ctype="application/json; charset=utf-8"):
        data = (
            body if isinstance(body, bytes) else json.dumps(body, default=str).encode()
        )
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(data)

    def allowed(self):
        # local host names only (DNS rebinding) and the session token on every API call
        host = self.headers.get("Host", "").split(":")[0]
        if host not in ("127.0.0.1", "localhost"):
            return False
        return not self.path.startswith("/api/") or (
            self.headers.get("X-Defenses-Token") == TOKEN
        )

    def body(self):
        size = int(self.headers.get("Content-Length") or 0)
        return json.loads(self.rfile.read(size) or b"{}")

    def do_GET(self):
        if not self.allowed():
            return self.send(403, {"error": "forbidden"})
        url = urlparse(self.path)
        q = {k: v[0] for k, v in parse_qs(url.query).items()}
        try:
            if url.path in ("/", "/index.html"):
                html = (STATIC / "index.html").read_text(encoding="utf-8")
                return self.send(
                    200, html.replace("{{TOKEN}}", TOKEN).encode(), TYPES[".html"]
                )
            if url.path.startswith("/static/"):
                path = (STATIC / url.path[len("/static/") :]).resolve()
                if STATIC.resolve() not in path.parents or not path.is_file():
                    return self.send(404, {"error": "not found"})
                return self.send(
                    200,
                    path.read_bytes(),
                    TYPES.get(path.suffix, "application/octet-stream"),
                )
            if url.path == "/api/state":
                text, cfg, error = read_config()
                countries = sorted((cfg.get("countries") or {}).keys())
                country = (
                    q.get("country")
                    or str(cfg.get("general", {}).get("country", ""))
                    or (countries[0] if countries else "")
                )
                return self.send(
                    200,
                    {
                        "config_path": str(SETTINGS["config_path"]),
                        "config_error": error,
                        "countries": countries,
                        "default_country": cfg.get("general", {}).get("country", ""),
                        "country": country,
                        "groups": [
                            {"id": g, "title": t, "about": a} for g, t, a in GROUPS
                        ],
                        "steps": step_states(country) if country else [],
                        "folders": (
                            [str(p) for p in folders(country)] if country else []
                        ),
                        "busy": JOBS.busy(),
                    },
                )
            if url.path == "/api/config":
                text, cfg, error = read_config()
                return self.send(200, {"text": text, "config": cfg, "error": error})
            if url.path == "/api/jobs":
                with JOBS.lock:
                    return self.send(200, [JOBS.summary(j) for j in JOBS.items])
            if url.path == "/api/job":
                job = JOBS.get(int(q.get("id", 0)))
                if job is None:
                    return self.send(404, {"error": "no such job"})
                since = int(q.get("since", 0))
                with JOBS.lock:
                    start = max(since - job["offset"], 0)
                    return self.send(
                        200,
                        {
                            "job": JOBS.summary(job),
                            "from": job["offset"] + start,
                            "lines": job["lines"][start : start + 4000],
                        },
                    )
            if url.path == "/api/results":
                return self.send(200, results(q.get("country", "")))
            if url.path == "/api/logs":
                out_dir, _ = folders(q.get("country", ""))
                logs = out_dir / "logs"
                files = (
                    sorted(
                        logs.glob("*.log"),
                        key=lambda p: p.stat().st_mtime,
                        reverse=True,
                    )
                    if logs.exists()
                    else []
                )
                return self.send(
                    200,
                    [
                        {
                            "name": p.name,
                            "size": p.stat().st_size,
                            "mtime": p.stat().st_mtime,
                        }
                        for p in files[:300]
                    ],
                )
            if url.path == "/api/log":
                out_dir, _ = folders(q.get("country", ""))
                path = out_dir / "logs" / Path(q.get("name", "")).name
                if not path.is_file():
                    return self.send(404, {"error": "no such log"})
                data = path.read_bytes()[-600_000:]
                return self.send(200, {"text": data.decode("utf-8", "replace")})
            if url.path == "/api/file":
                path = (ROOT / q.get("path", "")).resolve()
                if (
                    ROOT / "reports"
                ).resolve() not in path.parents or not path.is_file():
                    return self.send(404, {"error": "not found"})
                return self.send(
                    200,
                    path.read_bytes(),
                    TYPES.get(path.suffix, "application/octet-stream"),
                )
            if url.path == "/api/columns":
                path = q.get("path", "")
                if not Path(path).is_file():
                    return self.send(200, {"columns": []})
                info = pyogrio.read_info(path)
                return self.send(
                    200,
                    {
                        "columns": [
                            {"name": n, "type": str(t)}
                            for n, t in zip(info["fields"], info["dtypes"])
                        ],
                        "features": int(info.get("features", -1)),
                    },
                )
            return self.send(404, {"error": "not found"})
        except Exception as exc:  # noqa: BLE001 - reported to the page
            return self.send(500, {"error": f"{type(exc).__name__}: {exc}"})

    def do_POST(self):
        if not self.allowed():
            return self.send(403, {"error": "forbidden"})
        url = urlparse(self.path)
        try:
            data = self.body()
            if url.path == "/api/run":
                step = next((s for s in STEPS if s[0] == data.get("module")), None)
                if step is None:
                    return self.send(400, {"error": "unknown step"})
                variables = {
                    k: v
                    for k, v in (data.get("variables") or {}).items()
                    if re.fullmatch(r"[A-Z_]+", k)
                }
                for line in str(data.get("extra", "")).splitlines():
                    if "=" in line:
                        k, v = line.split("=", 1)
                        if re.fullmatch(r"[A-Z_]+", k.strip()):
                            variables[k.strip()] = v.strip()
                job_id = JOBS.add(
                    step[0], step[2], str(data.get("country", "")), variables
                )
                return self.send(200, {"id": job_id})
            if url.path == "/api/stop":
                return self.send(200, {"ok": JOBS.stop(int(data.get("id", 0)))})
            if url.path == "/api/clear":
                JOBS.clear()
                return self.send(200, {"ok": True})
            if url.path == "/api/config/raw":
                text = str(data.get("text", ""))
                tomllib.loads(text)
                SETTINGS["config_path"].write_text(text, encoding="utf-8")
                return self.send(200, {"ok": True})
            if url.path == "/api/config/set":
                for change in data.get("changes", []):
                    gis.Config.set_toml_value(
                        SETTINGS["config_path"],
                        change["table"],
                        change["key"],
                        change["value"],
                    )
                return self.send(200, {"ok": True})
            if url.path == "/api/paths":
                return self.send(
                    200,
                    {
                        p: Path(p).exists()
                        for p in data.get("paths", [])
                        if isinstance(p, str)
                    },
                )
            if url.path == "/api/browse":
                paths, error = browse(
                    data.get("kind", "file"), bool(data.get("multiple"))
                )
                return self.send(200, {"paths": paths, "error": error})
            if url.path == "/api/open":
                target = data.get("target", "")
                out_dir, work_dir = folders(data.get("country", ""))
                path = {
                    "output": out_dir,
                    "work": work_dir,
                    "reports": ROOT / "reports" / data.get("country", ""),
                    "config": SETTINGS["config_path"].parent,
                }.get(target)
                return self.send(200, {"ok": bool(path) and open_folder(path)})
            if url.path == "/api/quit":
                if JOBS.busy() and not data.get("force"):
                    return self.send(409, {"error": "a step is running"})
                # answer first; the server stops a moment later, so that the answer is not cut off
                self.send(200, {"ok": True})
                self.wfile.flush()
                threading.Timer(0.5, self.server.shutdown).start()
                return None
            return self.send(404, {"error": "not found"})
        except tomllib.TOMLDecodeError as exc:
            return self.send(400, {"error": f"TOML error: {exc}"})
        except Exception as exc:  # noqa: BLE001 - reported to the page
            return self.send(500, {"error": f"{type(exc).__name__}: {exc}"})


def open_window(url):
    # an app window of Edge or Chrome looks like a desktop application; the default browser otherwise
    list_exe = [
        shutil.which("msedge"),
        shutil.which("chrome"),
        shutil.which("google-chrome"),
        shutil.which("chromium"),
        r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
        r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
        r"C:\Program Files\Google\Chrome\Application\chrome.exe",
    ]
    for exe in [e for e in list_exe if e and Path(e).exists()]:
        try:
            subprocess.Popen([exe, f"--app={url}", "--window-size=1480,940"])
            return
        except OSError:
            continue
    webbrowser.open(url)


if __name__ == "__main__":
    run = gis.Run.env_overrides(
        CONFIG_PATH=CONFIG_PATH, PORT=PORT, OPEN_WINDOW=OPEN_WINDOW
    )
    config_path = Path(run["CONFIG_PATH"])
    SETTINGS = {
        "config_path": config_path if config_path.is_absolute() else ROOT / config_path
    }
    TOKEN = secrets.token_urlsafe(24)
    CACHE = {}
    BROWSE_LOCK = threading.Lock()
    JOBS = Jobs()
    server = ThreadingHTTPServer(("127.0.0.1", run["PORT"]), Handler)
    url = f"http://127.0.0.1:{server.server_address[1]}/"
    print(
        f"Defenses GUI at {url} (configuration {SETTINGS['config_path']}). Close with Ctrl+C or Quit."
    )
    if run["OPEN_WINDOW"]:
        open_window(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    for job in list(JOBS.items):
        if job["status"] in ("running", "queued"):
            JOBS.stop(job["id"])
    server.server_close()
