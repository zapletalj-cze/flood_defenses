# Author: Jakub Zapletal
# Date: 2026-11-12
# Checks of the Defenses GUI server (gui/app.py) on a small synthetic setup: security (token, host, paths), state and
# script defaults, results, logs, configuration editing (form and raw), path and column checks, runs with the queue,
# console lines, progress, failure, stop, cancel, clear, file dialog and folder opening, quit. No browser needed.
# Run: pixi run python -m tests.test_gui  (DEFENSES_TEST_DIR sets the working folder, default: temporary)

import ast
import json
import os
import re
import subprocess
import sys
import tempfile
import threading
import time
import tomllib
import types
import urllib.error
import urllib.request
from pathlib import Path

import geopandas as gpd
import pandas as pd
from shapely.geometry import Point, box

from gui.app import STEPS, Jobs

EPSG = 32633  # synthetic test data only
repo = Path(__file__).resolve().parents[1]
root = Path(
    os.environ.get("DEFENSES_TEST_DIR") or tempfile.mkdtemp(prefix="defenses_gui_")
)
root.mkdir(parents=True, exist_ok=True)
print(f"Test folder: {root}")
out = root / "output" / "TG"
work = root / "work" / "TG"
(out / "logs").mkdir(parents=True, exist_ok=True)
(work / "units" / "000_000").mkdir(parents=True, exist_ok=True)
(work / "units" / "000_001").mkdir(parents=True, exist_ok=True)

# --- Synthetic outputs ---
gpd.GeoDataFrame(
    {
        "area_m2": [40000.0, 60000.0],
        "protection_level_m": [2.0, 3.0],
        "truncated": [False, True],
    },
    geometry=[box(0, 0, 200, 200), box(300, 0, 600, 200)],
    crs=EPSG,
).to_file(out / "protected_areas.gpkg")
gpd.GeoDataFrame(
    {
        "leak_id": ["TG_L1", "TG_L2", "TG_L3"],
        "cause": ["line_end", "line_end", "river_crossing"],
        "x": [10.0, 20.0, 30.0],
        "y": [1.0, 2.0, 3.0],
        "levee_length_m": [400.0, 100.0, 250.0],
        "lost_margin_median_m": [1.5, 0.5, 2.0],
        "line_ids": ["1", "2", "3,4"],
    },
    geometry=[Point(10, 1), Point(20, 2), Point(30, 3)],
    crs=EPSG,
).to_file(out / "leaks.gpkg")
gpd.GeoDataFrame(
    {
        "no_ridge": [True, False, False],
        "ridge_share": [0.0, 1.0, 0.9],
        "shift_median_m": [0.0, 6.0, 0.0],
    },
    geometry=[box(0, 0, 1, 1).boundary] * 3,
    crs=EPSG,
).to_file(out / "levees_crest.gpkg")
(out / "performance.json").write_text(
    json.dumps({"seconds_per_domain_cell": 2.4e-6, "country_seconds": 5400})
)
pd.DataFrame({"unit_id": ["000_000", "000_001"]}).to_csv(out / "units.csv", index=False)
(work / "units" / "000_000" / "done.json").write_text("{}")
(out / "logs" / "protected_areas_20261112_120000.log").write_text(
    "2026-11-12 12:00:00,000 INFO protected_areas: start\n"
    "2026-11-12 12:00:01,000 WARNING protected_areas: something to look at\n"
    "2026-11-12 12:00:02,000 ERROR protected_areas: something wrong\n"
)
gpd.GeoDataFrame(
    {"name": ["main"], "uparea": [1200.5]},
    geometry=[box(0, 0, 10, 10).boundary],
    crs=EPSG,
).to_file(root / "rivers.gpkg")
config_text = f"""# GUI test configuration
[general]
country = "TG"
work_dir = {json.dumps(str(root / 'work'))}
output_dir = {json.dumps(str(root / 'output'))}
n_workers = 1                   # parallel unit processes
max_ram_gb_per_worker = 2.0

[countries.TG]
dtm = []
levees = []
rivers = [{json.dumps(str(root / 'rivers.gpkg'))}]
river_upstream_area_field = ""  # the NAME of the upstream area column

[countries.TG.rp]
max_damage_csv = ""
"""
config = root / "defenses.toml"
config.write_text(config_text)

# --- Server ---
env = dict(
    os.environ,
    DEFENSES_CONFIG_PATH=str(config),
    DEFENSES_PORT="0",
    DEFENSES_OPEN_WINDOW="0",
    PYTHONUNBUFFERED="1",
)
server = subprocess.Popen(
    [sys.executable, "-m", "gui.app"],
    cwd=repo,
    env=env,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    text=True,
)
first = server.stdout.readline()
m = re.search(r"http://127\.0\.0\.1:(\d+)/", first)
if not m:
    print(first + server.stdout.read())
    sys.exit(1)
port = int(m.group(1))
threading.Thread(target=lambda: [None for _ in server.stdout], daemon=True).start()
base = f"http://127.0.0.1:{port}"


def request(path, body=None, token=True, host=None):
    # (status, parsed JSON or bytes)
    req = urllib.request.Request(
        base + path,
        data=json.dumps(body).encode() if body is not None else None,
        method="POST" if body is not None else "GET",
    )
    req.add_header("Content-Type", "application/json")
    if token:
        req.add_header("X-Defenses-Token", TOKEN)
    if host:
        req.add_header("Host", host)
    try:
        with urllib.request.urlopen(req, timeout=30) as res:
            data = res.read()
            status = res.status
    except urllib.error.HTTPError as err:
        data, status = err.read(), err.code
    try:
        return status, json.loads(data)
    except ValueError:
        return status, data


def wait_job(job_id, states, seconds=120):
    end = time.time() + seconds
    while time.time() < end:
        _, data = request(f"/api/job?id={job_id}&since=0")
        if data["job"]["status"] in states:
            return data
        time.sleep(0.3)
    return data


list_res = []
TOKEN = ""
status, page = request("/", token=False)
TOKEN = (
    re.search(rb'DEFENSES_TOKEN = "([^"]+)"', page).group(1).decode()
    if status == 200
    else ""
)
list_res.append(
    ("page served with a session token", status == 200 and len(TOKEN) >= 20)
)

# --- Security ---
s_logo, logo = request("/static/logo.svg", token=False)
s_trav, _ = request("/static/../app.py", token=False)
s_trav2, _ = request("/static/%2e%2e/app.py", token=False)
s_none, _ = request("/api/state", token=False)
s_host, _ = request("/api/state", host="evil.example.com")
s_404, _ = request("/api/nothing")
list_res.append(
    (
        f"static files served, paths outside gui/static refused ({s_trav}, {s_trav2}), API without token {s_none}, "
        f"foreign Host {s_host}, unknown endpoint {s_404}",
        s_logo == 200
        and b"<svg" in logo
        and s_trav == 404
        and s_trav2 == 404
        and s_none == 403
        and s_host == 403
        and s_404 == 404,
    )
)

# --- State and the defaults of the run variables read from the scripts ---
_, state = request("/api/state")
pa = next(s for s in state["steps"] if s["module"] == "processing.protected_areas")
tree = ast.parse(
    (repo / "processing" / "protected_areas.py").read_text(encoding="utf-8")
)
force_default = next(
    ast.literal_eval(n.value)
    for n in tree.body
    if isinstance(n, ast.Assign) and n.targets[0].id == "FORCE"
)
rp = next(s for s in state["steps"] if s["module"] == "processing.rp_index")
list_res.append(
    (
        f"state: country TG, {len(state['steps'])} steps in {len(state['groups'])} groups, FORCE of "
        f"protected_areas = {force_default} as in the script, LAYER on rp_index, protected areas done",
        state["countries"] == ["TG"]
        and state["country"] == "TG"
        and len(state["steps"]) == len(STEPS)
        and {o["name"]: o["default"] for o in pa["options"]}["FORCE"] == force_default
        and any(o["name"] == "LAYER" and o["default"] == "A" for o in rp["options"])
        and pa["status"] == "done"
        and pa["log"] == "protected_areas_20261112_120000.log",
    )
)

# --- Results, logs, files ---
_, res = request("/api/results?country=TG")
causes = {c["cause"]: (c["places"], c["km"]) for c in res["leaks"]["by_cause"]}
list_res.append(
    (
        f"results: 2 protected areas {res['protected']['area_km2']} km2, leaks {causes}, units "
        f"{res['units']}, levees {res['levees']}, speed",
        res["protected"]["count"] == 2
        and abs(res["protected"]["area_km2"] - 0.1) < 1e-9
        and causes == {"line_end": (2, 0.5), "river_crossing": (1, 0.25)}
        and res["leaks"]["top"][0]["leak_id"] == "TG_L1"
        and res["units"] == {"total": 2, "done": 1}
        and res["levees"] == {"lines": 3, "no_ridge": 1, "moved": 1}
        and res["performance"]["country_seconds"] == 5400
        and any(f["name"] == "leaks.gpkg" for f in res["files"]),
    )
)
_, logs = request("/api/logs?country=TG")
_, log = request("/api/log?country=TG&name=protected_areas_20261112_120000.log")
s_logtrav, _ = request("/api/log?country=TG&name=../../defenses.toml")
s_file, _ = request("/api/file?path=../defenses.toml")
s_file2, _ = request("/api/file?path=config/defenses.toml")
list_res.append(
    (
        f"logs listed and read; log names and files outside the logs and reports folders refused "
        f"({s_logtrav}, {s_file}, {s_file2})",
        [entry["name"] for entry in logs] == ["protected_areas_20261112_120000.log"]
        and "WARNING" in log["text"]
        and s_logtrav == 404
        and s_file == 404
        and s_file2 == 404,
    )
)

# --- Configuration: paths, columns, form changes, raw text ---
_, found = request(
    "/api/paths", {"paths": [str(root / "rivers.gpkg"), str(root / "missing.tif")]}
)
_, cols = request(
    f"/api/columns?path={urllib.request.quote(str(root / 'rivers.gpkg'))}"
)
_, cols_missing = request("/api/columns?path=nothing.gpkg")
list_res.append(
    (
        f"path check (found {list(found.values())}) and columns of the rivers ({[c['name'] for c in cols['columns']]})",
        list(found.values()) == [True, False]
        and {"name", "uparea"} <= {c["name"] for c in cols["columns"]}
        and next(c for c in cols["columns"] if c["name"] == "uparea")[
            "type"
        ].startswith("float")
        and cols_missing["columns"] == [],
    )
)
s_set, _ = request(
    "/api/config/set",
    {
        "changes": [
            {"table": "general", "key": "n_workers", "value": 3},
            {
                "table": "countries.TG",
                "key": "river_upstream_area_field",
                "value": "uparea",
            },
            {
                "table": "countries.TG",
                "key": "levees",
                "value": ["D:/a b/levees.gpkg", "D:/x#1.gpkg"],
            },
            {"table": "countries.NEW", "key": "dtm", "value": []},
        ]
    },
)
text = config.read_text()
cfg = tomllib.loads(text)
list_res.append(
    (
        "form changes saved in place: values set, comments kept, a new country table added",
        s_set == 200
        and cfg["general"]["n_workers"] == 3
        and cfg["countries"]["TG"]["river_upstream_area_field"] == "uparea"
        and cfg["countries"]["TG"]["levees"] == ["D:/a b/levees.gpkg", "D:/x#1.gpkg"]
        and cfg["countries"]["NEW"] == {"dtm": []}
        and "n_workers = 3                   # parallel unit processes" in text
        and '"uparea"  # the NAME of the upstream area column' in text
        and text.startswith("# GUI test configuration"),
    )
)
s_bad, bad = request("/api/config/raw", {"text": "[general\ncountry = 1"})
unchanged = config.read_text() == text
s_raw, _ = request(
    "/api/config/raw",
    {
        "text": text.replace(
            "max_ram_gb_per_worker = 2.0", "max_ram_gb_per_worker = 4.0"
        )
    },
)
_, conf = request("/api/config")
list_res.append(
    (
        f"raw text: invalid TOML refused ({s_bad}: {str(bad.get('error', ''))[:40]}…) and the file kept; valid text saved",
        s_bad == 400
        and "TOML" in bad.get("error", "")
        and unchanged
        and s_raw == 200
        and conf["config"]["general"]["max_ram_gb_per_worker"] == 4.0
        and conf["error"] == "",
    )
)

# --- Runs: validation of the request, failure, variables ---
s_unknown, _ = request("/api/run", {"module": "os.system", "country": "TG"})
_, run1 = request(
    "/api/run",
    {
        "module": "processing.check_inputs",
        "country": "TG",
        "variables": {"FORCE": "1", "bad name": "x", "lower": "1"},
        "extra": "UNITS=000_000\nlowercase=1\nnot a variable",
    },
)
job1 = wait_job(run1["id"], {"done", "failed"})
list_res.append(
    (
        f"unknown modules refused ({s_unknown}); a step runs with only valid run variables "
        f"({job1['job']['variables']}), its console lines arrive, a failing check ends as failed "
        f"(exit {job1['job']['returncode']})",
        s_unknown == 400
        and job1["job"]["variables"] == {"FORCE": "1", "UNITS": "000_000"}
        and job1["job"]["status"] == "failed"
        and job1["job"]["returncode"] != 0
        and len(job1["lines"]) > 0
        and job1["job"]["n_lines"] == len(job1["lines"]),
    )
)
_, part = request(f"/api/job?id={run1['id']}&since={max(len(job1['lines']) - 2, 0)}")
list_res.append(
    (
        "console lines can be fetched from an offset",
        part["lines"] == job1["lines"][-2:] and part["from"] == len(job1["lines"]) - 2,
    )
)

# --- Queue, stop, cancel, clear ---
_, long_job = request("/api/run", {"module": "tests.test_kernels", "country": "TG"})
_, queued = request("/api/run", {"module": "processing.check_inputs", "country": "TG"})
_, cancel = request("/api/run", {"module": "processing.check_inputs", "country": "TG"})
wait_job(long_job["id"], {"running"}, 30)
_, jobs = request("/api/jobs")
states = {j["id"]: j["status"] for j in jobs}
_, stop_q = request("/api/stop", {"id": cancel["id"]})
time.sleep(1.0)
_, stop_r = request("/api/stop", {"id": long_job["id"]})
stopped = wait_job(long_job["id"], {"stopped", "done", "failed"}, 60)
second = wait_job(queued["id"], {"done", "failed"}, 120)
_, cancelled = request(f"/api/job?id={cancel['id']}&since=0")
_, stop_none = request("/api/stop", {"id": 999})
leftover = (
    subprocess.run(
        ["pgrep", "-f", "tests.test_kernels"], capture_output=True, text=True
    ).stdout.split()
    if os.name != "nt"
    else []
)
list_res.append(
    (
        f"queue runs one step at a time ({states[long_job['id']]}, {states[queued['id']]}); stop ends the running "
        f"step ({stopped['job']['status']}, exit {stopped['job']['returncode']}, processes left {len(leftover)}), "
        f"the next starts ({second['job']['status']}); a queued step is cancelled ({cancelled['job']['status']})",
        states[long_job["id"]] == "running"
        and states[queued["id"]] == "queued"
        and stop_q["ok"]
        and stop_r["ok"]
        and stopped["job"]["status"] == "stopped"
        and not leftover
        and second["job"]["status"] in ("done", "failed")
        and second["job"]["started"] >= stopped["job"]["ended"] - 0.5
        and cancelled["job"]["status"] == "cancelled"
        and cancelled["job"]["started"] is None
        and stop_none["ok"] is False,
    )
)
request("/api/clear", {})
_, after_clear = request("/api/jobs")
list_res.append(("clear removes finished runs", after_clear == []))

# --- Progress and ETA from the console lines (unit counter of protected_areas.py and others) ---
fake = types.SimpleNamespace(lock=threading.Condition())
job = {
    "lines": [],
    "offset": 0,
    "progress": None,
    "eta": "",
    "warnings": 0,
    "errors": 0,
}
for line in [
    "2026 INFO protected_areas: Units: 7 to process",
    "2026 INFO protected_areas: [3/7] unit 001_001: 9 protected areas | elapsed 0:00:15, ETA 0:01:02",
    "2026 WARNING gis: placeholder",
    "Traceback (most recent call last):",
]:
    Jobs.append(fake, job, line)
list_res.append(
    (
        f"progress {job['progress']}, ETA {job['eta']}, warnings {job['warnings']} and errors {job['errors']} "
        "read from the console",
        job["progress"] == [3, 7]
        and job["eta"] == "0:01:02"
        and job["warnings"] == 1
        and job["errors"] == 1,
    )
)

# --- Folders and file dialog: answers, never a server error ---
s_open, opened = request("/api/open", {"target": "output", "country": "TG"})
s_open2, opened2 = request("/api/open", {"target": "nothing", "country": "TG"})
checks = [s_open == 200, s_open2 == 200 and opened2["ok"] is False]
if os.name != "nt" and not os.environ.get("DISPLAY"):
    s_browse, browsed = request("/api/browse", {"kind": "file", "multiple": True})
    checks.append(s_browse == 200 and browsed["paths"] == [] and browsed["error"] != "")
list_res.append(
    (
        "open folder and file dialog answer without a server error (no display: a message instead of the dialog)",
        all(checks),
    )
)

# --- Quit: refused while a step runs, unless forced ---
_, busy_job = request("/api/run", {"module": "tests.test_kernels", "country": "TG"})
wait_job(busy_job["id"], {"running"}, 30)
s_quit, _ = request("/api/quit", {})
s_force, _ = request("/api/quit", {"force": True})
try:
    server.wait(timeout=30)
except subprocess.TimeoutExpired:
    server.kill()
list_res.append(
    (
        f"quit refused while a step runs ({s_quit}), forced quit stops the server and the step (exit {server.returncode})",
        s_quit == 409 and s_force == 200 and server.returncode == 0,
    )
)

for label, ok in list_res:
    print(f"{'PASS' if ok else 'FAIL'}  {label}")
sys.exit(0 if all(ok for _, ok in list_res) else 1)
