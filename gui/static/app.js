// Author: Jakub Zapletal
// Date: 2026-11-12
// Page logic of the Defenses GUI: pipeline with run options, console with queue and progress, results,
// configuration (form and raw TOML) and logs. Talks to gui/app.py through /api with the session token.
"use strict";

const TOKEN = window.DEFENSES_TOKEN;
const ICONS = {
  flow: '<path d="M4 6h6M4 12h10M4 18h16"/><circle cx="14" cy="6" r="2"/><circle cx="18" cy="12" r="2"/>',
  chart: '<path d="M4 20V10M10 20V4M16 20v-7M22 20H2"/>',
  sliders: '<path d="M4 6h10M18 6h2M4 12h4M12 12h8M4 18h12M20 18h0"/><circle cx="16" cy="6" r="2"/><circle cx="10" cy="12" r="2"/><circle cx="18" cy="18" r="2"/>',
  doc: '<path d="M7 3h7l5 5v13H7z"/><path d="M14 3v5h5M10 13h6M10 17h6"/>',
  power: '<path d="M12 3v8"/><path d="M6.3 7.5a8 8 0 1 0 11.4 0"/>',
  stop: '<rect x="6" y="6" width="12" height="12" rx="2"/>',
  chevron: '<path d="M6 9l6 6 6-6"/>',
  play: '<path d="M8 5.5v13l10-6.5z"/>',
  folder: '<path d="M3 7a2 2 0 0 1 2-2h4l2 2h8a2 2 0 0 1 2 2v8a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2z"/>',
  refresh: '<path d="M20 11a8 8 0 1 0-2.3 5.7"/><path d="M20 4v7h-7"/>',
  check: '<path d="M5 12.5l4.5 4.5L19 7.5"/>',
  alert: '<path d="M12 4l9 16H3z"/><path d="M12 10v4M12 17.5v.5"/>',
  plus: '<path d="M12 5v14M5 12h14"/>',
  trash: '<path d="M5 7h14M10 7V4h4v3M7 7l1 13h8l1-13"/>',
  search: '<circle cx="11" cy="11" r="6"/><path d="M20 20l-4.5-4.5"/>',
  layers: '<path d="M12 4l9 5-9 5-9-5z"/><path d="M3 14l9 5 9-5"/>',
  shield: '<path d="M12 3l8 3v6c0 5-3.5 8-8 9-4.5-1-8-4-8-9V6z"/>',
  drop: '<path d="M12 3c4 5 6 8 6 11a6 6 0 0 1-12 0c0-3 2-6 6-11z"/>',
  grid: '<rect x="4" y="4" width="7" height="7" rx="1.5"/><rect x="13" y="4" width="7" height="7" rx="1.5"/><rect x="4" y="13" width="7" height="7" rx="1.5"/><rect x="13" y="13" width="7" height="7" rx="1.5"/>',
  gauge: '<path d="M4 16a8 8 0 1 1 16 0"/><path d="M12 16l4-5"/>',
  copy: '<rect x="8" y="8" width="12" height="12" rx="2"/><path d="M16 8V5a1 1 0 0 0-1-1H5a1 1 0 0 0-1 1v10a1 1 0 0 0 1 1h3"/>',
  save: '<path d="M5 4h11l3 3v13H5z"/><path d="M8 4v5h7V4M8 20v-6h8v6"/>',
  wave: '<path d="M3 9c2.5-2 4.5-2 6 0s3.5 2 6 0 4.5-2 6 0M3 15c2.5-2 4.5-2 6 0s3.5 2 6 0 4.5-2 6 0"/>',
};
const OPTION_TEXT = {
  FORCE: ["switch", "Recompute (FORCE)"],
  RUN_BENCHMARK: ["switch", "Read benchmark"],
  HINTS: ["switch", "Missing-barrier hints"],
  UNITS: ["text", "Units", "all, or 012_034,012_035"],
  RUN_ID: ["text", "Run id", "base"],
  DATE_FROM: ["text", "From", "YYYY-MM-DD"],
  DATE_TO: ["text", "To", "YYYY-MM-DD"],
  LAYER: ["choice", "File", ["A", "B"]],
};
// help texts, only where a label alone is not clear
const HELP = {
  FORCE: "Recompute what this step has already finished. Off: units with a finished checkpoint are skipped, so a rerun only does what is missing. Turn it on after changing levees, rivers or parameters.",
  UNITS: "Unit ids separated by commas, e.g. 012_034,012_035; empty = all units. The ids are in units.csv and units_core.gpkg (Build units). Run two units first to check the time and the results.",
  LAYER: "A: protected areas with all rivers as water sources at once (protected_areas.gpkg). B: per river group, only the rivers of the group are water sources (protected_areas_by_group.gpkg); the outputs get the suffix _b.",
  RUN_ID: "Name of the folder of this run's plots and results. A new name keeps an earlier run for comparison.",
  HINTS: "Also look for places that stay dry in flood scenes behind no mapped levee (missing barriers). Reads the whole unit window per scene, so it is slower. File A only.",
  "processing.protected_areas": "Finds the areas behind levees that river water reaches only at a higher level than without the levees. Writes file A (all rivers at once), file B (per river group) and leaks.gpkg: where levees protect nothing, and why.",
  "processing.transport_embankments": "Optional. Needs a file of road and railway lines in transport_lines (Configuration, Raw TOML). The lines only show where to look; the DTM decides what is an embankment and where its openings are. Not yet part of the barriers.",
  "processing.rp_index": "Estimates the return period of protection of every protected area from the value exposed behind the levees and the size of the river. Needs the damage tables (Configuration) and the GHSL building volume.",
  "processing.rp_bounds": "Narrows the return period with observations: GFM flood scenes held or flooded, at the discharge of GloFAS or EFAS that day. Needs GFM evidence and the prepared discharge.",
  workers: "Units computed at the same time. Each needs about the RAM of one unit window (Build units log); workers × that must fit in the memory. More workers than CPU cores do not help.",
  upstream: "The NAME of the column of the river file with the upstream area in km² (e.g. uparea in MERIT), not a value. Needed for file B (river groups) and the return period index.",
  river_guides: "Optional lines of the real river channels (e.g. waterways). The river network (e.g. MERIT, about 90 m) is moved onto them and then into the lowest DTM cells. Without them the rivers stay where they are.",
  tile_b: "Areas of different river groups overlap: one place can be protected against the main river and against a tributary. This area counts every group; the report gives the union, every place once.",
  tile_leaks: "Places where water gets behind levees that protect nothing (leaks.gpkg): an open river crossing, a low crest, a free line end, a stream behind the levee, or terrain behind higher than the crest (no levee there, or its crest underestimated).",
  tile_speed: "Computing time per cell of the processing domain, measured on the finished units, and what it means for the whole country with the configured workers.",
};
const help = (key) => (HELP[key] ? `<button type="button" class="help" data-help="${esc(key)}" aria-label="Help">?</button>` : "");
function showHelp(btn) {
  let pop = $("#helpPop");
  if (!pop) {
    pop = document.createElement("div");
    pop.id = "helpPop";
    pop.className = "help-pop";
    document.body.appendChild(pop);
  }
  if (pop.dataset.key === btn.dataset.help && pop.classList.contains("show")) return hideHelp();
  pop.textContent = HELP[btn.dataset.help];
  pop.dataset.key = btn.dataset.help;
  pop.classList.add("show");
  S.helpOpened = Date.now();
  const r = btn.getBoundingClientRect();
  const w = Math.min(340, window.innerWidth - 32);
  pop.style.width = `${w}px`;
  pop.style.left = `${Math.max(16, Math.min(r.left + r.width / 2 - w / 2, window.innerWidth - w - 16))}px`;
  const below = r.bottom + 10 + pop.offsetHeight < window.innerHeight;
  pop.style.top = `${below ? r.bottom + 8 : r.top - pop.offsetHeight - 8}px`;
}
function hideHelp() {
  const pop = $("#helpPop");
  if (pop) {
    pop.classList.remove("show");
    pop.dataset.key = "";
  }
}
document.addEventListener("click", (e) => {
  const btn = e.target.closest(".help");
  if (btn) {
    e.preventDefault();
    e.stopPropagation();
    showHelp(btn);
  } else if (!e.target.closest("#helpPop")) hideHelp();
}, true);
document.addEventListener("keydown", (e) => { if (e.key === "Escape") hideHelp(); });
// scrolling closes the bubble, except the scroll that the click itself causes
window.addEventListener("scroll", () => { if (Date.now() - (S.helpOpened || 0) > 400) hideHelp(); }, true);

const PAGES = {
  pipeline: ["Pipeline", "Run the steps in order. Every card shows what the step has written."],
  results: ["Results", "Protected areas, leaks of the levees, units and performance of the country."],
  config: ["Configuration", "Inputs and folders of the runs. Every other value comes from config/advanced.toml."],
  logs: ["Logs", "Run logs of the scripts, newest first."],
};
const S = {
  state: null, country: "", page: "pipeline", jobs: [], job: null, from: 0, follow: true,
  opts: {}, results: null, cfg: null, cfgText: "", dirty: {}, cfgTab: "form", cfgCountry: "",
  logs: [], logName: null, logText: "", logFilter: "all", logSearch: "", columns: [],
};

// --- Helpers ---
const $ = (sel, root = document) => root.querySelector(sel);
const esc = (s) => String(s ?? "").replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[c]);
const icon = (name) => `<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round">${ICONS[name] || ""}</svg>`;
const ico = (name) => `<span class="ico">${icon(name)}</span>`;
function hydrate(root = document) {
  root.querySelectorAll(".ico[data-icon]").forEach((el) => { if (!el.innerHTML) el.innerHTML = icon(el.dataset.icon); });
}
async function api(path, body) {
  const res = await fetch(path, {
    method: body ? "POST" : "GET",
    headers: { "X-Defenses-Token": TOKEN, "Content-Type": "application/json" },
    body: body ? JSON.stringify(body) : undefined,
  });
  const data = await res.json().catch(() => ({}));
  if (!res.ok) throw new Error(data.error || res.statusText);
  return data;
}
function ago(ts) {
  if (!ts) return "never";
  const s = Date.now() / 1000 - ts;
  if (s < 60) return "just now";
  if (s < 3600) return `${Math.round(s / 60)} min ago`;
  if (s < 86400) return `${Math.round(s / 3600)} h ago`;
  return new Date(ts * 1000).toLocaleDateString(undefined, { day: "numeric", month: "short", year: "numeric" });
}
function duration(sec) {
  if (sec == null) return "";
  sec = Math.max(0, Math.round(sec));
  const h = Math.floor(sec / 3600), m = Math.floor((sec % 3600) / 60), s = sec % 60;
  return h ? `${h} h ${m} min` : m ? `${m} min ${s} s` : `${s} s`;
}
const bytes = (n) => (n > 1e9 ? `${(n / 1e9).toFixed(1)} GB` : n > 1e6 ? `${(n / 1e6).toFixed(1)} MB` : `${Math.max(1, Math.round(n / 1e3))} kB`);
const num = (v, d = 0) => (v == null || Number.isNaN(v) ? "–" : Number(v).toLocaleString(undefined, { maximumFractionDigits: d, minimumFractionDigits: d }));
const STATUS = { idle: "Not run", done: "Done", running: "Running", stopping: "Stopping", queued: "Queued", failed: "Failed", stopped: "Stopped", cancelled: "Cancelled" };
const chip = (status) => `<span class="chip ${esc(status)}">${esc(STATUS[status] || status)}</span>`;
let toastTimer;
function toast(msg, err = false) {
  const el = $("#toast");
  el.textContent = msg;
  el.className = `toast show${err ? " err" : ""}`;
  clearTimeout(toastTimer);
  toastTimer = setTimeout(() => (el.className = "toast"), 3600);
}

// --- Boot and navigation ---
async function boot() {
  hydrate();
  $("#country").addEventListener("change", async (e) => {
    S.country = e.target.value;
    S.results = null;
    S.logName = null;
    await loadState();
    render();
  });
  window.addEventListener("hashchange", route);
  $("#page").addEventListener("focusout", () => setTimeout(() => {
    if (S.pendingRender && !$("#page").contains(document.activeElement)) {
      S.pendingRender = false;
      if (S.page === "pipeline" || S.page === "results") render();
    }
  }, 0));
  $("#consoleHead").addEventListener("click", (e) => { if (!e.target.closest("button")) toggleConsole(); });
  $("#consoleToggle").addEventListener("click", () => toggleConsole());
  $("#stopBtn").addEventListener("click", stopJob);
  $("#quitBtn").addEventListener("click", quit);
  $("#consoleLog").addEventListener("scroll", (e) => {
    const el = e.target;
    S.follow = el.scrollTop + el.clientHeight >= el.scrollHeight - 30;
  });
  await loadState();
  route();
  pollJobs();
  setInterval(pollJobs, 1200);
}
async function loadState() {
  try {
    S.state = await api(`/api/state?country=${encodeURIComponent(S.country)}`);
  } catch (err) {
    toast(err.message, true);
    return;
  }
  if (!S.country) S.country = S.state.country;
  const sel = $("#country");
  const list = S.state.countries.length ? S.state.countries : [S.country].filter(Boolean);
  sel.innerHTML = list.map((c) => `<option value="${esc(c)}" ${c === S.country ? "selected" : ""}>${esc(c)}${c === S.state.default_country ? "  · default" : ""}</option>`).join("");
  $("#cfgPath").textContent = S.state.config_path;
  const done = S.state.steps.filter((s) => s.status === "done").length;
  $("#countryMeta").innerHTML = `<span class="chip plain">${done} of ${S.state.steps.length} steps done</span>`;
  $("#banner").innerHTML = S.state.config_error
    ? `<div class="banner err">${ico("alert")}<div><b>The configuration cannot be read.</b> ${esc(S.state.config_error)} Fix it under Configuration → Raw TOML.</div></div>`
    : "";
}
function route() {
  const page = (location.hash || "#pipeline").slice(1).split("?")[0];
  S.page = PAGES[page] ? page : "pipeline";
  document.querySelectorAll("#nav a").forEach((a) => a.classList.toggle("active", a.dataset.page === S.page));
  render();
}
function render() {
  const [title, sub] = PAGES[S.page];
  $("#pageTitle").textContent = title;
  $("#pageSub").textContent = sub;
  $("#topActions").innerHTML = "";
  ({ pipeline: renderPipeline, results: renderResults, config: renderConfig, logs: renderLogs })[S.page]();
}

// --- Pipeline ---
function optionValue(step, opt) {
  S.opts[step.module] = S.opts[step.module] || {};
  if (!(opt.name in S.opts[step.module])) {
    const d = opt.default;
    S.opts[step.module][opt.name] = Array.isArray(d) ? d.join(",") : d ?? (OPTION_TEXT[opt.name]?.[0] === "switch" ? false : "");
  }
  return S.opts[step.module][opt.name];
}
function optionsHtml(step) {
  return step.options.map((opt) => {
    const [kind, label, ph] = OPTION_TEXT[opt.name] || ["text", opt.name, ""];
    const value = optionValue(step, opt);
    const id = `${step.module}|${opt.name}`;
    if (kind === "choice") {
      const titles = { A: "A · all rivers", B: "B · per river group" };
      return `<span class="opt-pair"><label class="opt-text">${esc(label)}<div class="select-wrap"><select class="input small" data-opt="${esc(id)}">${ph.map((o) => `<option value="${esc(o)}" ${String(value) === o ? "selected" : ""}>${esc(titles[o] || o)}</option>`).join("")}</select></div></label>${help(opt.name)}</span>`;
    }
    if (kind === "switch") {
      return `<span class="opt-pair"><label class="switch" title="${esc(opt.name)}"><input type="checkbox" data-opt="${esc(id)}" ${value ? "checked" : ""}><span class="track"></span>${esc(label)}</label>${help(opt.name)}</span>`;
    }
    return `<span class="opt-pair"><label class="opt-text">${esc(label)}<input class="input small mono" data-opt="${esc(id)}" value="${esc(value)}" placeholder="${esc(ph)}"></label>${help(opt.name)}</span>`;
  }).join("");
}
function variablesOf(step) {
  const out = {};
  step.options.forEach((opt) => {
    const v = optionValue(step, opt);
    if (typeof v === "boolean") out[opt.name] = v ? "1" : "0";
    else if (String(v).trim()) out[opt.name] = String(v).trim();
  });
  return out;
}
async function runStep(module) {
  const step = S.state.steps.find((s) => s.module === module);
  try {
    const { id } = await api("/api/run", { module, country: S.country, variables: variablesOf(step) });
    S.job = id;
    S.from = 0;
    S.autoFollow = true;
    $("#consoleLog").innerHTML = "";
    openConsole();
    toast(`${step.title} queued for ${S.country}`);
    await pollJobs();
    await loadState();
    if (S.page === "pipeline") renderPipeline();
  } catch (err) {
    toast(err.message, true);
  }
}
async function runCore() {
  // the core steps without the optional transport embankments (they need transport_lines)
  const list = S.state.steps.filter((s) => s.group === "core" && s.module !== "processing.transport_embankments");
  if (!confirm(`Queue the core pipeline for ${S.country}?\n\n${list.map((s, i) => `${i + 1}. ${s.title}`).join("\n")}`)) return;
  for (const s of list) {
    await api("/api/run", { module: s.module, country: S.country, variables: variablesOf(s) });
  }
  S.job = null;
  S.autoFollow = true;
  openConsole();
  toast(`${list.length} steps queued`);
  await pollJobs();
  await loadState();
  renderPipeline();
}
function ring(done, total) {
  const r = 46, c = 2 * Math.PI * r, frac = total ? done / total : 0;
  return `<div class="ring"><svg viewBox="0 0 104 104"><circle cx="52" cy="52" r="${r}" fill="none" stroke="#e9eef5" stroke-width="9"/>
    <circle cx="52" cy="52" r="${r}" fill="none" stroke="url(#rg)" stroke-width="9" stroke-linecap="round" stroke-dasharray="${c}" stroke-dashoffset="${c * (1 - frac)}"/>
    <defs><linearGradient id="rg" x1="0" y1="0" x2="1" y2="1"><stop offset="0" stop-color="#5aa6ee"/><stop offset="1" stop-color="#3fbfae"/></linearGradient></defs></svg>
    <div class="val"><b>${total ? Math.round(100 * frac) : 0}%</b><span>${total ? `${done} / ${total} units` : "no units yet"}</span></div></div>`;
}
async function renderPipeline() {
  $("#topActions").innerHTML = `<button class="btn ghost" id="checkBtn">${ico("check")}Check inputs</button><button class="btn primary" id="coreBtn">${ico("play")}Run core pipeline</button>`;
  $("#checkBtn").onclick = () => runStep("processing.check_inputs");
  $("#coreBtn").onclick = runCore;
  const st = S.state;
  if (!st) return;
  let n = 0;
  const groups = st.groups.map((g) => {
    const steps = st.steps.filter((s) => s.group === g.id);
    const cards = steps.map((s) => {
      n += 1;
      const outs = s.outputs.map((o) => `<span class="out ${o.found ? "found" : ""}">${esc(o.name)}</span>`).join("");
      return `<article class="step ${s.status} ${s.status === "running" ? "is-running" : ""}" data-module="${esc(s.module)}">
        <div class="step-head"><div class="step-num">${n}</div><div><div class="step-title">${esc(s.title)}${HELP[s.module] ? "&nbsp;" + help(s.module) : ""}</div><div class="step-about">${esc(s.about)}</div></div>${chip(s.status)}</div>
        <div class="card-progress" hidden><div class="pbar"><span></span></div><span class="ptext"></span></div>
        ${outs ? `<div class="outputs">${outs}</div>` : ""}
        ${s.options.length ? `<div class="options">${optionsHtml(s)}</div>` : ""}
        <div class="step-foot"><span class="when">${s.last ? `Last written ${ago(s.last)}` : "Nothing written yet"}</span>
          <div class="btns">${s.log ? `<button class="btn ghost small" data-log="${esc(s.log)}">${ico("doc")}Log</button>` : ""}
          <button class="btn soft small" data-run="${esc(s.module)}">${ico("play")}Run</button></div></div>
      </article>`;
    }).join("");
    return `<section class="group"><div class="section-head"><div><h2>${esc(g.title)}</h2><p>${esc(g.about)}</p></div></div><div class="steps">${cards}</div></section>`;
  }).join("");
  const r = S.results || {};
  const units = r.units || { done: 0, total: 0 };
  $("#page").innerHTML = `
    <div class="hero">
      ${ring(units.done, units.total)}
      <div>
        <div class="eyebrow">Country</div>
        <h2>${esc(S.country || "–")}</h2>
        <div class="facts">
          <div class="fact"><b>${r.protected ? num(r.protected.count) : "–"}</b><span>protected areas (A)</span></div>
          <div class="fact"><b>${r.protected ? num(r.protected.area_km2, 1) : "–"}</b><span>km² protected</span></div>
          <div class="fact"><b>${r.by_group ? num(r.by_group.count) : "–"}</b><span>per river group (B)</span></div>
          <div class="fact"><b>${r.leaks ? num(r.leaks.count) : "–"}</b><span>leaking places</span></div>
        </div>
      </div>
      <div class="hero-actions">
        <button class="btn ghost" data-open="output">${ico("folder")}Output folder</button>
        <a class="btn ghost" href="#results">${ico("chart")}See results</a>
      </div>
    </div>${groups}`;
  bindPipeline();
  if (!S.results) {
    try {
      S.results = await api(`/api/results?country=${encodeURIComponent(S.country)}`);
      if (S.page === "pipeline") renderPipeline();
    } catch (err) {
      S.results = {};
    }
  }
}
function bindPipeline() {
  const page = $("#page");
  page.querySelectorAll("[data-run]").forEach((b) => (b.onclick = () => runStep(b.dataset.run)));
  page.querySelectorAll("[data-log]").forEach((b) => (b.onclick = () => { S.logName = b.dataset.log; location.hash = "#logs"; }));
  page.querySelectorAll("[data-open]").forEach((b) => (b.onclick = () => openFolder(b.dataset.open)));
  page.querySelectorAll("[data-opt]").forEach((el) => {
    el.onchange = el.oninput = () => {
      const [module, name] = el.dataset.opt.split("|");
      S.opts[module][name] = el.type === "checkbox" ? el.checked : el.value;
    };
  });
}
async function openFolder(target) {
  try {
    const { ok } = await api("/api/open", { target, country: S.country });
    if (!ok) toast("The folder does not exist yet or cannot be opened here", true);
  } catch (err) {
    toast(err.message, true);
  }
}

// --- Console ---
function toggleConsole(force) {
  const el = $("#console");
  const open = force ?? el.classList.contains("collapsed");
  el.classList.toggle("collapsed", !open);
  document.body.classList.toggle("console-open", open);
}
const openConsole = () => toggleConsole(true);
async function pollJobs() {
  try {
    S.jobs = await api("/api/jobs");
  } catch (err) {
    return;
  }
  const active = S.jobs.find((j) => j.status === "running" || j.status === "stopping");
  const current = S.jobs.find((j) => j.id === S.job);
  if (!current || (["done", "failed", "stopped", "cancelled"].includes(current.status) && active && active.id !== current.id && S.autoFollow !== false)) {
    const pick = active || S.jobs[S.jobs.length - 1];
    if (pick && pick.id !== S.job) {
      S.job = pick.id;
      S.from = 0;
      $("#consoleLog").innerHTML = "";
    }
  }
  renderQueue();
  if (S.job) await pollLines();
  else {
    $("#consoleStatus").style.display = "none";
    $("#stopBtn").disabled = true;
    $("#consoleMeta").textContent = "Runs of this session appear here";
  }
  const finished = S.jobs.filter((j) => ["done", "failed", "stopped"].includes(j.status) && !(S.seen || new Set()).has(j.id));
  S.seen = S.seen || new Set(S.jobs.filter((j) => j.ended).map((j) => j.id));
  for (const j of finished) {
    if (S.seen.has(j.id)) continue;
    S.seen.add(j.id);
    toast(`${j.title}: ${STATUS[j.status].toLowerCase()} in ${duration(j.ended - j.started)}`, j.status !== "done");
    S.results = null;
    await loadState();
    if (S.page === "pipeline" || S.page === "results") renderSoon();
  }
}
// a re-render waits while a field of the page has the focus (typing is not interrupted)
function renderSoon() {
  const el = document.activeElement;
  if (el && $("#page").contains(el) && el.matches("input, select, textarea")) {
    S.pendingRender = true;
    return;
  }
  render();
}
function renderQueue() {
  const html = S.jobs.slice().reverse().map((j) => `<div class="q-item ${j.id === S.job ? "active" : ""}" data-job="${j.id}">
      <b>${esc(j.title)}</b>${chip(j.status)}
      <span class="t">${esc(j.country)} · ${j.started ? `${duration((j.ended || Date.now() / 1000) - j.started)}` : "waiting"}${Object.keys(j.variables).length ? " · " + esc(Object.entries(j.variables).map(([k, v]) => `${k}=${v}`).join(" ")) : ""}</span></div>`).join("");
  $("#queue").innerHTML = html || `<div class="hint" style="padding:10px">No runs in this session yet.</div>`;
  $("#queue").querySelectorAll("[data-job]").forEach((el) => (el.onclick = () => {
    S.job = Number(el.dataset.job);
    S.from = 0;
    $("#consoleLog").innerHTML = "";
    S.autoFollow = false;
    pollJobs();
  }));
  const running = S.jobs.some((j) => j.status === "running");
  $("#consolePulse").classList.toggle("on", running);
  updateProgress();
}
// progress of the running step on its card and in the window title; a moving bar while no counts are known
function progressText(j) {
  const p = j.progress;
  const elapsed = j.started ? duration(Date.now() / 1000 - j.started) : "";
  return p ? `${Math.round((100 * p[0]) / p[1])} % · ${p[0]} / ${p[1]} units${j.eta ? ` · ETA ${j.eta}` : ""}` : `running ${elapsed}`;
}
function setBar(bar, j) {
  const p = j.progress;
  bar.classList.toggle("indet", !p);
  bar.firstElementChild.style.width = p ? `${(100 * p[0]) / p[1]}%` : "";
}
function updateProgress() {
  const run = S.jobs.find((j) => j.status === "running");
  document.querySelectorAll("article.step .card-progress").forEach((box) => {
    const card = box.closest("article.step");
    const mine = run && run.module === card.dataset.module && run.country === S.country;
    box.hidden = !mine;
    card.classList.toggle("is-running", !!mine);
    if (mine) {
      setBar(box.querySelector(".pbar"), run);
      box.querySelector(".ptext").textContent = progressText(run);
    }
  });
  const waiting = S.jobs.filter((j) => j.status === "queued").length;
  document.title = run ? `${run.progress ? Math.round((100 * run.progress[0]) / run.progress[1]) + " % · " : ""}${run.title} · Defenses` : "Defenses";
  S.waiting = waiting;
}
const lineClass = (line) => (/ ERROR |Traceback|Error:/.test(line) ? "ln-err" : / WARNING /.test(line) ? "ln-warn" : "ln-info");
const lineHtml = (line) => `<span class="${lineClass(line)}">${esc(line)}</span>`;
async function pollLines() {
  let data;
  try {
    data = await api(`/api/job?id=${S.job}&since=${S.from}`);
  } catch (err) {
    return;
  }
  const j = data.job;
  if (data.lines.length) {
    const el = $("#consoleLog");
    el.insertAdjacentHTML("beforeend", data.lines.map(lineHtml).join(""));
    while (el.childElementCount > 6000) el.firstElementChild.remove();
    S.from = data.from + data.lines.length;
    if (S.follow) el.scrollTop = el.scrollHeight;
  }
  $("#consoleJob").textContent = `${j.title} · ${j.country}`;
  $("#consoleStatus").style.display = "";
  $("#consoleStatus").className = `chip ${j.status}`;
  $("#consoleStatus").textContent = STATUS[j.status] || j.status;
  const elapsed = j.started ? (j.ended || Date.now() / 1000) - j.started : null;
  $("#consoleMeta").textContent = [elapsed != null ? duration(elapsed) : "", j.warnings ? `${j.warnings} warnings` : "", j.errors ? `${j.errors} errors` : "", S.waiting ? `${S.waiting} waiting in the queue` : ""].filter(Boolean).join(" · ");
  const p = j.progress;
  const frac = j.status === "done" ? 1 : p ? p[0] / p[1] : 0;
  const bar = $("#consoleBar").parentElement;
  bar.classList.toggle("indet", j.status === "running" && !p);
  $("#consoleBar").style.width = j.status === "running" && !p ? "" : `${Math.round(frac * 100)}%`;
  $("#consoleEta").textContent = j.status === "running" ? progressText(j) : j.status === "done" ? "finished" : "";
  $("#stopBtn").disabled = !["running", "queued"].includes(j.status);
}
async function stopJob() {
  if (!S.job) return;
  try {
    await api("/api/stop", { id: S.job });
    toast("Stopping the step and its workers");
  } catch (err) {
    toast(err.message, true);
  }
}
async function quit() {
  const running = S.jobs.some((j) => j.status === "running");
  if (!confirm(running ? "A step is still running. Stop it and quit the GUI?" : "Quit the Defenses GUI?")) return;
  try {
    await api("/api/quit", { force: running });
    document.body.innerHTML = `<div class="empty" style="margin-top:20vh"><img src="/static/logo.svg" width="64" alt=""><h3>Defenses is closed</h3><p>Start it again with <code>pixi run gui</code>.</p></div>`;
  } catch (err) {
    toast(err.message, true);
  }
}

// --- Results ---
async function renderResults() {
  $("#topActions").innerHTML = `<button class="btn ghost" id="refreshBtn">${ico("refresh")}Refresh</button><button class="btn soft" id="reportBtn">${ico("doc")}Make report</button>`;
  $("#refreshBtn").onclick = () => { S.results = null; renderResults(); };
  $("#reportBtn").onclick = () => runStep("processing.data_report");
  if (!S.results) {
    $("#page").innerHTML = `<div class="card empty">${ico("refresh")}<h3>Reading the outputs…</h3></div>`;
    try {
      S.results = await api(`/api/results?country=${encodeURIComponent(S.country)}`);
    } catch (err) {
      $("#page").innerHTML = `<div class="banner err">${ico("alert")}<div>${esc(err.message)}</div></div>`;
      return;
    }
  }
  const r = S.results;
  if (!r.exists) {
    $("#page").innerHTML = `<div class="card empty">${ico("layers")}<h3>No outputs for ${esc(S.country)} yet</h3><p>Run the core pipeline first. Outputs go to <span class="mono">${esc(r.output_dir)}</span>.</p></div>`;
    return;
  }
  const perf = r.performance || {};
  const tiles = [
    ["shield", "Protected areas · A", r.protected ? num(r.protected.count) : "–", r.protected ? `${num(r.protected.area_km2, 1)} km² · median level ${num(r.protected.level_median_m, 2)} m` : "not computed yet", "#e2f5f1"],
    ["wave", "Per river group · B", r.by_group ? num(r.by_group.count) : "–", r.by_group ? `${num(r.by_group.groups)} groups · ${num(r.by_group.area_km2, 1)} km²` : "needs the upstream area column", "#e7f1fc", "tile_b"],
    ["alert", "Leaking places", r.leaks ? num(r.leaks.count) : "–", r.leaks ? `${num(r.leaks.by_cause.reduce((a, b) => a + b.km, 0), 1)} km of levee protect nothing` : "not computed yet", "#fcf3e2", "tile_leaks"],
    ["layers", "Levee lines", r.levees ? num(r.levees.lines) : "–", r.levees ? `${num(r.levees.moved)} moved onto the crest · ${num(r.levees.no_ridge)} without a ridge` : "not computed yet", "#efecfc"],
    ["grid", "Units", r.units ? `${num(r.units.done)}<small>/ ${num(r.units.total)}</small>` : "–", r.units ? (r.units.done === r.units.total ? "all finished" : `${num(r.units.total - r.units.done)} remaining`) : "run Build units", "#e7f1fc"],
    ["gauge", "Speed", perf.seconds_per_domain_cell ? `${num(perf.seconds_per_domain_cell * 1e6, 2)}<small>µs/cell</small>` : "–", perf.country_seconds ? `whole country about ${duration(perf.country_seconds)}` : "after the first units", "#e2f5f1", "tile_speed"],
  ].map(([ic, k, v, d, tone, hk]) => `<div class="tile" style="--tone:${tone}"><div class="k">${esc(k)}${help(hk)}</div><div class="v">${v}</div><div class="d">${esc(d)}</div></div>`).join("");
  const maxKm = r.leaks ? Math.max(...r.leaks.by_cause.map((c) => c.km), 0.001) : 1;
  const bars = r.leaks && r.leaks.by_cause.length
    ? r.leaks.by_cause.map((c) => `<div class="bar-row"><div class="name"><b>${esc(c.text)}</b><span class="mono">${esc(c.cause)}</span></div><div class="track-h"><span style="width:${(100 * c.km) / maxKm}%"></span></div><div class="num">${num(c.km, 1)} km · ${num(c.places)} places</div></div>`).join("")
    : `<div class="hint">No leaks: every sampled levee holds water back or stands on higher ground.</div>`;
  const top = r.leaks && r.leaks.top.length
    ? `<div class="table-wrap"><table class="table"><thead><tr><th>Leak</th><th>Cause</th><th class="num">x</th><th class="num">y</th><th class="num">Levee m</th><th class="num">Lost margin m</th><th></th></tr></thead><tbody>
      ${r.leaks.top.map((l) => `<tr><td class="mono">${esc(l.leak_id)}</td><td><span class="chip plain">${esc(l.cause)}</span></td><td class="num mono">${l.x?.toFixed(0)}</td><td class="num mono">${l.y?.toFixed(0)}</td><td class="num">${num(l.levee_length_m)}</td><td class="num">${num(l.lost_margin_median_m, 2)}</td>
      <td><button class="icon-btn" title="Copy coordinates" data-copy="${l.x?.toFixed(1)}, ${l.y?.toFixed(1)}">${ico("copy")}</button></td></tr>`).join("")}</tbody></table></div>`
    : "";
  const files = r.files.length
    ? `<div class="table-wrap"><table class="table"><thead><tr><th>File</th><th class="num">Size</th><th class="num">Written</th></tr></thead><tbody>${r.files.map((f) => `<tr><td class="mono">${esc(f.name)}</td><td class="num">${bytes(f.size)}</td><td class="num">${ago(f.mtime)}</td></tr>`).join("")}</tbody></table></div>`
    : `<div class="hint">No files yet.</div>`;
  $("#page").innerHTML = `
    <div class="tiles">${tiles}</div>
    <section class="card card-pad"><div class="section-head"><div><h2>Where levees protect nothing, and why</h2><p>From leaks.gpkg: levee length per cause; start with the largest places.</p></div></div>
      <div class="bars">${bars}</div>${top ? `<div style="margin-top:20px">${top}</div>` : ""}</section>
    <div class="split">
      <section class="card card-pad"><div class="section-head"><div><h2>Overview</h2><p>From the latest report (reports/${esc(S.country)}).</p></div><button class="btn ghost small" data-open="reports">${ico("folder")}Reports</button></div>
        ${r.overview ? `<img class="overview" data-src="/api/file?path=${encodeURIComponent(r.overview)}&t=${Date.now()}" alt="Overview of ${esc(S.country)}">` : `<div class="empty">${ico("chart")}<h3>No overview yet</h3><p>Make a report to draw the footprint, levees and protected areas.</p></div>`}</section>
      <section class="card card-pad"><div class="section-head"><div><h2>Output files</h2><p class="mono">${esc(r.output_dir)}</p></div><div style="display:flex;gap:8px"><button class="btn ghost small" data-open="output">${ico("folder")}Output</button><button class="btn ghost small" data-open="work">${ico("folder")}Work</button></div></div>${files}</section>
    </div>`;
  if (r.overview) {
    // images need the token: loaded through fetch
    const img = $(".overview");
    fetch(img.dataset.src, { headers: { "X-Defenses-Token": TOKEN } }).then((res) => res.blob()).then((b) => (img.src = URL.createObjectURL(b)));
  }
  $("#page").querySelectorAll("[data-open]").forEach((b) => (b.onclick = () => openFolder(b.dataset.open)));
  $("#page").querySelectorAll("[data-copy]").forEach((b) => (b.onclick = async () => {
    try {
      await navigator.clipboard.writeText(b.dataset.copy);
      toast(`Copied ${b.dataset.copy}`);
    } catch (err) {
      toast(b.dataset.copy);
    }
  }));
}

// --- Configuration ---
const PATH_FIELDS = [
  ["dtm", "DTM", "file", true, "One or more GeoTIFF files of the terrain"],
  ["levees", "Levees", "file", true, "Levee lines you prepare (one layer per file)"],
  ["rivers", "Rivers", "file", true, "River network, e.g. MERIT reaches with an upstream area column"],
  ["river_guides", "River guide lines", "file", true, "Optional: lines of the river channels for snapping, e.g. waterways"],
];
function cfgGet(table, key) {
  const id = `${table}|${key}`;
  if (id in S.dirty) return S.dirty[id].value;
  let node = S.cfg || {};
  for (const part of table.split(".").filter(Boolean)) node = (node || {})[part];
  return (node || {})[key];
}
function setDirty(table, key, value) {
  S.dirty[`${table}|${key}`] = { table, key, value };
  renderSaveBar();
}
function renderSaveBar() {
  const n = Object.keys(S.dirty).length;
  const bar = $("#saveBar");
  if (!bar) return;
  bar.style.display = n ? "flex" : "none";
  $("#saveCount").textContent = `${n} unsaved change${n === 1 ? "" : "s"}`;
}
function pathList(table, key, label, kind, multiple, hint) {
  let value = cfgGet(table, key);
  const list = Array.isArray(value) ? value : value ? [value] : [];
  const rows = (list.length ? list : [""]).map((p, i) => `<div class="path-row"><span class="dot" data-path="${esc(p)}"></span>
    <input class="input mono" value="${esc(p)}" data-pl="${esc(table)}|${esc(key)}|${i}" placeholder="D:/data/…">
    ${multiple ? `<button class="icon-btn" title="Remove" data-del="${esc(table)}|${esc(key)}|${i}">${ico("trash")}</button>` : ""}</div>`).join("");
  return `<div class="field"><label>${esc(label)}${help(key)}</label><div class="path-list">${rows}</div>
    <div class="path-actions"><button class="btn ghost small" data-browse="${esc(table)}|${esc(key)}|${kind}|${multiple ? 1 : 0}">${ico("folder")}Browse…</button>
    ${multiple ? `<button class="btn ghost small" data-add="${esc(table)}|${esc(key)}">${ico("plus")}Add path</button>` : ""}</div>
    <span class="hint">${esc(hint)}</span></div>`;
}
function readPathList(table, key, multiple) {
  const inputs = [...document.querySelectorAll(`[data-pl^="${CSS.escape(`${table}|${key}|`)}"]`)];
  const list = inputs.map((i) => i.value.trim()).filter(Boolean);
  return multiple ? list : list[0] || "";
}
async function renderConfig() {
  $("#topActions").innerHTML = `<div class="tabs"><button data-tab="form" class="${S.cfgTab === "form" ? "active" : ""}">Form</button><button data-tab="raw" class="${S.cfgTab === "raw" ? "active" : ""}">Raw TOML</button></div>
    <button class="btn ghost" data-open="config">${ico("folder")}Folder</button><button class="btn soft" id="cfgCheck">${ico("check")}Check inputs</button>`;
  $("#topActions").querySelectorAll("[data-tab]").forEach((b) => (b.onclick = () => { S.cfgTab = b.dataset.tab; renderConfig(); }));
  $("#topActions").querySelector("[data-open]").onclick = () => openFolder("config");
  $("#cfgCheck").onclick = () => runStep("processing.check_inputs");
  if (!S.cfg || !Object.keys(S.dirty).length) {
    try {
      const data = await api("/api/config");
      S.cfg = data.config;
      S.cfgText = data.text;
      S.cfgError = data.error;
    } catch (err) {
      toast(err.message, true);
      return;
    }
  }
  if (S.cfgTab === "raw") return renderRaw();
  const countries = Object.keys((S.cfg || {}).countries || {}).sort();
  S.cfgCountry = countries.includes(S.cfgCountry) ? S.cfgCountry : countries.includes(S.country) ? S.country : countries[0] || "";
  const c = S.cfgCountry;
  const t = `countries.${c}`;
  const g = (k) => cfgGet("general", k);
  $("#page").innerHTML = `
    ${S.cfgError ? `<div class="banner err">${ico("alert")}<div>${esc(S.cfgError)} Use the Raw TOML tab to fix it.</div></div>` : ""}
    <section class="card card-pad"><div class="section-head"><div><h2>General</h2><p>Folders, the default country and the parallel workers.</p></div></div>
      <div class="grid-3">
        <div class="field"><label>Default country</label><div class="select-wrap"><select data-gen="country">${countries.map((x) => `<option ${x === g("country") ? "selected" : ""}>${esc(x)}</option>`).join("")}</select></div><span class="hint">Used when DEFENSES_COUNTRY is not set</span></div>
        <div class="field"><label>Workers${help("workers")}</label><input class="input" type="number" min="1" data-gen="n_workers" value="${esc(g("n_workers") ?? "")}"><span class="hint">RAM need = workers × unit RAM</span></div>
        <div class="field"><label>RAM per worker (GB)</label><input class="input" type="number" min="1" step="0.5" data-gen="max_ram_gb_per_worker" value="${esc(g("max_ram_gb_per_worker") ?? "")}"><span class="hint">build_units stops above it</span></div>
      </div>
      <div class="grid-2" style="margin-top:18px">
        ${pathList("general", "work_dir", "Work folder", "folder", false, "Intermediates, on a fast local disk")}
        ${pathList("general", "output_dir", "Output folder", "folder", false, "Results and logs")}
      </div></section>
    <section class="card card-pad"><div class="section-head"><div><h2>Country inputs</h2><p>Files you prepare for the country.</p></div>
      <div style="display:flex;gap:8px;align-items:center"><div class="select-wrap" style="width:150px"><select id="cfgCountry">${countries.map((x) => `<option ${x === c ? "selected" : ""}>${esc(x)}</option>`).join("")}</select></div>
      <button class="btn ghost small" id="addCountry">${ico("plus")}Add country</button></div></div>
      ${c ? `<div class="grid-2">${PATH_FIELDS.map(([k, label, kind, multi, hint]) => pathList(t, k, label, kind, multi, hint)).join("")}</div>
      <div class="grid-2" style="margin-top:18px">
        <div class="field"><label>Upstream area column of the rivers${help("upstream")}</label>
          <div class="select-wrap"><select id="upField"><option value="">– none –</option></select></div>
          <span class="hint" id="upHint">The NAME of the column with the upstream area in km² (needed for the river groups and the return period)</span></div>
      </div>` : `<div class="empty">${ico("plus")}<h3>No country yet</h3><p>Add the first country of your runs.</p></div>`}
    </section>
    ${c ? `<section class="card card-pad"><div class="section-head"><div><h2>Damage tables</h2><p>Used by the return period index.</p></div></div>
      <div class="grid-2">${pathList(`${t}.rp`, "max_damage_csv", "Maximum damage per category", "file", false, "CSV table")}${pathList(`${t}.rp`, "depth_damage_csv", "Depth-damage functions", "file", false, "CSV: depth_m and one damage fraction column per category")}</div></section>` : ""}
    <div class="save-bar" id="saveBar" style="display:none"><span id="saveCount"></span><div style="display:flex;gap:8px"><button class="btn ghost small" id="discardBtn">Discard</button><button class="btn primary small" id="saveBtn">${ico("save")}Save</button></div></div>`;
  bindConfig(c, t);
  renderSaveBar();
  checkPaths();
  if (c) loadColumns(t);
}
function bindConfig(c, t) {
  const page = $("#page");
  page.querySelectorAll("[data-gen]").forEach((el) => (el.onchange = () => {
    const k = el.dataset.gen;
    setDirty("general", k, el.type === "number" ? Number(el.value) : el.value);
  }));
  page.querySelectorAll("[data-pl]").forEach((el) => (el.onchange = () => {
    const [table, key] = el.dataset.pl.split("|");
    const multi = PATH_FIELDS.some((f) => f[0] === key && f[3]);
    setDirty(table, key, readPathList(table, key, multi));
    checkPaths();
    if (key === "rivers") loadColumns(table);
  }));
  page.querySelectorAll("[data-del]").forEach((b) => (b.onclick = () => {
    const [table, key, i] = b.dataset.del.split("|");
    const list = readPathList(table, key, true).filter((_, k) => k !== Number(i));
    setDirty(table, key, list);
    renderConfig();
  }));
  page.querySelectorAll("[data-add]").forEach((b) => (b.onclick = () => {
    const [table, key] = b.dataset.add.split("|");
    setDirty(table, key, [...readPathList(table, key, true), ""].filter((p, i, a) => p || i === a.length - 1));
    renderConfig();
  }));
  page.querySelectorAll("[data-browse]").forEach((b) => (b.onclick = async () => {
    const [table, key, kind, multi] = b.dataset.browse.split("|");
    try {
      const { paths, error } = await api("/api/browse", { kind, multiple: multi === "1" });
      if (error) toast(error, true);
      if (!paths.length) return;
      const value = multi === "1" ? [...readPathList(table, key, true), ...paths] : paths[0];
      setDirty(table, key, value);
      renderConfig();
      if (key === "rivers") loadColumns(table);
    } catch (err) {
      toast(err.message, true);
    }
  }));
  const cc = $("#cfgCountry");
  if (cc) cc.onchange = () => { S.cfgCountry = cc.value; renderConfig(); };
  $("#addCountry").onclick = () => {
    const code = (prompt("Code of the new country (e.g. CZE):") || "").trim().toUpperCase();
    if (!/^[A-Z0-9_]{2,8}$/.test(code)) return code && toast("Use 2–8 letters, digits or _", true);
    ["dtm", "levees", "rivers", "river_guides"].forEach((k) => setDirty(`countries.${code}`, k, []));
    setDirty(`countries.${code}`, "river_upstream_area_field", "");
    setDirty(`countries.${code}.rp`, "max_damage_csv", "");
    setDirty(`countries.${code}.rp`, "depth_damage_csv", "");
    S.cfgCountry = code;
    S.cfg.countries = S.cfg.countries || {};
    S.cfg.countries[code] = {};
    renderConfig();
  };
  const up = $("#upField");
  if (up) up.onchange = () => setDirty(t, "river_upstream_area_field", up.value);
  $("#saveBtn").onclick = saveConfig;
  $("#discardBtn").onclick = () => { S.dirty = {}; S.cfg = null; renderConfig(); };
}
async function checkPaths() {
  const dots = [...document.querySelectorAll(".dot[data-path]")];
  const inputs = [...document.querySelectorAll("[data-pl]")];
  const paths = inputs.map((i) => i.value.trim()).filter(Boolean);
  if (!paths.length) return;
  try {
    const found = await api("/api/paths", { paths });
    inputs.forEach((inp, k) => {
      const dot = dots[k];
      if (!dot) return;
      const p = inp.value.trim();
      dot.className = `dot ${p ? (found[p] ? "ok" : "no") : ""}`;
      dot.title = p ? (found[p] ? "Found" : "Not found on this computer") : "";
    });
  } catch (err) {
    /* the dots stay grey */
  }
}
async function loadColumns(table) {
  const sel = $("#upField");
  if (!sel) return;
  const current = cfgGet(table, "river_upstream_area_field") || "";
  const rivers = [].concat(cfgGet(table, "rivers") || []).filter(Boolean);
  let cols = [];
  if (rivers.length) {
    try {
      cols = (await api(`/api/columns?path=${encodeURIComponent(rivers[0])}`)).columns;
    } catch (err) {
      cols = [];
    }
  }
  const numeric = cols.filter((c) => /int|float|double|real/i.test(c.type));
  const names = numeric.map((c) => c.name);
  if (current && !names.includes(current)) names.unshift(current);
  sel.innerHTML = `<option value="">– none –</option>` + names.map((n) => `<option ${n === current ? "selected" : ""}>${esc(n)}</option>`).join("");
  const hint = $("#upHint");
  if (current && cols.length && !cols.some((c) => c.name === current)) {
    hint.innerHTML = `<span style="color:var(--coral)">“${esc(current)}” is not a column of ${esc(rivers[0].split("/").pop())}.</span> Pick the column with the upstream area in km².`;
  } else if (!rivers.length) {
    hint.textContent = "Set the rivers first; their numeric columns are offered here.";
  } else if (!cols.length) {
    hint.textContent = "The rivers file cannot be read on this computer.";
  }
}
async function saveConfig() {
  const changes = Object.values(S.dirty);
  try {
    await api("/api/config/set", { changes });
    S.dirty = {};
    S.cfg = null;
    toast(`Saved ${changes.length} change${changes.length === 1 ? "" : "s"} to the configuration`);
    await loadState();
    renderConfig();
  } catch (err) {
    toast(err.message, true);
  }
}
function renderRaw() {
  $("#page").innerHTML = `<section class="card card-pad"><div class="section-head"><div><h2>${esc(S.state.config_path.split(/[\\/]/).pop())}</h2><p>Edited as text; comments stay. It is checked before saving.</p></div>
      <div style="display:flex;gap:8px"><button class="btn ghost small" id="rawRevert">${ico("refresh")}Revert</button><button class="btn primary small" id="rawSave">${ico("save")}Save</button></div></div>
      <textarea class="mono" id="rawText" rows="30" spellcheck="false">${esc(S.cfgText)}</textarea></section>`;
  $("#rawRevert").onclick = () => { S.cfg = null; S.dirty = {}; renderConfig(); };
  $("#rawSave").onclick = async () => {
    try {
      await api("/api/config/raw", { text: $("#rawText").value });
      S.cfg = null;
      toast("Configuration saved");
      await loadState();
      renderConfig();
    } catch (err) {
      toast(err.message, true);
    }
  };
}

// --- Logs ---
async function renderLogs() {
  $("#topActions").innerHTML = `<button class="btn ghost" id="logsRefresh">${ico("refresh")}Refresh</button>`;
  $("#logsRefresh").onclick = renderLogs;
  try {
    S.logs = await api(`/api/logs?country=${encodeURIComponent(S.country)}`);
  } catch (err) {
    S.logs = [];
  }
  if (!S.logs.length) {
    $("#page").innerHTML = `<div class="card empty">${ico("doc")}<h3>No logs for ${esc(S.country)} yet</h3><p>Every script writes its log to the output folder.</p></div>`;
    return;
  }
  if (!S.logs.some((l) => l.name === S.logName)) S.logName = S.logs[0].name;
  $("#page").innerHTML = `<div class="logs">
    <div class="card log-list">${S.logs.map((l) => {
      const m = l.name.match(/^(.*)_(\d{8})_(\d{6})\.log$/);
      const title = m ? m[1].replace(/_/g, " ") : l.name;
      const clock = m ? `${m[3].slice(0, 2)}:${m[3].slice(2, 4)} · ` : "";
      return `<div class="log-item ${l.name === S.logName ? "active" : ""}" data-name="${esc(l.name)}"><b>${esc(title)}</b><span>${clock}${ago(l.mtime)} · ${bytes(l.size)}</span></div>`;
    }).join("")}</div>
    <div class="card viewer"><div class="viewer-head"><strong class="mono" style="font-size:12.5px">${esc(S.logName)}</strong>
      <div class="tabs" id="logFilter">${["all", "warnings", "errors"].map((f) => `<button data-f="${f}" class="${S.logFilter === f ? "active" : ""}">${f[0].toUpperCase() + f.slice(1)}</button>`).join("")}</div>
      <input class="input small" id="logSearch" placeholder="Search" value="${esc(S.logSearch)}"></div><pre id="logView">Loading…</pre></div></div>`;
  $("#page").querySelectorAll("[data-name]").forEach((el) => (el.onclick = () => { S.logName = el.dataset.name; renderLogs(); }));
  $("#logFilter").querySelectorAll("button").forEach((b) => (b.onclick = () => { S.logFilter = b.dataset.f; renderLogs(); }));
  $("#logSearch").oninput = (e) => { S.logSearch = e.target.value; showLog(); };
  try {
    S.logText = (await api(`/api/log?country=${encodeURIComponent(S.country)}&name=${encodeURIComponent(S.logName)}`)).text;
  } catch (err) {
    S.logText = err.message;
  }
  showLog();
}
function interleave(line, pattern) {
  let out = "", last = 0;
  line.replace(pattern, (m, offset) => {
    out += esc(line.slice(last, offset)) + `<mark>${esc(m)}</mark>`;
    last = offset + m.length;
    return m;
  });
  return out + esc(line.slice(last));
}
function showLog() {
  const q = S.logSearch.toLowerCase();
  const lines = S.logText.split("\n").filter((l) => (S.logFilter === "warnings" ? / WARNING | ERROR /.test(l) : S.logFilter === "errors" ? / ERROR |Traceback/.test(l) : true) && (!q || l.toLowerCase().includes(q)));
  const view = $("#logView");
  if (!view) return;
  const pattern = q ? new RegExp(q.replace(/[.*+?^${}()|[\]\\]/g, "\\$&"), "gi") : null;
  view.innerHTML = lines.map((l) => {
    // the search is marked in the text, then escaped piece by piece
    const marked = pattern ? interleave(l, pattern) : esc(l);
    return `<span class="${lineClass(l)}">${marked}</span>`;
  }).join("") || `<span class="ln-info">No matching lines.</span>`;
}

boot();
