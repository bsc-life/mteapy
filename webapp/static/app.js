// mteapy task visualizer -- frontend
//
// Two halves: (1) task table + upload/sample controls, driving the live
// /api/datasets/.../scores and /api/datasets/.../network endpoints; (2) the
// route network renderer, which is mteapy's own decorator callbacks
// (evidence/score-driven coloring, flux-driven edge width, the GPR/complex
// detail panel) wired into netgraph-viz.js's generic layout/draw/zoom/click
// layer (this file's own copy, in this same directory -- see its header for
// why it's not a shared cross-project package).

let datasetId = null;
let currentSample = null;
let currentModel = null;     // model key, e.g. "Human-GEM-2.0.1"
let currentTaskList = null;  // task list name within that model
let modelsInfo = [];         // [{key, name, version, available, task_lists: [{name, n_tasks, ...}]}]
let currentRun = null;       // {run_id, config, samples, byTask: {task_id: {scores, complete, tied}}} once a run is done
let datasetColumns = [];     // numeric columns of the loaded table (samples, or a DE table's columns)
let datasetSuggested = {};   // {lfc_col, pvalue_col}: best guess at a DE table's columns
let checkBlocked = false;    // the pre-run check found an error: Run stays disabled
let methodsInfo = [];        // [{id, label, description, params: [{name, label, type, default, help, choices}]}]
let allTasks = [];      // [{task_id, description, n_routes, score, is_tied, is_complete, winning_route_ids}]
let sortKey = "task_id";
let sortAsc = true;
let selectedKey = null; // task_id within the current model/task list
let currentData = null; // {task_id, task_description, tissue, has_data, panels: [...]}
let currentIndex = 0;
let orientation = "TB";

function rowKey(t) { return t.task_id; }

const statusMsg = document.getElementById("status-msg");
const sampleSelect = document.getElementById("sample-select");
const taskTbody = document.getElementById("task-tbody");
const taskSearch = document.getElementById("task-search");
const modelSelect = document.getElementById("model-select");
const taskListSelect = document.getElementById("tasklist-select");

function setStatus(text) { statusMsg.textContent = text; }

// ---------- bootstrapping ----------

async function loadModels() {
  const res = await fetch("/api/models");
  const data = await res.json();
  modelsInfo = data.models;
  window.GTEX_DATASET_ID = data.gtex_dataset_id;
  modelSelect.innerHTML = "";
  modelsInfo.forEach(m => {
    const opt = document.createElement("option");
    opt.value = m.key;
    opt.textContent = (m.version ? `${m.name} ${m.version}` : m.name) + (m.available ? "" : " (unusable)");
    opt.title = m.available ? (m.description || "") : m.problems.join("; ");
    opt.disabled = !m.available;
    modelSelect.appendChild(opt);
  });
  const first = modelsInfo.find(m => m.available);
  if (!first) {
    const why = modelsInfo.flatMap(m => m.problems).join("; ");
    setStatus("No usable model found" + (why ? `: ${why}` : " (see the server log)."));
    return;
  }
  modelSelect.disabled = false;
  modelSelect.value = first.key;
  await onModelChange();
}

async function onModelChange(preferredTaskList = null) {
  currentModel = modelSelect.value;
  const info = modelsInfo.find(m => m.key === currentModel);
  taskListSelect.innerHTML = "";
  info.task_lists.forEach(tl => {
    const opt = document.createElement("option");
    opt.value = tl.name;
    opt.textContent = `${tl.name} (${tl.n_tasks_with_routes}/${tl.n_tasks} tasks with routes)`;
    taskListSelect.appendChild(opt);
  });
  taskListSelect.disabled = info.task_lists.length === 0;
  if (!info.task_lists.length) { allTasks = []; renderTaskTable(); setStatus("This model has no task lists."); return; }
  // Prefer the first list that actually has routes.
  const preferred = info.task_lists.find(tl => tl.name === preferredTaskList)
    || info.task_lists.find(tl => tl.n_tasks_with_routes > 0) || info.task_lists[0];
  taskListSelect.value = preferred.name;
  await onTaskListChange();
}

async function onTaskListChange() {
  currentTaskList = taskListSelect.value;
  clearRun();
  selectedKey = null;
  currentData = null;
  document.getElementById("panel-title-text").textContent = "Select a task to render its network";
  document.getElementById("task-desc").textContent = "";
  document.getElementById("panel-host").innerHTML = "";
  setStatus(`Loading ${currentModel} (the first use of a model takes a few seconds)...`);
  modelSelect.disabled = taskListSelect.disabled = true;
  try {
    const res = await fetch(`/api/models/${encodeURIComponent(currentModel)}/task_lists/${encodeURIComponent(currentTaskList)}/tasks`);
    if (!res.ok) { setStatus(`Could not load tasks: ${(await res.json()).detail || res.statusText}`); return; }
    const data = await res.json();
    allTasks = data.tasks.map(t => ({ ...t, score: null, is_tied: false, is_complete: false, winning_route_ids: [] }));
    renderTaskTable();
    setStatus(`${allTasks.length} tasks in ${currentTaskList}.`);
  } finally {
    modelSelect.disabled = false;
    taskListSelect.disabled = false;
  }
  updateRunButton();
}

modelSelect.addEventListener("change", () => onModelChange());
taskListSelect.addEventListener("change", onTaskListChange);

document.getElementById("load-gtex-btn").addEventListener("click", async () => {
  if (!window.GTEX_DATASET_ID) { setStatus("No bundled GTEx dataset available on this server."); return; }
  const res = await fetch(`/api/datasets/${window.GTEX_DATASET_ID}`);
  const info = await res.json();
  onDatasetLoaded(info);
  setStatus(`Loaded GTEx example: ${info.samples.length} tissues.`);
});

document.getElementById("file-input").addEventListener("change", async (ev) => {
  const file = ev.target.files[0];
  if (!file) return;
  setStatus(`Uploading ${file.name}...`);
  const form = new FormData();
  form.append("file", file);
  const res = await fetch("/api/datasets", { method: "POST", body: form });
  if (!res.ok) { setStatus(`Upload failed: ${(await res.json()).detail || res.statusText}`); return; }
  const info = await res.json();
  onDatasetLoaded(info);
  setStatus(`Loaded ${file.name}: ${info.n_genes} genes, ${info.samples.length} sample(s).`);
});

function populateSampleSelect(samples) {
  sampleSelect.innerHTML = "";
  samples.forEach(s => {
    const opt = document.createElement("option");
    opt.value = s; opt.textContent = s;
    sampleSelect.appendChild(opt);
  });
  sampleSelect.disabled = false;
  sampleSelect.dispatchEvent(new Event("change"));
  updateRunButton();
}

sampleSelect.addEventListener("change", () => {
  if (!sampleSelect.value) return;
  currentSample = sampleSelect.value;
  applyRunScores();
  if (selectedKey) selectTask(selectedKey);
});

// ---------- analysis: methods, parameters, runs ----------

const methodSelect = document.getElementById("method-select");
const paramsHost = document.getElementById("method-params");
const runBtn = document.getElementById("run-btn");
const cancelBtn = document.getElementById("cancel-btn");
const runProgress = document.getElementById("run-progress");
const runStatus = document.getElementById("run-status");
let activeRunId = null;

function setRunStatus(text, isError = false) {
  runStatus.textContent = text;
  runStatus.classList.toggle("err", isError);
}

async function loadMethods() {
  const res = await fetch("/api/methods");
  methodsInfo = (await res.json()).methods;
  methodSelect.innerHTML = "";
  methodsInfo.forEach(m => {
    const opt = document.createElement("option");
    opt.value = m.id; opt.textContent = m.label; opt.title = m.description;
    methodSelect.appendChild(opt);
  });
  methodSelect.disabled = methodsInfo.length === 0;
  renderParams();
}

function currentMethod() { return methodsInfo.find(m => m.id === methodSelect.value); }

const INPUT_KIND_NOTE = {
  expression: "Input: an expression matrix (genes x samples, non-negative).",
  log_fold_change: "Input: a differential-expression table (log fold changes, so negative values; optional p-values).",
};

// A "column" parameter picks one of the loaded table's columns; the table's best
// guess (a log-FC / p-value looking name) is preselected, an optional one can be left empty.
function fillColumnSelect(select, p) {
  const keep = select.value;
  select.innerHTML = "";
  if (p.optional) { const o = document.createElement("option"); o.value = ""; o.textContent = "(none)"; select.appendChild(o); }
  if (!datasetColumns.length && !p.optional) { const o = document.createElement("option"); o.value = ""; o.textContent = "(load data)"; select.appendChild(o); }
  datasetColumns.forEach(c => { const o = document.createElement("option"); o.value = o.textContent = c; select.appendChild(o); });
  const guess = datasetSuggested[p.name];
  select.value = [...select.options].some(o => o.value === keep && keep !== "") ? keep
    : (guess && datasetColumns.includes(guess) ? guess : (p.optional ? "" : (datasetColumns[0] || "")));
}

function refreshColumnParams() {
  const m = currentMethod();
  if (!m) return;
  m.params.filter(p => p.type === "column").forEach(p => {
    const el = paramsHost.querySelector(`[data-param="${p.name}"]`);
    if (el) fillColumnSelect(el, p);
  });
}

function renderParams() {
  const m = currentMethod();
  paramsHost.innerHTML = "";
  document.getElementById("method-desc").textContent = m ? `${m.description} ${INPUT_KIND_NOTE[m.input_kind] || ""}` : "";
  if (!m) return;
  m.params.forEach(p => {
    const wrap = document.createElement("label");
    wrap.className = "param"; wrap.title = p.help;
    wrap.append(`${p.label}: `);
    let input;
    if (p.type === "column") {
      input = document.createElement("select");
      fillColumnSelect(input, p);
    } else if (p.type === "choice") {
      input = document.createElement("select");
      p.choices.forEach(c => { const o = document.createElement("option"); o.value = o.textContent = c; input.appendChild(o); });
      input.value = p.default;
    } else if (p.type === "bool") {
      input = document.createElement("input"); input.type = "checkbox"; input.checked = !!p.default;
    } else {
      input = document.createElement("input"); input.type = "number"; input.value = p.default ?? "";
      input.step = p.type === "int" ? "1" : "any";
      if (p.min !== null) input.min = p.min;
      if (p.max !== null) input.max = p.max;
      if (p.optional) input.placeholder = "auto";
    }
    input.dataset.param = p.name; input.dataset.type = p.type;
    wrap.dataset.when = JSON.stringify(p.when || []);
    wrap.appendChild(input);
    paramsHost.appendChild(wrap);
  });
  updateParamVisibility();
  scheduleCheck();
}

function updateParamVisibility() {
  const current = collectParams();
  paramsHost.querySelectorAll(".param").forEach(wrap => {
    const when = JSON.parse(wrap.dataset.when || "[]");
    wrap.hidden = !when.every(([name, value]) => current[name] === value);
  });
}

function collectParams() {
  const out = {};
  paramsHost.querySelectorAll("[data-param]").forEach(el => {
    const numeric = el.dataset.type === "int" || el.dataset.type === "float";
    out[el.dataset.param] = el.dataset.type === "bool" ? el.checked
      : numeric ? (el.value === "" ? null : Number(el.value))
      : el.dataset.type === "column" ? (el.value === "" ? null : el.value) : el.value;
  });
  return out;
}

function updateRunButton() {
  document.getElementById("save-btn").disabled = document.getElementById("tsv-btn").disabled = !currentRun;
  runBtn.disabled = !(datasetId && currentModel && currentTaskList && currentMethod() && !activeRunId) || checkBlocked;
  runBtn.title = datasetId ? "Score every task of the selected list against every sample"
                           : "Load an expression dataset first";
}

methodSelect.addEventListener("change", renderParams);
paramsHost.addEventListener("change", () => { updateParamVisibility(); scheduleCheck(); });

// ---------- pre-run input check ----------
// The server looks at the loaded data against the chosen method and parameters
// before anything runs: errors keep Run disabled, warnings (e.g. a "log-FC" column
// with no negative values, which looks like expression) are shown but never block.

const warnBox = document.getElementById("run-warnings");
let checkTimer = null, checkSeq = 0;

function showChecks(errors, warnings) {
  warnBox.innerHTML = "";
  errors.forEach(t => { const d = document.createElement("div"); d.className = "err"; d.textContent = t; warnBox.appendChild(d); });
  warnings.forEach(t => { const d = document.createElement("div"); d.className = "warn"; d.textContent = `Warning: ${t}`; warnBox.appendChild(d); });
  warnBox.hidden = errors.length + warnings.length === 0;
  checkBlocked = errors.length > 0;
  updateRunButton();
}

function scheduleCheck() {
  clearTimeout(checkTimer);
  checkTimer = setTimeout(runCheck, 250);
}

async function runCheck() {
  const m = currentMethod();
  const seq = ++checkSeq;
  if (activeRunId) return;       // a run in flight keeps showing the warnings it started with
  if (!(datasetId && currentModel && currentTaskList && m)) { showChecks([], []); return; }
  const body = { model: currentModel, task_list: currentTaskList, dataset_id: datasetId, method: m.id, params: collectParams() };
  try {
    const res = await fetch("/api/runs/check", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) });
    const data = await res.json();
    if (seq !== checkSeq) return;
    if (!res.ok) showChecks([typeof data.detail === "string" ? data.detail : JSON.stringify(data.detail)], []);
    else showChecks(data.errors, data.warnings);
  } catch (e) { if (seq === checkSeq) showChecks([], []); }
}

function onDatasetLoaded(info) {
  datasetId = info.dataset_id;
  datasetColumns = info.samples;
  datasetSuggested = info.suggested || {};
  clearRun();
  populateSampleSelect(info.samples);
  refreshColumnParams();
  scheduleCheck();
}

function clearRun(note) {
  currentRun = null;
  allTasks = allTasks.map(t => ({ ...t, score: null, pvalue: null, is_tied: false, is_complete: false, winning_route_ids: [] }));
  renderTaskTable();
  if (note) setRunStatus(note);
  else if (!activeRunId) setRunStatus("");
  updateRunButton();
}

function applyRunScores() {
  if (!currentRun || !currentSample) return;
  const idx = currentRun.samples.indexOf(currentSample);
  if (idx < 0) return;
  allTasks = allTasks.map(t => {
    const r = currentRun.byTask[t.task_id];
    return r ? { ...t, score: r.scores[idx], is_tied: r.tied[idx], is_complete: r.complete[idx],
                 pvalue: r.pvalues ? r.pvalues[idx] : null }
             : { ...t, score: null, is_tied: false, is_complete: false, pvalue: null };
  });
  renderTaskTable();
}

async function startRun() {
  const m = currentMethod();
  const body = { model: currentModel, task_list: currentTaskList, dataset_id: datasetId, method: m.id, params: collectParams() };
  runBtn.disabled = true;
  setRunStatus("Starting...");
  const res = await fetch("/api/runs", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) });
  if (!res.ok) {
    const detail = (await res.json()).detail;
    setRunStatus(`Could not start: ${typeof detail === "string" ? detail : JSON.stringify(detail)}`, true);
    updateRunButton();
    return;
  }
  const info = await res.json();
  activeRunId = info.run_id;
  cancelBtn.hidden = false;
  runProgress.hidden = false;
  currentRun = null;
  await pollRun(info);
}

async function pollRun(info) {
  const runId = info.run_id;
  for (;;) {
    runProgress.max = Math.max(info.total, 1);
    runProgress.value = info.done;
    setRunStatus(info.status === "queued" ? "Queued..." : info.total === 0 ? "Starting..."
      : info.unit === "permutations" ? `Permutation test ${info.done}/${info.total}...` : `Scoring ${info.done}/${info.total} samples...`);
    if (info.status !== "queued" && info.status !== "running") break;
    await new Promise(r => setTimeout(r, 350));
    const res = await fetch(`/api/runs/${runId}`);
    if (!res.ok) { setRunStatus("Lost track of the run.", true); break; }
    info = await res.json();
  }
  activeRunId = null;
  cancelBtn.hidden = true;
  runProgress.hidden = true;
  if (info.status === "done") {
    const res = await fetch(`/api/runs/${runId}/results`);
    const data = await res.json();
    currentRun = {
      run_id: runId, config: info.config, samples: data.samples,
      byTask: Object.fromEntries(data.tasks.map(t => [t.task_id, t])),
    };
    const p = Object.entries(info.config.params).map(([k, v]) => `${k}=${v}`).join(", ");
    setRunStatus(`${info.config.method} (${p}) done in ${info.seconds}s: ${data.tasks.length} tasks x ${data.samples.length} samples ` +
                 `(${info.genes_matched} model genes found in the data).`);
    showChecks([], info.warnings || []);
    if (data.samples.join("\u0000") !== [...sampleSelect.options].map(o => o.value).join("\u0000")) populateSampleSelect(data.samples);
    applyRunScores();
    if (selectedKey) selectTask(selectedKey);
  } else if (info.status === "cancelled") {
    setRunStatus("Run cancelled.");
  } else {
    setRunStatus(`Run failed: ${info.error}`, true);
  }
  updateRunButton();
}

// ---------- saving and loading results ----------

function download(url) {
  const a = document.createElement("a");
  a.href = url; a.download = "";
  document.body.appendChild(a); a.click(); a.remove();
}
document.getElementById("save-btn").addEventListener("click", () => currentRun && download(`/api/runs/${currentRun.run_id}/export`));
document.getElementById("tsv-btn").addEventListener("click", () => currentRun && download(`/api/runs/${currentRun.run_id}/export?format=tsv`));

const loadInput = document.getElementById("load-input");
document.getElementById("load-btn").addEventListener("click", () => loadInput.click());
loadInput.addEventListener("change", async () => {
  const file = loadInput.files[0];
  loadInput.value = "";
  if (!file) return;
  setRunStatus(`Loading ${file.name}...`);
  let bundle;
  try { bundle = JSON.parse(await file.text()); }
  catch { setRunStatus(`${file.name} is not a results file (invalid JSON).`, true); return; }
  const res = await fetch("/api/runs/import", { method: "POST", headers: { "Content-Type": "application/json" },
                                                body: JSON.stringify(bundle) });
  if (!res.ok) {
    const detail = (await res.json()).detail;
    setRunStatus(`Could not load ${file.name}: ${typeof detail === "string" ? detail : JSON.stringify(detail)}`, true);
    return;
  }
  const info = await res.json();
  const cfg = info.config;
  // Switch the whole view to what the results were computed on.
  modelSelect.value = cfg.model;
  await onModelChange(cfg.task_list);
  methodSelect.value = cfg.method;
  renderParams();
  paramsHost.querySelectorAll("[data-param]").forEach(el => {
    const v = cfg.params[el.dataset.param];
    if (el.tagName === "SELECT" && v && ![...el.options].some(o => o.value === v)) {
      const o = document.createElement("option"); o.value = o.textContent = v; el.appendChild(o);   // saved column, table not loaded
    }
    if (el.dataset.type === "bool") el.checked = !!v; else el.value = v ?? "";
  });
  updateParamVisibility();
  datasetId = info.dataset_id;
  datasetColumns = info.dataset_id ? info.samples : [];
  datasetSuggested = {};
  const data = await (await fetch(`/api/runs/${info.run_id}/results`)).json();
  currentRun = { run_id: info.run_id, config: cfg, samples: data.samples,
                 byTask: Object.fromEntries(data.tasks.map(t => [t.task_id, t])) };
  populateSampleSelect(data.samples);
  const p = Object.entries(cfg.params).map(([k, v]) => `${k}=${v}`).join(", ");
  setRunStatus(`Loaded ${cfg.method} (${p}), saved ${cfg.saved_at || "earlier"}: ${data.tasks.length} tasks x ` +
               `${data.samples.length} samples.` + (datasetId ? "" : " (Nothing to re-run: these results are loaded as scores only.)"));
  applyRunScores();
  showChecks([], info.warnings || []);
});

runBtn.addEventListener("click", startRun);
cancelBtn.addEventListener("click", async () => {
  if (activeRunId) await fetch(`/api/runs/${activeRunId}`, { method: "DELETE" });
});

// ---------- task table ----------

document.querySelectorAll("table.tasks th[data-sort]").forEach(th => {
  th.addEventListener("click", () => {
    const key = th.dataset.sort;
    sortAsc = sortKey === key ? !sortAsc : true;
    sortKey = key;
    renderTaskTable();
  });
});

taskSearch.addEventListener("input", renderTaskTable);

function renderTaskTable() {
  const q = taskSearch.value.trim().toLowerCase();
  let rows = allTasks.filter(t =>
    !q || t.task_id.toLowerCase().includes(q) || (t.description || "").toLowerCase().includes(q));

  rows.sort((a, b) => {
    let av = a[sortKey], bv = b[sortKey];
    if (sortKey === "task_id") { av = +av; bv = +bv; }
    if (av === null || av === undefined) av = -Infinity;
    if (bv === null || bv === undefined) bv = -Infinity;
    if (av < bv) return sortAsc ? -1 : 1;
    if (av > bv) return sortAsc ? 1 : -1;
    return 0;
  });

  const showP = !!currentRun && Object.values(currentRun.byTask).some(t => t.pvalues);
  document.getElementById("pvalue-th").hidden = !showP;
  document.querySelector(".wrap").classList.toggle("has-pvalue", showP);
  document.querySelector(".wrap").classList.toggle("legend-signed", !!currentRun && currentRun.config.method === "TIDE");
  taskTbody.innerHTML = "";
  if (!rows.length) {
    taskTbody.innerHTML = `<tr><td colspan="${showP ? 5 : 4}" class="placeholder">No matching tasks.</td></tr>`;
    return;
  }
  for (const t of rows) {
    const tr = document.createElement("tr");
    tr.className = "task-row" + (rowKey(t) === selectedKey ? " selected" : "");
    const scoreCell = t.score === null
      ? `<span class="placeholder">&mdash;</span>`
      : `<span class="score-pill" style="background:${scorePillColor(t)};">${t.score.toFixed(2)}</span>${t.is_tied ? " tied" : ""}`;
    tr.innerHTML = `
      <td>${t.task_id}</td>
      <td>${t.description || ""}</td>
      <td>${t.n_routes}</td>
      <td>${scoreCell}</td>
      ${showP ? `<td>${formatPvalue(t.pvalue)}</td>` : ""}
    `;
    tr.addEventListener("click", () => selectTask(t.task_id));
    taskTbody.appendChild(tr);
  }
}

// An empirical p-value is a count out of the permutations, so 0 means "below 1/N", not zero.
function formatPvalue(p) {
  if (p === null || p === undefined) return `<span class="placeholder">&mdash;</span>`;
  if (p === 0) {
    const n = currentRun && currentRun.config.params.n_permutations;
    return n ? `&lt; ${(1 / n).toPrecision(1)}` : "0";
  }
  return p < 0.001 ? p.toExponential(1) : p.toFixed(3);
}

function scorePillColor(t) {
  if (t.score === null) return "var(--score-empty)";
  if (t.score === 0) return "var(--evidence-none-bg)";
  return t.is_complete ? "#d2f0e3" : "#fbe9c9";
}

async function selectTask(taskId) {
  selectedKey = taskId;
  renderTaskTable();
  document.getElementById("panel-title-text").textContent = "Loading...";
  const hasSample = currentRun && currentSample;
  const qs = hasSample ? `?run_id=${encodeURIComponent(currentRun.run_id)}&sample=${encodeURIComponent(currentSample)}` : "";
  const res = await fetch(`/api/models/${encodeURIComponent(currentModel)}/task_lists/${encodeURIComponent(currentTaskList)}` +
                          `/tasks/${encodeURIComponent(taskId)}/network${qs}`);
  if (!res.ok) {
    document.getElementById("panel-title-text").textContent = `Error: ${(await res.json()).detail || res.statusText}`;
    return;
  }
  currentData = await res.json();
  currentIndex = 0;
  document.getElementById("task-desc").textContent = currentData.task_description || "";
  setStatus(hasSample ? "" : "No results yet -- showing all route variants, unscored (run an analysis to score them).");
  showPanel(0);
}

// ==========================================================================
// Route network renderer (ported from route_network_template_v2.html)
// ==========================================================================

// Colour by the size of a reaction's score on the log scale, accent for positive
// scores and a second hue for negative ones (TIDE: down-regulated; the sign is
// kept by the signed OR function). Expression scores are never negative.
const DOWN_COLOR = "#8e44ad";
function scaleFor(panel) {
  const mags = panel.nodes
    .filter(n => n.type === "reaction" && (n.evidence === "supported" || n.evidence === "ambiguous"))
    .map(n => Math.log1p(Math.abs(n.score || 0)));
  const max = Math.max(1e-6, ...mags);
  return score => {
    const t = Math.min(1, Math.log1p(Math.abs(score)) / max);
    return d3.interpolate(getComputedColor("--score-empty"), score < 0 ? DOWN_COLOR : getComputedColor("--accent"))(t);
  };
}
function getComputedColor(varName) {
  return getComputedStyle(document.documentElement).getPropertyValue(varName).trim() || "#999";
}

function nodeFillColor(d, colorScale) {
  if (d.type !== "reaction") {
    return d.io === "input" ? getComputedColor("--met-input")
      : d.io === "output" ? getComputedColor("--met-output")
      : getComputedColor("--surface");
  }
  if (d.evidence === "no_evidence") return getComputedColor("--evidence-none");
  if (d.evidence === "no_gpr") return getComputedColor("--score-empty");
  return colorScale(d.score || 0);
}

function widthScaleFor(panel) {
  const fluxes = panel.edges.map(d => Math.abs(d.flux || 0));
  const max = Math.max(1e-6, ...fluxes);
  return d3.scaleSqrt().domain([0, max]).range([1.2, 6]).clamp(true);
}

function nodeSize(d) { return d.type === "reaction" ? 22 : 18; }

// Topology-only mode (no sample loaded): every reaction naturally comes
// back as no_gpr/no_evidence against an empty signal, which would
// otherwise paint the whole route red/gray as if evidence were actually
// absent. Render plain neutral fills instead so "no data yet" reads as
// "no data yet", not as a negative finding.
function nodeFillColorNeutral(d) {
  if (d.type !== "reaction") {
    return d.io === "input" ? getComputedColor("--met-input")
      : d.io === "output" ? getComputedColor("--met-output")
      : getComputedColor("--surface");
  }
  return getComputedColor("--score-empty");
}

function renderPanel(container, panel, viewWidth, viewHeight, hasData) {
  const colorScale = scaleFor(panel);
  const widthScale = widthScaleFor(panel);
  const fill = d => hasData ? nodeFillColor(d, colorScale) : nodeFillColorNeutral(d);

  NetGraphViz.render(container, panel, {
    orientation,
    viewWidth,
    viewHeight,
    fitToView: true,
    nodeSize,
    nodeShape: d => d.type === "reaction" ? "rect" : "circle",
    nodeFill: fill,
    nodeStroke: d => d.type === "reaction" ? "rgba(11,11,11,0.25)" : "var(--met-node)",
    edgeColor: (edge, sourceNode, targetNode) => {
      const rxnNode = sourceNode.type === "reaction" ? sourceNode : targetNode;
      return fill(rxnNode);
    },
    edgeWidth: d => widthScale(Math.abs(d.flux || 0)),
    nodeLabel: d => d.type === "reaction" ? NetGraphViz.truncateLabel(d.label, 26) : d.label,
    nodeLabelStyle: d => d.type === "reaction" ? { fontWeight: 600 } : {},
    nodeTooltip: d => d.type === "reaction"
      ? `${d.label}${d.ec_code ? " (EC " + d.ec_code + ")" : ""} — flux ${d.flux.toFixed(3)}, score ${d.score.toFixed(2)}`
      : `${d.label} [${d.compartment}]`,
    onNodeClick: d => renderDetail(d),
  });
}

// ==========================================================================
// Export: the network picture, its data, and the context-aware evidence
// ==========================================================================

function downloadBlob(filename, blob) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url; a.download = filename;
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 2000);
}

function fileStem(...parts) {
  return parts.filter(p => p !== null && p !== undefined && p !== "").join("_").replace(/[^A-Za-z0-9._-]+/g, "-");
}

const STYLE_PROPS = ["fill", "fill-opacity", "stroke", "stroke-width", "stroke-opacity", "stroke-dasharray", "opacity",
  "font-family", "font-size", "font-weight", "font-style", "text-anchor", "dominant-baseline", "paint-order"];

// A standalone copy of an <svg>: every style the page resolves for it (CSS variables
// and theme included) written inline, the zoom/pan reset, and the viewBox set to
// the drawing itself, so the file looks the same anywhere and shows all of it.
function standaloneSvg(svgEl, { fit = true, background = null } = {}) {
  const clone = svgEl.cloneNode(true);
  const src = [svgEl, ...svgEl.querySelectorAll("*")];
  const dst = [clone, ...clone.querySelectorAll("*")];
  src.forEach((el, i) => {
    const cs = getComputedStyle(el);
    const css = STYLE_PROPS.map(p => `${p}:${cs.getPropertyValue(p)}`).join(";");
    dst[i].setAttribute("style", css);
    // the resolved values are in the style now; attributes holding var(--x) would mean nothing elsewhere
    for (const a of ["fill", "stroke", "opacity"]) dst[i].removeAttribute(a);
  });
  clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  let width, height;
  if (fit) {
    const g = clone.querySelector(":scope > g");
    const live = svgEl.querySelector(":scope > g");
    const box = live.getBBox();
    g.removeAttribute("transform");
    const pad = 20;
    const vb = [box.x - pad, box.y - pad, box.width + 2 * pad, box.height + 2 * pad];
    clone.setAttribute("viewBox", vb.join(" "));
    [width, height] = [vb[2], vb[3]];
  } else {
    const vb = svgEl.viewBox.baseVal;
    [width, height] = vb && vb.width ? [vb.width, vb.height] : [svgEl.clientWidth, svgEl.clientHeight];
  }
  clone.setAttribute("width", width);
  clone.setAttribute("height", height);
  if (background) {
    const [x, y, w, h] = clone.getAttribute("viewBox").split(" ");
    const rect = document.createElementNS("http://www.w3.org/2000/svg", "rect");
    for (const [k, v] of Object.entries({ x, y, width: w, height: h, fill: background })) rect.setAttribute(k, v);
    clone.insertBefore(rect, clone.firstChild);
  }
  return { text: new XMLSerializer().serializeToString(clone), width, height };
}

function svgToPng(svg, scale) {
  return new Promise((resolve, reject) => {
    const img = new Image();
    img.onload = () => {
      const canvas = document.createElement("canvas");
      canvas.width = Math.round(svg.width * scale);
      canvas.height = Math.round(svg.height * scale);
      canvas.getContext("2d").drawImage(img, 0, 0, canvas.width, canvas.height);
      canvas.toBlob(b => b ? resolve(b) : reject(new Error("PNG encoding failed")), "image/png");
    };
    img.onerror = () => reject(new Error("could not rasterise the SVG"));
    img.src = "data:image/svg+xml;charset=utf-8," + encodeURIComponent(svg.text);
  });
}

function exportName(ext, ...extra) {
  const panel = currentData.panels[currentIndex];
  return `${fileStem(currentTaskList, "task" + currentData.task_id, "route" + panel.route_id,
                     currentData.has_data ? currentData.tissue : null, ...extra)}.${ext}`;
}

const TSV_COLUMNS = ["task_list", "task_id", "route_id", "sample", "reaction_id", "reaction_name", "equation", "flux",
  "evidence", "reaction_score", "complex", "complex_score", "complex_in_winning", "gene", "gene_value"];

// One row per (reaction, candidate complex, gene): the evidence the score was built
// from. Reactions without a gene association keep a single row with empty complex fields.
function evidenceRows(reactionNodes) {
  const panel = currentData.panels[currentIndex];
  const hasData = currentData.has_data;
  const rows = [];
  for (const n of reactionNodes) {
    const winning = new Set((n.winning_complexes || []).map(c => JSON.stringify([...c].sort())));
    const base = [currentTaskList, currentData.task_id, panel.route_id, hasData ? currentData.tissue : "",
      n.reaction_id, n.label, n.equation, n.flux, hasData ? n.evidence : "", hasData ? n.score : ""];
    const complexes = n.complexes || [];
    if (!complexes.length) { rows.push([...base, "", "", "", "", ""]); continue; }
    complexes.forEach((c, ci) => c.forEach(g => rows.push([...base, ci + 1,
      hasData ? (n.complex_scores || [])[ci] : "", hasData ? winning.has(JSON.stringify([...c].sort())) : "",
      g, hasData ? (n.gene_values || {})[g] : ""])));
  }
  return rows;
}

function evidenceTsv(reactionNodes) {
  const cell = v => String(v ?? "").replace(/[\t\r\n]+/g, " ");
  return [TSV_COLUMNS, ...evidenceRows(reactionNodes)].map(r => r.map(cell).join("\t")).join("\n") + "\n";
}

async function exportNetwork(kind) {
  if (!currentData || !currentData.panels.length) return;
  const panel = currentData.panels[currentIndex];
  const reactions = panel.nodes.filter(n => n.type === "reaction");
  try {
    if (kind === "json") {
      downloadBlob(exportName("json"), new Blob([JSON.stringify({
        model: currentData.model, task_list: currentData.task_list, task_id: currentData.task_id,
        task_description: currentData.task_description, sample: currentData.has_data ? currentData.tissue : null,
        route_id: panel.route_id, nodes: panel.nodes, edges: panel.edges }, null, 1)], { type: "application/json" }));
    } else if (kind === "evidence") {
      downloadBlob(exportName("evidence.tsv"), new Blob([evidenceTsv(reactions)], { type: "text/tab-separated-values" }));
    } else {
      const svg = standaloneSvg(host.querySelector("svg"), { background: kind === "png" ? getComputedColor("--surface") || "#fff" : null });
      if (kind === "svg") downloadBlob(exportName("svg"), new Blob([svg.text], { type: "image/svg+xml" }));
      else downloadBlob(exportName("png"), await svgToPng(svg, 3));
    }
  } catch (err) {
    setStatus(`Export failed: ${err.message}`);
  }
}

const exportMenu = document.getElementById("export-menu");
exportMenu.querySelectorAll("button[data-export]").forEach(b => b.addEventListener("click", () => {
  exportMenu.open = false;
  exportNetwork(b.dataset.export);
}));
document.addEventListener("click", ev => { if (!exportMenu.contains(ev.target)) exportMenu.open = false; });

function exportReaction(d, kind) {
  const base = fileStem(currentTaskList, "task" + currentData.task_id, d.reaction_id,
                        currentData.has_data ? currentData.tissue : null);
  if (kind === "tsv") {
    downloadBlob(`${base}.evidence.tsv`, new Blob([evidenceTsv([d])], { type: "text/tab-separated-values" }));
  } else {
    const svg = standaloneSvg(document.getElementById("bipartite-svg"), { fit: false, background: getComputedColor("--surface") || "#fff" });
    downloadBlob(`${base}.complexes.svg`, new Blob([svg.text], { type: "image/svg+xml" }));
  }
}

const EVIDENCE_LABEL = {
  supported: ["Supported", "#0ca30c", "#d2f0e3"],
  ambiguous: ["Ambiguous", "#c98500", "#fbe9c9"],
  no_evidence: ["No evidence", "var(--evidence-none)", "var(--evidence-none-bg)"],
  no_gpr: ["No GPR", "#898781", "#e1e0d9"],
};
const IO_LABEL = {
  input: ["Task input", "var(--met-input)", "var(--met-input-bg)"],
  output: ["Task output", "var(--met-output)", "var(--met-output-bg)"],
};

function renderDetail(d) {
  if (d.type === "reaction") renderReactionDetail(d);
  else renderMetaboliteDetail(d);
}

function renderMetaboliteDetail(d) {
  document.getElementById("detail-title").textContent = "Metabolite detail";
  const panel = document.getElementById("detail-panel-inner");
  const idBadges = [
    d.bigg_id ? `<span class="ec-badge">BiGG: ${d.bigg_id}</span>` : "",
    d.currency ? `<span class="ec-badge">currency (drawn per-reaction)</span>` : "",
  ].join("");
  const ioRow = IO_LABEL[d.io]
    ? `<div style="margin:0.5rem 0;"><span class="evidence-badge" style="background:${IO_LABEL[d.io][2]};color:${IO_LABEL[d.io][1]};">${IO_LABEL[d.io][0]}</span></div>`
    : "";
  panel.innerHTML = `
    <div style="font-weight:600;">${d.label}${idBadges}</div>
    <div style="color:var(--muted);font-size:0.75rem;margin:0.2rem 0 0.6rem;">${d.full_id} &middot; compartment [${d.compartment}]</div>
    <div class="field-label">Formula</div>
    <div>${d.formula || "<span class=\"placeholder\">not annotated</span>"}</div>
    ${ioRow}
    <div class="placeholder" style="margin-top:0.9rem;">${d.currency
      ? "This metabolite is a highly-connected \"currency\" carrier (e.g. ATP/ADP, NAD(H), H⁺, water, CO₂, Pi) and is drawn as a separate copy at every reaction that uses it in this route, to keep the diagram readable."
      : "Click a reaction (square) to see its equation, EC/BiGG ids, expression evidence, and gene/complex diagram."}</div>
  `;
}

function renderReactionDetail(d) {
  document.getElementById("detail-title").textContent = "Reaction detail";
  const panel = document.getElementById("detail-panel-inner");
  const hasData = !!(currentData && currentData.has_data);
  const idBadges = [
    d.ec_code ? `<span class="ec-badge">EC ${d.ec_code}</span>` : "",
    d.bigg_id ? `<span class="ec-badge">BiGG: ${d.bigg_id}</span>` : "",
  ].join("");
  const evidenceRow = hasData
    ? (() => {
        const [evLabel, evColor, evBg] = EVIDENCE_LABEL[d.evidence] || ["Unknown", "#898781", "#e1e0d9"];
        return `<div><span class="evidence-badge" style="background:${evBg};color:${evColor};">${evLabel}</span>
          <span style="color:var(--muted);font-size:0.78rem;margin-left:0.5rem;">score ${d.score.toFixed(2)}</span></div>`;
      })()
    : `<div class="placeholder">No sample loaded -- load expression data to see evidence/score here.</div>`;
  panel.innerHTML = `
    <div style="font-weight:600;">${d.label}${idBadges}</div>
    <div style="color:var(--muted);font-size:0.75rem;margin:0.2rem 0 0.6rem;">${d.reaction_id}</div>
    <div style="margin-bottom:0.6rem;"><code>${d.equation}</code></div>
    <div style="color:var(--muted);font-size:0.78rem;margin-bottom:0.4rem;">flux (this route) ${d.flux.toFixed(3)} mmol/gDW/h</div>
    ${evidenceRow}
    <div class="detail-actions">
      <button id="export-rxn-tsv" title="Evidence of this reaction: complexes, genes, values and scores">Export evidence (TSV)</button>
      ${(d.complexes||[]).length ? '<button id="export-rxn-svg" title="The complex/gene diagram below">Export diagram (SVG)</button>' : ""}
    </div>
    <div class="field-label">${(d.complexes||[]).length} candidate complex${(d.complexes||[]).length===1?'':'es'} (OR alternatives; genes within one AND-required)</div>
    <div id="bipartite-wrap"><svg id="bipartite-svg"></svg></div>
  `;
  document.getElementById("export-rxn-tsv").addEventListener("click", () => exportReaction(d, "tsv"));
  const svgBtn = document.getElementById("export-rxn-svg");
  if (svgBtn) svgBtn.addEventListener("click", () => exportReaction(d, "svg"));
  drawBipartite(d.complexes || [], hasData ? (d.complex_scores || []) : [], hasData ? (d.gene_values || {}) : {}, hasData ? (d.winning_complexes || []) : []);
}

function drawBipartite(complexes, complexScores, geneValues, winningComplexes) {
  if (!complexes.length) {
    document.getElementById("bipartite-wrap").outerHTML = `<div class="placeholder">No gene association (transport by diffusion or spontaneous).</div>`;
    return;
  }
  const genes = [];
  const geneSeen = new Set();
  complexes.forEach(c => c.forEach(g => { if (!geneSeen.has(g)) { geneSeen.add(g); genes.push(g); } }));

  const winningSet = new Set(winningComplexes.map(c => JSON.stringify([...c].sort())));
  const isWinning = c => winningSet.has(JSON.stringify([...c].sort()));

  const maxVal = Math.max(1e-6, ...genes.map(g => geneValues[g] || 0));
  const geneScale = d3.scaleSequential().domain([0, Math.log1p(maxVal)])
    .interpolator(t => d3.interpolate(getComputedColor("--score-empty"), "#0ca30c")(t));

  const rowH = 24;
  const height = Math.max(complexes.length, genes.length) * rowH + 20;
  const complexX = 140, geneX = 300, width = 460;

  const svg = d3.select("#bipartite-svg")
    .attr("viewBox", [0, 0, width, height])
    .attr("width", width).attr("height", height);
  const complexY = complexes.map((_, i) => 16 + i * rowH);
  const geneY = genes.map((_, i) => 16 + i * rowH);
  const geneIndex = Object.fromEntries(genes.map((g, i) => [g, i]));

  complexes.forEach((c, ci) => c.forEach(g => {
    svg.append("line")
      .attr("x1", complexX + 8).attr("y1", complexY[ci])
      .attr("x2", geneX - 8).attr("y2", geneY[geneIndex[g]])
      .attr("stroke", isWinning(c) ? "#0ca30c" : "var(--hairline)")
      .attr("stroke-width", isWinning(c) ? 2 : 1.2);
  }));

  complexes.forEach((c, ci) => {
    svg.append("circle").attr("cx", complexX).attr("cy", complexY[ci]).attr("r", 7)
      .style("fill", "var(--complex-node)")
      .attr("stroke", isWinning(c) ? "#0ca30c" : "none").attr("stroke-width", 2.5);
    const scoreText = complexScores[ci] !== undefined ? ` — ${complexScores[ci].toFixed(2)}` : "";
    svg.append("text").attr("x", complexX - 14).attr("y", complexY[ci] + 3).attr("text-anchor", "end")
      .attr("font-size", 10).style("fill", "var(--ink)").text(`Complex ${ci + 1}${scoreText}${isWinning(c) ? " ✓" : ""}`);
  });

  genes.forEach((g, gi) => {
    const val = geneValues[g] || 0;
    svg.append("circle").attr("cx", geneX).attr("cy", geneY[gi]).attr("r", 6)
      .style("fill", geneScale(Math.log1p(val)));
    svg.append("text").attr("x", geneX + 12).attr("y", geneY[gi] + 3)
      .attr("font-size", 9).style("fill", "var(--ink)").text(`${g} (${val.toFixed(1)})`);
  });
}

const host = document.getElementById("panel-host");
const titleText = document.getElementById("panel-title-text");
const prevBtn = document.getElementById("prev-btn");
const nextBtn = document.getElementById("next-btn");
const orientBtn = document.getElementById("orient-btn");

// The network fills the width of its panel (so hiding the side panels gives it
// room) and the rest of the window's height, but never less than a usable minimum.
function availableNetworkHeight() {
  const top = host.getBoundingClientRect().top + window.scrollY;
  return Math.max(500, Math.round(window.innerHeight - top - 24));
}

let resizeTimer = null;
window.addEventListener("resize", () => {
  clearTimeout(resizeTimer);
  resizeTimer = setTimeout(() => { if (currentData) showPanel(currentIndex); }, 150);
});

function showPanel(index) {
  if (!currentData || !currentData.panels.length) return;
  currentIndex = (index + currentData.panels.length) % currentData.panels.length;
  const panel = currentData.panels[currentIndex];
  host.innerHTML = "";
  const nRxns = panel.nodes.filter(n => n.type === "reaction").length;
  const truncNote = currentData.truncated ? `, showing first ${currentData.panels.length} of ${currentData.n_routes_total}` : "";
  const countLabel = currentData.has_data
    ? `(${currentIndex + 1} of ${currentData.panels.length} tied)`
    : `(variant ${currentIndex + 1} of ${currentData.panels.length}${truncNote})`;
  titleText.textContent = currentData.panels.length > 1
    ? `Task ${currentData.task_id} — Route ${panel.route_id} ${countLabel} — ${nRxns} reactions`
    : `Task ${currentData.task_id} — Route ${panel.route_id} — ${nRxns} reactions`;
  renderPanel(host, panel, host.clientWidth || 1000, availableNetworkHeight(), currentData.has_data);
  exportMenu.setAttribute("aria-disabled", "false");
  prevBtn.style.visibility = currentData.panels.length > 1 ? "visible" : "hidden";
  nextBtn.style.visibility = currentData.panels.length > 1 ? "visible" : "hidden";
}

prevBtn.addEventListener("click", () => showPanel(currentIndex - 1));
nextBtn.addEventListener("click", () => showPanel(currentIndex + 1));
orientBtn.addEventListener("click", () => {
  orientation = orientation === "TB" ? "LR" : "TB";
  orientBtn.textContent = orientation === "TB" ? "↔ Horizontal" : "↕ Vertical";
  if (currentData) showPanel(currentIndex);
});

// ---------- show / hide panels ----------

const VIEW_PARTS = ["controls", "analysis", "legend", "tasks", "detail"];
const VIEW_KEYS = { c: "controls", a: "analysis", l: "legend", t: "tasks", d: "detail" };
const hiddenParts = new Set();

try { (JSON.parse(localStorage.getItem("mteapy.hidden") || "[]")).forEach(p => VIEW_PARTS.includes(p) && hiddenParts.add(p)); }
catch (e) { /* storage unavailable or corrupt: start with everything shown */ }

function applyView() {
  VIEW_PARTS.forEach(p => document.querySelector(".wrap").classList.toggle(`hide-${p}`, hiddenParts.has(p)));
  document.querySelectorAll(".view-toggles button").forEach(b => b.classList.toggle("on", !hiddenParts.has(b.dataset.view)));
  try { localStorage.setItem("mteapy.hidden", JSON.stringify([...hiddenParts])); } catch (e) { /* ignore */ }
  if (currentData) showPanel(currentIndex);   // re-fit the network to the new width
}

function toggleView(part) {
  if (hiddenParts.has(part)) hiddenParts.delete(part); else hiddenParts.add(part);
  applyView();
}

document.querySelectorAll("[data-view]").forEach(b => b.addEventListener("click", () => toggleView(b.dataset.view)));
document.addEventListener("keydown", ev => {
  if (ev.ctrlKey || ev.metaKey || ev.altKey) return;
  if (/^(INPUT|SELECT|TEXTAREA)$/.test(document.activeElement.tagName)) return;
  const part = VIEW_KEYS[ev.key.toLowerCase()];
  if (part) toggleView(part);
});
applyView();

loadMethods().then(updateRunButton);
loadModels();
