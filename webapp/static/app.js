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
  datasetId = window.GTEX_DATASET_ID;
  const res = await fetch(`/api/datasets/${datasetId}`);
  const info = await res.json();
  clearRun();
  populateSampleSelect(info.samples);
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
  datasetId = info.dataset_id;
  clearRun();
  populateSampleSelect(info.samples);
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

function renderParams() {
  const m = currentMethod();
  paramsHost.innerHTML = "";
  document.getElementById("method-desc").textContent = m ? m.description : "";
  if (!m) return;
  m.params.forEach(p => {
    const wrap = document.createElement("label");
    wrap.className = "param"; wrap.title = p.help;
    wrap.append(`${p.label}: `);
    let input;
    if (p.type === "choice") {
      input = document.createElement("select");
      p.choices.forEach(c => { const o = document.createElement("option"); o.value = o.textContent = c; input.appendChild(o); });
      input.value = p.default;
    } else if (p.type === "bool") {
      input = document.createElement("input"); input.type = "checkbox"; input.checked = !!p.default;
    } else {
      input = document.createElement("input"); input.type = "number"; input.value = p.default;
      input.step = p.type === "int" ? "1" : "any";
      if (p.min !== null) input.min = p.min;
      if (p.max !== null) input.max = p.max;
    }
    input.dataset.param = p.name; input.dataset.type = p.type;
    wrap.dataset.when = JSON.stringify(p.when || []);
    wrap.appendChild(input);
    paramsHost.appendChild(wrap);
  });
  updateParamVisibility();
}

// A parameter that only matters for some settings of others (e.g. the global
// threshold value) is shown only while those hold.
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
    out[el.dataset.param] = el.dataset.type === "bool" ? el.checked
      : (el.dataset.type === "int" || el.dataset.type === "float") ? Number(el.value) : el.value;
  });
  return out;
}

function updateRunButton() {
  document.getElementById("save-btn").disabled = document.getElementById("tsv-btn").disabled = !currentRun;
  runBtn.disabled = !(datasetId && currentModel && currentTaskList && currentMethod() && !activeRunId);
  runBtn.title = datasetId ? "Score every task of the selected list against every sample"
                           : "Load an expression dataset first";
}

methodSelect.addEventListener("change", renderParams);
paramsHost.addEventListener("change", updateParamVisibility);

function clearRun(note) {
  currentRun = null;
  allTasks = allTasks.map(t => ({ ...t, score: null, is_tied: false, is_complete: false, winning_route_ids: [] }));
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
    return r ? { ...t, score: r.scores[idx], is_tied: r.tied[idx], is_complete: r.complete[idx] }
             : { ...t, score: null, is_tied: false, is_complete: false };
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
    setRunStatus(info.status === "queued" ? "Queued..." : `Scoring ${info.done}/${info.total} samples...`);
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
    if (el.dataset.type === "bool") el.checked = !!v; else el.value = v;
  });
  updateParamVisibility();
  datasetId = info.dataset_id;
  const data = await (await fetch(`/api/runs/${info.run_id}/results`)).json();
  currentRun = { run_id: info.run_id, config: cfg, samples: data.samples,
                 byTask: Object.fromEntries(data.tasks.map(t => [t.task_id, t])) };
  populateSampleSelect(data.samples);
  const p = Object.entries(cfg.params).map(([k, v]) => `${k}=${v}`).join(", ");
  setRunStatus(`Loaded ${cfg.method} (${p}), saved ${cfg.saved_at || "earlier"}: ${data.tasks.length} tasks x ` +
               `${data.samples.length} samples.` + (datasetId ? "" : " (No expression to re-run: CellFie results are loaded as scores only.)"));
  applyRunScores();
  updateRunButton();
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

  taskTbody.innerHTML = "";
  if (!rows.length) {
    taskTbody.innerHTML = `<tr><td colspan="4" class="placeholder">No matching tasks.</td></tr>`;
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
    `;
    tr.addEventListener("click", () => selectTask(t.task_id));
    taskTbody.appendChild(tr);
  }
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

function scaleFor(panel) {
  const scores = panel.nodes
    .filter(n => n.type === "reaction" && (n.evidence === "supported" || n.evidence === "ambiguous"))
    .map(n => Math.log1p(n.score || 0));
  const max = Math.max(1e-6, ...scores);
  return d3.scaleSequential().domain([0, max]).interpolator(t => d3.interpolate(
    getComputedColor("--score-empty"), getComputedColor("--accent")
  )(t));
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
  return colorScale(Math.log1p(d.score || 0));
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
    <div class="field-label">${(d.complexes||[]).length} candidate complex${(d.complexes||[]).length===1?'':'es'} (OR alternatives; genes within one AND-required)</div>
    <div id="bipartite-wrap"><svg id="bipartite-svg"></svg></div>
  `;
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
