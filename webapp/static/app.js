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
let allTasks = [];      // [{source, task_id, description, n_routes, score, is_tied, is_complete, winning_route_ids}]
let sortKey = "task_id";
let sortAsc = true;
let selectedKey = null; // `${source}::${task_id}`
let currentData = null; // {task_id, task_description, tissue, has_data, panels: [...]}
let currentIndex = 0;
let orientation = "TB";

function rowKey(t) { return `${t.source}::${t.task_id}`; }

const statusMsg = document.getElementById("status-msg");
const sampleSelect = document.getElementById("sample-select");
const taskTbody = document.getElementById("task-tbody");
const taskSearch = document.getElementById("task-search");

function setStatus(text) { statusMsg.textContent = text; }

// ---------- bootstrapping ----------

async function loadTaskList() {
  const res = await fetch("/api/tasks");
  const data = await res.json();
  allTasks = data.tasks.map(t => ({ ...t, score: null, is_tied: false, is_complete: false, winning_route_ids: [] }));
  window.GTEX_DATASET_ID = data.gtex_dataset_id;
  renderTaskTable();
}

document.getElementById("load-gtex-btn").addEventListener("click", async () => {
  if (!window.GTEX_DATASET_ID) { setStatus("No bundled GTEx dataset available on this server."); return; }
  datasetId = window.GTEX_DATASET_ID;
  const res = await fetch(`/api/datasets/${datasetId}`);
  const info = await res.json();
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
}

sampleSelect.addEventListener("change", async () => {
  if (!datasetId || !sampleSelect.value) return;
  currentSample = sampleSelect.value;
  setStatus(`Scoring all tasks against ${currentSample}...`);
  const res = await fetch(`/api/datasets/${datasetId}/scores?sample=${encodeURIComponent(currentSample)}`);
  if (!res.ok) { setStatus(`Scoring failed: ${(await res.json()).detail || res.statusText}`); return; }
  const data = await res.json();
  const byKey = Object.fromEntries(data.results.map(r => [rowKey(r), r]));
  allTasks = allTasks.map(t => ({ ...t, ...(byKey[rowKey(t)] || {}) }));
  renderTaskTable();
  setStatus(`Scored ${data.results.length} tasks against ${currentSample}.`);
  if (selectedKey) selectTask(selectedKey);
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
    taskTbody.innerHTML = `<tr><td colspan="5" class="placeholder">No matching tasks.</td></tr>`;
    return;
  }
  for (const t of rows) {
    const tr = document.createElement("tr");
    tr.className = "task-row" + (rowKey(t) === selectedKey ? " selected" : "");
    const scoreCell = t.score === null
      ? `<span class="placeholder">&mdash;</span>`
      : `<span class="score-pill" style="background:${scorePillColor(t)};">${t.score.toFixed(2)}</span>${t.is_tied ? " tied" : ""}`;
    tr.innerHTML = `
      <td>${t.source}</td>
      <td>${t.task_id}</td>
      <td>${t.description || ""}</td>
      <td>${t.n_routes}</td>
      <td>${scoreCell}</td>
    `;
    tr.addEventListener("click", () => selectTask(t.source, t.task_id));
    taskTbody.appendChild(tr);
  }
}

function scorePillColor(t) {
  if (t.score === null) return "var(--score-empty)";
  if (t.score === 0) return "var(--evidence-none-bg)";
  return t.is_complete ? "#d2f0e3" : "#fbe9c9";
}

async function selectTask(source, taskId) {
  selectedKey = `${source}::${taskId}`;
  renderTaskTable();
  document.getElementById("panel-title-text").textContent = "Loading...";
  const hasSample = datasetId && currentSample;
  const qs = hasSample ? `?dataset_id=${encodeURIComponent(datasetId)}&sample=${encodeURIComponent(currentSample)}` : "";
  const res = await fetch(`/api/tasks/${encodeURIComponent(source)}/${encodeURIComponent(taskId)}/network${qs}`);
  if (!res.ok) {
    document.getElementById("panel-title-text").textContent = `Error: ${(await res.json()).detail || res.statusText}`;
    return;
  }
  currentData = await res.json();
  currentIndex = 0;
  document.getElementById("task-desc").textContent = currentData.task_description || "";
  setStatus(hasSample ? "" : "No sample loaded -- showing all route variants, unscored.");
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
  renderPanel(host, panel, 1000, 800, currentData.has_data);
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

loadTaskList();
