/**
 * netgraph-viz -- this webapp's generic dagre+D3 renderer for {nodes, edges}
 * bipartite graphs, with no assumptions about what a node or edge *means*.
 * It knows how to lay out, zoom/pan, draw, and dispatch clicks on a graph;
 * app.js supplies plain callbacks ("decorators", in the same spirit as
 * cobra-netgraph's reaction_decorator/metabolite_decorator on the Python
 * side) that turn each node/edge into a shape, a color, a label, and a click
 * handler. This file never reads a domain-specific field itself -- if you
 * need a new field to affect rendering, that logic belongs in app.js's
 * callbacks, not here.
 *
 * Kept as its own file (rather than folded into app.js) purely to keep the
 * generic layout/draw mechanics visually separate from mteapy's own
 * evidence/score/flux-driven styling -- not because it's meant to be reused
 * outside this webapp. An earlier version of this project spun this out as
 * a standalone cross-project package; it was folded back in because a
 * future gap-filling/reconstruction visualization tool would need a
 * different-enough interaction model (curating a large draft network vs.
 * inspecting one small enumerated route) that a shared renderer would end
 * up either bloated with options to please both, or fought against by one
 * of them. cobra-netgraph (the graph-construction layer this consumes
 * server-side) stayed a shared Python package because building the same
 * bipartite graph from a COBRA model genuinely is one shared problem.
 *
 * Graph shape expected: {"nodes": [{"id": ..., ...}], "edges": [{"source":
 * id, "target": id, ...}]} -- the same shape cobra_netgraph.bipartite's
 * build_bipartite_graph produces.
 *
 * Depends on d3 (v7) and dagre (@dagrejs/dagre) being loaded globally before
 * this script -- see index.html.
 */
(function (global) {
  "use strict";

  const DEFAULT_LAYOUT_CONFIG = {
    TB: { nodesep: 46, ranksep: 90, marginx: 30, marginy: 20 },
    LR: { nodesep: 46, ranksep: 80, marginx: 20, marginy: 30 },
  };

  /** Truncates a label to at most n characters, adding an ellipsis. */
  function truncateLabel(s, n) {
    return s.length > n ? s.slice(0, n - 1) + "…" : s;
  }

  /**
   * Runs a dagre hierarchical layout over a {nodes, edges} graph.
   *
   * @param graph {nodes: [{id, ...}], edges: [{source, target, ...}]}
   * @param options.orientation "TB" | "LR" (default "TB")
   * @param options.nodeSize (node) -> number -- the node's layout box side
   * @param options.layoutConfig per-orientation dagre graph options, keyed
   *   like DEFAULT_LAYOUT_CONFIG (merged over the default for the chosen
   *   orientation)
   * @returns {nodes: [node with x/y], links: [edge with resolved source/
   *   target node objects and dagre-computed `points`], size: {width,height}}
   */
  function layout(graph, options) {
    const { orientation = "TB", nodeSize, layoutConfig } = options || {};
    const config = { ...DEFAULT_LAYOUT_CONFIG[orientation], ...(layoutConfig && layoutConfig[orientation]) };

    const g = new dagre.graphlib.Graph();
    g.setGraph({ rankdir: orientation, ...config });
    g.setDefaultEdgeLabel(() => ({}));
    graph.nodes.forEach((d) => {
      const size = nodeSize(d) + 6; // small fixed breathing room around the drawn shape
      g.setNode(d.id, { width: size, height: size });
    });
    graph.edges.forEach((d) => g.setEdge(d.source, d.target));
    dagre.layout(g);

    const nodeById = Object.fromEntries(graph.nodes.map((d) => [d.id, { ...d }]));
    g.nodes().forEach((id) => {
      const gn = g.node(id);
      nodeById[id].x = gn.x;
      nodeById[id].y = gn.y;
    });
    const links = graph.edges.map((d) => ({
      ...d,
      source: nodeById[d.source],
      target: nodeById[d.target],
      points: g.edge(d.source, d.target).points,
    }));
    return { nodes: Object.values(nodeById), links, size: g.graph() };
  }

  /** Trims a dagre point path by `startOffset`/`endOffset` px at each end,
   * so an edge stops short of the node shapes it connects (room for arrowheads). */
  function clipPoints(points, startOffset, endOffset) {
    const pts = points.map((p) => ({ x: p.x, y: p.y }));
    const trimEnd = (arr, offset) => {
      let remaining = offset;
      while (arr.length > 1 && remaining > 0) {
        const a = arr[arr.length - 1];
        const b = arr[arr.length - 2];
        const dx = a.x - b.x;
        const dy = a.y - b.y;
        const segLen = Math.sqrt(dx * dx + dy * dy) || 1;
        if (segLen > remaining) {
          arr[arr.length - 1] = { x: a.x - (dx / segLen) * remaining, y: a.y - (dy / segLen) * remaining };
          remaining = 0;
        } else {
          arr.pop();
          remaining -= segLen;
        }
      }
      return arr;
    };
    return trimEnd(trimEnd(pts, endOffset).reverse(), startOffset).reverse();
  }

  /**
   * Lays out and renders a {nodes, edges} graph as an interactive
   * (zoom/pan, click-to-detail) SVG inside `container`.
   *
   * Every visual decision -- shape, size, fill, stroke, label, tooltip,
   * click behavior, edge color/width -- is delegated to callbacks in
   * `options`; netgraph-viz reads no domain-specific field itself. A
   * callback that needs graph-wide context (e.g. a color scale spanning all
   * node scores) should be built by the caller from `graph` before calling
   * render(), then closed over.
   *
   * Required options: nodeSize(node), nodeShape(node) -> "rect"|"circle",
   * nodeFill(node), edgeColor(edge, sourceNode, targetNode), edgeWidth(edge).
   *
   * Optional options: orientation ("TB", default), viewWidth (1000),
   * viewHeight (800), layoutConfig, nodeStroke(node) (default "none"),
   * nodeLabel(node), labelOffset(node) (default nodeSize(node)/2 + 13),
   * nodeLabelStyle(node) -> {fontSize, fontWeight, fill}, nodeTooltip(node),
   * onNodeClick(node), zoomExtent ([0.3, 5]).
   *
   * @returns the same {nodes, links, size} layout() produced, in case the
   *   caller wants it (e.g. to report the rendered graph's reaction count).
   */
  function render(container, graph, options) {
    const {
      orientation = "TB",
      viewWidth = 1000,
      viewHeight = 800,
      layoutConfig,
      nodeSize,
      nodeShape,
      nodeFill,
      nodeStroke = () => "none",
      edgeColor,
      edgeWidth,
      nodeLabel,
      labelOffset,
      nodeLabelStyle = () => ({}),
      nodeTooltip,
      onNodeClick,
      zoomExtent = [0.3, 5],
    } = options;

    if (!nodeSize || !nodeShape || !nodeFill || !edgeColor || !edgeWidth) {
      throw new Error("NetGraphViz.render requires nodeSize, nodeShape, nodeFill, edgeColor and edgeWidth callbacks");
    }

    const { nodes, links, size } = layout(graph, { orientation, nodeSize, layoutConfig });
    const width = Math.max(viewWidth, size.width + 40);
    const height = Math.max(viewHeight, size.height + 40);

    container.innerHTML = "";
    const svg = d3
      .select(container)
      .append("svg")
      .attr("viewBox", [0, 0, width, height])
      .style("width", "100%")
      .style("height", viewHeight + "px");

    const defs = svg.append("defs");
    const g = svg.append("g");
    svg.call(d3.zoom().scaleExtent(zoomExtent).on("zoom", (ev) => g.attr("transform", ev.transform)));

    const lineGen = d3.line().x((p) => p.x).y((p) => p.y).curve(d3.curveLinear);
    const uid = `ngv-${Math.random().toString(36).slice(2)}`;

    g.append("g")
      .attr("fill", "none")
      .attr("stroke-opacity", 0.85)
      .selectAll("path")
      .data(links)
      .join("path")
      .attr("stroke-width", (d) => edgeWidth(d))
      .attr("d", (d) =>
        lineGen(clipPoints(d.points, nodeSize(d.source) / 2 + 3, nodeSize(d.target) / 2 + 5))
      )
      .each(function (d, i) {
        const color = edgeColor(d, d.source, d.target);
        const markerId = `${uid}-arrow-${i}`;
        defs
          .append("marker")
          .attr("id", markerId)
          .attr("viewBox", "0 0 8 8")
          .attr("refX", 7)
          .attr("refY", 4)
          .attr("markerWidth", 6)
          .attr("markerHeight", 6)
          .attr("markerUnits", "userSpaceOnUse")
          .attr("orient", "auto-start-reverse")
          .append("path")
          .attr("d", "M0,0 L8,4 L0,8 Z")
          .attr("fill", color);
        d3.select(this).attr("stroke", color).attr("marker-end", `url(#${markerId})`);
      });

    const node = g
      .append("g")
      .selectAll("g")
      .data(nodes)
      .join("g")
      .style("cursor", onNodeClick ? "pointer" : "default")
      .attr("transform", (d) => `translate(${d.x},${d.y})`);

    node.each(function (d) {
      const sel = d3.select(this);
      const shape = nodeShape(d);
      const s = nodeSize(d);
      if (shape === "circle") {
        sel
          .append("circle")
          .attr("r", s / 2)
          .style("fill", nodeFill(d))
          .attr("stroke", nodeStroke(d))
          .attr("stroke-width", 1.6);
      } else {
        sel
          .append("rect")
          .attr("width", s)
          .attr("height", s)
          .attr("x", -s / 2)
          .attr("y", -s / 2)
          .attr("rx", 5)
          .style("fill", nodeFill(d))
          .attr("stroke", nodeStroke(d))
          .attr("stroke-width", 1);
      }
    });

    if (nodeLabel) {
      const label = node
        .append("text")
        .text((d) => nodeLabel(d))
        .attr("text-anchor", "middle")
        .attr("dy", (d) => (labelOffset ? labelOffset(d) : nodeSize(d) / 2 + 13))
        .attr("font-size", (d) => nodeLabelStyle(d).fontSize || 9.5)
        .attr("font-weight", (d) => nodeLabelStyle(d).fontWeight || 400)
        .style("fill", (d) => nodeLabelStyle(d).fill || "var(--ink)");

      label.each(function () {
        const bbox = this.getBBox();
        d3.select(this.parentNode)
          .insert("rect", "text")
          .attr("x", bbox.x - 2)
          .attr("y", bbox.y - 1)
          .attr("width", bbox.width + 4)
          .attr("height", bbox.height + 2)
          .attr("fill", "var(--surface)")
          .attr("opacity", 0.88);
      });
    }

    if (nodeTooltip) {
      node.append("title").text((d) => nodeTooltip(d));
    }

    if (onNodeClick) {
      node.on("click", (event, d) => onNodeClick(d));
    }

    return { nodes, links, size };
  }

  global.NetGraphViz = { render, layout, clipPoints, truncateLabel };
})(typeof window !== "undefined" ? window : globalThis);
