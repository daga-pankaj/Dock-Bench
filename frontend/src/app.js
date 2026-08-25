/* Docking Bench front end.
 *
 * State lives in `S`. The viewer is a single Mol* plugin holding the
 * receptor structure plus, once a result is selected, one structure per
 * docked pose. This file is bundled by `npm run build` (see build.mjs)
 * into ../app.js, which is what index.html actually loads.
 */

import { createPluginUI } from "molstar/lib/mol-plugin-ui/index";
import { renderReact18 } from "molstar/lib/mol-plugin-ui/react18";
import { DefaultPluginUISpec } from "molstar/lib/mol-plugin-ui/spec";
import { PluginCommands } from "molstar/lib/mol-plugin/commands";
import { PluginConfig } from "molstar/lib/mol-plugin/config";
import { StructureElement, StructureProperties } from "molstar/lib/mol-model/structure";
import { StructureSelectionQueries as Queries } from "molstar/lib/mol-plugin-state/helpers/structure-selection-query";
import { MolScriptBuilder as MS } from "molstar/lib/mol-script/language/builder";
import { Color } from "molstar/lib/mol-util/color";
import { Vec3 } from "molstar/lib/mol-math/linear-algebra";
import "molstar/lib/mol-util/polyfill";

const $ = (id) => document.getElementById(id);

const S = {
  receptor: null,      // {id, name, het_ligands, n_protein_atoms}
  pdbText: "",
  ligands: [],         // parsed ligand specs from the backend
  pickMode: "ligand",
  picks: [],           // [x,y,z] click points in "pick atoms" mode
  activeHet: null,
  job: null,
  poll: null,
  selectedLigand: null,
  poseReps: [],        // Mol* representation refs, one per pose currently shown
  activePose: 0,
};

let plugin = null;
let receptorStructure = null;   // StateObjectSelector for the loaded receptor Structure
let receptorComponentRefs = []; // components (+ their reps) built by styleReceptor, for teardown
let boxRefs = [];                // state refs for the docking-box wireframe structure
let picksRefs = [];              // state refs for the pick-point markers structure
let poseStructureRefs = [];      // top-level structure refs, one per pose (for teardown)
let spinning = false;

// ---------------------------------------------------------------------------
// small helpers
// ---------------------------------------------------------------------------

function toast(msg, bad = false) {
  const el = $("toast");
  el.textContent = msg;
  el.className = "toast" + (bad ? " bad" : "");
  clearTimeout(toast._t);
  toast._t = setTimeout(() => el.classList.add("hidden"), 4200);
}

async function api(path, opts = {}) {
  const res = await fetch(path, opts);
  if (!res.ok) {
    let detail = res.statusText;
    try { detail = (await res.json()).detail || detail; } catch (_) {}
    throw new Error(detail);
  }
  const type = res.headers.get("content-type") || "";
  return type.includes("json") ? res.json() : res.text();
}

const num = (id) => parseFloat($(id).value);

// Delete a batch of state objects (and everything built on top of them).
async function deleteRefs(refs) {
  if (!refs.length) return;
  const b = plugin.state.data.build();
  refs.forEach((r) => b.delete(r));
  await b.commit();
}

// ---------------------------------------------------------------------------
// viewer bootstrap
// ---------------------------------------------------------------------------

async function initViewer() {
  const defaultSpec = DefaultPluginUISpec();
  const spec = {
    ...defaultSpec,
    layout: { initial: { isExpanded: false, showControls: false, controlsDisplay: "reactive" } },
    components: {
      ...defaultSpec.components,
      controls: { ...defaultSpec.components?.controls, top: "none", bottom: "none", left: "none", right: "none" },
      remoteState: "none",
    },
    config: [
      [PluginConfig.Viewport.ShowExpand, false],
      [PluginConfig.Viewport.ShowControls, false],
      [PluginConfig.Viewport.ShowSettings, false],
      [PluginConfig.Viewport.ShowSelectionMode, false],
      [PluginConfig.Viewport.ShowAnimation, false],
    ],
  };

  plugin = await createPluginUI({ target: $("viewer"), spec, render: renderReact18 });

  await PluginCommands.Canvas3D.SetSettings(plugin, {
    settings: {
      renderer: { ...plugin.canvas3d.props.renderer, backgroundColor: Color(0x000000) },
    },
  });

  plugin.behaviors.interaction.click.subscribe(onPluginClick);
}

function onPluginClick(event) {
  try {
    const loci = event.current.loci;
    if (!loci || !StructureElement.Loci.is(loci) || StructureElement.Loci.size(loci) === 0) return;
    const loc = StructureElement.Location.create(loci.structure);
    StructureElement.Loci.getFirstLocation(loci, loc);
    if (!loc.unit) return;

    const p = Vec3.zero();
    loc.unit.conformation.position(loc.element, p);

    onAtomClick({
      resn: StructureProperties.atom.label_comp_id(loc),
      hetflag: StructureProperties.residue.group_PDB(loc) === "HETATM",
      chain: StructureProperties.chain.auth_asym_id(loc),
      resi: StructureProperties.residue.auth_seq_id(loc),
      x: p[0], y: p[1], z: p[2],
      atom: StructureProperties.atom.label_atom_id(loc),
    });
  } catch (err) {
    console.warn("click handling failed", err);
  }
}

// ---------------------------------------------------------------------------
// representations
// ---------------------------------------------------------------------------

// MolScript expression selecting one specific het group by chain/resname/resno
// (auth_* fields — the receptor is parsed straight off the PDB text, which
// backend/prep.py also reads by raw author-numbering columns).
function hetExpr(het) {
  const tests = {
    "residue-test": MS.core.logic.and([
      MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_comp_id(), het.resname]),
      MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_seq_id(), parseInt(het.resseq, 10)]),
    ]),
  };
  if (het.chain && het.chain !== "_") {
    tests["chain-test"] = MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_asym_id(), het.chain]);
  }
  return MS.struct.generator.atomGroups(tests);
}

function zoomToExpression(expr, durationMs = 300) {
  if (!plugin || !receptorStructure?.data) return;
  try {
    const loci = StructureElement.Loci.fromExpression(receptorStructure.data, expr);
    plugin.managers.camera.focusLoci(loci, { durationMs });
  } catch (err) {
    console.warn("zoom failed", err);
  }
}

async function addRep(componentRef, type, color, colorParams, typeParams) {
  if (!componentRef) return null;
  return plugin.builders.structure.representation.addRepresentation(componentRef, {
    type, color, colorParams, typeParams,
  });
}

async function styleReceptor() {
  if (!plugin || !receptorStructure) return;
  await deleteRefs(receptorComponentRefs);
  receptorComponentRefs = [];

  const mode = $("styleSel").value;
  const track = (ref) => { if (ref) receptorComponentRefs.push(ref); };

  const polymerComp = await plugin.builders.structure.tryCreateComponentFromExpression(
    receptorStructure, Queries.polymer.expression, "polymer");
  track(polymerComp);

  if (mode === "sticks") {
    track(await addRep(polymerComp, "ball-and-stick", "element-symbol", undefined, { sizeFactor: 0.2 }));
  } else if (mode === "line") {
    track(await addRep(polymerComp, "line", "element-symbol"));
  } else {
    track(await addRep(polymerComp, "cartoon", "chain-id", undefined, { opacity: 0.92 }));
  }

  if (mode === "cartoon+sticks") {
    const c = boxCenter();
    if (c) {
      const r = Math.max(...boxSize()) / 2;
      const pocketExpr = MS.struct.modifier.exceptBy({
        0: MS.struct.modifier.includeSurroundings({
          0: Queries.ligand.expression, radius: r, "as-whole-residues": true,
        }),
        by: Queries.ligand.expression,
      });
      const pocketComp = await plugin.builders.structure.tryCreateComponentFromExpression(
        receptorStructure, pocketExpr, "pocket");
      track(pocketComp);
      track(await addRep(pocketComp, "ball-and-stick", "element-symbol", undefined, { sizeFactor: 0.2 }));
    }
  }

  // bound ligands / cofactors — always sticks so they stay clickable
  const ligandExpr = MS.struct.combinator.merge([Queries.ligand.expression, Queries.ion.expression]);
  const ligandComp = await plugin.builders.structure.tryCreateComponentFromExpression(
    receptorStructure, ligandExpr, "ligands");
  track(ligandComp);
  track(await addRep(ligandComp, "ball-and-stick", "uniform", { value: Color(0x4caf50) }, { sizeFactor: 0.25 }));

  if ($("showWaters").checked) {
    const waterComp = await plugin.builders.structure.tryCreateComponentFromExpression(
      receptorStructure, Queries.water.expression, "waters");
    track(waterComp);
    track(await addRep(waterComp, "spacefill", "uniform", { value: Color(0x5aa9ff) }, { sizeFactor: 0.25 }));
  }

  if ($("showSurface").checked) {
    track(await addRep(polymerComp, "molecular-surface", "uniform", { value: Color(0x3a4552) }, { alpha: 0.65 }));
  }

  if (S.activeHet) {
    const activeComp = await plugin.builders.structure.tryCreateComponentFromExpression(
      receptorStructure, hetExpr(S.activeHet), "active-site");
    track(activeComp);
    track(await addRep(activeComp, "ball-and-stick", "uniform", { value: Color(0xf5a623) }, { sizeFactor: 0.28 }));
  }
}

// ---------------------------------------------------------------------------
// docking box + pick markers (drawn as synthetic tiny structures — dummy
// carbon atoms so Mol*'s normal distance-based bonding renders the edges of
// the box as a wireframe; each lives in its own structure so it can never
// bond across into the receptor or into the other synthetic structure)
// ---------------------------------------------------------------------------

function pdbAtomLine(serial, name, resName, chain, resSeq, x, y, z, element) {
  const rj = (v, n) => String(v).padStart(n);
  const lj = (v, n) => String(v).padEnd(n);
  return (
    "HETATM" + rj(serial, 5) + " " + lj(name, 4) + " " + lj(resName, 3) + " " +
    (chain || "Z") + rj(resSeq, 4) + "    " +
    rj(x.toFixed(3), 8) + rj(y.toFixed(3), 8) + rj(z.toFixed(3), 8) +
    rj("1.00", 6) + rj("0.00", 6) + "          " + rj(element, 2)
  );
}

// Points spaced ~1.2 A apart along the segment so Mol*'s default bonding
// (distance based) chains them into a visible line.
function edgePoints(a, b, spacing = 1.2) {
  const d = Math.hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
  const n = Math.max(1, Math.round(d / spacing));
  const pts = [];
  for (let i = 0; i <= n; i++) {
    const t = i / n;
    pts.push([a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t]);
  }
  return pts;
}

async function loadSyntheticPdb(pdbText, resName) {
  const data = await plugin.builders.data.rawData({ data: pdbText, label: resName });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, "pdb");
  const model = await plugin.builders.structure.createModel(trajectory);
  const structure = await plugin.builders.structure.createStructure(model);
  return { structure, refs: [data, trajectory, model, structure] };
}

async function drawBox() {
  if (!plugin) return;
  await deleteRefs(boxRefs);
  boxRefs = [];

  const c = boxCenter();
  const s = boxSize();
  if (!c || !s.every((v) => Number.isFinite(v) && v > 0)) return;

  const half = s.map((v) => v / 2);
  const corner = (bit) => [0, 1, 2].map((i) => c[i] + (bit & (1 << i) ? half[i] : -half[i]));
  const lines = [];
  let serial = 1;
  for (let i = 0; i < 8; i++) {
    for (let bit = 0; bit < 3; bit++) {
      const j = i ^ (1 << bit);
      if (j <= i) continue;
      for (const p of edgePoints(corner(i), corner(j))) {
        lines.push(pdbAtomLine(serial, "C", "BOX", "Z", 1, p[0], p[1], p[2], "C"));
        serial++;
      }
    }
  }
  const { structure, refs } = await loadSyntheticPdb(lines.join("\n") + "\nEND\n", "dockbox");
  const rep = await addRep(structure, "line", "uniform", { value: Color(0x4da3ff) }, { lineSizeAttenuation: false });
  boxRefs = [...refs, rep].filter(Boolean);

  const vol = (s[0] * s[1] * s[2]).toFixed(0);
  $("siteInfo").textContent = `box ${s.map((v) => v.toFixed(1)).join(" x ")} A  (${vol} A^3)`;
}

async function rebuildPicks() {
  if (!plugin) return;
  await deleteRefs(picksRefs);
  picksRefs = [];
  if (!S.picks.length) return;

  const lines = S.picks.map((p, i) => pdbAtomLine(i + 1, "C", "PIK", "Z", i + 1, p[0], p[1], p[2], "C"));
  const { structure, refs } = await loadSyntheticPdb(lines.join("\n") + "\nEND\n", "picks");
  const rep = await addRep(structure, "spacefill", "uniform", { value: Color(0xf0b528) }, { sizeFactor: 0.4 });
  picksRefs = [...refs, rep].filter(Boolean);
}

function setBox(center, size) {
  ["cx", "cy", "cz"].forEach((id, i) => ($(id).value = center[i].toFixed(2)));
  ["sx", "sy", "sz"].forEach((id, i) => ($(id).value = Math.max(8, size[i]).toFixed(1)));
  drawBox();
}

function boxCenter() {
  const c = [num("cx"), num("cy"), num("cz")];
  return c.every((v) => Number.isFinite(v)) ? c : null;
}

function boxSize() {
  return [num("sx"), num("sy"), num("sz")];
}

function clearPicks() {
  S.picks = [];
  S.activeHet = null;
  rebuildPicks();
  document.querySelectorAll(".het").forEach((e) => e.classList.remove("active"));
  styleReceptor();
  drawBox();
}

// ---------------------------------------------------------------------------
// site picking
// ---------------------------------------------------------------------------

function boxFromPoints(points) {
  const pad = num("padding") || 0;
  const mins = [0, 1, 2].map((i) => Math.min(...points.map((p) => p[i])));
  const maxs = [0, 1, 2].map((i) => Math.max(...points.map((p) => p[i])));
  const center = [0, 1, 2].map((i) => (mins[i] + maxs[i]) / 2);
  const size = [0, 1, 2].map((i) => maxs[i] - mins[i] + 2 * pad);
  return [center, size];
}

function useHet(het) {
  S.activeHet = het;
  S.picks = [];
  rebuildPicks();
  const pad = num("padding") || 0;
  setBox(het.center, het.size.map((v) => v + 2 * pad));
  document.querySelectorAll(".het").forEach((e) =>
    e.classList.toggle("active", e.dataset.key === het.key));
  styleReceptor();
  zoomToExpression(hetExpr(het));
  toast(`site set from ${het.resname} ${het.chain}${het.resseq}`);
}

function onAtomClick(atom) {
  if (!atom) return;
  const isWater = ["HOH", "WAT", "DOD"].includes((atom.resn || "").toUpperCase());

  if (S.pickMode === "ligand" && atom.hetflag && !isWater) {
    const key = `${(atom.resn || "").toUpperCase()}:${atom.chain || "_"}:${atom.resi}`;
    const het = (S.receptor.het_ligands || []).find((h) => h.key === key);
    if (het) { useHet(het); return; }
  }

  // Anything else becomes a pick point; three or four picks define a site well.
  S.picks.push([atom.x, atom.y, atom.z]);
  rebuildPicks();
  S.activeHet = null;
  const [center, size] = boxFromPoints(S.picks);
  const minEdge = S.picks.length === 1 ? 20 : 14;
  setBox(center, size.map((v) => Math.max(v, minEdge)));
  toast(`${S.picks.length} pick(s) - ${atom.resn}${atom.resi} ${atom.atom}`);
}

// ---------------------------------------------------------------------------
// receptor loading
// ---------------------------------------------------------------------------

function renderHetList() {
  const box = $("hetList");
  box.innerHTML = "";
  const hets = S.receptor.het_ligands || [];
  if (!hets.length) {
    box.innerHTML = '<span class="muted small">no bound ligands found - use "pick atoms"</span>';
    return;
  }
  hets.forEach((h) => {
    const el = document.createElement("div");
    el.className = "het";
    el.dataset.key = h.key;
    el.innerHTML = `${h.resname} ${h.chain}${h.resseq}<span class="n">${h.n_atoms}</span>`;
    el.onclick = () => useHet(h);
    box.appendChild(el);
  });
}

async function loadReceptor(meta) {
  S.receptor = meta;
  S.pdbText = await api(`/api/receptor/${meta.id}/pdb`);
  $("receptorInfo").textContent =
    `${meta.name} - ${meta.n_protein_atoms} protein atoms - ${(meta.het_ligands || []).length} het groups`;

  clearPoses();
  await deleteRefs(boxRefs); boxRefs = [];
  await deleteRefs(picksRefs); picksRefs = [];
  if (receptorStructure) { await deleteRefs([receptorStructure]); receptorStructure = null; receptorComponentRefs = []; }
  S.picks = [];
  S.activeHet = null;

  const data = await plugin.builders.data.rawData({ data: S.pdbText, label: meta.name });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, "pdb");
  const model = await plugin.builders.structure.createModel(trajectory);
  receptorStructure = await plugin.builders.structure.createStructure(model);

  await styleReceptor();
  plugin.canvas3d?.requestCameraReset();
  renderHetList();

  // Biggest non-ion HET group is almost always the co-crystallised ligand.
  const best = (meta.het_ligands || []).find((h) => !h.is_ion);
  if (best) useHet(best);
  else drawBox();
}

// ---------------------------------------------------------------------------
// ligands
// ---------------------------------------------------------------------------

function showLigands(data) {
  S.ligands = data.ligands.filter((l) => l.valid);
  const bad = data.ligands.length - data.n_valid;
  $("ligInfo").textContent = `${data.n_valid} ligand(s) ready${bad ? `, ${bad} rejected` : ""}`;
  const box = $("ligPreview");
  box.innerHTML = "";
  data.ligands.slice(0, 500).forEach((l) => {
    const d = document.createElement("div");
    d.className = l.valid ? "" : "bad";
    d.textContent = `${l.name}  ${l.smiles || "(3D block)"}  - ${l.note}`;
    box.appendChild(d);
  });
  if (data.ligands.length > 500) {
    const d = document.createElement("div");
    d.className = "muted";
    d.textContent = `... ${data.ligands.length - 500} more`;
    box.appendChild(d);
  }
}

// ---------------------------------------------------------------------------
// jobs
// ---------------------------------------------------------------------------

async function runDocking() {
  if (!S.receptor) return toast("load a receptor first", true);
  if (!S.ligands.length) return toast("no ligands loaded - click 'Load from box'", true);
  const center = boxCenter();
  if (!center) return toast("define the active site first", true);

  const body = {
    receptor_id: S.receptor.id,
    center,
    size: boxSize(),
    ligands: S.ligands.map((l) => ({ name: l.name, smiles: l.smiles, molblock: l.molblock })),
    exhaustiveness: parseInt($("exh").value, 10),
    num_modes: parseInt($("modes").value, 10),
    energy_range: num("erange"),
    seed: parseInt($("seed").value, 10),
    scoring: $("scoring").value,
    engine: $("engineSel").value,
    use_gpu: $("deviceSel").value !== "cpu",
    keep_waters: $("keepWaters").checked,
    trim_distant_chains: $("trimChains").checked,
  };

  $("runBtn").disabled = true;
  try {
    const job = await api("/api/jobs", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
    });
    S.job = job;
    S.selectedLigand = null;
    clearPoses();
    $("cancelBtn").disabled = false;
    $("dlZip").href = `/api/jobs/${job.id}/download`;
    $("dlCsv").href = `/api/jobs/${job.id}/results.csv`;
    $("dlZip").classList.remove("disabled");
    $("dlCsv").classList.remove("disabled");
    renderJob(job);
    startPolling(job.id);
    toast(`job ${job.id} started - ${job.n_ligands} ligand(s)`);
  } catch (err) {
    toast(err.message, true);
    $("runBtn").disabled = false;
  }
}

function startPolling(jobId) {
  clearInterval(S.poll);
  S.poll = setInterval(async () => {
    try {
      const job = await api(`/api/jobs/${jobId}`);
      S.job = job;
      renderJob(job);
      if (job.status === "done" || job.status === "cancelled") {
        clearInterval(S.poll);
        $("runBtn").disabled = false;
        $("cancelBtn").disabled = true;
        const ok = job.ligands.filter((l) => l.status === "done").length;
        toast(`job ${job.status}: ${ok}/${job.n_ligands} docked`);
      }
    } catch (err) {
      clearInterval(S.poll);
      toast(err.message, true);
    }
  }, 1500);
}

function renderJob(job) {
  const pct = job.n_ligands ? (100 * job.n_finished) / job.n_ligands : 0;
  $("bar").style.width = `${pct}%`;
  $("runInfo").textContent = `${job.status} - ${job.n_finished}/${job.n_ligands}`;
  $("resultsEmpty").classList.add("hidden");

  const rows = [...job.ligands].sort((a, b) => {
    if (a.best_score === null && b.best_score === null) return a.name.localeCompare(b.name);
    if (a.best_score === null) return 1;
    if (b.best_score === null) return -1;
    return a.best_score - b.best_score;
  });

  const tbody = $("results").querySelector("tbody");
  tbody.innerHTML = "";
  rows.forEach((lig) => {
    const tr = document.createElement("tr");
    if (lig.name === S.selectedLigand) tr.className = "sel";
    const score = lig.best_score === null ? "-" : lig.best_score.toFixed(2);
    tr.innerHTML = `
      <td title="${lig.smiles || ""}">${lig.name}</td>
      <td><span class="tag ${lig.status}">${lig.status}</span></td>
      <td class="score">${score}</td>
      <td>${lig.n_poses || ""}</td>
      <td></td>`;
    if (lig.error) tr.querySelector("td:nth-child(2)").title = lig.error;

    const actions = tr.lastElementChild;
    if (lig.status === "done") {
      const view = document.createElement("button");
      view.className = "ghost";
      view.textContent = "view";
      view.onclick = () => showPoses(lig);
      const dl = document.createElement("a");
      dl.className = "btnlink ghost";
      dl.textContent = "sdf";
      dl.href = `/api/jobs/${job.id}/ligands/${encodeURIComponent(lig.name)}/poses.sdf`;
      actions.append(view, " ", dl);
      tr.style.cursor = "pointer";
      tr.onclick = (e) => { if (e.target.tagName !== "A" && e.target.tagName !== "BUTTON") showPoses(lig); };
    }
    tbody.appendChild(tr);
  });
}

// ---------------------------------------------------------------------------
// pose display
// ---------------------------------------------------------------------------

async function clearPoses() {
  await deleteRefs(poseStructureRefs);
  poseStructureRefs = [];
  S.poseReps = [];
  $("poseBar").classList.add("hidden");
}

async function showPoses(lig) {
  if (!S.job) return;
  try {
    const sdf = await api(`/api/jobs/${S.job.id}/ligands/${encodeURIComponent(lig.name)}/poses.sdf`);
    await clearPoses();
    S.selectedLigand = lig.name;

    const blocks = sdf.split("$$$$").map((b) => b.trim()).filter((b) => b.length > 10);
    for (const block of blocks) {
      const data = await plugin.builders.data.rawData({ data: block + "\n$$$$\n", label: lig.name });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, "sdf");
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      const rep = await addRep(structure, "ball-and-stick", "uniform", { value: Color(0x4dd0e1) }, { sizeFactor: 0.16 });
      poseStructureRefs.push(data, trajectory, model, structure);
      S.poseReps.push(rep);
    }

    const bar = $("poseBar");
    bar.innerHTML = `<span class="title">${lig.name} - poses:</span>`;
    S.poseReps.forEach((_, i) => {
      const b = document.createElement("button");
      const score = lig.scores[i];
      b.className = "pose";
      b.textContent = `${i + 1}${score !== undefined ? `  ${score.toFixed(2)}` : ""}`;
      b.onclick = () => selectPose(i);
      bar.appendChild(b);
    });
    const all = document.createElement("button");
    all.className = "pose";
    all.textContent = "all";
    all.onclick = () => selectPose(-1);
    bar.appendChild(all);
    bar.classList.remove("hidden");

    selectPose(0);
    renderJob(S.job);
  } catch (err) {
    toast(err.message, true);
  }
}

function selectPose(idx) {
  S.activePose = idx;
  S.poseReps.forEach((rep, i) => {
    if (!rep) return;
    plugin.state.data.updateCellState(rep.ref, { isHidden: !(idx === -1 || i === idx) });
  });
  document.querySelectorAll("#poseBar .pose").forEach((b, i) => {
    b.classList.toggle("active", i === (idx === -1 ? S.poseReps.length : idx));
  });
}

// ---------------------------------------------------------------------------
// camera helpers
// ---------------------------------------------------------------------------

function zoomBy(factor) {
  if (!plugin?.canvas3d) return;
  try {
    const cam = plugin.canvas3d.camera;
    const delta = Vec3.zero();
    Vec3.sub(delta, cam.state.position, cam.state.target);
    Vec3.scale(delta, delta, factor);
    const position = Vec3.zero();
    Vec3.add(position, cam.state.target, delta);
    cam.setState({ position }, 0);
    plugin.canvas3d.requestDraw();
  } catch (err) {
    console.warn("zoom failed", err);
  }
}

function toggleSpin() {
  if (!plugin?.canvas3d) return;
  spinning = !spinning;
  plugin.canvas3d.setProps({
    trackball: {
      animate: spinning
        ? { name: "spin", params: { speed: 0.4, axis: Vec3.create(0, 1, 0) } }
        : { name: "off", params: {} },
    },
  });
  $("spinBtn").classList.toggle("active", spinning);
}

function takeSnapshot() {
  try {
    const canvas = plugin.canvas3dContext.canvas;
    const a = document.createElement("a");
    a.href = canvas.toDataURL("image/png");
    a.download = `${(S.receptor && S.receptor.id) || "view"}.png`;
    a.click();
  } catch (err) {
    toast(err.message, true);
  }
}

// ---------------------------------------------------------------------------
// wiring
// ---------------------------------------------------------------------------

function wire() {
  $("fetchBtn").onclick = async () => {
    const id = $("pdbId").value.trim();
    if (!id) return;
    $("fetchBtn").disabled = true;
    try {
      await loadReceptor(await api("/api/receptor/fetch", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ pdb_id: id }),
      }));
    } catch (err) { toast(err.message, true); }
    $("fetchBtn").disabled = false;
  };
  $("pdbId").onkeydown = (e) => { if (e.key === "Enter") $("fetchBtn").click(); };

  $("pdbFile").onchange = async (e) => {
    const file = e.target.files[0];
    if (!file) return;
    const fd = new FormData();
    fd.append("file", file);
    try {
      await loadReceptor(await api("/api/receptor/upload", { method: "POST", body: fd }));
    } catch (err) { toast(err.message, true); }
    e.target.value = "";
  };

  $("parseText").onclick = async () => {
    const text = $("smilesBox").value.trim();
    if (!text) return toast("paste SMILES first", true);
    try {
      showLigands(await api("/api/ligands/parse_text", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ text }),
      }));
    } catch (err) { toast(err.message, true); }
  };

  $("ligFile").onchange = async (e) => {
    const file = e.target.files[0];
    if (!file) return;
    const fd = new FormData();
    fd.append("file", file);
    try {
      showLigands(await api("/api/ligands/parse_file", { method: "POST", body: fd }));
    } catch (err) { toast(err.message, true); }
    e.target.value = "";
  };

  document.querySelectorAll(".segbtn").forEach((b) => {
    b.onclick = () => {
      document.querySelectorAll(".segbtn").forEach((x) => x.classList.remove("active"));
      b.classList.add("active");
      S.pickMode = b.dataset.mode;
    };
  });

  $("clearPicks").onclick = clearPicks;
  $("cubeBtn").onclick = () => {
    const s = boxSize().filter(Number.isFinite);
    if (!s.length) return;
    const edge = Math.max(...s);
    ["sx", "sy", "sz"].forEach((id) => ($(id).value = edge.toFixed(1)));
    drawBox();
  };

  ["cx", "cy", "cz", "sx", "sy", "sz"].forEach((id) => ($(id).oninput = drawBox));
  $("styleSel").onchange = styleReceptor;
  $("showSurface").onchange = styleReceptor;
  $("showWaters").onchange = styleReceptor;
  $("resetView").onclick = () => plugin.canvas3d?.requestCameraReset();
  $("zoomInBtn").onclick = () => zoomBy(0.8);
  $("zoomOutBtn").onclick = () => zoomBy(1.25);
  $("zoomSiteBtn").onclick = () => {
    if (S.activeHet) zoomToExpression(hetExpr(S.activeHet));
    else if (S.picks.length) toast("zoom not available for free-hand picks yet - drag/scroll to frame the box", true);
    else toast("set a docking box first", true);
  };
  $("spinBtn").onclick = toggleSpin;
  $("fullscreenBtn").onclick = () => {
    const el = $("viewer");
    if (document.fullscreenElement) document.exitFullscreen();
    else el.requestFullscreen?.();
  };
  $("snapshotBtn").onclick = takeSnapshot;

  $("runBtn").onclick = runDocking;
  $("cancelBtn").onclick = async () => {
    if (!S.job) return;
    try {
      await api(`/api/jobs/${S.job.id}/cancel`, { method: "POST" });
      toast("cancelling - running ligands finish first");
    } catch (err) { toast(err.message, true); }
  };

  ["dlZip", "dlCsv"].forEach((id) => {
    $(id).onclick = (e) => { if ($(id).classList.contains("disabled")) e.preventDefault(); };
  });

  $("engineSel").onchange = () => {
    $("gpuRow").style.display = $("engineSel").value === "gnina" ? "flex" : "none";
  };
}

async function boot() {
  await initViewer();
  wire();
  try {
    const st = await api("/api/status");
    const sel = $("scoring");
    st.scoring_functions.forEach((s) => {
      const o = document.createElement("option");
      o.value = o.textContent = s;
      sel.appendChild(o);
    });
    $("modes").value = st.default_num_modes;

    const pill = $("engine");
    const parts = [];
    parts.push(st.smina ? (st.smina_version || "smina") : "smina not found");
    parts.push(st.gnina ? (st.gnina_version || "gnina") : "gnina not found");
    pill.className = "pill " + (st.smina || st.gnina ? "ok" : "bad");
    pill.textContent = `${parts.join(" - ")} - ${st.workers} workers${st.gpu_available ? " - GPU detected" : ""}`;

    $("gpuInfo").textContent = st.gpu_available
      ? "NVIDIA GPU detected on this machine"
      : "no NVIDIA GPU detected - gnina will run on CPU";
  } catch (err) {
    $("engine").className = "pill bad";
    $("engine").textContent = "backend unreachable";
  }
}

boot();
