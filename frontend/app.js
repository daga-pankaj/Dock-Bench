/* Docking Bench front end.
 *
 * State lives in `S`. The viewer is a single 3Dmol instance holding the
 * receptor model plus, once a result is selected, one model per docked pose.
 */

const $ = (id) => document.getElementById(id);

const S = {
  receptor: null,      // {id, name, het_ligands, n_protein_atoms}
  pdbText: "",
  ligands: [],         // parsed ligand specs from the backend
  pickMode: "ligand",
  picks: [],           // {x,y,z} click points in "pick atoms" mode
  activeHet: null,
  job: null,
  poll: null,
  selectedLigand: null,
  poseModels: [],      // 3Dmol model handles for the currently shown poses
  activePose: 0,
};

let viewer = null;
let boxShape = null;
let pickSpheres = [];

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

// ---------------------------------------------------------------------------
// viewer
// ---------------------------------------------------------------------------

function initViewer() {
  viewer = $3Dmol.createViewer($("viewer"), { backgroundColor: "#0b0e12" });
}

function styleReceptor() {
  if (!viewer || !S.receptor) return;
  const mode = $("styleSel").value;
  viewer.setStyle({}, {});

  const protein = { hetflag: false };
  if (mode === "sticks") viewer.setStyle(protein, { stick: { radius: 0.12, colorscheme: "grayCarbon" } });
  else if (mode === "line") viewer.setStyle(protein, { line: { colorscheme: "grayCarbon" } });
  else viewer.setStyle(protein, { cartoon: { color: "spectrum", opacity: 0.92 } });

  if (mode === "cartoon+sticks") {
    const c = boxCenter();
    if (c) {
      // Residues lining the box get sticks so the pocket is readable.
      viewer.addStyle(
        { hetflag: false, within: { distance: Math.max(...boxSize()) / 2, sel: { hetflag: true } } },
        { stick: { radius: 0.1, colorscheme: "grayCarbon" } }
      );
    }
  }

  // Bound ligands / cofactors always drawn as sticks so they can be clicked.
  viewer.addStyle({ hetflag: true, resn: ["HOH", "WAT", "DOD"], invert: true },
    { stick: { radius: 0.2, colorscheme: "greenCarbon" } });

  if ($("showWaters").checked) {
    viewer.addStyle({ resn: ["HOH", "WAT", "DOD"] }, { sphere: { radius: 0.25, color: "#5aa9ff" } });
  }

  viewer.removeAllSurfaces();
  if ($("showSurface").checked) {
    viewer.addSurface($3Dmol.SurfaceType.VDW,
      { opacity: 0.65, color: "#3a4552" }, { hetflag: false });
  }

  if (S.activeHet) {
    viewer.addStyle(hetSel(S.activeHet), { stick: { radius: 0.28, colorscheme: "orangeCarbon" } });
  }
  viewer.render();
}

function hetSel(het) {
  const sel = { resn: het.resname, resi: parseInt(het.resseq, 10) };
  if (het.chain && het.chain !== "_") sel.chain = het.chain;
  return sel;
}

function boxCenter() {
  const c = [num("cx"), num("cy"), num("cz")];
  return c.every((v) => Number.isFinite(v)) ? c : null;
}

function boxSize() {
  return [num("sx"), num("sy"), num("sz")];
}

function drawBox() {
  if (!viewer) return;
  if (boxShape) { viewer.removeShape(boxShape); boxShape = null; }
  const c = boxCenter();
  const s = boxSize();
  if (!c || !s.every((v) => Number.isFinite(v) && v > 0)) { viewer.render(); return; }
  boxShape = viewer.addBox({
    center: { x: c[0], y: c[1], z: c[2] },
    dimensions: { w: s[0], h: s[1], d: s[2] },
    wireframe: true,
    color: "#4da3ff",
  });
  const vol = (s[0] * s[1] * s[2]).toFixed(0);
  $("siteInfo").textContent = `box ${s.map((v) => v.toFixed(1)).join(" x ")} A  (${vol} A^3)`;
  viewer.render();
}

function setBox(center, size) {
  ["cx", "cy", "cz"].forEach((id, i) => ($(id).value = center[i].toFixed(2)));
  ["sx", "sy", "sz"].forEach((id, i) => ($(id).value = Math.max(8, size[i]).toFixed(1)));
  drawBox();
}

function clearPicks() {
  pickSpheres.forEach((s) => viewer.removeShape(s));
  pickSpheres = [];
  S.picks = [];
  S.activeHet = null;
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
  pickSpheres.forEach((s) => viewer.removeShape(s));
  pickSpheres = [];
  const pad = num("padding") || 0;
  setBox(het.center, het.size.map((v) => v + 2 * pad));
  document.querySelectorAll(".het").forEach((e) =>
    e.classList.toggle("active", e.dataset.key === het.key));
  styleReceptor();
  viewer.zoomTo(hetSel(het));
  viewer.zoom(0.35);
  viewer.render();
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
  pickSpheres.push(viewer.addSphere({
    center: { x: atom.x, y: atom.y, z: atom.z },
    radius: 0.55, color: "#f0b429", opacity: 0.9,
  }));
  S.activeHet = null;
  const [center, size] = boxFromPoints(S.picks);
  const minEdge = S.picks.length === 1 ? 20 : 14;
  setBox(center, size.map((v) => Math.max(v, minEdge)));
  toast(`${S.picks.length} pick(s) - ${atom.resn}${atom.resi} ${atom.atom}`);
  viewer.render();
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

  viewer.clear();
  S.poseModels = [];
  boxShape = null;
  pickSpheres = [];
  S.picks = [];
  S.activeHet = null;
  viewer.addModel(S.pdbText, "pdb");
  styleReceptor();
  viewer.setClickable({}, true, (atom) => onAtomClick(atom));
  viewer.zoomTo();
  viewer.render();
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
    keep_waters: $("keepWaters").checked,
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

function clearPoses() {
  S.poseModels.forEach((m) => viewer.removeModel(m));
  S.poseModels = [];
  $("poseBar").classList.add("hidden");
}

async function showPoses(lig) {
  if (!S.job) return;
  try {
    const sdf = await api(`/api/jobs/${S.job.id}/ligands/${encodeURIComponent(lig.name)}/poses.sdf`);
    clearPoses();
    S.selectedLigand = lig.name;

    const blocks = sdf.split("$$$$").map((b) => b.trim()).filter((b) => b.length > 10);
    blocks.forEach((block) => {
      const model = viewer.addModel(block + "\n$$$$\n", "sdf");
      S.poseModels.push(model);
    });

    const bar = $("poseBar");
    bar.innerHTML = `<span class="title">${lig.name} - poses:</span>`;
    S.poseModels.forEach((_, i) => {
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
  S.poseModels.forEach((m, i) => {
    const show = idx === -1 || i === idx;
    if (!show) {
      m.setStyle({}, {});   // empty style hides the model
      return;
    }
    m.setStyle({}, {
      stick: { radius: idx === -1 ? 0.12 : 0.16, colorscheme: "cyanCarbon" },
    });
  });
  document.querySelectorAll("#poseBar .pose").forEach((b, i) => {
    b.classList.toggle("active", i === (idx === -1 ? S.poseModels.length : idx));
  });
  viewer.render();
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
  $("resetView").onclick = () => { viewer.zoomTo(); viewer.render(); };

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
}

async function boot() {
  initViewer();
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
    if (st.smina) {
      pill.className = "pill ok";
      pill.textContent = `${st.smina_version || "smina"} - ${st.workers} workers`;
    } else {
      pill.className = "pill bad";
      pill.textContent = "smina not found - run ./setup.sh";
    }
  } catch (err) {
    $("engine").className = "pill bad";
    $("engine").textContent = "backend unreachable";
  }
}

boot();
