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
import { Structure, StructureElement, StructureProperties } from "molstar/lib/mol-model/structure";
import { StructureSelectionQueries as Queries } from "molstar/lib/mol-plugin-state/helpers/structure-selection-query";
import { PluginStateObject as PSO, PluginStateTransform } from "molstar/lib/mol-plugin-state/objects";
import { InteractionsRepresentationProvider } from "molstar/lib/mol-model-props/computed/representations/interactions";
import { InteractionTypeColorThemeProvider } from "molstar/lib/mol-model-props/computed/themes/interaction-type";
import { MolScriptBuilder as MS } from "molstar/lib/mol-script/language/builder";
import { Color } from "molstar/lib/mol-util/color";
import { Vec3 } from "molstar/lib/mol-math/linear-algebra";
import { StateObject } from "molstar/lib/mol-state";
import { Task } from "molstar/lib/mol-task";
import { ParamDefinition as PD } from "molstar/lib/mol-util/param-definition";
import "molstar/lib/mol-util/polyfill";

// Combines two already-loaded structures (receptor + one docked pose) into a
// single Structure so Mol*'s non-covalent-interactions computation - which
// only looks within one Structure - can see contacts between them. Adapted
// verbatim from Mol*'s own src/apps/docking-viewer/index.ts.
const MergeStructures = PluginStateTransform.BuiltIn({
  name: "docking-bench-merge-structures",
  display: { name: "Merge Structures", description: "Merge Structure" },
  from: PSO.Root,
  to: PSO.Molecule.Structure,
  params: {
    structures: PD.ObjectList({ ref: PD.Text("") }, ({ ref }) => ref, { isHidden: true }),
  },
})({
  apply({ params, dependencies }) {
    return Task.create("Merge Structures", async () => {
      if (params.structures.length === 0) return StateObject.Null;
      const first = dependencies[params.structures[0].ref].data;
      const builder = Structure.Builder({ masterModel: first.models[0] });
      for (const { ref } of params.structures) {
        const s = dependencies[ref].data;
        for (const unit of s.units) {
          builder.addUnit(unit.kind, unit.model, unit.conformation.operator, unit.elements, unit.traits);
        }
      }
      return new PSO.Molecule.Structure(builder.getStructure(), { label: "Merged Structure" });
    });
  },
});

const $ = (id) => document.getElementById(id);

// Docking-parameter defaults. Peptides/peptidomimetics are floppier and need
// a much wider conformer search than a typical small-molecule ligand.
const STANDARD_DOCKING_PARAMS = { exhaustiveness: 8, energy_range: 4 };
const PEPTIDE_DOCKING_PARAMS = { exhaustiveness: 64, energy_range: 10 };

const S = {
  receptor: null,      // {id, name, het_ligands, n_protein_atoms}
  pdbText: "",
  ligands: [],         // parsed ligand specs from the backend
  pickMode: "ligand",
  picks: [],           // [x,y,z] click points in "pick atoms" mode
  activeHet: null,
  excludeChains: new Set(), // chain IDs to delete from the receptor before docking (e.g. a peptide ligand)
  constraints: [],     // pharmacophore constraints: {id, kind, x, y, z, radius, required, label}
  pendingPickId: null, // constraint row currently waiting for a viewer click to set its point
  covalentAtomPickActive: false, // waiting for user to click nucleophilic atom in viewer
  job: null,
  poll: null,
  selectedLigand: null,
  selectedLigandObj: null, // the ligand row object last passed to showPoses()
  poseReps: [],        // Mol* representation refs, one per pose currently shown
  activePose: 0,
  posePlayTimer: null,   // setInterval handle while the pose bar's play button is cycling poses
  hideOrigLigand: false, // hide co-crystallized ligand from the complex
  hideProtein: false,    // hide the receptor's polymer cartoon from the complex
  diagramHidden: false,   // 2D interaction diagram panel toggled off
  diagramUrl: null,        // current object URL behind #diagramImg, revoked on change
  depictionUrl: null,      // object URL behind the viewer's 2D ligand overlay
  receptorPanelHidden: false,     // "1 · Receptor" card minimized
  activeSitePanelHidden: false,   // "2 · Active site" card minimized
  pharmacophorePanelHidden: false, // "3 · Pharmacophore constraints" card minimized
  covalentPanelHidden: false,      // "3b · Covalent docking" card minimized
  ligandsPanelHidden: false,       // "4 · Ligands" card minimized
  expandedResultRows: new Set(),   // ligand names whose per-pose score/strain breakdown is open
};

let plugin = null;
let receptorStructure = null;   // StateObjectSelector for the loaded receptor Structure
let receptorComponentRefs = []; // components (+ their reps) built by styleReceptor, for teardown
let boxRefs = [];                // state refs for the docking-box wireframe structure
let picksRefs = [];              // state refs for the pick-point markers structure
let constraintMarkerRefs = [];   // state refs for the pharmacophore-constraint markers structure
let poseStructureRefs = [];      // top-level structure refs, one per pose (for teardown)
let poseStructures = [];         // Structure selectors, one per pose (for the interactions merge below)
let interactionRefs = [];        // merged receptor+pose structure showing the active pose's contacts
let receptorLigandRep = null;    // the receptor's own bound-ligand(s) rep, for the "hide" checkbox

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
  // The real, full Mol* viewer app (sequence panel, its own controls/log,
  // native selection + representation tools) - not a stripped-down custom
  // canvas. We only layer the active-site click hook on top of it, plus one
  // override: Mol*'s default spec wires a "remote state" feature that saves/
  // loads session state to a public server (webchem.ncbr.muni.cz) run by the
  // Mol* project - out of place in a tool that's otherwise entirely local
  // (see README: "no external services beyond fetching PDB files from RCSB").
  // Drop the panel AND blank the server URLs, so there's no leftover path to
  // it (a config value, not just a hidden button).
  const defaultSpec = DefaultPluginUISpec();
  const spec = {
    ...defaultSpec,
    components: { ...defaultSpec.components, remoteState: "none" },
    config: [
      ...(defaultSpec.config || []),
      [PluginConfig.State.DefaultServer, ""],
      [PluginConfig.State.CurrentServer, ""],
    ],
  };
  plugin = await createPluginUI({ target: $("viewer"), spec, render: renderReact18 });

  await PluginCommands.Canvas3D.SetSettings(plugin, {
    settings: {
      renderer: { ...plugin.canvas3d.props.renderer, backgroundColor: Color(0x000000) },
    },
  });

  await plugin.managers.structure.component.setOptions({
    ...plugin.managers.structure.component.state.options,
    hydrogens: "only-polar",
  });

  plugin.behaviors.interaction.click.subscribe(onPluginClick);
}

// Residue name used by the synthetic box-wireframe / pick-marker atoms
// (see drawBox/rebuildPicks below) — clicks landing on our own overlay
// geometry must be ignored, the same way NGL's shape/buffer primitives were
// never atom-pickable in the first place.
const SYNTHETIC_RESNAMES = new Set(["BOX", "PIK", "PHC"]);

// Side-chain aromatic ring atom names per residue, keyed by residue name -
// used to snap an "aromatic" constraint pick to the ring's centroid (the
// real pi-stacking anchor) instead of whichever single ring atom got
// clicked. TRP has two rings; whichever one contains the clicked atom wins.
const AROMATIC_RINGS = {
  PHE: [["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]],
  TYR: [["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]],
  HIS: [["CG", "ND1", "CD2", "CE1", "NE2"]],
  TRP: [
    ["CG", "CD1", "NE1", "CE2", "CD2"],
    ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
  ],
};

// Given the clicked location, resolve the aromatic ring (if any) the clicked
// atom belongs to and return {centroid, normal} - null if the residue/atom isn't
// part of a recognized aromatic ring, or the ring came back incomplete.
// Normal vector computed via cross product of first two edges, normalized.
function ringCentroidAt(loc, resn, atomName) {
  const rings = AROMATIC_RINGS[resn];
  if (!rings) return null;
  const ring = rings.find((names) => names.includes(atomName));
  if (!ring) return null;
  const chain = StructureProperties.chain.auth_asym_id(loc);
  const resi = StructureProperties.residue.auth_seq_id(loc);
  const probe = StructureElement.Location.create(loc.structure);
  probe.unit = loc.unit;
  const p = Vec3.zero();
  const pts = [];
  for (const el of loc.unit.elements) {
    probe.element = el;
    if (StructureProperties.chain.auth_asym_id(probe) !== chain) continue;
    if (StructureProperties.residue.auth_seq_id(probe) !== resi) continue;
    if (!ring.includes(StructureProperties.atom.label_atom_id(probe))) continue;
    loc.unit.conformation.position(el, p);
    pts.push([p[0], p[1], p[2]]);
  }
  if (pts.length < ring.length) return null; // incomplete ring - don't trust a partial centroid
  const n = pts.length;
  const sum = pts.reduce((a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]], [0, 0, 0]);
  const centroid = [sum[0] / n, sum[1] / n, sum[2] / n];

  // Compute plane normal: cross product of two edge vectors
  if (pts.length < 3) return { centroid, normal: null };
  const v1 = [pts[1][0] - pts[0][0], pts[1][1] - pts[0][1], pts[1][2] - pts[0][2]];
  const v2 = [pts[2][0] - pts[0][0], pts[2][1] - pts[0][1], pts[2][2] - pts[0][2]];
  const cross = [
    v1[1] * v2[2] - v1[2] * v2[1],
    v1[2] * v2[0] - v1[0] * v2[2],
    v1[0] * v2[1] - v1[1] * v2[0]
  ];
  const mag = Math.sqrt(cross[0] ** 2 + cross[1] ** 2 + cross[2] ** 2);
  const normal = mag > 0 ? [cross[0] / mag, cross[1] / mag, cross[2] / mag] : null;

  return { centroid, normal };
}

function onPluginClick(event) {
  try {
    const loci = event.current.loci;
    if (!loci || !StructureElement.Loci.is(loci) || StructureElement.Loci.size(loci) === 0) return;
    const loc = StructureElement.Location.create(loci.structure);
    StructureElement.Loci.getFirstLocation(loci, loc);
    if (!loc.unit) return;

    const resn = StructureProperties.atom.label_comp_id(loc);
    if (SYNTHETIC_RESNAMES.has(resn)) return;

    const p = Vec3.zero();
    loc.unit.conformation.position(loc.element, p);

    const atomName = StructureProperties.atom.label_atom_id(loc);
    onAtomClick({
      resn,
      hetflag: StructureProperties.residue.group_PDB(loc) === "HETATM",
      chain: StructureProperties.chain.auth_asym_id(loc),
      resi: StructureProperties.residue.auth_seq_id(loc),
      x: p[0], y: p[1], z: p[2],
      atom: atomName,
      ring: ringCentroidAt(loc, resn, atomName), // [x,y,z] centroid, or null - see AROMATIC_RINGS
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

// Raw addRepresentation ignores the viewer's "Hydrogens" option (only presets
// apply it), so mirror it here for the built-in representation types.
function hydrogenTypeParams() {
  const h = plugin.managers.structure.component.state.options.hydrogens;
  return { ignoreHydrogens: h !== "all", ignoreHydrogensVariant: h === "only-polar" ? "non-polar" : "all" };
}

async function addRep(componentRef, type, color, colorParams, typeParams) {
  if (!componentRef) return null;
  if (typeof type === "string") typeParams = { ...hydrogenTypeParams(), ...typeParams };
  return plugin.builders.structure.representation.addRepresentation(componentRef, {
    type, color, colorParams, typeParams,
  });
}

// Applies a sane default look once on load (and whenever the active site
// changes, to (re)draw the orange highlight). Past that, representations are
// the user's to manage from Mol*'s own right-hand panel - we don't fight it
// with a competing style dropdown.
async function styleReceptor() {
  if (!plugin || !receptorStructure) return;
  await deleteRefs(receptorComponentRefs);
  receptorComponentRefs = [];
  const track = (ref) => { if (ref) receptorComponentRefs.push(ref); };
  const hideOrig = S.hideOrigLigand;
  const hideIfOrig = (rep) => {
    if (rep && hideOrig) plugin.state.data.updateCellState(rep.ref, { isHidden: true });
  };
  const hideIfProtein = (rep) => {
    if (rep && S.hideProtein) plugin.state.data.updateCellState(rep.ref, { isHidden: true });
  };

  const polymerComp = await plugin.builders.structure.tryCreateComponentFromExpression(
    receptorStructure, Queries.polymer.expression, "polymer");
  track(polymerComp);
  const polymerRep = await addRep(polymerComp, "cartoon", "chain-id", undefined, { opacity: 0.92 });
  track(polymerRep);
  hideIfProtein(polymerRep);

  // bound ligands / cofactors — always sticks so they stay clickable (unless
  // the "hide existing ligand" checkbox is on, e.g. to declutter a pocket
  // that already has a co-crystallised ligand while looking at a docked one)
  const ligandExpr = MS.struct.combinator.merge([Queries.ligand.expression, Queries.ion.expression]);
  const ligandComp = await plugin.builders.structure.tryCreateComponentFromExpression(
    receptorStructure, ligandExpr, "ligands");
  track(ligandComp);
  receptorLigandRep = await addRep(ligandComp, "ball-and-stick", "uniform", { value: Color(0x4caf50) }, { sizeFactor: 0.25 });
  track(receptorLigandRep);
  hideIfOrig(receptorLigandRep);

  if (S.activeHet) {
    const activeComp = await plugin.builders.structure.tryCreateComponentFromExpression(
      receptorStructure, hetExpr(S.activeHet), "active-site");
    track(activeComp);
    // the orange active-site highlight sits on the same atoms as the ligand
    // rep above (when the active site *is* a bound ligand) - "hide existing
    // ligand" should take that with it too.
    const activeRep = await addRep(activeComp, "ball-and-stick", "uniform", { value: Color(0xf5a623) }, { sizeFactor: 0.28 });
    track(activeRep);
    hideIfOrig(activeRep);
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
        lines.push(pdbAtomLine(serial, "C", "BOX", "Z", serial, p[0], p[1], p[2], "C"));
        serial++;
      }
    }
  }
  const { structure, refs } = await loadSyntheticPdb(lines.join("\n") + "\nEND\n", "dockbox");
  const rep = await addRep(structure, "line", "uniform", { value: Color(0x4da3ff) }, { lineSizeAttenuation: false });
  boxRefs = [...refs, rep].filter(Boolean);

  const vol = (s[0] * s[1] * s[2]).toFixed(0);
  const siteInfo = $("siteInfo");
  if (siteInfo) {
    siteInfo.textContent = `box ${s.map((v) => v.toFixed(1)).join(" x ")} A  (${vol} A^3)`;
  }
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

async function rebuildConstraintMarkers() {
  if (!plugin) return;
  await deleteRefs(constraintMarkerRefs);
  constraintMarkerRefs = [];
  const placed = S.constraints.filter((c) => [c.x, c.y, c.z].every(Number.isFinite));
  if (!placed.length) return;

  // Build one structure per constraint type for separate coloring
  for (const c of placed) {
    const color = CONSTRAINT_COLORS[c.kind] || 0xff6bd6;
    const lines = [];
    let serial = 1;

    // Center point
    lines.push(pdbAtomLine(serial++, "C", "PHC", "Z", 1, c.x, c.y, c.z, "C"));

    // Aromatic: disk in plane with normal
    if (c.kind === "aromatic" && [c.nx, c.ny, c.nz].every(Number.isFinite)) {
      const diskRadius = 1.5;
      const nSteps = 16;
      // Find two orthogonal vectors in the plane
      const n = [c.nx, c.ny, c.nz];
      const up = Math.abs(n[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0];
      const u = [
        up[1] * n[2] - up[2] * n[1],
        up[2] * n[0] - up[0] * n[2],
        up[0] * n[1] - up[1] * n[0]
      ];
      const uMag = Math.sqrt(u[0] ** 2 + u[1] ** 2 + u[2] ** 2);
      const u1 = [u[0] / uMag, u[1] / uMag, u[2] / uMag];
      const v = [
        n[1] * u1[2] - n[2] * u1[1],
        n[2] * u1[0] - n[0] * u1[2],
        n[0] * u1[1] - n[1] * u1[0]
      ];
      for (let i = 0; i <= nSteps; i++) {
        const theta = (2 * Math.PI * i) / nSteps;
        const px = c.x + diskRadius * (Math.cos(theta) * u1[0] + Math.sin(theta) * v[0]);
        const py = c.y + diskRadius * (Math.cos(theta) * u1[1] + Math.sin(theta) * v[1]);
        const pz = c.z + diskRadius * (Math.cos(theta) * u1[2] + Math.sin(theta) * v[2]);
        lines.push(pdbAtomLine(serial++, "C", "PHC", "Z", serial, px, py, pz, "C"));
      }
    }

    // H-bond: arrow showing direction
    else if (["hbond_donor", "hbond_acceptor"].includes(c.kind) && [c.dx, c.dy, c.dz].every(Number.isFinite)) {
      const arrowLen = 2.0;
      const tipLen = 0.5;
      // Shaft
      for (let t = 0; t <= 1; t += 0.2) {
        const px = c.x + t * arrowLen * c.dx;
        const py = c.y + t * arrowLen * c.dy;
        const pz = c.z + t * arrowLen * c.dz;
        lines.push(pdbAtomLine(serial++, "C", "PHC", "Z", serial, px, py, pz, "C"));
      }
      // Arrowhead (simple cone approximation)
      const tip = [c.x + arrowLen * c.dx, c.y + arrowLen * c.dy, c.z + arrowLen * c.dz];
      const base = [c.x + (arrowLen - tipLen) * c.dx, c.y + (arrowLen - tipLen) * c.dy, c.z + (arrowLen - tipLen) * c.dz];
      const perp = Math.abs(c.dz) < 0.9 ? [0, 0, 1] : [1, 0, 0];
      const p1 = [
        perp[1] * c.dz - perp[2] * c.dy,
        perp[2] * c.dx - perp[0] * c.dz,
        perp[0] * c.dy - perp[1] * c.dx
      ];
      const p1Mag = Math.sqrt(p1[0] ** 2 + p1[1] ** 2 + p1[2] ** 2);
      const p1n = [p1[0] / p1Mag, p1[1] / p1Mag, p1[2] / p1Mag];
      for (let i = 0; i < 4; i++) {
        const theta = (Math.PI * i) / 2;
        const px = base[0] + 0.3 * Math.cos(theta) * p1n[0];
        const py = base[1] + 0.3 * Math.cos(theta) * p1n[1];
        const pz = base[2] + 0.3 * Math.cos(theta) * p1n[2];
        for (const pt of edgePoints([px, py, pz], tip, 0.5)) {
          lines.push(pdbAtomLine(serial++, "C", "PHC", "Z", serial, pt[0], pt[1], pt[2], "C"));
        }
      }
    }

    const { structure, refs } = await loadSyntheticPdb(lines.join("\n") + "\nEND\n", `constraint-${c.id}`);
    const rep = await addRep(structure, "ball-and-stick", "uniform", { value: Color(color) }, { sizeFactor: 0.25 });
    constraintMarkerRefs.push(...refs, rep);
  }
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

  // Covalent atom picking mode - user selecting nucleophilic residue atom
  if (S.covalentAtomPickActive) {
    const nucleophileMap = {
      'CYS': 'SG',
      'SER': 'OG',
      'THR': 'OG1',
      'TYR': 'OH',
      'LYS': 'NZ'
    };

    const resname = (atom.resn || "").toUpperCase();
    const expectedAtom = nucleophileMap[resname];

    // Validate nucleophilic residue type
    if (!expectedAtom) {
      toast(`⚠ ${resname} not a typical nucleophile. Supported: CYS, SER, THR, TYR, LYS`, true);
      S.covalentAtomPickActive = false;
      $("pickCovalentAtomBtn").classList.remove("armed");
      return;
    }

    // Warn if wrong atom selected (but allow it)
    if (atom.atom !== expectedAtom) {
      toast(`⚠ Picked ${atom.atom}, but ${resname} nucleophile is typically ${expectedAtom}. Proceeding anyway.`, true);
    }

    // Fill residue field: chain:resname:resnum
    const chain = atom.chain || "";
    const spec = `${chain}:${resname}:${atom.resi}`;
    $("covalentResidue").value = spec;

    S.covalentAtomPickActive = false;
    $("pickCovalentAtomBtn").classList.remove("armed");
    toast(`✓ Covalent residue set: ${spec} (${atom.atom})`);
    return;
  }

  // A constraint row is waiting for its point - this click sets it, no
  // matter what the box-picking mode ("click ligand" / "pick atoms") is.
  if (S.pendingPickId != null) {
    const c = S.constraints.find((c) => c.id === S.pendingPickId);
    if (!c) { S.pendingPickId = null; return; }

    const needsDirection = ["hbond_donor", "hbond_acceptor"].includes(c.kind);
    const hasPoint = [c.x, c.y, c.z].every(Number.isFinite);

    // H-bond: two-click mode - first click = point, second click = direction
    if (needsDirection && hasPoint) {
      // Second click: compute direction vector from point to this atom
      const dx = atom.x - c.x;
      const dy = atom.y - c.y;
      const dz = atom.z - c.z;
      const mag = Math.sqrt(dx ** 2 + dy ** 2 + dz ** 2);
      if (mag < 0.1) {
        toast("direction atom too close to constraint point - pick a different atom", true);
        return;
      }
      c.dx = dx / mag; c.dy = dy / mag; c.dz = dz / mag;
      S.pendingPickId = null;
      renderConstraints();
      rebuildConstraintMarkers();
      toast(`direction set from ${atom.resn}${atom.resi} ${atom.atom}`);
      return;
    }

    // First click (or single-click for non-H-bond types)
    S.pendingPickId = needsDirection && hasPoint ? null : S.pendingPickId;

    // Aromatic: use ring centroid + plane normal
    if (c.kind === "aromatic" && atom.ring) {
      c.x = atom.ring.centroid[0];
      c.y = atom.ring.centroid[1];
      c.z = atom.ring.centroid[2];
      if (atom.ring.normal) {
        c.nx = atom.ring.normal[0];
        c.ny = atom.ring.normal[1];
        c.nz = atom.ring.normal[2];
      }
      S.pendingPickId = null;
      renderConstraints();
      rebuildConstraintMarkers();
      toast(`aromatic plane set from ${atom.resn}${atom.resi} ring`);
      return;
    }

    // Everything else: single click sets point
    c.x = atom.x; c.y = atom.y; c.z = atom.z;
    if (needsDirection) {
      toast(`point set - now click a second atom to define ${c.kind === "hbond_donor" ? "donor" : "acceptor"} direction`);
      // Keep S.pendingPickId active for second click
    } else {
      S.pendingPickId = null;
      renderConstraints();
      rebuildConstraintMarkers();
      toast(`constraint point set from ${atom.resn}${atom.resi} ${atom.atom}`);
    }
    renderConstraints();
    rebuildConstraintMarkers();
    return;
  }

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
  if (!box) return;
  box.innerHTML = "";
  const hets = S.receptor.het_ligands || [];
  if (!hets.length) {
    box.innerHTML = '<span class="muted small">no bound ligands found - use "pick atoms" or "find sites"</span>';
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

// Every chain, so a co-crystallized peptide/peptidomimetic ligand (its own
// ATOM chain, invisible to list_het_ligands which only sees HETATM residues)
// can be ticked for deletion from the receptor before docking.
function renderChainList() {
  const box = $("chainList");
  if (!box) return;
  box.innerHTML = "";
  const chains = S.receptor.chains || [];
  if (chains.length < 2) {
    box.innerHTML = '<span class="muted small">single chain - nothing extra to exclude</span>';
    return;
  }
  chains.forEach((c) => {
    const row = document.createElement("label");
    row.className = "chain-row";
    const checked = S.excludeChains.has(c.chain) ? "checked" : "";
    row.innerHTML =
      `<input type="checkbox" ${checked}> delete chain ${c.chain}` +
      `<span class="n">${c.n_atoms} atom(s), ${c.n_hetatm} het</span>`;
    row.querySelector("input").onchange = (e) => {
      if (e.target.checked) S.excludeChains.add(c.chain);
      else S.excludeChains.delete(c.chain);
    };
    box.appendChild(row);
  });
}

async function loadReceptor(meta) {
  S.receptor = meta;
  S.pdbText = await api(`/api/receptor/${meta.id}/pdb`);
  const receptorInfo = $("receptorInfo");
  if (receptorInfo) {
    receptorInfo.textContent =
      `${meta.name} - ${meta.n_protein_atoms} protein atoms - ${(meta.het_ligands || []).length} het groups`;
  }

  clearPoses();
  await deleteRefs(boxRefs); boxRefs = [];
  await deleteRefs(picksRefs); picksRefs = [];
  await deleteRefs(constraintMarkerRefs); constraintMarkerRefs = [];
  if (receptorStructure) { await deleteRefs([receptorStructure]); receptorStructure = null; receptorComponentRefs = []; }
  S.picks = [];
  S.activeHet = null;
  S.excludeChains = new Set();
  S.constraints = [];
  S.pendingPickId = null;
  renderConstraints();

  const data = await plugin.builders.data.rawData({ data: S.pdbText, label: meta.name });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, "pdb");
  const model = await plugin.builders.structure.createModel(trajectory);
  receptorStructure = await plugin.builders.structure.createStructure(model);

  await styleReceptor();
  plugin.canvas3d?.requestCameraReset();
  renderHetList();
  renderChainList();

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
// pharmacophore constraints
// ---------------------------------------------------------------------------

let constraintSeq = 0;

const CONSTRAINT_KIND_LABELS = {
  hbond_donor: "H-bond donor",
  hbond_acceptor: "H-bond acceptor",
  positional: "positional (any atom)",
  aromatic: "aromatic ring (pi-stacking)",
  metal_coordination: "metal coordination",
  smarts: "SMARTS match",
};

// Color per constraint type for visualization
const CONSTRAINT_COLORS = {
  hbond_donor: 0x4a90e2,      // blue
  hbond_acceptor: 0xe24a4a,   // red
  positional: 0xf5a623,       // yellow
  aromatic: 0x9b59d6,         // purple
  metal_coordination: 0x27ae60, // green
  smarts: 0x1abc9c,           // cyan
};

// Typical contact distance per kind, shown as a hint next to the radius
// field - not enforced, just a starting point for the ones that aren't 3.5A.
const CONSTRAINT_RADIUS_HINT = {
  hbond_donor: "typical 2.5-3.5 Å",
  hbond_acceptor: "typical 2.5-3.5 Å",
  positional: "",
  aromatic: "ring-centroid distance; typical 4-5.5 Å",
  metal_coordination: "metal-ligand bond; typical 2.0-2.6 Å",
  smarts: "",
};

function addConstraintRow() {
  S.constraints.push({
    id: ++constraintSeq,
    kind: "hbond_donor",
    x: undefined, y: undefined, z: undefined,
    // Direction vector for H-bond donor/acceptor (dx, dy, dz)
    dx: undefined, dy: undefined, dz: undefined,
    // Normal vector for aromatic plane (nx, ny, nz)
    nx: undefined, ny: undefined, nz: undefined,
    radius: 3.5,
    required: true,
    label: "",
    smarts: "",
  });
  renderConstraints();
}

function removeConstraintRow(id) {
  S.constraints = S.constraints.filter((c) => c.id !== id);
  if (S.pendingPickId === id) S.pendingPickId = null;
  renderConstraints();
  rebuildConstraintMarkers();
}

function renderConstraints() {
  const box = $("constraintList");
  if (!box) return;
  box.innerHTML = "";
  S.constraints.forEach((c) => {
    const row = document.createElement("div");
    row.className = "constraint-row";
    const kindOpts = Object.entries(CONSTRAINT_KIND_LABELS)
      .map(([v, label]) => `<option value="${v}" ${v === c.kind ? "selected" : ""}>${label}</option>`)
      .join("");
    const pt = (v) => (Number.isFinite(v) ? v.toFixed(2) : "");
    const smartsRow = c.kind === "smarts"
      ? `<input class="c-smarts" type="text" placeholder="SMARTS pattern, e.g. [OX1]=[C,S]" value="${c.smarts || ""}">`
      : "";

    // Status badge for direction/normal
    const hasPoint = [c.x, c.y, c.z].every(Number.isFinite);
    const hasDir = ["hbond_donor", "hbond_acceptor"].includes(c.kind) && [c.dx, c.dy, c.dz].every(Number.isFinite);
    const hasNormal = c.kind === "aromatic" && [c.nx, c.ny, c.nz].every(Number.isFinite);
    let statusBadge = "";
    if (["hbond_donor", "hbond_acceptor"].includes(c.kind)) {
      statusBadge = hasPoint && hasDir ? '<span class="tag done">✓ directional</span>' :
                    hasPoint ? '<span class="tag warn">needs direction</span>' : '';
    } else if (c.kind === "aromatic") {
      statusBadge = hasNormal ? '<span class="tag done">✓ planar</span>' : '';
    }

    row.innerHTML = `
      <select class="c-kind" title="${CONSTRAINT_RADIUS_HINT[c.kind] || ""}">${kindOpts}</select>
      <button class="ghost c-pick" type="button">pick point</button>
      ${statusBadge}
      <input class="c-x" type="number" step="0.1" placeholder="x" value="${pt(c.x)}">
      <input class="c-y" type="number" step="0.1" placeholder="y" value="${pt(c.y)}">
      <input class="c-z" type="number" step="0.1" placeholder="z" value="${pt(c.z)}">
      <label class="inline" title="${CONSTRAINT_RADIUS_HINT[c.kind] || ""}">radius <input class="c-radius" type="number" step="0.1" min="0.5" max="15" value="${c.radius}"> Å</label>
      <label class="inline"><input class="c-required" type="checkbox" ${c.required ? "checked" : ""}> required</label>
      <input class="c-label" type="text" placeholder="label (optional)" value="${c.label}">
      ${smartsRow}
      <button class="ghost c-remove" type="button" title="remove">✕</button>`;

    row.querySelector(".c-kind").onchange = (e) => { c.kind = e.target.value; renderConstraints(); };
    row.querySelector(".c-x").oninput = (e) => { c.x = parseFloat(e.target.value); rebuildConstraintMarkers(); };
    row.querySelector(".c-y").oninput = (e) => { c.y = parseFloat(e.target.value); rebuildConstraintMarkers(); };
    row.querySelector(".c-z").oninput = (e) => { c.z = parseFloat(e.target.value); rebuildConstraintMarkers(); };
    row.querySelector(".c-radius").oninput = (e) => { c.radius = parseFloat(e.target.value); };
    row.querySelector(".c-required").onchange = (e) => { c.required = e.target.checked; };
    row.querySelector(".c-label").oninput = (e) => { c.label = e.target.value; };
    row.querySelector(".c-smarts")?.addEventListener("input", (e) => { c.smarts = e.target.value; });
    row.querySelector(".c-pick").classList.toggle("armed", S.pendingPickId === c.id);
    row.querySelector(".c-pick").onclick = () => {
      S.pendingPickId = S.pendingPickId === c.id ? null : c.id;
      renderConstraints();
      if (S.pendingPickId != null) toast("click an atom in the viewer to set this constraint's point");
    };
    row.querySelector(".c-remove").onclick = () => removeConstraintRow(c.id);

    box.appendChild(row);
  });
}

// Built job-submission payload for the constraints currently defined, or
// null (with a toast) if any row is missing its point - such a constraint
// can't be checked against anything.
function buildConstraintsPayload() {
  const incomplete = S.constraints.some((c) => ![c.x, c.y, c.z].every(Number.isFinite));
  if (incomplete) {
    toast("a pharmacophore constraint is missing its point - pick one or remove the row", true);
    return null;
  }
  const badSmarts = S.constraints.some((c) => c.kind === "smarts" && !c.smarts?.trim());
  if (badSmarts) {
    toast("a SMARTS constraint is missing its pattern", true);
    return null;
  }
  return S.constraints.map((c) => ({
    kind: c.kind,
    point: [c.x, c.y, c.z],
    radius: c.radius,
    required: c.required,
    label: c.label,
    smarts: c.kind === "smarts" ? c.smarts.trim() : undefined,
  }));
}

// Loads `smiles` as the ligand to dock (named after the selected bound
// ligand) and turns on RMSD validation against it - shared by both
// self-dock buttons below (RCSB-looked-up SMILES, or one pasted by hand).
// The RMSD itself isn't computed here - ticking validateRmsd just makes the
// next runDocking() send reference_ligand_key, and backend/jobs.py computes
// it per pose once that job actually runs (see backend/rmsd.py).
async function loadSelfDockLigand(smiles, label) {
  const data = await api("/api/ligands/parse_text", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ text: `${smiles} ${S.activeHet.resname}_bound` }),
  });
  if (!data.n_valid) throw new Error("that SMILES didn't parse");
  showLigands(data);
  $("validateRmsd").checked = true;
  toast(`loaded ${S.activeHet.resname}${label ? ` (${label})` : ""} for self-docking - run docking to get RMSD`);
}

// Convenience for repeated redocking tests: fetch the selected bound
// ligand's own SMILES from RCSB's Chemical Component Dictionary (keyed by
// its 3-letter residue code) instead of pasting it in by hand every time.
async function selfDockFromBoundLigand() {
  if (!S.activeHet) return toast("pick a bound ligand above first", true);
  const btn = $("selfDockBtn");
  btn.disabled = true;
  try {
    const q = new URLSearchParams({ receptor_id: S.receptor.id, het_key: S.activeHet.key });
    const { smiles, name, covalent_note } = await api(
      `/api/pdb_components/${encodeURIComponent(S.activeHet.resname)}/smiles?${q}`);
    await loadSelfDockLigand(smiles, name);
    if (covalent_note) toast(covalent_note);
  } catch (err) {
    toast(err.message, true);
  }
  btn.disabled = false;
}

// Same, but for a ligand that isn't on RCSB at all (a novel/in-house
// compound crystallized in your own receptor) - you supply its SMILES
// yourself instead of it being looked up.
async function selfDockWithManualSmiles() {
  if (!S.activeHet) return toast("pick a bound ligand above first", true);
  const smiles = $("selfDockSmiles").value.trim();
  if (!smiles) return toast("paste the bound ligand's SMILES first", true);
  try {
    await loadSelfDockLigand(smiles);
  } catch (err) {
    toast(err.message, true);
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
  const constraints = buildConstraintsPayload();
  if (constraints === null) return;
  if ($("validateRmsd").checked && !S.activeHet) {
    return toast("RMSD validation needs a bound ligand selected above as the site", true);
  }

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
    minimize_after_docking: $("minimizeAfterDocking").checked,
    keep_waters: $("keepWaters").checked,
    remove_ligand: $("removeLigand").checked,
    exclude_chains: [...S.excludeChains],
    trim_distant_chains: $("trimChains").checked,
    constraints,
    reference_ligand_key: $("validateRmsd").checked ? S.activeHet.key : null,
    moe_placement: $("moePlacement").value,
    moe_refinement: $("moeRefinement").value,
    ph4_id: S.ph4_id || null,
    covalent_enabled: $("covalentEnabled").checked,
    covalent_reaction: $("covalentReaction").value || null,
    covalent_residue: $("covalentResidue").value.trim() || null,
    custom_rxn_id: S.custom_rxn_id || null,
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
        // A structure mismatch (docked ligand isn't really the same molecule
        // as the bound reference - wrong SMILES, wrong bound ligand picked,
        // etc) makes obrms report "inf" for every pose with no error
        // anywhere - silent enough that it's worth calling out explicitly
        // rather than just leaving the RMSD column blank. See backend/rmsd.py.
        const wantedRmsd = job.params.reference_ligand_pdb;
        const doneLigs = job.ligands.filter((l) => l.status === "done");
        const gotNoRmsd = doneLigs.length > 0 && doneLigs.every(
          (l) => !l.rmsd_to_reference || l.rmsd_to_reference.every((r) => r === null));
        if (wantedRmsd && gotNoRmsd) {
          toast("RMSD validation was on, but no docked pose matched the bound ligand's heavy-atom structure - double check the SMILES you docked is really the same molecule as the crystal ligand", true);
        }
      }
    } catch (err) {
      clearInterval(S.poll);
      toast(err.message, true);
    }
  }, 1500);
}

// The badge summarizes the top-ranked pose's pharmacophore check (see
// backend/jobs.py:_run_one); the title lists each constraint's distance so
// hovering explains "why" without needing the log.
function constraintBadge(lig) {
  if (lig.constraint_satisfied === null || lig.constraint_satisfied === undefined) return "";
  const good = lig.constraint_satisfied;
  const detail = (lig.constraint_detail || [])
    .map((d) => `${d.label}: ${d.distance === null ? "no matching atom" : d.distance.toFixed(2) + " Å"}${d.required ? "" : " (optional)"}`)
    .join("\n");
  return `<span class="tag ${good ? "done" : "unset"}" title="${detail.replace(/"/g, "&quot;")}">${good ? "✓ met" : "✗ unmet"}</span>`;
}


const escapeHtml = (text) => String(text).replace(/[&<>"']/g,
  (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[c]);

// MOE's per-pose SD fields, in pipeline order (mirrors docking.MOE_SCORE_FIELDS).
const MOE_SCORE_COLUMNS = [
  ["S", "final score (rescoring function)"],
  ["E_place", "placement score"],
  ["E_score1", "London dG after placement"],
  ["E_refine", "refinement energy"],
  ["E_score2", "rescoring function after refinement"],
  ["E_conf", "ligand conformer energy"],
  ["rmsd_refine", "RMSD moved during refinement (Å)"],
];

// Per-pose table for a MOE result: every MOE score, plus strain, RMSD to the
// bound ligand and (covalent) the group each pose is bonded through.
function moePoseTable(lig) {
  const fmt = (v, d = 2) => (v === null || v === undefined ? "" : v.toFixed(d));
  const hasStrain = (lig.strain_energies || []).some((v) => v !== null && v !== undefined);
  const hasRmsd = (lig.rmsd_to_reference || []).some((v) => v !== null && v !== undefined);
  const hasAttach = (lig.covalent_attachments || []).some(Boolean);
  const head = ["<th>pose</th>",
    ...MOE_SCORE_COLUMNS.map(([k, help]) => `<th title="${help}">${k}</th>`),
    hasStrain ? '<th title="ligand-only MMFF strain (kcal/mol)">strain</th>' : "",
    hasRmsd ? '<th title="RMSD to the bound ligand (Å)">RMSD</th>' : "",
    hasAttach ? '<th title="ligand atoms around the bond to the residue (* = residue atom)">bonded via</th>' : "",
  ].join("");
  const body = lig.pose_scores.map((ps, i) => `<tr>
      <td>#${i + 1}</td>
      ${MOE_SCORE_COLUMNS.map(([k]) => `<td>${fmt(ps[k])}</td>`).join("")}
      ${hasStrain ? `<td>${fmt(lig.strain_energies[i], 1)}</td>` : ""}
      ${hasRmsd ? `<td>${fmt(lig.rmsd_to_reference[i])}</td>` : ""}
      ${hasAttach ? `<td class="mono">${escapeHtml(lig.covalent_attachments[i] || "")}</td>` : ""}
    </tr>`).join("");
  return `<table class="moescores"><thead><tr>${head}</tr></thead><tbody>${body}</tbody></table>`;
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
    const isMoe = lig.pose_scores && lig.pose_scores.length > 0;
    const warnings = lig.warnings || [];
    const hasBreakdown = (lig.scores && lig.scores.length > 1) || isMoe || warnings.length > 0;
    const expanded = S.expandedResultRows.has(lig.name);
    // rmsd_to_reference[0] is the top-ranked pose's RMSD to the receptor's
    // own bound-ligand position (backend/rmsd.py) - only present for a
    // validation run, and only non-null for the ligand that matched it.
    const topRmsd = lig.rmsd_to_reference && lig.rmsd_to_reference[0];
    const rmsdText = topRmsd === null || topRmsd === undefined ? "" :
      `<span class="tag ${topRmsd <= 2.0 ? "done" : "unset"}">${topRmsd.toFixed(2)}</span>`;
    tr.innerHTML = `
      <td title="${lig.smiles || ""}">${lig.name}</td>
      <td><span class="tag ${lig.status}">${lig.status}</span>${warnings.length ? ` <span class="tag unset" title="${escapeHtml(warnings.join("\n"))}">⚠</span>` : ""}</td>
      <td class="score">${score}${hasBreakdown ? ` <button class="ghost mini posetoggle" title="${isMoe ? "show every pose's MOE scores" : "show every pose's score/strain"}">${expanded ? "▴" : "▾"}</button>` : ""}</td>
      <td>${lig.n_poses || ""}</td>
      <td>${constraintBadge(lig)}</td>
      <td>${rmsdText}</td>
      <td></td>`;
    if (lig.error) tr.querySelector("td:nth-child(2)").title = lig.error;

    if (hasBreakdown) {
      tr.querySelector(".posetoggle").onclick = (e) => {
        e.stopPropagation();
        if (S.expandedResultRows.has(lig.name)) S.expandedResultRows.delete(lig.name);
        else S.expandedResultRows.add(lig.name);
        renderJob(job);
      };
    }

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

    if (hasBreakdown && expanded) {
      const detailTr = document.createElement("tr");
      detailTr.className = "posedetail";
      const warnHtml = warnings.map((w) => `<div class="warnline">⚠ ${escapeHtml(w)}</div>`).join("");
      const cells = isMoe ? moePoseTable(lig) : lig.scores.map((s, i) => {
        const strain = lig.strain_energies && lig.strain_energies[i];
        const strainText = strain === null || strain === undefined ? "" : ` <span class="muted">strain ${strain.toFixed(1)}</span>`;
        const r = lig.rmsd_to_reference && lig.rmsd_to_reference[i];
        const rmsdText = r === null || r === undefined ? "" : ` <span class="muted">rmsd ${r.toFixed(2)} Å</span>`;
        return `<span class="posecell">#${i + 1} ${s.toFixed(2)} kcal/mol${strainText}${rmsdText}</span>`;
      }).join("");
      detailTr.innerHTML = `<td colspan="7">${warnHtml}${cells}</td>`;
      tbody.appendChild(detailTr);
    }
  });
}

// ---------------------------------------------------------------------------
// pose display
// ---------------------------------------------------------------------------

async function clearPoses() {
  stopPosePlayback(false);
  await deleteRefs(interactionRefs);
  interactionRefs = [];
  await deleteRefs(poseStructureRefs);
  poseStructureRefs = [];
  poseStructures = [];
  S.poseReps = [];
  $("poseBar").classList.add("hidden");
  hideLigandDepiction();
  if (S.diagramUrl) { URL.revokeObjectURL(S.diagramUrl); S.diagramUrl = null; }
  $("diagramImg").classList.add("hidden");
  $("diagramEmpty").textContent = "select a pose to see its interaction diagram";
  $("diagramEmpty").classList.remove("hidden");
}

async function showPoses(lig) {
  if (!S.job) return;
  try {
    const sdf = await api(`/api/jobs/${S.job.id}/ligands/${encodeURIComponent(lig.name)}/poses.sdf`);
    await clearPoses();
    S.selectedLigand = lig.name;
    S.selectedLigandObj = lig;

    const blocks = sdf.split("$$$$").map((b) => b.trim()).filter((b) => b.length > 10);
    for (const block of blocks) {
      const data = await plugin.builders.data.rawData({ data: block + "\n$$$$\n", label: lig.name });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, "sdf");
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      // atom-type (element-symbol / "CPK") coloring, same scheme molecular
      // viewers usually mean by "color by atom type"
      const rep = await addRep(structure, "ball-and-stick", "element-symbol", undefined, { sizeFactor: 0.16 });
      poseStructureRefs.push(data, trajectory, model, structure);
      poseStructures.push(structure);
      S.poseReps.push(rep);
    }

    const bar = $("poseBar");
    bar.innerHTML = `<span class="title">${lig.name} - poses:</span>`;
    const stepBtn = (text, title, onclick, id) => {
      const b = document.createElement("button");
      if (id) b.id = id;
      b.className = "pose-play";
      b.textContent = text;
      b.title = title;
      b.onclick = onclick;
      bar.appendChild(b);
    };
    stepBtn("\u23EE", "previous pose", () => { stopPosePlayback(false); stepPose(-1); });
    stepBtn("\u25B6", "play through the docked poses", togglePosePlayback, "posePlayBtn");
    stepBtn("\u23ED", "next pose", () => { stopPosePlayback(false); stepPose(1); });
    S.poseReps.forEach((_, i) => {
      const b = document.createElement("button");
      const score = lig.scores[i];
      b.className = "pose";
      b.textContent = `${i + 1}${score !== undefined ? `  ${score.toFixed(2)}` : ""}`;
      b.onclick = () => { stopPosePlayback(false); selectPose(i); };
      bar.appendChild(b);
    });
    const all = document.createElement("button");
    all.className = "pose";
    all.textContent = "all";
    all.onclick = () => { stopPosePlayback(false); selectPose(-1); };
    bar.appendChild(all);
    bar.classList.remove("hidden");

    selectPose(0);
    showLigandDepiction(lig);
    renderJob(S.job);
  } catch (err) {
    toast(err.message, true);
  }
}

// `light` skips the per-pose contacts + 2D diagram (both hit the backend /
// rebuild structures) so playback can step quickly; they're refreshed once
// playback stops.
function selectPose(idx, { light = false } = {}) {
  S.activePose = idx;
  S.poseReps.forEach((rep, i) => {
    if (!rep) return;
    plugin.state.data.updateCellState(rep.ref, { isHidden: !(idx === -1 || i === idx) });
  });
  document.querySelectorAll("#poseBar .pose").forEach((b, i) => {
    b.classList.toggle("active", i === (idx === -1 ? S.poseReps.length : idx));
  });
  if (light) return;
  // "all" (idx -1) has no single pose to show contacts for; a specific pose does.
  showInteractions(idx).catch((err) => console.warn("interactions failed", err));
  showInteractionDiagram(idx).catch((err) => console.warn("2d diagram failed", err));
}

// Moves the active pose by `delta`, wrapping; from "all" it lands on the
// first (forward) or last (back) pose.
function stepPose(delta, opts) {
  const n = S.poseReps.length;
  if (!n) return;
  const next = S.activePose < 0 ? (delta > 0 ? 0 : n - 1) : (S.activePose + delta + n) % n;
  selectPose(next, opts);
}

const POSE_PLAY_INTERVAL_MS = 1200;

function togglePosePlayback() {
  if (S.posePlayTimer) { stopPosePlayback(true); return; }
  const n = S.poseReps.length;
  if (n < 2) return;
  const step = () => stepPose(1, { light: true });
  step();
  S.posePlayTimer = setInterval(step, POSE_PLAY_INTERVAL_MS);
  const btn = $("posePlayBtn");
  if (btn) { btn.textContent = "\u275A\u275A"; btn.title = "pause"; btn.classList.add("active"); }
}

// `refresh` re-selects the pose playback stopped on so its contacts and 2D
// diagram (skipped while cycling) get drawn.
function stopPosePlayback(refresh) {
  if (!S.posePlayTimer) return;
  clearInterval(S.posePlayTimer);
  S.posePlayTimer = null;
  const btn = $("posePlayBtn");
  if (btn) { btn.textContent = "\u25B6"; btn.title = "play through the docked poses"; btn.classList.remove("active"); }
  if (refresh) selectPose(S.activePose);
}

// Small 2D structure of the docked ligand pinned to the top-left of Mol*'s 3D
// canvas. Mounted inside .msp-viewport (not #viewer, which also holds Mol*'s
// side panels) and re-attached each time, since Mol* can remount the
// viewport on layout changes. Same <img>-not-innerHTML reasoning as the
// interaction diagram below.
function ligandDepictionEl() {
  let el = $("ligandDepiction");
  if (!el) {
    el = document.createElement("div");
    el.id = "ligandDepiction";
    el.className = "hidden";
    el.innerHTML = '<div class="depname"></div><img alt="2D structure of the docked ligand">';
  }
  const viewport = $("viewer").querySelector(".msp-viewport") || $("viewer");
  if (el.parentElement !== viewport) viewport.appendChild(el);
  return el;
}

async function showLigandDepiction(lig) {
  const el = ligandDepictionEl();
  const img = el.querySelector("img");
  try {
    const res = await fetch(`/api/jobs/${S.job.id}/ligands/${encodeURIComponent(lig.name)}/depiction.svg`);
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    const blob = await res.blob();
    if (S.selectedLigand !== lig.name) return; // another ligand was picked meanwhile
    if (S.depictionUrl) URL.revokeObjectURL(S.depictionUrl);
    S.depictionUrl = URL.createObjectURL(blob);
    img.src = S.depictionUrl;
    el.querySelector(".depname").textContent = lig.name;
    el.classList.remove("hidden");
  } catch (err) {
    console.warn("2d depiction failed", err);
    hideLigandDepiction();
  }
}

function hideLigandDepiction() {
  const el = $("ligandDepiction");
  if (el) el.classList.add("hidden");
  if (S.depictionUrl) { URL.revokeObjectURL(S.depictionUrl); S.depictionUrl = null; }
}

// Fetches the backend's ProLIF-detected, LigPlot-style 2D diagram for one pose
// and swaps it into #diagramImg. Rendered via a blob URL on an <img>, not
// innerHTML'd as inline SVG - the backend escapes PDB-derived text before it
// reaches the SVG, but an <img>-rendered SVG can't execute embedded scripts
// even if that escaping ever had a gap, which innerHTML wouldn't guarantee.
async function showInteractionDiagram(idx) {
  const img = $("diagramImg"), empty = $("diagramEmpty");
  if (S.diagramUrl) { URL.revokeObjectURL(S.diagramUrl); S.diagramUrl = null; }
  if (idx < 0 || S.diagramHidden || !S.job || !S.selectedLigandObj) {
    img.classList.add("hidden");
    empty.textContent = "select a pose to see its interaction diagram";
    empty.classList.remove("hidden");
    return;
  }
  empty.textContent = "loading diagram…";
  empty.classList.remove("hidden");
  img.classList.add("hidden");
  try {
    const name = encodeURIComponent(S.selectedLigandObj.name);
    const res = await fetch(`/api/jobs/${S.job.id}/ligands/${name}/poses/${idx + 1}/interactions2d.svg`);
    if (!res.ok) {
      let detail = res.statusText;
      try { detail = (await res.json()).detail || detail; } catch (_) {}
      throw new Error(detail);
    }
    const svgText = await res.text();
    const blob = new Blob([svgText], { type: "image/svg+xml" });
    S.diagramUrl = URL.createObjectURL(blob);
    img.src = S.diagramUrl;
    img.classList.remove("hidden");
    empty.classList.add("hidden");
  } catch (err) {
    empty.textContent = "diagram unavailable";
    empty.classList.remove("hidden");
    toast(`interaction diagram: ${err.message}`, true);
  }
}

// Non-covalent contacts (H-bonds, hydrophobic, ionic, ...) between the given
// pose and the receptor - drawn as dashed lines, colored by interaction
// type. Computed on a merged copy of (receptor + that one pose) since Mol*'s
// interaction analysis only looks within a single Structure.
async function showInteractions(idx) {
  await deleteRefs(interactionRefs);
  interactionRefs = [];
  if (idx < 0 || !receptorStructure || !poseStructures[idx]) return;
  try {
    // With "hide ligand" on, merge in a ligand/ion-free copy of the receptor
    // instead: otherwise the co-crystallised ligand (which usually overlaps
    // the pose) comes back through the merged structure - drawn as a
    // "contact residue" and picked up by the interaction lines.
    let receptorSrc = receptorStructure;
    if (S.hideOrigLigand) {
      const strippedExpr = MS.struct.modifier.exceptBy({
        0: MS.struct.generator.all(),
        by: MS.struct.combinator.merge([Queries.ligand.expression, Queries.ion.expression]),
      });
      const stripped = await plugin.builders.structure.tryCreateComponentFromExpression(
        receptorStructure, strippedExpr, "receptor-no-ligand");
      if (stripped) { interactionRefs.push(stripped); receptorSrc = stripped; }
    }
    const structures = [{ ref: receptorSrc.ref }, { ref: poseStructures[idx].ref }];
    const merged = await plugin.state.data.build().toRoot()
      .apply(MergeStructures, { structures }, { dependsOn: structures.map((s) => s.ref) })
      .commit();

    // Queries.ligand also catches the receptor's own bound ligand(s), if any
    // (e.g. the co-crystallised BTN in 1STP) - exclude those so this is
    // specifically the docked pose we just merged in.
    const originalLigandExprs = (S.receptor?.het_ligands || []).map(hetExpr);
    const poseExpr = originalLigandExprs.length
      ? MS.struct.modifier.exceptBy({ 0: Queries.ligand.expression, by: MS.struct.combinator.merge(originalLigandExprs) })
      : Queries.ligand.expression;

    const poseComp = await plugin.builders.structure.tryCreateComponentFromExpression(merged, poseExpr, "pose-ligand");
    await addRep(poseComp, InteractionsRepresentationProvider, InteractionTypeColorThemeProvider,
      undefined, { includeParent: true, parentDisplay: "between" });

    // Residues close enough to plausibly interact, as sticks - the base
    // cartoon representation doesn't draw side chains, so without this the
    // dashed contact lines above have nothing visible to point at.
    const nearbyExpr = MS.struct.modifier.exceptBy({
      0: MS.struct.modifier.includeSurroundings({ 0: poseExpr, radius: 4.5, "as-whole-residues": true }),
      by: poseExpr,
    });
    const nearbyComp = await plugin.builders.structure.tryCreateComponentFromExpression(merged, nearbyExpr, "pose-contact-residues");
    await addRep(nearbyComp, "ball-and-stick", "element-symbol", undefined, { sizeFactor: 0.18 });

    interactionRefs.push(merged);
  } catch (err) {
    console.warn("interactions failed", err);
  }
}

// ---------------------------------------------------------------------------
// scoring dropdown helper
// ---------------------------------------------------------------------------

function updateScoringOptions(engine) {
  const scoringSelect = $("scoring");
  scoringSelect.innerHTML = "";

  // For MOE this picks the rescoring function applied after refinement
  // (placement is always London dG, as in MOE's own Dock panel).
  $("scoringLabel").textContent = engine === "moe" ? "rescoring" : "scoring";
  $("scoring").title = engine === "moe" ? "MOE rescoring after refinement - placement always uses London dG" : "";
  if (engine === "moe" && window.APP_STATUS?.moe_scoring_functions) {
    window.APP_STATUS.moe_scoring_functions.forEach(fn => {
      const opt = document.createElement("option");
      opt.value = fn;
      opt.textContent = fn;
      scoringSelect.appendChild(opt);
    });
  } else if (window.APP_STATUS?.scoring_functions) {
    window.APP_STATUS.scoring_functions.forEach(fn => {
      const opt = document.createElement("option");
      opt.value = fn;
      opt.textContent = fn;
      scoringSelect.appendChild(opt);
    });
    scoringSelect.value = "vina";
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

  $("peptideMode").onchange = (e) => {
    const params = e.target.checked ? PEPTIDE_DOCKING_PARAMS : STANDARD_DOCKING_PARAMS;
    $("exh").value = params.exhaustiveness;
    $("erange").value = params.energy_range;
  };

  $("clearPicks").onclick = clearPicks;
  $("selfDockBtn").onclick = selfDockFromBoundLigand;
  $("selfDockManualBtn").onclick = selfDockWithManualSmiles;
  $("selfDockSmiles").onkeydown = (e) => { if (e.key === "Enter") $("selfDockManualBtn").click(); };
  $("addConstraint").onclick = addConstraintRow;
  $("hideLigandBtn").onclick = () => {
    S.hideOrigLigand = !S.hideOrigLigand;
    $("hideLigandBtn").classList.toggle("active", S.hideOrigLigand);
    styleReceptor();
    showInteractions(S.activePose).catch((err) => console.warn("interactions failed", err));
  };
  // viewer-header shortcut for the receptor card's "delete existing ligand(s)
  // before docking" checkbox - the two stay in sync either way
  const syncDeleteLigandBtn = () => $("deleteLigandBtn").classList.toggle("active", $("removeLigand").checked);
  $("deleteLigandBtn").onclick = () => {
    $("removeLigand").checked = !$("removeLigand").checked;
    syncDeleteLigandBtn();
  };
  $("removeLigand").addEventListener("change", syncDeleteLigandBtn);
  syncDeleteLigandBtn();
  $("hideProteinBtn").onclick = () => {
    S.hideProtein = !S.hideProtein;
    $("hideProteinBtn").classList.toggle("active", S.hideProtein);
    styleReceptor();
  };
  $("toggleDiagramBtn").onclick = () => {
    S.diagramHidden = !S.diagramHidden;
    $("toggleDiagramBtn").textContent = S.diagramHidden ? "+" : "_";
    $("diagramBody").classList.toggle("hidden", S.diagramHidden);
    if (!S.diagramHidden) showInteractionDiagram(S.activePose).catch((err) => console.warn("2d diagram failed", err));
  };
  const panelToggle = (btnId, bodyId, stateKey) => {
    $(btnId).onclick = () => {
      S[stateKey] = !S[stateKey];
      $(btnId).textContent = S[stateKey] ? "+" : "_";
      $(bodyId).classList.toggle("hidden", S[stateKey]);
    };
  };
  panelToggle("toggleReceptorBtn", "receptorBody", "receptorPanelHidden");
  panelToggle("toggleActiveSiteBtn", "activeSiteBody", "activeSitePanelHidden");
  panelToggle("togglePharmacophoreBtn", "pharmacophoreBody", "pharmacophorePanelHidden");
  panelToggle("toggleCovalentBtn", "covalentBody", "covalentPanelHidden");
  panelToggle("toggleLigandsBtn", "ligandsBody", "ligandsPanelHidden");

  $("cubeBtn").onclick = () => {
    const s = boxSize().filter(Number.isFinite);
    if (!s.length) return;
    const edge = Math.max(...s);
    ["sx", "sy", "sz"].forEach((id) => ($(id).value = edge.toFixed(1)));
    drawBox();
  };

  ["cx", "cy", "cz", "sx", "sy", "sz"].forEach((id) => ($(id).oninput = drawBox));
  $("zoomSiteBtn").onclick = () => {
    if (S.activeHet) zoomToExpression(hetExpr(S.activeHet));
    else if (S.picks.length) toast("zoom not available for free-hand picks yet - drag/scroll to frame the box", true);
    else toast("set a docking box first", true);
  };

  $("moePh4File").onchange = async (e) => {
    const file = e.target.files[0];
    if (!file) return;
    const fd = new FormData();
    fd.append("file", file);
    try {
      const data = await api("/api/pharmacophore/upload", { method: "POST", body: fd });
      S.ph4_id = data.ph4_id;
      $("moePh4Name").textContent = data.filename;
      toast(`uploaded ${data.filename}`);
    } catch (err) {
      toast(err.message || "failed to upload .ph4 file", true);
    }
  };

  $("customRxnFile").onchange = async (e) => {
    const file = e.target.files[0];
    if (!file) return;
    const fd = new FormData();
    fd.append("file", file);
    try {
      const data = await api("/api/reactions/upload", { method: "POST", body: fd });
      S.custom_rxn_id = data.rxn_id;
      $("customRxnName").textContent = data.filename;
      // Add to dropdown as custom option
      const sel = $("covalentReaction");
      // Remove old custom option if exists
      const oldCustom = sel.querySelector('option[data-custom="true"]');
      if (oldCustom) oldCustom.remove();
      const opt = document.createElement("option");
      opt.value = "";  // Empty - custom_rxn_id used instead
      opt.textContent = `Custom: ${data.filename}`;
      opt.dataset.custom = "true";
      sel.insertBefore(opt, sel.firstChild.nextSibling);
      sel.value = "";
      toast(`uploaded ${data.filename}`);
    } catch (err) {
      toast(err.message || "failed to upload .rxn file", true);
    }
  };

  $("pickCovalentAtomBtn").onclick = () => {
    S.covalentAtomPickActive = !S.covalentAtomPickActive;
    $("pickCovalentAtomBtn").classList.toggle("armed", S.covalentAtomPickActive);
    if (S.covalentAtomPickActive) {
      toast("Click a nucleophilic atom (CYS SG, SER OG, LYS NZ, etc.) in the viewer");
    } else {
      toast("Atom picking cancelled");
    }
  };

  $("findResiduesBtn").onclick = async () => {
    if (!S.receptor) return toast("load receptor first", true);
    const resname = $("residueSearchInput").value.trim().toUpperCase();
    if (!resname || !/^[A-Z]{3}$/.test(resname)) {
      return toast("enter 3-letter residue code (e.g., CYS, SER, LYS)", true);
    }
    const btn = $("findResiduesBtn");
    btn.disabled = true;
    try {
      const data = await api(`/api/receptor/${S.receptor.id}/residues/${resname}`);
      const residueList = $("residueList");
      residueList.innerHTML = "";
      if (!data.residues.length) {
        $("residueInfo").textContent = `no ${resname} residues found`;
        return;
      }
      $("residueInfo").textContent = `${data.residues.length} ${resname} residue(s) - click to use`;
      data.residues.forEach((spec) => {
        const div = document.createElement("div");
        div.className = "residue-item";
        div.textContent = spec;
        div.onclick = () => {
          $("covalentResidue").value = spec;
          toast(`selected ${spec}`);
        };
        residueList.appendChild(div);
      });
    } catch (err) {
      toast(err.message, true);
    }
    btn.disabled = false;
  };

  $("residueSearchInput").onkeydown = (e) => {
    if (e.key === "Enter") $("findResiduesBtn").click();
  };

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
    const engine = $("engineSel").value;
    const isMoe = engine === "moe";
    const isGnina = engine === "gnina";

    // Update engine info text
    const engineInfoTexts = {
      smina: "CPU-based Vina fork - good for general docking",
      gnina: "CNN scoring + GPU support - better pose prediction",
      moe: "Commercial engine - supports covalent docking, pharmacophore constraints"
    };
    $("engineInfo").textContent = engineInfoTexts[engine] || "";

    // Show/hide GPU row for gnina
    $("gpuRow").style.display = isGnina ? "flex" : "none";

    // Show/hide MOE-specific params
    $("moeParamsRow").style.display = isMoe ? "grid" : "none";

    // Show/hide covalent docking panel (MOE only)
    $("covalentCard").style.display = isMoe ? "block" : "none";
    if (!isMoe) {
      $("covalentEnabled").checked = false;
      $("covalentParams").style.display = "none";
    }

    // Hide smina/gnina params for MOE
    $("sminaGninaParamsRow").style.display = isMoe ? "none" : "flex";
    $("sminaGninaMinimizeRow").style.display = isMoe ? "none" : "flex";

    // Show/hide pharmacophore constraints based on engine
    $("sminaGninaConstraints").style.display = isMoe ? "none" : "block";
    $("moeConstraints").style.display = isMoe ? "block" : "none";

    // Update scoring dropdown based on engine
    updateScoringOptions(engine);
  };

  // Covalent checkbox toggles params
  $("covalentEnabled").onchange = () => {
    $("covalentParams").style.display = $("covalentEnabled").checked ? "block" : "none";
  };
}

async function boot() {
  await initViewer();
  wire();
  try {
    const st = await api("/api/status");
    window.APP_STATUS = st;  // Store for updateScoringOptions

    // Populate engine selector based on available engines
    const engineSel = $("engineSel");
    engineSel.innerHTML = "";
    if (st.smina) {
      const opt = document.createElement("option");
      opt.value = "smina";
      opt.textContent = "smina (CPU)";
      engineSel.appendChild(opt);
    }
    // Always show Gnina option
    const optGnina = document.createElement("option");
    optGnina.value = "gnina";
    optGnina.textContent = st.gnina ? "gnina (CPU/GPU)" : "gnina (not installed)";
    optGnina.disabled = !st.gnina;
    engineSel.appendChild(optGnina);

    if (st.moe_available) {
      const opt = document.createElement("option");
      opt.value = "moe";
      opt.textContent = "MOE";
      engineSel.appendChild(opt);
    }

    // Initial scoring options for default engine
    updateScoringOptions(engineSel.value);

    $("modes").value = st.default_num_modes;

    const pill = $("engine");
    const parts = [];
    parts.push(st.smina ? (st.smina_version || "smina") : "smina not found");
    parts.push(st.gnina ? (st.gnina_version || "gnina") : "gnina not found");
    if (st.moe_available) {
      parts.push("MOE available");
    }
    pill.className = "pill " + (st.smina || st.gnina || st.moe_available ? "ok" : "bad");
    pill.textContent = `${parts.join(" - ")} - ${st.workers} workers${st.gpu_available ? " - GPU detected" : ""}`;

    $("gpuInfo").textContent = st.gpu_available
      ? "NVIDIA GPU detected on this machine"
      : "no NVIDIA GPU detected - gnina will run on CPU";

    // Populate covalent reactions dropdown
    const covalentReactionSel = $("covalentReaction");
    covalentReactionSel.innerHTML = '<option value="">Select reaction...</option>';
    if (st.covalent_reactions && st.covalent_reactions.length > 0) {
      st.covalent_reactions.forEach((rxn) => {
        const opt = document.createElement("option");
        opt.value = rxn;
        opt.textContent = rxn.replace(/_/g, " ");
        covalentReactionSel.appendChild(opt);
      });
    }

    // Trigger initial engine UI update
    $("engineSel").dispatchEvent(new Event("change"));
  } catch (err) {
    $("engine").className = "pill bad";
    $("engine").textContent = "backend unreachable";
  }
}

// Read-only debug hook (mirrors the pattern Mol*'s own apps use, e.g.
// window.DockingViewer) - lets the browser console or a test script inspect
// live state without poking at internals.
window.__docking = { get S() { return S; }, get plugin() { return plugin; } };

boot();
