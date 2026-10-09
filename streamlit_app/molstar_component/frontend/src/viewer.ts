/* ViewerController: owns the single Mol* plugin and every state ref the
 * component builds on top of it (receptor + its styling components, docking
 * box, pick/constraint markers, docked poses, per-pose interactions).
 *
 * Ported from the original frontend/src/app.js viewer code. All mutations go
 * through `enqueue` so overlapping Streamlit renders can never interleave
 * half-built Mol* state.
 */

import { createPluginUI } from "molstar/lib/mol-plugin-ui/index";
import { renderReact18 } from "molstar/lib/mol-plugin-ui/react18";
import { DefaultPluginUISpec } from "molstar/lib/mol-plugin-ui/spec";
import type { PluginUIContext } from "molstar/lib/mol-plugin-ui/context";
import { PluginCommands } from "molstar/lib/mol-plugin/commands";
import { PluginConfig } from "molstar/lib/mol-plugin/config";
import { Structure, StructureElement, StructureProperties } from "molstar/lib/mol-model/structure";
import { StructureSelectionQueries as Queries } from "molstar/lib/mol-plugin-state/helpers/structure-selection-query";
import { PluginStateObject as PSO, PluginStateTransform } from "molstar/lib/mol-plugin-state/objects";
import { InteractionsRepresentationProvider } from "molstar/lib/mol-model-props/computed/representations/interactions";
import { InteractionTypeColorThemeProvider } from "molstar/lib/mol-model-props/computed/themes/interaction-type";
import { MolScriptBuilder as MS } from "molstar/lib/mol-script/language/builder";
import type { Expression } from "molstar/lib/mol-script/language/expression";
import { Color } from "molstar/lib/mol-util/color";
import { Vec3 } from "molstar/lib/mol-math/linear-algebra";
import { StateObject } from "molstar/lib/mol-state";
import { Task } from "molstar/lib/mol-task";
import { ParamDefinition as PD } from "molstar/lib/mol-util/param-definition";
import "molstar/lib/mol-util/polyfill";

import {
  SYNTHETIC_RESNAMES, boxPdb, constraintColor, constraintPdb, isPlaced, isValidBox,
  picksPdb, ringFor, ringGeometry,
} from "./synthetic";
import type {
  BoxArg, ClickedAtom, ConstraintMarker, FocusArg, HetLigand, HetRef, PosesArg, ReceptorArg,
  StyleArg, Vec3Tuple,
} from "./types";

// Any state-tree ref/selector we later delete; Mol* selectors and cells both
// expose `.ref`.
type Ref = { ref: string };

// Combines two already-loaded structures (receptor + one docked pose) into a
// single Structure so Mol*'s non-covalent-interactions computation - which
// only looks within one Structure - can see contacts between them. Adapted
// from Mol*'s own src/apps/docking-viewer/index.ts.
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
      const first = dependencies![params.structures[0].ref].data as Structure;
      const builder = Structure.Builder({ masterModel: first.models[0] });
      for (const { ref } of params.structures) {
        const s = dependencies![ref].data as Structure;
        for (const unit of s.units) {
          builder.addUnit(unit.kind, unit.model, unit.conformation.operator, unit.elements, unit.traits);
        }
      }
      return new PSO.Molecule.Structure(builder.getStructure(), { label: "Merged Structure" });
    });
  },
});

// MolScript expression selecting one het group by chain/resname/resno
// (auth_* fields, matching backend/prep.py's raw PDB columns).
export function hetExpr(het: HetRef): Expression {
  const tests: Record<string, Expression> = {
    "residue-test": MS.core.logic.and([
      MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_comp_id(), het.resname]),
      MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_seq_id(), parseInt(String(het.resseq), 10)]),
    ]),
  };
  if (het.chain && het.chain !== "_") {
    tests["chain-test"] = MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_asym_id(), het.chain]);
  }
  return MS.struct.generator.atomGroups(tests);
}

function residueExpr(chain: string, resseq: number): Expression {
  return MS.struct.generator.atomGroups({
    "chain-test": MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_asym_id(), chain]),
    "residue-test": MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_seq_id(), resseq]),
  });
}

const ligandsAndIons = () => MS.struct.combinator.merge([Queries.ligand.expression, Queries.ion.expression]);

export class ViewerController {
  plugin!: PluginUIContext;
  private queue: Promise<unknown> = Promise.resolve();

  private receptorStructure: any = null;  // StateObjectSelector of the receptor Structure
  private receptorComponentRefs: Ref[] = [];
  private boxRefs: Ref[] = [];
  private picksRefs: Ref[] = [];
  private constraintRefs: Ref[] = [];
  private poseStructureRefs: Ref[] = [];
  private poseStructures: any[] = [];
  private poseReps: (Ref | null)[] = [];
  private interactionRefs: Ref[] = [];

  private style: StyleArg = { hide_ligand: false, hide_protein: false, active_het: null, prepared: false };
  private hetLigands: HetLigand[] = [];

  constructor(private onAtomClick: (atom: ClickedAtom) => void) {}

  // -- lifecycle ------------------------------------------------------------

  async init(target: HTMLElement): Promise<void> {
    // The full Mol* app, minus the "remote state" feature that would save
    // sessions to a public Mol* server - everything else here stays local.
    const defaultSpec = DefaultPluginUISpec();
    const spec = {
      ...defaultSpec,
      components: { ...defaultSpec.components, remoteState: "none" as const },
      config: [
        ...(defaultSpec.config || []),
        [PluginConfig.State.DefaultServer, ""],
        [PluginConfig.State.CurrentServer, ""],
      ] as [any, unknown][],
    };
    this.plugin = await createPluginUI({ target, spec, render: renderReact18 });

    await PluginCommands.Canvas3D.SetSettings(this.plugin, {
      settings: {
        renderer: { ...this.plugin.canvas3d!.props.renderer, backgroundColor: Color(0x000000) },
      },
    });
    await this.plugin.managers.structure.component.setOptions({
      ...this.plugin.managers.structure.component.state.options,
      hydrogens: "only-polar",
    });
    this.plugin.behaviors.interaction.click.subscribe((e) => this.handleClick(e));
  }

  /** Run `fn` after every previously queued mutation finished. */
  enqueue<T>(fn: () => Promise<T>): Promise<T> {
    const next = this.queue.then(fn, fn);
    this.queue = next.catch((err) => console.warn("viewer update failed", err));
    return next;
  }

  /** The element the 2D depiction overlay is pinned into. */
  viewportElement(fallback: HTMLElement): HTMLElement {
    return (fallback.querySelector(".msp-viewport") as HTMLElement | null) || fallback;
  }

  get poseCount(): number {
    return this.poseReps.length;
  }

  // -- helpers --------------------------------------------------------------

  private async deleteRefs(refs: Ref[]): Promise<void> {
    const live = refs.filter(Boolean);
    if (!live.length) return;
    const b = this.plugin.state.data.build();
    live.forEach((r) => b.delete(r.ref));
    await b.commit();
  }

  // Raw addRepresentation ignores the viewer's "Hydrogens" option (only
  // presets apply it), so mirror it for the built-in representation types.
  private hydrogenTypeParams() {
    const h = this.plugin.managers.structure.component.state.options.hydrogens;
    return { ignoreHydrogens: h !== "all", ignoreHydrogensVariant: h === "only-polar" ? "non-polar" : "all" };
  }

  private async addRep(componentRef: any, type: any, color: any, colorParams?: any, typeParams?: any): Promise<Ref | null> {
    if (!componentRef) return null;
    if (typeof type === "string") typeParams = { ...this.hydrogenTypeParams(), ...typeParams };
    return this.plugin.builders.structure.representation.addRepresentation(componentRef, {
      type, color, colorParams, typeParams,
    } as any);
  }

  private async loadStructure(text: string, format: "pdb" | "sdf", label: string) {
    const data = await this.plugin.builders.data.rawData({ data: text, label });
    const trajectory = await this.plugin.builders.structure.parseTrajectory(data, format);
    const model = await this.plugin.builders.structure.createModel(trajectory);
    const structure = await this.plugin.builders.structure.createStructure(model);
    return { structure, refs: [data, trajectory, model, structure] as Ref[] };
  }

  private setHidden(rep: Ref | null | undefined, hidden: boolean): void {
    if (rep) this.plugin.state.data.updateCellState(rep.ref, { isHidden: hidden });
  }

  // -- clicks ---------------------------------------------------------------

  private ringCentroidAt(loc: StructureElement.Location, resn: string, atomName: string) {
    const ring = ringFor(resn, atomName);
    if (!ring) return null;
    const chain = StructureProperties.chain.auth_asym_id(loc);
    const resi = StructureProperties.residue.auth_seq_id(loc);
    const probe = StructureElement.Location.create(loc.structure);
    probe.unit = loc.unit;
    const p = Vec3.zero();
    const pts: Vec3Tuple[] = [];
    const elements = loc.unit.elements;
    for (let i = 0; i < elements.length; i++) {
      const el = elements[i];
      probe.element = el;
      if (StructureProperties.chain.auth_asym_id(probe) !== chain) continue;
      if (StructureProperties.residue.auth_seq_id(probe) !== resi) continue;
      if (!ring.includes(StructureProperties.atom.label_atom_id(probe))) continue;
      loc.unit.conformation.position(el, p);
      pts.push([p[0], p[1], p[2]]);
    }
    return ringGeometry(pts, ring.length);
  }

  private handleClick(event: any): void {
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
      this.onAtomClick({
        resn,
        hetflag: StructureProperties.residue.group_PDB(loc) === "HETATM",
        chain: StructureProperties.chain.auth_asym_id(loc),
        resi: StructureProperties.residue.auth_seq_id(loc),
        x: p[0], y: p[1], z: p[2],
        atom: atomName,
        ring: this.ringCentroidAt(loc, resn, atomName),
      });
    } catch (err) {
      console.warn("click handling failed", err);
    }
  }

  // -- receptor -------------------------------------------------------------

  async setReceptor(receptor: ReceptorArg | null): Promise<void> {
    await this.clearPoses();
    await this.deleteRefs(this.boxRefs); this.boxRefs = [];
    await this.deleteRefs(this.picksRefs); this.picksRefs = [];
    await this.deleteRefs(this.constraintRefs); this.constraintRefs = [];
    if (this.receptorStructure) {
      await this.deleteRefs([this.receptorStructure]);
      this.receptorStructure = null;
      this.receptorComponentRefs = [];
    }
    if (!receptor) return;
    const { structure } = await this.loadStructure(receptor.pdb, "pdb", receptor.name);
    this.receptorStructure = structure;
    await this.styleReceptor();
    this.plugin.canvas3d?.requestCameraReset();
  }

  setHetLigands(hets: HetLigand[]): void {
    this.hetLigands = hets || [];
  }

  async setStyle(style: StyleArg): Promise<void> {
    this.style = style;
    await this.styleReceptor();
  }

  // A sane default look; past that, representations are the user's to
  // manage from Mol*'s own right-hand panel.
  private async styleReceptor(): Promise<void> {
    if (!this.receptorStructure) return;
    await this.deleteRefs(this.receptorComponentRefs);
    this.receptorComponentRefs = [];
    const track = (ref: Ref | null | undefined) => { if (ref) this.receptorComponentRefs.push(ref); };
    const { hide_ligand: hideLigand, hide_protein: hideProtein, active_het: activeHet, prepared } = this.style;
    const builders = this.plugin.builders.structure;
    const rs = this.receptorStructure;

    const polymerComp = await builders.tryCreateComponentFromExpression(rs, Queries.polymer.expression, "polymer");
    track(polymerComp);
    const polymerRep = await this.addRep(polymerComp, "cartoon", "chain-id", undefined, { opacity: 0.92 });
    track(polymerRep);
    if (hideProtein) this.setHidden(polymerRep, true);

    // bound ligands / cofactors - sticks, so they stay clickable
    const ligandComp = await builders.tryCreateComponentFromExpression(rs, ligandsAndIons(), "ligands");
    track(ligandComp);
    const ligandRep = await this.addRep(ligandComp, "ball-and-stick", "uniform", { value: Color(0x4caf50) }, { sizeFactor: 0.25 });
    track(ligandRep);
    if (hideLigand) this.setHidden(ligandRep, true);

    if (activeHet) {
      const activeComp = await builders.tryCreateComponentFromExpression(rs, hetExpr(activeHet), "active-site");
      track(activeComp);
      const activeRep = await this.addRep(activeComp, "ball-and-stick", "uniform", { value: Color(0xf5a623) }, { sizeFactor: 0.28 });
      track(activeRep);
      if (hideLigand) this.setHidden(activeRep, true);

      // After QuickPrep: site residues as sticks with polar hydrogens, so
      // His tautomers and charged/neutral acids/bases read off directly.
      if (prepared) {
        const pocketExpr = MS.struct.modifier.exceptBy({
          0: MS.struct.modifier.includeSurroundings({ 0: hetExpr(activeHet), radius: 5, "as-whole-residues": true }),
          by: hetExpr(activeHet),
        });
        const pocketComp = await builders.tryCreateComponentFromExpression(rs, pocketExpr, "quickprep-pocket");
        track(pocketComp);
        const pocketRep = await this.addRep(pocketComp, "ball-and-stick", "element-symbol", undefined,
          { sizeFactor: 0.2, ignoreHydrogens: true, ignoreHydrogensVariant: "non-polar" });
        track(pocketRep);
        if (hideProtein) this.setHidden(pocketRep, true);
      }
    }
  }

  focus(focus: FocusArg): void {
    if (!this.receptorStructure?.data) return;
    const expr = focus.kind === "het" ? hetExpr(focus.het) : residueExpr(focus.chain, focus.resseq);
    try {
      const loci = StructureElement.Loci.fromExpression(this.receptorStructure.data, expr);
      this.plugin.managers.camera.focusLoci(loci, { durationMs: 300 });
    } catch (err) {
      console.warn("zoom failed", err);
    }
  }

  // -- overlays -------------------------------------------------------------

  async setBox(box: BoxArg | null): Promise<void> {
    await this.deleteRefs(this.boxRefs);
    this.boxRefs = [];
    if (!box || !isValidBox(box.center, box.size)) return;
    const { structure, refs } = await this.loadStructure(boxPdb(box.center, box.size), "pdb", "dockbox");
    const rep = await this.addRep(structure, "line", "uniform", { value: Color(0x4da3ff) }, { lineSizeAttenuation: false });
    this.boxRefs = [...refs, rep].filter(Boolean) as Ref[];
  }

  async setPicks(picks: Vec3Tuple[]): Promise<void> {
    await this.deleteRefs(this.picksRefs);
    this.picksRefs = [];
    if (!picks.length) return;
    const { structure, refs } = await this.loadStructure(picksPdb(picks), "pdb", "picks");
    const rep = await this.addRep(structure, "spacefill", "uniform", { value: Color(0xf0b528) }, { sizeFactor: 0.4 });
    this.picksRefs = [...refs, rep].filter(Boolean) as Ref[];
  }

  async setConstraints(constraints: ConstraintMarker[]): Promise<void> {
    await this.deleteRefs(this.constraintRefs);
    this.constraintRefs = [];
    for (const c of constraints.filter(isPlaced)) {
      const { structure, refs } = await this.loadStructure(constraintPdb(c), "pdb", `constraint-${c.id}`);
      const rep = await this.addRep(structure, "ball-and-stick", "uniform", { value: Color(constraintColor(c.kind)) }, { sizeFactor: 0.25 });
      this.constraintRefs.push(...refs, ...(rep ? [rep] : []));
    }
  }

  // -- poses ----------------------------------------------------------------

  async clearPoses(): Promise<void> {
    await this.deleteRefs(this.interactionRefs);
    this.interactionRefs = [];
    await this.deleteRefs(this.poseStructureRefs);
    this.poseStructureRefs = [];
    this.poseStructures = [];
    this.poseReps = [];
  }

  async setPoses(poses: PosesArg | null): Promise<void> {
    await this.clearPoses();
    if (!poses) return;
    for (const block of poses.blocks) {
      const { structure, refs } = await this.loadStructure(block + "\n$$$$\n", "sdf", poses.name);
      // element-symbol ("CPK") coloring
      const rep = await this.addRep(structure, "ball-and-stick", "element-symbol", undefined, { sizeFactor: 0.16 });
      this.poseStructureRefs.push(...refs);
      this.poseStructures.push(structure);
      this.poseReps.push(rep);
    }
  }

  /** Show only pose `idx` (-1 = all). `light` skips the contacts rebuild. */
  async selectPose(idx: number, light: boolean): Promise<void> {
    this.poseReps.forEach((rep, i) => this.setHidden(rep, !(idx === -1 || i === idx)));
    if (!light) await this.showInteractions(idx);
  }

  // Non-covalent contacts between the given pose and the receptor, computed
  // on a merged copy of (receptor + that pose) since Mol*'s interaction
  // analysis only looks within a single Structure.
  async showInteractions(idx: number): Promise<void> {
    await this.deleteRefs(this.interactionRefs);
    this.interactionRefs = [];
    if (idx < 0 || !this.receptorStructure || !this.poseStructures[idx]) return;
    try {
      const builders = this.plugin.builders.structure;
      // With "hide ligand" on, merge a ligand/ion-free copy of the receptor,
      // so the co-crystallised ligand doesn't come back through the merge.
      let receptorSrc = this.receptorStructure;
      if (this.style.hide_ligand) {
        const strippedExpr = MS.struct.modifier.exceptBy({ 0: MS.struct.generator.all(), by: ligandsAndIons() });
        const stripped = await builders.tryCreateComponentFromExpression(this.receptorStructure, strippedExpr, "receptor-no-ligand");
        if (stripped) { this.interactionRefs.push(stripped); receptorSrc = stripped; }
      }
      const structures = [{ ref: receptorSrc.ref }, { ref: this.poseStructures[idx].ref }];
      const merged = await this.plugin.state.data.build().toRoot()
        .apply(MergeStructures, { structures }, { dependsOn: structures.map((s) => s.ref) })
        .commit();

      // Queries.ligand also catches the receptor's own bound ligand(s) -
      // exclude those so this is specifically the docked pose.
      const originalLigandExprs = this.hetLigands.map(hetExpr);
      const poseExpr = originalLigandExprs.length
        ? MS.struct.modifier.exceptBy({ 0: Queries.ligand.expression, by: MS.struct.combinator.merge(originalLigandExprs) })
        : Queries.ligand.expression;

      const poseComp = await builders.tryCreateComponentFromExpression(merged, poseExpr, "pose-ligand");
      await this.addRep(poseComp, InteractionsRepresentationProvider, InteractionTypeColorThemeProvider,
        undefined, { includeParent: true, parentDisplay: "between" });

      // Nearby residues as sticks, so the contact lines point at something.
      const nearbyExpr = MS.struct.modifier.exceptBy({
        0: MS.struct.modifier.includeSurroundings({ 0: poseExpr, radius: 4.5, "as-whole-residues": true }),
        by: poseExpr,
      });
      const nearbyComp = await builders.tryCreateComponentFromExpression(merged, nearbyExpr, "pose-contact-residues");
      await this.addRep(nearbyComp, "ball-and-stick", "element-symbol", undefined, { sizeFactor: 0.18 });

      this.interactionRefs.push(merged as unknown as Ref);
    } catch (err) {
      console.warn("interactions failed", err);
    }
  }
}
