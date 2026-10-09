/* Args the Python wrapper sends (see ../../__init__.py) and events sent back. */

export type Vec3Tuple = [number, number, number];

export interface ReceptorArg {
  id: string;
  name: string;
  pdb: string;
}

export interface HetRef {
  resname: string;
  chain: string;
  resseq: string | number;
}

export interface HetLigand extends HetRef {
  key: string;
  n_atoms: number;
  center: Vec3Tuple;
  size: Vec3Tuple;
  is_ion: boolean;
}

export interface StyleArg {
  hide_ligand: boolean;
  hide_protein: boolean;
  active_het: HetRef | null;
  prepared: boolean;
}

export interface BoxArg {
  center: Vec3Tuple;
  size: Vec3Tuple;
}

export interface ConstraintMarker {
  id: number | string;
  kind: string;
  x?: number | null;
  y?: number | null;
  z?: number | null;
  nx?: number | null;
  ny?: number | null;
  nz?: number | null;
  dx?: number | null;
  dy?: number | null;
  dz?: number | null;
}

export interface PosesArg {
  key: string;
  name: string;
  blocks: string[];
  scores: (number | null)[];
}

export interface DepictionArg {
  name: string;
  svg: string;
}

export type FocusArg =
  | { nonce: number; kind: "het"; het: HetRef }
  | { nonce: number; kind: "residue"; chain: string; resseq: number };

export interface ComponentArgs {
  receptor: ReceptorArg | null;
  het_ligands: HetLigand[];
  style: StyleArg | null;
  box: BoxArg | null;
  picks: Vec3Tuple[];
  constraints: ConstraintMarker[];
  poses: PosesArg | null;
  depiction: DepictionArg | null;
  focus: FocusArg | null;
  height: number;
}

export interface ClickedAtom {
  resn: string;
  hetflag: boolean;
  chain: string;
  resi: number;
  atom: string;
  x: number;
  y: number;
  z: number;
  ring: { centroid: Vec3Tuple; normal: Vec3Tuple | null } | null;
}

export type ViewerEventBody =
  | { type: "click"; atom: ClickedAtom }
  | { type: "pose"; index: number };

export type ViewerEvent = ViewerEventBody & { nonce: number };
