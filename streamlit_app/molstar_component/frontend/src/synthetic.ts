/* Synthetic overlay geometry: the docking box, pick markers and
 * pharmacophore-constraint markers are drawn as tiny PDB structures of dummy
 * carbon atoms, so Mol*'s normal distance-based bonding renders their edges.
 * Everything here is pure (no Mol* imports) so it can be unit tested.
 */

import type { ConstraintMarker, Vec3Tuple } from "./types";

// Residue names used by the synthetic atoms. Clicks landing on our own
// overlay geometry must be ignored.
export const BOX_RESNAME = "BOX";
export const PICK_RESNAME = "PIK";
export const CONSTRAINT_RESNAME = "PHC";
export const SYNTHETIC_RESNAMES = new Set([BOX_RESNAME, PICK_RESNAME, CONSTRAINT_RESNAME]);

// Color per constraint kind (mirrors the original app).
export const CONSTRAINT_COLORS: Record<string, number> = {
  hbond_donor: 0x4a90e2,
  hbond_acceptor: 0xe24a4a,
  positional: 0xf5a623,
  aromatic: 0x9b59d6,
  metal_coordination: 0x27ae60,
  smarts: 0x1abc9c,
};
export const DEFAULT_CONSTRAINT_COLOR = 0xff6bd6;

export function pdbAtomLine(
  serial: number, name: string, resName: string, chain: string, resSeq: number,
  x: number, y: number, z: number, element: string,
): string {
  const rj = (v: string | number, n: number) => String(v).padStart(n);
  const lj = (v: string | number, n: number) => String(v).padEnd(n);
  return (
    "HETATM" + rj(serial, 5) + " " + lj(name, 4) + " " + lj(resName, 3) + " " +
    (chain || "Z") + rj(resSeq, 4) + "    " +
    rj(x.toFixed(3), 8) + rj(y.toFixed(3), 8) + rj(z.toFixed(3), 8) +
    rj("1.00", 6) + rj("0.00", 6) + "          " + rj(element, 2)
  );
}

// Points spaced ~`spacing` Å apart along a segment so Mol*'s distance-based
// bonding chains them into a visible line. Always includes both endpoints.
export function edgePoints(a: Vec3Tuple, b: Vec3Tuple, spacing = 1.2): Vec3Tuple[] {
  const d = Math.hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
  const n = Math.max(1, Math.round(d / spacing));
  const pts: Vec3Tuple[] = [];
  for (let i = 0; i <= n; i++) {
    const t = i / n;
    pts.push([a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t]);
  }
  return pts;
}

const isFiniteVec = (v: (number | null | undefined)[]): v is number[] =>
  v.every((x) => typeof x === "number" && Number.isFinite(x));

export function isValidBox(center: Vec3Tuple, size: Vec3Tuple): boolean {
  return isFiniteVec(center) && isFiniteVec(size) && size.every((v) => v > 0);
}

// Wireframe of an axis-aligned box: the 12 edges, each sampled by edgePoints.
export function boxPoints(center: Vec3Tuple, size: Vec3Tuple): Vec3Tuple[] {
  const half = size.map((v) => v / 2);
  const corner = (bit: number): Vec3Tuple =>
    [0, 1, 2].map((i) => center[i] + (bit & (1 << i) ? half[i] : -half[i])) as Vec3Tuple;
  const pts: Vec3Tuple[] = [];
  for (let i = 0; i < 8; i++) {
    for (let bit = 0; bit < 3; bit++) {
      const j = i ^ (1 << bit);
      if (j <= i) continue;
      pts.push(...edgePoints(corner(i), corner(j)));
    }
  }
  return pts;
}

export function boxPdb(center: Vec3Tuple, size: Vec3Tuple): string {
  const lines = boxPoints(center, size).map((p, i) =>
    pdbAtomLine(i + 1, "C", BOX_RESNAME, "Z", i + 1, p[0], p[1], p[2], "C"));
  return lines.join("\n") + "\nEND\n";
}

export function picksPdb(picks: Vec3Tuple[]): string {
  const lines = picks.map((p, i) => pdbAtomLine(i + 1, "C", PICK_RESNAME, "Z", i + 1, p[0], p[1], p[2], "C"));
  return lines.join("\n") + "\nEND\n";
}

const cross = (a: number[], b: number[]): Vec3Tuple =>
  [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const normalize = (v: number[]): Vec3Tuple => {
  const m = Math.hypot(v[0], v[1], v[2]);
  return [v[0] / m, v[1] / m, v[2] / m];
};

export function isPlaced(c: ConstraintMarker): boolean {
  return isFiniteVec([c.x, c.y, c.z]);
}

// Marker points for one constraint: the center point, plus a disk in the
// ring plane for an aromatic constraint with a normal, or a direction arrow
// for an H-bond constraint with a direction.
export function constraintPoints(c: ConstraintMarker): Vec3Tuple[] {
  const x = c.x as number, y = c.y as number, z = c.z as number;
  const pts: Vec3Tuple[] = [[x, y, z]];

  if (c.kind === "aromatic" && isFiniteVec([c.nx, c.ny, c.nz])) {
    const diskRadius = 1.5;
    const nSteps = 16;
    const n = [c.nx as number, c.ny as number, c.nz as number];
    const up = Math.abs(n[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0];
    const u1 = normalize(cross(up, n));
    const v = cross(n, u1);
    for (let i = 0; i <= nSteps; i++) {
      const theta = (2 * Math.PI * i) / nSteps;
      pts.push([0, 1, 2].map((k) =>
        [x, y, z][k] + diskRadius * (Math.cos(theta) * u1[k] + Math.sin(theta) * v[k])) as Vec3Tuple);
    }
  } else if ((c.kind === "hbond_donor" || c.kind === "hbond_acceptor") && isFiniteVec([c.dx, c.dy, c.dz])) {
    const d = [c.dx as number, c.dy as number, c.dz as number];
    const arrowLen = 2.0;
    const tipLen = 0.5;
    // Shaft: t = 0, 0.2, ..., 1.0 (same float stepping as the original).
    for (let t = 0; t <= 1; t += 0.2) {
      pts.push([x + t * arrowLen * d[0], y + t * arrowLen * d[1], z + t * arrowLen * d[2]]);
    }
    const tip: Vec3Tuple = [x + arrowLen * d[0], y + arrowLen * d[1], z + arrowLen * d[2]];
    const base = [x + (arrowLen - tipLen) * d[0], y + (arrowLen - tipLen) * d[1], z + (arrowLen - tipLen) * d[2]];
    const perp = Math.abs(d[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0];
    const p1n = normalize(cross(perp, d));
    for (let i = 0; i < 4; i++) {
      const cosT = Math.cos((Math.PI * i) / 2);
      const start: Vec3Tuple = [base[0] + 0.3 * cosT * p1n[0], base[1] + 0.3 * cosT * p1n[1], base[2] + 0.3 * cosT * p1n[2]];
      pts.push(...edgePoints(start, tip, 0.5));
    }
  }
  return pts;
}

export function constraintPdb(c: ConstraintMarker): string {
  const lines = constraintPoints(c).map((p, i) =>
    // first atom keeps resSeq 1, the rest use serial+1 (as the original did)
    pdbAtomLine(i + 1, "C", CONSTRAINT_RESNAME, "Z", i === 0 ? 1 : i + 2, p[0], p[1], p[2], "C"));
  return lines.join("\n") + "\nEND\n";
}

export function constraintColor(kind: string): number {
  return CONSTRAINT_COLORS[kind] ?? DEFAULT_CONSTRAINT_COLOR;
}

// Side-chain aromatic ring atom names per residue - used to snap a click to
// the ring's centroid. TRP has two rings; whichever contains the atom wins.
export const AROMATIC_RINGS: Record<string, string[][]> = {
  PHE: [["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]],
  TYR: [["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]],
  HIS: [["CG", "ND1", "CD2", "CE1", "NE2"]],
  TRP: [
    ["CG", "CD1", "NE1", "CE2", "CD2"],
    ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
  ],
};

export function ringFor(resn: string, atomName: string): string[] | null {
  return AROMATIC_RINGS[resn]?.find((names) => names.includes(atomName)) ?? null;
}

// Centroid + plane normal (cross product of the first two edges) of ring
// points; null when the ring came back incomplete.
export function ringGeometry(pts: Vec3Tuple[], ringSize: number):
  { centroid: Vec3Tuple; normal: Vec3Tuple | null } | null {
  if (pts.length < ringSize) return null;
  const n = pts.length;
  const sum = pts.reduce((a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]], [0, 0, 0]);
  const centroid: Vec3Tuple = [sum[0] / n, sum[1] / n, sum[2] / n];
  if (pts.length < 3) return { centroid, normal: null };
  const v1 = [pts[1][0] - pts[0][0], pts[1][1] - pts[0][1], pts[1][2] - pts[0][2]];
  const v2 = [pts[2][0] - pts[0][0], pts[2][1] - pts[0][1], pts[2][2] - pts[0][2]];
  const c = cross(v1, v2);
  const mag = Math.hypot(c[0], c[1], c[2]);
  return { centroid, normal: mag > 0 ? [c[0] / mag, c[1] / mag, c[2] / mag] : null };
}

// Stable string fingerprint of any JSON-able value, for arg diffing.
export function fingerprint(value: unknown): string {
  return JSON.stringify(value ?? null);
}
