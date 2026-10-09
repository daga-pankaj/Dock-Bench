import { describe, expect, it } from "vitest";

import { stepIndex } from "./poseBar";
import {
  boxPdb, boxPoints, constraintColor, constraintPdb, constraintPoints, edgePoints, fingerprint,
  isPlaced, isValidBox, pdbAtomLine, picksPdb, ringFor, ringGeometry,
} from "./synthetic";

describe("pdbAtomLine", () => {
  it("writes fixed PDB columns", () => {
    const line = pdbAtomLine(1, "C", "BOX", "Z", 1, 1.5, -2.25, 10, "C");
    expect(line.slice(0, 6)).toBe("HETATM");
    expect(line.slice(6, 11)).toBe("    1");
    expect(line.slice(17, 20)).toBe("BOX");
    expect(line[21]).toBe("Z");
    expect(parseFloat(line.slice(30, 38))).toBe(1.5);
    expect(parseFloat(line.slice(38, 46))).toBe(-2.25);
    expect(parseFloat(line.slice(46, 54))).toBe(10);
    expect(line.slice(76, 78)).toBe(" C");
  });

  it("defaults an empty chain to Z", () => {
    expect(pdbAtomLine(1, "C", "PIK", "", 1, 0, 0, 0, "C")[21]).toBe("Z");
  });
});

describe("edgePoints", () => {
  it("includes both endpoints and spaces points ~1.2 Å", () => {
    const pts = edgePoints([0, 0, 0], [12, 0, 0]);
    expect(pts).toHaveLength(11);
    expect(pts[0]).toEqual([0, 0, 0]);
    expect(pts[10]).toEqual([12, 0, 0]);
  });

  it("returns two points for a zero-length segment", () => {
    expect(edgePoints([1, 1, 1], [1, 1, 1])).toHaveLength(2);
  });
});

describe("box", () => {
  it("validates center and size", () => {
    expect(isValidBox([0, 0, 0], [10, 10, 10])).toBe(true);
    expect(isValidBox([0, NaN, 0], [10, 10, 10])).toBe(false);
    expect(isValidBox([0, 0, 0], [10, 0, 10])).toBe(false);
  });

  it("samples all 12 edges within the box bounds", () => {
    const pts = boxPoints([1, 2, 3], [12, 12, 12]);
    // 12 edges x 11 points each
    expect(pts).toHaveLength(132);
    for (const p of pts) {
      expect(p[0]).toBeGreaterThanOrEqual(-5 - 1e-9);
      expect(p[0]).toBeLessThanOrEqual(7 + 1e-9);
    }
    const corners = pts.filter((p) => Math.abs(Math.abs(p[0] - 1) - 6) < 1e-9
      && Math.abs(Math.abs(p[1] - 2) - 6) < 1e-9 && Math.abs(Math.abs(p[2] - 3) - 6) < 1e-9);
    expect(new Set(corners.map((c) => c.join())).size).toBe(8);
  });

  it("box and picks PDB end with END", () => {
    expect(boxPdb([0, 0, 0], [5, 5, 5]).endsWith("\nEND\n")).toBe(true);
    expect(picksPdb([[1, 2, 3]]).split("\n")[0]).toContain("PIK");
  });
});

describe("constraint markers", () => {
  it("needs a finite point to be placed", () => {
    expect(isPlaced({ id: 1, kind: "positional", x: 1, y: 2, z: 3 })).toBe(true);
    expect(isPlaced({ id: 1, kind: "positional", x: 1, y: null, z: 3 })).toBe(false);
  });

  it("is a single point without direction/normal", () => {
    expect(constraintPoints({ id: 1, kind: "positional", x: 1, y: 2, z: 3 })).toEqual([[1, 2, 3]]);
    expect(constraintPoints({ id: 1, kind: "hbond_donor", x: 1, y: 2, z: 3 })).toHaveLength(1);
  });

  it("adds a 17-point disk of radius 1.5 in the ring plane for aromatic", () => {
    const pts = constraintPoints({ id: 1, kind: "aromatic", x: 0, y: 0, z: 0, nx: 0, ny: 0, nz: 1 });
    expect(pts).toHaveLength(18);
    for (const p of pts.slice(1)) {
      expect(Math.hypot(p[0], p[1], p[2])).toBeCloseTo(1.5);
      expect(p[2]).toBeCloseTo(0);
    }
  });

  it("adds a shaft and arrowhead reaching 2 Å along the H-bond direction", () => {
    const pts = constraintPoints({ id: 1, kind: "hbond_acceptor", x: 0, y: 0, z: 0, dx: 1, dy: 0, dz: 0 });
    expect(pts.length).toBeGreaterThan(6);
    const maxX = Math.max(...pts.map((p) => p[0]));
    expect(maxX).toBeCloseTo(2);
  });

  it("colors by kind with a fallback", () => {
    expect(constraintColor("aromatic")).toBe(0x9b59d6);
    expect(constraintColor("unknown")).toBe(0xff6bd6);
    expect(constraintPdb({ id: 1, kind: "positional", x: 0, y: 0, z: 0 })).toContain("PHC");
  });
});

describe("aromatic rings", () => {
  it("finds the ring containing the clicked atom", () => {
    expect(ringFor("PHE", "CZ")).toHaveLength(6);
    expect(ringFor("TRP", "NE1")).toHaveLength(5);
    expect(ringFor("TRP", "CZ2")).toHaveLength(6);
    expect(ringFor("ALA", "CB")).toBeNull();
    expect(ringFor("PHE", "CA")).toBeNull();
  });

  it("returns centroid and unit normal, or null for an incomplete ring", () => {
    const hexagon = Array.from({ length: 6 }, (_, i) =>
      [Math.cos((i * Math.PI) / 3), Math.sin((i * Math.PI) / 3), 2] as [number, number, number]);
    const g = ringGeometry(hexagon, 6)!;
    expect(g.centroid[0]).toBeCloseTo(0);
    expect(g.centroid[2]).toBeCloseTo(2);
    expect(Math.abs(g.normal![2])).toBeCloseTo(1);
    expect(ringGeometry(hexagon.slice(0, 4), 6)).toBeNull();
  });
});

describe("pose stepping", () => {
  it("wraps and enters from 'all'", () => {
    expect(stepIndex(0, 1, 3)).toBe(1);
    expect(stepIndex(2, 1, 3)).toBe(0);
    expect(stepIndex(0, -1, 3)).toBe(2);
    expect(stepIndex(-1, 1, 3)).toBe(0);
    expect(stepIndex(-1, -1, 3)).toBe(2);
  });
});

describe("fingerprint", () => {
  it("is stable for equal values and treats undefined as null", () => {
    expect(fingerprint({ a: [1, 2] })).toBe(fingerprint({ a: [1, 2] }));
    expect(fingerprint(undefined)).toBe(fingerprint(null));
  });
});
