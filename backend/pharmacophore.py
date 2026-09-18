"""Pharmacophoric constraint checking: post-dock pose filter/rerank.

smina/gnina have no concept of "this ligand atom must land near that point"
during their own search - there's no anchor-point bias to hook into. The only
way to make a pharmacophore requirement (an H-bond to a specific residue, a
group forced into a specific sub-pocket) actually steer the results with this
stack is to check each already-docked pose's geometry and rerank/filter on
it: poses meeting every required constraint sort first, and the rest are
penalised in proportion to how far they miss by.
"""

from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem

KINDS = (
    "hbond_donor",
    "hbond_acceptor",
    "positional",
    "aromatic",          # ring-centroid proximity, as a pi-stacking proxy
    "metal_coordination",  # Lewis-base atom near a metal ion position
    "smarts",             # any atom matching a user-supplied SMARTS pattern
)

# RDKit's own Lipinski donor/acceptor atom-matching SMARTS (the same patterns
# rdkit.Chem.Lipinski uses for NumHDonors/NumHAcceptors) - used here to find
# *which* ligand atom a pose would have to place near the constraint point.
_DONOR_SMARTS = Chem.MolFromSmarts(
    "[$([N;!H0;v3]),$([N;!H0;+1;v4]),$([O,S;H1;+0]),$([n;H1;+0])]"
)
_ACCEPTOR_SMARTS = Chem.MolFromSmarts(
    "[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),"
    "$([N;v3;!$(N-*=!@[O,N,P,S])]),$([nH0,o,s;+0])]"
)


@dataclass
class Constraint:
    kind: str
    point: tuple[float, float, float]
    radius: float = 3.5
    required: bool = True
    weight: float = 1.0
    label: str = ""
    smarts: str | None = None  # only used (and required) when kind == "smarts"

    def __post_init__(self) -> None:
        if self.kind not in KINDS:
            raise ValueError(f"unknown constraint kind: {self.kind}")
        if self.kind == "smarts":
            if not self.smarts:
                raise ValueError("a 'smarts' constraint needs a `smarts` pattern")
            if Chem.MolFromSmarts(self.smarts) is None:
                raise ValueError(f"invalid SMARTS pattern: {self.smarts!r}")

    @classmethod
    def from_dict(cls, d: dict) -> "Constraint":
        return cls(
            kind=d["kind"],
            point=tuple(d["point"]),
            radius=float(d.get("radius", 3.5)),
            required=bool(d.get("required", True)),
            weight=float(d.get("weight", 1.0)),
            label=d.get("label", ""),
            smarts=d.get("smarts"),
        )


@dataclass
class ConstraintResult:
    label: str
    kind: str
    required: bool
    distance: float | None  # None => the ligand has no atom of the needed kind at all
    satisfied: bool

    def as_dict(self) -> dict:
        return {
            "label": self.label,
            "kind": self.kind,
            "required": self.required,
            "distance": self.distance,
            "satisfied": self.satisfied,
        }


def _candidate_atoms(mol: Chem.Mol, kind: str) -> set[int]:
    """Ligand atom indices to check against a constraint's point, for the
    atom-index-based kinds. "aromatic" (ring-centroid) and "smarts"
    (per-constraint pattern) go through `_candidate_points` instead."""
    if kind == "positional":
        return {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1}
    if kind == "hbond_donor":
        patt = _DONOR_SMARTS
    elif kind in ("hbond_acceptor", "metal_coordination"):
        # No separate RDKit pattern exists for "chelates a metal" - the atoms
        # that do (carbonyl/hydroxyl/carboxylate O, amine/aromatic N,
        # thiol/thioether S) are the same lone-pair-bearing set that accepts
        # H-bonds, so it's reused here as a chemically-reasonable proxy.
        patt = _ACCEPTOR_SMARTS
    else:
        return set()
    return {m[0] for m in mol.GetSubstructMatches(patt)}


def _ring_centroids(mol: Chem.Mol, conf) -> list[tuple[float, float, float]]:
    """Centroid of each aromatic 5- or 6-membered ring - pi-stacking is a
    ring-to-ring (or ring-to-point) alignment, not a single-atom contact, so
    the "aromatic" kind checks proximity of ring centers rather than atoms."""
    centroids = []
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) not in (5, 6):
            continue
        if not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            continue
        pts = [conf.GetAtomPosition(i) for i in ring]
        n = len(pts)
        centroids.append((sum(p.x for p in pts) / n, sum(p.y for p in pts) / n, sum(p.z for p in pts) / n))
    return centroids


def _candidate_points(mol: Chem.Mol, conf, c: Constraint) -> list[tuple[float, float, float]]:
    """3D points on the ligand to check against constraint `c`'s point."""
    if c.kind == "aromatic":
        return _ring_centroids(mol, conf)
    if c.kind == "smarts":
        patt = Chem.MolFromSmarts(c.smarts) if c.smarts else None
        if patt is None:
            return []
        idxs = {m[0] for m in mol.GetSubstructMatches(patt)}
    else:
        idxs = _candidate_atoms(mol, c.kind)
    points = []
    for i in idxs:
        p = conf.GetAtomPosition(i)
        points.append((p.x, p.y, p.z))
    return points


def evaluate(mol: Chem.Mol, constraints: list[Constraint]) -> list[ConstraintResult]:
    """Check one already-posed RDKit mol (needs a 3D conformer) against every constraint."""
    conf = mol.GetConformer()
    results = []
    for c in constraints:
        best = None
        for p in _candidate_points(mol, conf, c):
            d = ((p[0] - c.point[0]) ** 2 + (p[1] - c.point[1]) ** 2 + (p[2] - c.point[2]) ** 2) ** 0.5
            if best is None or d < best:
                best = d
        results.append(
            ConstraintResult(
                label=c.label or c.kind,
                kind=c.kind,
                required=c.required,
                distance=best,
                satisfied=best is not None and best <= c.radius,
            )
        )
    return results


def penalty(results: list[ConstraintResult], constraints: list[Constraint]) -> float:
    """Soft score added on top of affinity for reranking - 0 when everything is met."""
    total = 0.0
    for r, c in zip(results, constraints):
        if r.distance is None:
            total += c.weight * c.radius * 2  # no matching atom on the ligand at all
        elif not r.satisfied:
            total += c.weight * (r.distance - c.radius)
    return total


def all_required_met(results: list[ConstraintResult]) -> bool:
    return all(r.satisfied for r in results if r.required)
