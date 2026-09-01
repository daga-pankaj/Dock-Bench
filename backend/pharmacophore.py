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

KINDS = ("hbond_donor", "hbond_acceptor", "positional")

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

    def __post_init__(self) -> None:
        if self.kind not in KINDS:
            raise ValueError(f"unknown constraint kind: {self.kind}")

    @classmethod
    def from_dict(cls, d: dict) -> "Constraint":
        return cls(
            kind=d["kind"],
            point=tuple(d["point"]),
            radius=float(d.get("radius", 3.5)),
            required=bool(d.get("required", True)),
            weight=float(d.get("weight", 1.0)),
            label=d.get("label", ""),
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
    if kind == "positional":
        return {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1}
    patt = _DONOR_SMARTS if kind == "hbond_donor" else _ACCEPTOR_SMARTS
    return {m[0] for m in mol.GetSubstructMatches(patt)}


def evaluate(mol: Chem.Mol, constraints: list[Constraint]) -> list[ConstraintResult]:
    """Check one already-posed RDKit mol (needs a 3D conformer) against every constraint."""
    conf = mol.GetConformer()
    results = []
    for c in constraints:
        best = None
        for idx in _candidate_atoms(mol, c.kind):
            p = conf.GetAtomPosition(idx)
            d = ((p.x - c.point[0]) ** 2 + (p.y - c.point[1]) ** 2 + (p.z - c.point[2]) ** 2) ** 0.5
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
