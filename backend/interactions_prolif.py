"""ProLIF-based interaction detection for one already-docked pose.

Detection only - no rendering (see backend/interaction_render.py). Keeping
the two separate means the rendering step stays reusable independent of how
interactions are found - it only ever consumes the plain `Interaction` list
below.

Unlike PLIP (which needs a single combined receptor+ligand PDB text block,
with the ligand's atom index recovered from PLIP's own read order), ProLIF
works directly on RDKit mols for both sides and matches its own binding-site
residues by distance (`Fingerprint`'s `vicinity_cutoff`) - no complex-PDB
assembly or synthetic chain/resnum bookkeeping needed.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import prolif as plf
from rdkit import Chem

# ProLIF interaction names -> our own coarser kind buckets (the same ones
# backend/interaction_render.py knows how to style/label). PiStacking's two
# sub-flavors (FaceToFace/EdgeToFace) aren't in the default interaction set
# but are included here in case a caller ever asks for "all".
_KIND_BY_PROLIF_NAME = {
    "HBDonor": "hbond",
    "HBAcceptor": "hbond",
    "Hydrophobic": "hydrophobic",
    "PiStacking": "pistacking",
    "FaceToFace": "pistacking",
    "EdgeToFace": "pistacking",
    "PiCation": "pication",
    "CationPi": "pication",
    "Anionic": "saltbridge",
    "Cationic": "saltbridge",
    "XBDonor": "halogen",
    "XBAcceptor": "halogen",
    "MetalDonor": "metal",
    "MetalAcceptor": "metal",
}


@dataclass
class Interaction:
    kind: str  # hbond | hydrophobic | pistacking | pication | saltbridge | halogen | metal
    restype: str
    resnr: int
    reschain: str
    distance: float
    ligand_atom_idx: int | None = None  # None when ProLIF reports a ring/no single anchor atom

    def as_dict(self) -> dict:
        return {
            "kind": self.kind,
            "restype": self.restype,
            "resnr": self.resnr,
            "reschain": self.reschain,
            "distance": self.distance,
            "ligand_atom_idx": self.ligand_atom_idx,
        }


def extract_pose(sdf_path: Path, pose_index: int) -> Chem.Mol:
    """pose_index is 1-based, matching docking.Pose.rank and the pose bar's labels."""
    if pose_index < 1:
        raise ValueError("pose_index is 1-based")
    try:
        supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=False)
    except OSError as exc:
        raise ValueError(f"could not read poses file: {exc}") from exc
    for i, mol in enumerate(supplier, start=1):
        if i != pose_index:
            continue
        if mol is None:
            raise ValueError(f"pose {pose_index} could not be parsed")
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            pass  # ProLIF re-derives its own perception from the 3D coords anyway
        return mol
    raise ValueError(f"pose_index {pose_index} out of range")


def _protein_molecule(receptor_pdb_text: str) -> "plf.Molecule":
    mol = Chem.MolFromPDBBlock(receptor_pdb_text, sanitize=False, removeHs=False)
    if mol is None:
        raise ValueError("receptor PDB text could not be parsed")
    Chem.SanitizeMol(mol, catchErrors=True)
    return plf.Molecule.from_rdkit(mol)


def _ligand_molecule(pose_mol: Chem.Mol) -> "plf.Molecule":
    mol = Chem.Mol(pose_mol)
    try:
        Chem.AssignStereochemistryFrom3D(mol)  # ProLIF/RDKit's aromaticity+H-bond
        # patterns key off stereo/hybridization perceived from the 3D pose, not
        # just 2D graph connectivity from the docking SDF.
    except Exception:
        pass
    return plf.Molecule.from_rdkit(mol)


def detect_interactions(receptor_pdb_text: str, pose_mol: Chem.Mol) -> list[Interaction]:
    protein = _protein_molecule(receptor_pdb_text)
    ligand = _ligand_molecule(pose_mol)

    fp = plf.Fingerprint()
    fp.run_from_iterable([ligand], protein, progress=False)
    frame_ifp = fp.ifp.get(0, {})

    hits: list[Interaction] = []
    for (_lig_res, prot_res), by_name in frame_ifp.items():
        for name, occurrences in by_name.items():
            kind = _KIND_BY_PROLIF_NAME.get(name)
            if kind is None:
                continue
            for occ in occurrences:
                lig_indices = occ["indices"]["ligand"]
                ligand_atom_idx = lig_indices[0] if lig_indices else None
                hits.append(
                    Interaction(
                        kind,
                        prot_res.name,
                        prot_res.number,
                        prot_res.chain or "",
                        occ["distance"],
                        ligand_atom_idx,
                    )
                )
    return hits


def analyze_pose(receptor_pdb_text: str, sdf_path: Path, pose_index: int) -> tuple[Chem.Mol, list[Interaction]]:
    """Top-level entrypoint: extract one pose and run ProLIF against the
    receptor. An empty interaction list is a valid result (the pose just
    doesn't contact anything notable), not an error."""
    mol = extract_pose(sdf_path, pose_index)
    hits = detect_interactions(receptor_pdb_text, mol)
    return mol, hits
