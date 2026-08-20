"""Structure preparation: receptor cleanup, site detection, ligand 3D building.

smina links Open Babel, so it reads a plain PDB receptor and SDF ligands and
writes SDF poses. That means no PDBQT conversion step is needed anywhere here.
"""

from __future__ import annotations

import io
import re
from dataclasses import dataclass
from typing import Iterable

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

WATER = {"HOH", "DOD", "WAT", "H2O", "SOL"}

# Ions and common crystallisation additives are not druglike "site markers".
IONS = {
    "NA", "K", "MG", "CA", "ZN", "FE", "FE2", "MN", "CU", "CU1", "CO", "NI",
    "CD", "HG", "CL", "BR", "IOD", "F", "SO4", "PO4", "NO3", "ACT", "EDO",
    "GOL", "PEG", "PGE", "MPD", "DMS", "TRS", "IMD", "FMT", "ACY", "CIT",
    "EPE", "MES", "BME", "IPA", "SCN", "AZI", "NH4",
}

# Metals worth keeping in the receptor file even though they are HETATM.
KEEP_METALS = {"ZN", "MG", "MN", "FE", "FE2", "CA", "CU", "CU1", "NI", "CO", "K", "NA"}


@dataclass
class HetLigand:
    """A non-water HETATM residue that can seed the docking box."""

    key: str
    resname: str
    chain: str
    resseq: str
    n_atoms: int
    center: list[float]
    size: list[float]
    is_ion: bool

    def as_dict(self) -> dict:
        return {
            "key": self.key,
            "resname": self.resname,
            "chain": self.chain,
            "resseq": self.resseq,
            "n_atoms": self.n_atoms,
            "center": [round(v, 3) for v in self.center],
            "size": [round(v, 3) for v in self.size],
            "is_ion": self.is_ion,
        }


def _atom_records(pdb_text: str) -> Iterable[tuple[str, str]]:
    for line in pdb_text.splitlines():
        tag = line[:6].strip()
        if tag in ("ATOM", "HETATM"):
            yield tag, line


def _xyz(line: str) -> tuple[float, float, float]:
    return (float(line[30:38]), float(line[38:46]), float(line[46:54]))


def _resname(line: str) -> str:
    return line[17:20].strip().upper()


def list_het_ligands(pdb_text: str, min_atoms: int = 4) -> list[HetLigand]:
    """Group HETATM residues (waters excluded) so the UI can offer them as sites."""
    groups: dict[str, list[str]] = {}
    for tag, line in _atom_records(pdb_text):
        if tag != "HETATM":
            continue
        resname = _resname(line)
        if resname in WATER:
            continue
        chain = line[21].strip() or "_"
        resseq = line[22:27].strip()
        groups.setdefault(f"{resname}:{chain}:{resseq}", []).append(line)

    out: list[HetLigand] = []
    for key, lines in groups.items():
        resname, chain, resseq = key.split(":")
        is_ion = resname in IONS
        if len(lines) < min_atoms and not is_ion:
            continue
        coords = [_xyz(ln) for ln in lines]
        mins = [min(c[i] for c in coords) for i in range(3)]
        maxs = [max(c[i] for c in coords) for i in range(3)]
        center = [(mins[i] + maxs[i]) / 2 for i in range(3)]
        size = [maxs[i] - mins[i] for i in range(3)]
        out.append(HetLigand(key, resname, chain, resseq, len(lines), center, size, is_ion))
    # Biggest residues first: the real ligand usually tops the list.
    out.sort(key=lambda h: (h.is_ion, -h.n_atoms))
    return out


def clean_receptor(
    pdb_text: str,
    keep_waters: bool = False,
    keep_metals: bool = True,
    keep_het: Iterable[str] = (),
) -> str:
    """Strip everything smina should not treat as the rigid receptor.

    Keeps all ATOM records, drops alternate locations other than A, and keeps
    only the HETATM groups explicitly requested (plus metals, by default).
    """
    keep_het = {k.upper() for k in keep_het}
    kept: list[str] = []
    for line in pdb_text.splitlines():
        tag = line[:6].strip()
        if tag in ("TER", "END"):
            continue
        if tag not in ("ATOM", "HETATM"):
            continue
        altloc = line[16]
        if altloc not in (" ", "A"):
            continue
        resname = _resname(line)
        if tag == "HETATM":
            if resname in WATER:
                if not keep_waters:
                    continue
            elif resname in keep_het:
                pass
            elif keep_metals and resname in KEEP_METALS:
                pass
            else:
                continue
        # Hydrogens from the crystal structure are unreliable; smina adds its own.
        element = line[76:78].strip().upper()
        if element == "H" or (not element and line[12:16].strip().startswith("H")):
            continue
        kept.append(line[:66].ljust(66))
    kept.append("END")
    return "\n".join(kept) + "\n"


# --------------------------------------------------------------------------
# Ligand input
# --------------------------------------------------------------------------


@dataclass
class LigandSpec:
    name: str
    smiles: str | None = None
    molblock: str | None = None


def _dedupe_names(specs: list[LigandSpec]) -> list[LigandSpec]:
    seen: dict[str, int] = {}
    for spec in specs:
        base = re.sub(r"[^A-Za-z0-9._+-]", "_", spec.name.strip()) or "ligand"
        n = seen.get(base, 0)
        seen[base] = n + 1
        spec.name = base if n == 0 else f"{base}_{n + 1}"
    return specs


def parse_smiles_block(text: str, prefix: str = "lig") -> list[LigandSpec]:
    """Parse one SMILES per line, with an optional whitespace-separated name."""
    specs: list[LigandSpec] = []
    for i, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(None, 1)
        smiles = parts[0]
        name = parts[1].strip() if len(parts) > 1 else f"{prefix}_{i}"
        specs.append(LigandSpec(name=name, smiles=smiles))
    return _dedupe_names(specs)


def parse_sdf_block(text: str) -> list[LigandSpec]:
    """Read an SDF/MOL file; 3D coordinates already present are reused as-is."""
    specs: list[LigandSpec] = []
    supplier = Chem.ForwardSDMolSupplier(io.BytesIO(text.encode()), removeHs=False)
    for i, mol in enumerate(supplier, start=1):
        if mol is None:
            continue
        name = mol.GetProp("_Name").strip() if mol.HasProp("_Name") else ""
        specs.append(
            LigandSpec(
                name=name or f"sdf_{i}",
                smiles=Chem.MolToSmiles(Chem.RemoveHs(mol)) if mol.GetNumAtoms() else None,
                molblock=Chem.MolToMolBlock(mol) if mol.GetNumConformers() else None,
            )
        )
    return _dedupe_names(specs)


def parse_ligand_file(filename: str, text: str) -> list[LigandSpec]:
    """Dispatch on extension: .sdf/.mol keep coordinates, everything else is SMILES."""
    lower = filename.lower()
    if lower.endswith((".sdf", ".mol", ".sd")):
        return parse_sdf_block(text)
    if lower.endswith(".csv") or lower.endswith(".tsv"):
        sep = "\t" if lower.endswith(".tsv") else ","
        specs: list[LigandSpec] = []
        lines = [ln for ln in text.splitlines() if ln.strip()]
        if not lines:
            return []
        header = [c.strip().lower() for c in lines[0].split(sep)]
        smi_col = next((i for i, c in enumerate(header) if c in ("smiles", "smi")), None)
        name_col = next((i for i, c in enumerate(header) if c in ("name", "id", "title")), None)
        body = lines[1:] if smi_col is not None else lines
        if smi_col is None:
            smi_col = 0
        for i, line in enumerate(body, start=1):
            cells = line.split(sep)
            if smi_col >= len(cells):
                continue
            name = cells[name_col].strip() if name_col is not None and name_col < len(cells) else f"row_{i}"
            specs.append(LigandSpec(name=name, smiles=cells[smi_col].strip()))
        return _dedupe_names(specs)
    return parse_smiles_block(text, prefix="file")


def build_3d(spec: LigandSpec, seed: int = 0xF00D) -> str:
    """Return an SDF molblock with 3D coordinates, ready for smina.

    Raises ValueError with a human-readable reason when the input cannot be
    embedded, so the job can mark that one ligand failed and carry on.
    """
    if spec.molblock:
        mol = Chem.MolFromMolBlock(spec.molblock, removeHs=False)
        if mol is not None and mol.GetNumConformers():
            mol.SetProp("_Name", spec.name)
            return Chem.MolToMolBlock(mol)

    if not spec.smiles:
        raise ValueError("no SMILES and no usable 3D block")
    mol = Chem.MolFromSmiles(spec.smiles)
    if mol is None:
        raise ValueError(f"unparseable SMILES: {spec.smiles}")

    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.useSmallRingTorsions = True
    if AllChem.EmbedMolecule(mol, params) != 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            raise ValueError("3D embedding failed")
    try:
        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        else:
            AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    except Exception:  # optimisation is a nicety, the embedding is what matters
        pass
    mol.SetProp("_Name", spec.name)
    return Chem.MolToMolBlock(mol)
