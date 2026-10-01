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


def _element(line: str) -> str:
    """Element symbol for an ATOM/HETATM line: the dedicated column when
    present, else a guess from the atom-name field (good enough for the
    common single-letter elements ligands are mostly made of)."""
    e = line[76:78].strip().upper()
    if e:
        return e
    return re.sub(r"[^A-Za-z]", "", line[12:16])[:1].upper()


def _atom_name(line: str) -> str:
    """PDB atom-name field (columns 13-16), stripped/uppercased - e.g. "SG"
    for a cysteine sulfur. Needed for chain:resnum:atom_name identifiers
    (gnina's --covalent_rec_atom, and backend.covalent.receptor_atom_coords)."""
    return line[12:16].strip().upper()


def het_ligand_pdb_block(pdb_text: str, key: str) -> str | None:
    """Raw HETATM lines for one list_het_ligands() group, keyed the same way
    (resname:chain:resseq). Used to pull a bound ligand's real coordinates
    back out of the ORIGINAL uploaded receptor for RMSD validation, since
    clean_receptor may have already stripped it out of the docking receptor.
    None if the key doesn't match anything.

    Alternate-location copies other than "A" are dropped, same as
    clean_receptor - otherwise a ligand modelled with alt confs would come
    back with 2-3x its real atom count, which makes every RMSD come out None
    (the docked pose can never match that composition)."""
    lines = [
        line
        for tag, line in _atom_records(pdb_text)
        if tag == "HETATM"
        and line[16] in (" ", "A")
        and f"{_resname(line)}:{line[21].strip() or '_'}:{line[22:27].strip()}" == key
    ]
    if not lines:
        return None
    return "\n".join(lines) + "\nEND\n"


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


@dataclass
class ChainInfo:
    """Per-chain atom counts, so the UI can offer a co-crystallized peptide/
    peptidomimetic ligand (its own ATOM chain, not a HETATM group) for
    exclusion from the receptor - list_het_ligands only sees HETATM residues."""

    chain: str
    n_atoms: int
    n_hetatm: int
    n_residues: int

    def as_dict(self) -> dict:
        return {
            "chain": self.chain,
            "n_atoms": self.n_atoms,
            "n_hetatm": self.n_hetatm,
            "n_residues": self.n_residues,
        }


def list_chains(pdb_text: str) -> list[ChainInfo]:
    """Atom/residue counts per chain, biggest first (the biggest is almost
    always the target protein; small ATOM chains are candidate peptide
    ligands)."""
    stats: dict[str, dict] = {}
    for tag, line in _atom_records(pdb_text):
        chain = line[21].strip() or "_"
        d = stats.setdefault(chain, {"n_atoms": 0, "n_hetatm": 0, "residues": set()})
        if tag == "HETATM":
            d["n_hetatm"] += 1
        else:
            d["n_atoms"] += 1
        d["residues"].add(line[22:27])
    out = [
        ChainInfo(chain, d["n_atoms"], d["n_hetatm"], len(d["residues"]))
        for chain, d in stats.items()
    ]
    out.sort(key=lambda c: -(c.n_atoms + c.n_hetatm))
    return out


def chains_near_point(
    pdb_text: str, center: tuple[float, float, float], cutoff: float = 15.0
) -> set[str]:
    """Chain IDs with at least one atom within `cutoff` A of `center`.

    Crystal asymmetric units often hold several copies of the same protein
    (or unrelated neighbour chains); this finds the ones actually touching
    the chosen site so the rest can be dropped as "extra monomers".
    """
    cx, cy, cz = center
    cutoff2 = cutoff * cutoff
    near: set[str] = set()
    for tag, line in _atom_records(pdb_text):
        chain = line[21].strip() or "_"
        if chain in near:
            continue
        x, y, z = _xyz(line)
        if (x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2 <= cutoff2:
            near.add(chain)
    return near


def clean_receptor(
    pdb_text: str,
    keep_waters: bool = False,
    keep_metals: bool = True,
    keep_het: Iterable[str] = (),
    strip_ligand: bool = True,
    exclude_chains: Iterable[str] = (),
    site_center: tuple[float, float, float] | None = None,
    chain_cutoff: float = 15.0,
) -> str:
    """Strip everything smina should not treat as the rigid receptor.

    Keeps all ATOM records, drops alternate locations other than A. HETATM
    groups (the usual co-crystallized small-molecule ligand) are dropped by
    default (`strip_ligand=True`), except metals and anything named in
    `keep_het`. `exclude_chains` additionally drops whole chains outright -
    for a co-crystallized peptide/peptidomimetic ligand, which is its own
    ATOM chain rather than a HETATM group and so isn't touched by
    `strip_ligand`.

    When `site_center` is given, whole chains with no atom within
    `chain_cutoff` A of it are dropped too - extra copies of the protein (or
    other neighbour chains) elsewhere in the asymmetric unit that cannot
    reach the docking box anyway.
    """
    keep_het = {k.upper() for k in keep_het}
    exclude_chains = {c.upper() for c in exclude_chains}
    keep_chains = chains_near_point(pdb_text, site_center, chain_cutoff) if site_center else None
    kept: list[str] = []
    for line in pdb_text.splitlines():
        tag = line[:6].strip()
        if tag in ("TER", "END"):
            continue
        if tag not in ("ATOM", "HETATM"):
            continue
        chain = line[21].strip() or "_"
        if chain.upper() in exclude_chains:
            continue
        if keep_chains is not None and chain not in keep_chains:
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
            elif not strip_ligand:
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
