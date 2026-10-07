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
        # A QuickPrep'd receptor carries ligand hydrogens; site size and atom
        # count stay on heavy atoms, as for a crystal structure.
        if _element(line) in ("H", "D"):
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
    keep_hydrogens: bool = False,
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

    `keep_hydrogens` is for a receptor that was already protonated on
    purpose (MOE QuickPrep): its hydrogens and the formal-charge columns are
    kept instead of being stripped as crystal-structure noise.
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
        if keep_hydrogens:
            kept.append(line[:80])
            continue
        # Hydrogens from the crystal structure are unreliable; smina adds its own.
        element = line[76:78].strip().upper()
        if element == "H" or (not element and line[12:16].strip().startswith("H")):
            continue
        kept.append(line[:66].ljust(66))
    kept.append("END")
    return "\n".join(kept) + "\n"


# --------------------------------------------------------------------------
# QuickPrep report
# --------------------------------------------------------------------------


def _formal_charge(line: str) -> int:
    """PDB columns 79-80, e.g. "1+" / "1-"; blank means neutral."""
    field = line[78:80].strip()
    if len(field) == 2 and field[0].isdigit() and field[1] in "+-":
        return int(field[0]) * (1 if field[1] == "+" else -1)
    return 0


def _protonation_state(resname: str, atoms: dict[str, str]) -> str | None:
    """Protonation state worth reporting for one residue (atom name -> line),
    or None when it is the usual one at pH 7. HIS is always reported, since
    which nitrogen carries the H is exactly what Protonate3D decides."""
    names = set(atoms)
    charge = sum(_formal_charge(line) for line in atoms.values())
    if resname == "HIS":
        d1, e2 = "HD1" in names, "HE2" in names
        return "HIP (+1)" if d1 and e2 else "HID" if d1 else "HIE" if e2 else None
    if resname == "ASP" and charge == 0 and names & {"HD1", "HD2"}:
        return "ASH (neutral)"
    if resname == "GLU" and charge == 0 and names & {"HE1", "HE2"}:
        return "GLH (neutral)"
    if resname == "LYS" and charge == 0 and "HZ1" in names:
        return "LYN (neutral)"
    if resname == "ARG" and charge == 0 and "HH11" in names:
        return "ARN (neutral)"
    if resname == "CYS" and charge < 0:
        return "CYM (-1)"
    if resname == "TYR" and charge < 0:
        return "TYM (-1)"
    return None


def quickprep_report(
    original_pdb: str,
    prepared_pdb: str,
    site_coords: list[tuple[float, float, float]] = (),
    site_radius: float = 6.0,
) -> dict:
    """What MOE QuickPrep changed: hydrogens added, waters removed, residues
    built, how far heavy atoms moved in the minimization, and the
    protonation states it chose. Residues within `site_radius` of any
    `site_coords` atom (the site ligand) are flagged and listed first."""
    def records(text):
        out = {}
        for tag, line in _atom_records(text):
            if line[16] not in (" ", "A"):
                continue
            key = (line[21], line[22:27].strip(), _resname(line), _atom_name(line))
            out[key] = line
        return out

    orig, prepd = records(original_pdb), records(prepared_pdb)
    is_h = lambda line: _element(line) in ("H", "D")
    n_water = lambda recs: len({k[:3] for k in recs if k[2] in WATER})

    shifts = []
    for key, line in prepd.items():
        if is_h(line) or key not in orig or is_h(orig[key]):
            continue
        a, b = _xyz(orig[key]), _xyz(line)
        shifts.append(sum((a[i] - b[i]) ** 2 for i in range(3)) ** 0.5)

    residues: dict[tuple, dict[str, str]] = {}
    for key, line in prepd.items():
        if line.startswith("ATOM"):
            residues.setdefault(key[:3], {})[key[3]] = line

    states = []
    r2 = site_radius ** 2
    for (chain, resseq, resname), atoms in residues.items():
        state = _protonation_state(resname, atoms)
        if state is None:
            continue
        near = any(
            sum((p - c) ** 2 for p, c in zip(_xyz(line), site)) <= r2
            for line in atoms.values() for site in site_coords
        )
        states.append({"residue": f"{resname} {chain.strip() or '_'}{resseq}", "state": state, "near_site": near})
    states.sort(key=lambda s: not s["near_site"])

    heavy = lambda recs: sum(1 for line in recs.values() if not is_h(line))
    return {
        "hydrogens_added": sum(1 for l in prepd.values() if is_h(l)) - sum(1 for l in orig.values() if is_h(l)),
        "heavy_atoms_added": max(0, heavy(prepd) - heavy(orig)),
        "waters_removed": n_water(orig) - n_water(prepd),
        "rmsd_heavy": round((sum(d * d for d in shifts) / len(shifts)) ** 0.5, 3) if shifts else 0.0,
        "max_shift": round(max(shifts), 3) if shifts else 0.0,
        "n_moved": sum(1 for d in shifts if d > 0.5),
        "protonation": states,
    }


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


def parse_residue_spec(spec: str) -> tuple[str, str, int] | None:
    """Parse 'chain:resname:resnum' or ':resname:resnum' format.

    Returns:
        (chain, resname, resnum) or None if invalid

    Examples:
        "A:CYS:189" → ("A", "CYS", 189)
        ":SER:25"   → ("", "SER", 25)
    """
    match = re.match(r'^([A-Z]?):([A-Z]{3}):(\d+)$', spec.strip())
    if not match:
        return None

    chain, resname, resnum_str = match.groups()
    return (chain, resname, int(resnum_str))


def validate_residue_exists(receptor_pdb_text: str, spec: str) -> bool:
    """Check if specified residue exists in PDB text."""
    parsed = parse_residue_spec(spec)
    if not parsed:
        return False

    chain, resname, resnum = parsed

    # Look for ATOM/HETATM lines matching the residue
    # PDB format: columns 22 (chain), 18-20 (resname), 23-26 (resnum)
    for line in receptor_pdb_text.splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        if len(line) < 27:
            continue

        line_chain = line[21:22].strip()
        line_resname = line[17:20].strip()
        try:
            line_resnum = int(line[22:26].strip())
        except ValueError:
            continue

        # Match chain (if specified), resname, and resnum
        if (not chain or line_chain == chain) and line_resname == resname and line_resnum == resnum:
            return True

    return False


def list_residues_by_name(receptor_pdb_text: str, resname: str) -> list[str]:
    """List all occurrences of residue name in PDB text.

    Returns list of specs like ["A:CYS:25", "A:CYS:51", "B:CYS:102"]
    """
    seen = set()
    result = []

    for line in receptor_pdb_text.splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        if len(line) < 27:
            continue

        line_chain = line[21:22].strip()
        line_resname = line[17:20].strip()
        try:
            line_resnum = int(line[22:26].strip())
        except ValueError:
            continue

        if line_resname == resname:
            spec = f"{line_chain or '_'}:{resname}:{line_resnum}"
            if spec not in seen:
                seen.add(spec)
                result.append(spec)

    return sorted(result)
