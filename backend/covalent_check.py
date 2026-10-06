"""Ligand-side sanity checks for MOE covalent docking.

MOE's covalent placement reacts *every* group in the ligand that matches the
reaction's ligand reactant, and silently docks whatever it finds. Two ways
that goes wrong, both seen on 3PBS (aztreonam acyl-enzyme):

* the ligand isn't in its pre-reaction form - the PDB Chemical Component
  Dictionary stores a covalently bound ligand as the cut adduct, so an
  acyl-enzyme's carbonyl comes back as an aldehyde. An ester-forming
  reaction then can't touch it and acylates some other carboxylic acid on
  the ligand instead;
* the ligand has several matching groups - MOE docks one product per group,
  and the top-scored pose may be bonded through the wrong one.

`reactant_warning` flags both before docking, `acid_form_from_link` rebuilds
the pre-reaction acid for a self-dock of a LINKed acyl adduct, and
`pose_attachment` reports which ligand atom each docked pose is bonded
through.
"""

from __future__ import annotations

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem


def _rxn_molblocks(text: str) -> tuple[int, list[str]]:
    """(number of reactants, reactant+product MOL blocks) from an MDL RXN file."""
    parts = text.split("$MOL")
    header = parts[0].splitlines()
    n_react = int(header[4][:3]) if len(header) > 4 and header[4][:3].strip().isdigit() else 0
    # each block starts right after the "$MOL" line; its first line (the
    # molecule name) is often blank, so drop only that one line break
    blocks = [p.split("\n", 1)[1] if "\n" in p else "" for p in parts[1:]]
    return n_react, blocks


def _normalize_als(block: str) -> str:
    """Rewrite "M  ALS" atom-list lines into the fixed columns RDKit reads.
    MOESketch pads the element list one column off from the MDL spec, which
    RDKit reads as an empty element."""
    out = []
    for line in block.splitlines():
        if line.startswith("M  ALS"):
            tok = line[6:].split()
            idx, n, flag, elems = int(tok[0]), int(tok[1]), tok[2], tok[3:]
            line = f"M  ALS {idx:3d}{n:3d} {flag} " + "".join(f"{e:<4}" for e in elems)
        out.append(line)
    return "\n".join(out) + "\n"


def ligand_reactant_query(rxn_path: Path) -> Chem.Mol | None:
    """Query molecule for the ligand-side reactant of a MOE covalent .rxn.

    MOE's convention is reactant 1 = receptor residue, reactant 2 = ligand.
    None if the file can't be read that way (the check is then skipped).
    """
    try:
        n_react, blocks = _rxn_molblocks(Path(rxn_path).read_text())
    except (OSError, ValueError):
        return None
    if n_react < 2 or len(blocks) < 2:
        return None
    try:
        query = Chem.MolFromMolBlock(_normalize_als(blocks[1]), sanitize=False, removeHs=False)
    except Exception:
        return None
    if query is None or query.GetNumAtoms() == 0:
        return None
    query.UpdatePropertyCache(strict=False)
    return query


def count_reactive_sites(mol: Chem.Mol, query: Chem.Mol) -> int:
    """Number of distinct ligand groups matching the reactant query.

    Reactant drawings mix explicit-H atom lists (e.g. [H,C]) with "as drawn"
    substitution counts that only hold on a heavy-atom graph, so match both
    representations and keep whichever finds more. Matches are deduplicated
    by the heavy atoms involved (AddHs appends hydrogens, so heavy-atom
    indices are shared between the two graphs).
    """
    heavy = Chem.RemoveHs(mol, sanitize=False)
    heavy.UpdatePropertyCache(strict=False)
    with_h = Chem.AddHs(heavy)
    best = 0
    for target in (heavy, with_h):
        sites = {
            frozenset(i for i in match if target.GetAtomWithIdx(i).GetAtomicNum() > 1)
            for match in target.GetSubstructMatches(query, useChirality=False, uniquify=True, maxMatches=1000)
        }
        best = max(best, len(sites))
    return best


def reactant_warning(molblock: str, rxn_path: Path | None) -> str | None:
    """Human-readable warning when the ligand matches the reaction's ligand
    reactant zero times or more than once; None when it's exactly one (or the
    check can't run)."""
    if rxn_path is None:
        return None
    query = ligand_reactant_query(rxn_path)
    mol = Chem.MolFromMolBlock(molblock, removeHs=False)
    if query is None or mol is None:
        return None
    n = count_reactive_sites(mol, query)
    name = Path(rxn_path).stem
    if n == 0:
        return (f"ligand has no group matching the ligand reactant of {name} - "
                "is it drawn in its pre-reaction form? (e.g. a carboxylic acid, not the "
                "aldehyde a PDB acyl adduct comes back as)")
    if n > 1:
        return (f"ligand has {n} groups matching {name}; MOE docks one adduct per group - "
                "check which atom each pose is bonded through")
    return None


# ---------------------------------------------------------------------------
# self-dock of a LINKed adduct
# ---------------------------------------------------------------------------


def _link_atoms(pdb_text: str, het_key: str) -> list[str]:
    """Atom names on the het side of every LINK record touching het_key
    (resname:chain:resseq) whose other side is not another het group."""
    if het_key.count(":") != 2:
        return []
    resname, chain, resseq = het_key.split(":")
    chain = "" if chain == "_" else chain
    names = []
    for line in pdb_text.splitlines():
        if not line.startswith("LINK"):
            continue
        line = line.ljust(80)
        sides = [
            (line[12:16].strip(), line[17:20].strip(), line[21].strip(), line[22:27].strip()),
            (line[42:46].strip(), line[47:50].strip(), line[51].strip(), line[52:57].strip()),
        ]
        for mine, other in (sides, sides[::-1]):
            if (mine[1], mine[2], mine[3]) == (resname, chain, resseq) and other[1] != resname:
                names.append(mine[0])
    return names


def acid_form_from_link(pdb_text: str, het_key: str, ccd_smiles: str) -> str | None:
    """Pre-reaction SMILES for a bound acyl adduct, else None.

    When het_key is LINKed to the protein through a carbonyl carbon that has
    only one other heavy neighbour (an aldehyde once the bond to the residue
    is cut), add the OH back so an ester/thioester-forming reaction reacts
    there. Bond orders come from the CCD SMILES, stereo from the crystal
    coordinates.
    """
    link_names = _link_atoms(pdb_text, het_key)
    if not link_names:
        return None
    from . import prep

    block = prep.het_ligand_pdb_block(pdb_text, het_key)
    template = Chem.MolFromSmiles(ccd_smiles) if ccd_smiles else None
    if block is None or template is None:
        return None
    xtal = Chem.MolFromPDBBlock(block, removeHs=True, sanitize=False)
    if xtal is None:
        return None
    try:
        mol = AllChem.AssignBondOrdersFromTemplate(template, xtal)
    except Exception:
        return None
    # template atom -> crystal atom, so the edit lands on the CCD molecule
    # and keeps its chiral tags (crystal geometry can be too flat to
    # re-derive them)
    match = mol.GetSubstructMatch(template)
    if not match:
        return None

    rw = Chem.RWMol(template)
    changed = False
    for t_idx, x_idx in enumerate(match):
        info = mol.GetAtomWithIdx(x_idx).GetPDBResidueInfo()
        atom = template.GetAtomWithIdx(t_idx)
        if info is None or info.GetName().strip() not in link_names or atom.GetSymbol() != "C":
            continue
        has_carbonyl = any(
            b.GetBondType() == Chem.BondType.DOUBLE and b.GetOtherAtom(atom).GetSymbol() == "O"
            for b in atom.GetBonds()
        )
        if has_carbonyl and atom.GetDegree() == 2:
            o = rw.AddAtom(Chem.Atom(8))
            rw.AddBond(t_idx, o, Chem.BondType.SINGLE)
            rw.GetAtomWithIdx(t_idx).SetNumExplicitHs(0)
            changed = True
    if not changed:
        return None
    out = rw.GetMol()
    try:
        Chem.SanitizeMol(out)
    except Exception:
        return None
    return Chem.MolToSmiles(out)


# ---------------------------------------------------------------------------
# docked pose
# ---------------------------------------------------------------------------


def pose_attachment(mol: Chem.Mol) -> str | None:
    """SMILES of the ligand around the atom a MOE covalent pose is bonded
    through (that atom and everything within two bonds of it, residue atoms
    excluded, the residue atom written as "*"), e.g. "*C(=O)C(C)N" - enough to tell which of several
    reactive groups was used. None for non-covalent poses.

    MOE writes the residue fragment first and flags its atoms 0 at the start
    of the useAtoms field; that leading run of zeros is the residue.
    """
    if not mol.HasProp("useAtoms"):
        return None
    flags = mol.GetProp("useAtoms").split()
    n_res = 0
    for f in flags:
        if f != "0":
            break
        n_res += 1
    if n_res == 0 or n_res >= mol.GetNumAtoms():
        return None
    residue = set(range(n_res))
    attach = None
    for i in residue:
        for nb in mol.GetAtomWithIdx(i).GetNeighbors():
            if nb.GetIdx() not in residue:
                attach, nucleophile = nb.GetIdx(), i
                break
        if attach is not None:
            break
    if attach is None:
        return None
    env = {attach}
    frontier = {attach}
    for _ in range(2):
        frontier = {
            nb.GetIdx()
            for i in frontier
            for nb in mol.GetAtomWithIdx(i).GetNeighbors()
            if nb.GetIdx() not in residue and nb.GetAtomicNum() > 1
        } - env
        env |= frontier
    try:
        symbols = [a.GetSymbol() for a in mol.GetAtoms()]
        symbols[nucleophile] = "*"
        return Chem.MolFragmentToSmiles(mol, atomsToUse=sorted(env | {nucleophile}), canonical=True,
                                        atomSymbols=symbols)
    except Exception:
        return None
