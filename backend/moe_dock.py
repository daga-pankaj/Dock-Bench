#!/usr/bin/env python3
"""MOE Docking command-line wrapper."""

import subprocess
import shutil
from pathlib import Path
from dataclasses import dataclass
from typing import Tuple, Optional


class MoeBatchNotFound(RuntimeError):
    pass


class MoeDockingFailed(RuntimeError):
    pass


@dataclass
class MoePose:
    """MOE docking pose result."""
    rank: int
    score: float
    rmsd: Optional[float] = None


def find_moebatch() -> Optional[str]:
    """Locate moebatch executable."""
    # Check common install paths
    candidates = [
        Path("/Applications/moe2024.0604/bin/moebatch"),
        Path("/Applications/moe/bin/moebatch"),
    ]

    # Also check versioned installs
    apps = Path("/Applications")
    if apps.is_dir():
        for moe_dir in sorted(apps.glob("moe*"), reverse=True):
            candidates.append(moe_dir / "bin" / "moebatch")

    for candidate in candidates:
        if candidate.exists():
            return str(candidate)

    return shutil.which("moebatch")


def _get_nucleophile_atom_name(resname: str) -> str:
    """Get nucleophilic atom name for residue type.

    Args:
        resname: 3-letter residue code (CYS, SER, LYS, etc.)

    Returns:
        Atom name for nucleophilic atom (SG, OG, NZ, etc.)
    """
    nucleophile_map = {
        'CYS': 'SG',   # Cysteine sulfur
        'SER': 'OG',   # Serine oxygen
        'THR': 'OG1',  # Threonine oxygen
        'TYR': 'OH',   # Tyrosine hydroxyl
        'LYS': 'NZ',   # Lysine nitrogen
    }

    resname_upper = resname.upper()
    if resname_upper not in nucleophile_map:
        raise ValueError(
            f"Unsupported nucleophilic residue: {resname}. "
            f"Supported: {', '.join(nucleophile_map.keys())}"
        )

    return nucleophile_map[resname_upper]


def _build_moe_atom_selection(residue_spec: str) -> tuple[str, int, str]:
    """Convert residue spec to MOE atom selection components.

    Args:
        residue_spec: "chain:resname:resnum" (e.g., "A:CYS:189")

    Returns:
        Tuple of (chain, resnum, atom_name) for SVL filtering
    """
    from . import prep

    parsed = prep.parse_residue_spec(residue_spec)
    if not parsed:
        raise ValueError(f"Invalid residue spec: {residue_spec}")

    chain, resname, resnum = parsed
    atom_name = _get_nucleophile_atom_name(resname)

    return (chain, resnum, atom_name)


def _build_simple_atom_filter(chain: str, resnum: int, atom_name: str) -> str:
    """SVL filter for the nucleophile atom; rUID matches the PDB residue number.

    Chain is ignored: MOE renames chains on PDB import (tag becomes the file name).
    """
    return f"(rUID aResidue protein_atoms == {resnum} and aName protein_atoms == '{atom_name}')"


def build_svl_script(
    receptor: Path,
    ligand: Path,
    out_mdb: Path,
    out_sdf: Path,
    center: Tuple[float, float, float],
    size: Tuple[float, float, float],
    num_poses: int = 5,
    scoring: str = "London dG",
    placement: str = "Triangle Matcher",
    refinement: str = "Rigid Receptor",
    ph4_file: Optional[Path] = None,
    template_path: Optional[Path] = None,
    covalent_reaction_file: Optional[Path] = None,
    covalent_residue: Optional[str] = None,
) -> str:
    """Generate SVL script from template.

    Args:
        covalent_reaction_file: Path to .rxn file (e.g., .../michael_acceptor_s.rxn)
        covalent_residue: Target residue spec "chain:resname:resnum" (e.g., "A:CYS:189")
    """
    if template_path is None:
        template_path = Path(__file__).parent / "moe_dock_template.svl"

    template = template_path.read_text()

    # Calculate site radius from box size (sphere enclosing box)
    import math
    site_radius = math.sqrt((size[0]/2)**2 + (size[1]/2)**2 + (size[2]/2)**2)

    # Build covalent options
    covalent_options = ""
    if covalent_reaction_file and covalent_residue:
        # Convert residue spec to components
        chain, resnum, atom_name = _build_moe_atom_selection(covalent_residue)
        filter_condition = _build_simple_atom_filter(chain, resnum, atom_name)

        covalent_options = f"""
    // Covalent docking - residue {covalent_residue}
    local rxsite_atoms = protein_atoms | {filter_condition};
    if length rxsite_atoms == 0 then
        exit 'Error: {covalent_residue} atom {atom_name} not found (rUID={resnum})';
    endif;

    opt.dockmode = 'Covalent';
    opt.placement = 'Covalent Docking';
    opt.placement_opt = [
        rxnFile: '{covalent_reaction_file.absolute()}',
        rxsite: 'asel',
        mode: 'Covalent',
        outrmsd: 0,
        silent: 0,
        debug: 0
    ];
    opt.refine_opt = [
        fixrec: 'Fix',
        cutoff: 6,
        wholeres: 1,
        mmgbvi: 1,
        packsidechains: 1,
        rigidlig: 0,
        rigidPlanar: 0,
        tether: 10,
        gtest: 0.01,
        maxit: 500,
        OverrideSetup: 1,
        k_potl: 100,
        roffset: 0.4,
        silent: 0,
        debug: 0
    ];
    opt.rxsite = rxsite_atoms(1);
"""

    script = template.format(
        RECEPTOR_FILE=str(receptor.absolute()),
        LIGAND_FILE=str(ligand.absolute()),
        OUTPUT_MDB=str(out_mdb.absolute()),
        OUTPUT_SDF=str(out_sdf.absolute()),
        CENTER_X=center[0],
        CENTER_Y=center[1],
        CENTER_Z=center[2],
        SITE_RADIUS=site_radius,
        NUM_POSES=num_poses,
        SCORING=scoring,
        PLACEMENT=placement,
        REFINEMENT=refinement,
        PH4_FILE=str(ph4_file.absolute()) if ph4_file else '',
        COVALENT_OPTIONS=covalent_options,
    )

    return script


def dock(
    receptor: Path,
    ligand: Path,
    out_mdb: Path,
    out_sdf: Path,
    center: Tuple[float, float, float],
    size: Tuple[float, float, float],
    num_poses: int = 5,
    scoring: str = "London dG",
    placement: str = "Triangle Matcher",
    refinement: str = "Rigid Receptor",
    log_path: Optional[Path] = None,
    timeout: int = 3600,
    auto_detect_site: bool = True,
    ph4_file: Optional[Path] = None,
    covalent_reaction_file: Optional[Path] = None,
    covalent_residue: Optional[str] = None,
) -> Path:
    """Run MOE docking via moebatch.

    Args:
        receptor: Receptor PDB/MOE file
        ligand: Ligand SDF/MOE file
        out_mdb: Output MOE database (temp, deleted after SDF export)
        out_sdf: Output SDF file with poses and scores
        center: Box center (x, y, z)
        size: Box size (x, y, z)
        num_poses: Number of poses to generate
        scoring: Scoring function (London dG, GBVI/WSA dG, Affinity dG, ASE)
        placement: Placement method (Triangle Matcher, Alpha PMI, Pharmacophore)
        refinement: Refinement (None, Rigid Receptor, Induced Fit, Rigid Body)
        log_path: Optional log file path
        timeout: Max runtime in seconds
        auto_detect_site: If True, use HETATM center with radius 10Å if found
        ph4_file: Optional MOE pharmacophore (.ph4) file for constraint-based docking
        covalent_reaction_file: Path to .rxn file for covalent docking
        covalent_residue: Target residue spec "chain:resname:resnum"

    Returns:
        Path to output SDF file
    """
    moebatch = find_moebatch()
    if not moebatch:
        raise MoeBatchNotFound("moebatch not found - install MOE")

    # Auto-detect binding site from HETATM if requested
    if auto_detect_site and receptor.suffix.lower() == '.pdb':
        from . import prep
        pdb_text = receptor.read_text()
        het_ligands = prep.list_het_ligands(pdb_text, min_atoms=4)
        # Use first non-ion ligand
        site_ligand = next((lig for lig in het_ligands if not lig.is_ion), None)
        if site_ligand:
            center = tuple(site_ligand.center)
            # Use radius of 10 Å, convert to box size (diameter = 20)
            size = (20.0, 20.0, 20.0)
            if log_path:
                log_msg = f"Auto-detected site from HETATM {site_ligand.resname} (chain {site_ligand.chain}): center={center}, radius=10Å\n"
                log_path.write_text(log_msg)

    # Convert PDB to MOE format if needed
    receptor_for_moe = receptor
    if receptor.suffix.lower() == '.pdb':
        receptor_moe = receptor.with_suffix('.moe')

        # Inline SVL for PDB to MOE conversion
        # Note: Uses ReadPDB/WriteMOE/Close - simple one-liner via -exec
        svl_cmd = f"ReadPDB '{receptor.absolute()}'; WriteMOE '{receptor_moe.absolute()}'; Close [];"

        try:
            result = subprocess.run(
                [moebatch, '-exec', svl_cmd],
                capture_output=True,
                timeout=60,
                text=True
            )

            if not receptor_moe.exists():
                stderr = (result.stderr or result.stdout or "").strip()
                raise MoeDockingFailed(f"PDB to MOE conversion failed: {stderr or 'no output file created'}")

            receptor_for_moe = receptor_moe

        except subprocess.TimeoutExpired:
            raise MoeDockingFailed("PDB to MOE conversion timed out after 60s")
        except subprocess.CalledProcessError as e:
            stderr = (e.stderr or e.stdout or "").strip()
            raise MoeDockingFailed(f"moebatch conversion error: {stderr or 'unknown error'}")

    # Generate SVL script that exports to SDF directly
    script = build_svl_script(
        receptor_for_moe, ligand, out_mdb, out_sdf, center, size,
        num_poses, scoring, placement, refinement, ph4_file,
        covalent_reaction_file=covalent_reaction_file,
        covalent_residue=covalent_residue,
    )

    # Write script to file and run via moebatch (needs file for Open/atom selection)
    script_path = out_mdb.parent / f"{out_mdb.stem}_dock.svl"
    script_path.write_text(script)

    cmd = [moebatch, "-load", str(script_path), "-exec", "main []"]

    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=timeout,
    )

    if log_path:
        log_path.write_text(
            f"$ {' '.join(cmd)}\n\n{proc.stdout}\n{proc.stderr}"
        )

    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()
        raise MoeDockingFailed(f"{tail[-1] if tail else f'moebatch exited {proc.returncode}'} (script: {script_path})")

    if not out_sdf.exists():
        # Keep script on failure for debugging
        raise MoeDockingFailed(f"moebatch produced no SDF output (script: {script_path})")

    # Clean up script on success
    if script_path.exists():
        script_path.unlink()
    if out_mdb.exists():
        out_mdb.unlink()
    # Clean up converted receptor if created
    if receptor_for_moe != receptor and receptor_for_moe.exists():
        receptor_for_moe.unlink()

    return out_sdf


def parse_mdb_results(mdb_path: Path) -> list[MoePose]:
    """Parse MOE database results (requires MOE Python API or db parser)."""
    # Placeholder - would need MOE's Python API or custom MDB parser
    # For now return empty list
    return []


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="MOE Docking CLI")
    parser.add_argument("--receptor", required=True, type=Path, help="Receptor PDB")
    parser.add_argument("--ligand", required=True, type=Path, help="Ligand SDF")
    parser.add_argument("--out", required=True, type=Path, help="Output MDB")
    parser.add_argument("--out_sdf", type=Path, help="Output SDF")
    parser.add_argument("--center_x", type=float, required=True)
    parser.add_argument("--center_y", type=float, required=True)
    parser.add_argument("--center_z", type=float, required=True)
    parser.add_argument("--size_x", type=float, required=True)
    parser.add_argument("--size_y", type=float, required=True)
    parser.add_argument("--size_z", type=float, required=True)
    parser.add_argument("--num_poses", type=int, default=5)
    parser.add_argument("--scoring", default="London dG",
                       choices=["London dG", "GBVI/WSA dG", "Affinity dG", "ASE", "Alpha HB"])
    parser.add_argument("--placement", default="Triangle Matcher",
                       choices=["Triangle Matcher", "Alpha PMI", "Pharmacophore",
                               "Proxy Triangle", "Alpha Triangle"])
    parser.add_argument("--refinement", default="Rigid Receptor",
                       choices=["None", "Rigid Receptor", "Induced Fit", "Rigid Body"])
    parser.add_argument("--log", type=Path, help="Log file path")
    parser.add_argument("--timeout", type=int, default=3600, help="Timeout in seconds")

    args = parser.parse_args()

    center = (args.center_x, args.center_y, args.center_z)
    size = (args.size_x, args.size_y, args.size_z)
    out_sdf = args.out_sdf or args.out.with_suffix(".sdf")

    dock(
        args.receptor,
        args.ligand,
        args.out,
        out_sdf,
        center,
        size,
        num_poses=args.num_poses,
        scoring=args.scoring,
        placement=args.placement,
        refinement=args.refinement,
        log_path=args.log,
        timeout=args.timeout,
    )

    print(f"Docking complete: {out_sdf}")
