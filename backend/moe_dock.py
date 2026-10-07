#!/usr/bin/env python3
"""MOE Docking command-line wrapper."""

import re
import subprocess
import threading
import time
from pathlib import Path
from typing import Tuple, Optional

from . import config


class MoeBatchNotFound(RuntimeError):
    pass


class MoeDockingFailed(RuntimeError):
    pass


# Server-wide cap on concurrent moebatch processes (see
# config.MOE_MAX_CONCURRENT); every moebatch call goes through _run_moebatch.
_LICENSE_SLOTS = threading.BoundedSemaphore(config.MOE_MAX_CONCURRENT)
_LICENSE_ERROR = "Can't check out license"


def _run_moebatch(cmd: list[str], timeout: int, cwd: Path) -> subprocess.CompletedProcess:
    """Run moebatch once a local slot is free, retrying while the license
    server has no tokens left (other users hold them) for up to
    config.MOE_LICENSE_WAIT seconds. Runs in `cwd`, where Dock leaves its
    dock_lig_*_bound_*.mdb scratch files. TimeoutExpired propagates."""
    deadline = time.monotonic() + config.MOE_LICENSE_WAIT
    while True:
        with _LICENSE_SLOTS:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, cwd=cwd)
        if _LICENSE_ERROR not in (proc.stdout or "") + (proc.stderr or ""):
            return proc
        if time.monotonic() + config.MOE_LICENSE_RETRY > deadline:
            return proc
        time.sleep(config.MOE_LICENSE_RETRY)


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


# MOE's Ligand Interactions diagram (the GUI's, from prolig2d.svl) with its
# legend underneath, written as SVG. Shared by the docking script and
# render_interactions.
LID_FUNCTION = """
function DrawLigandInteractions;

local function WriteLigandInteractions [lig, rec, file]
    local key = first cat DrawLigandInteractions [[lig], [rec], []];
    local leg = call ['_DrawLigandInteractionsLegend', [[], []], 'prolig2d.svl'];
    local s1 = (gr_header key).size, s2 = (gr_header leg).size;
    local w = max [s1(1), s2(1)];
    local both = gr_create [size: [w, s1(2) + s2(2)]];
    gr_draw [both, key, [(w - s1(1)) / 2, s2(2), 1, 1]];
    gr_draw [both, leg, [(w - s2(1)) / 2, 0, 1, 1]];
    gr_fwrite [file, both, 'image/svg+xml', [0, 0, 1, 1]];
    gr_destroy [key, leg, both];
endfunction
"""

# One diagram per pose, drawn inside the docking run while the receptor is
# still loaded: no second moebatch run, so no second license checkout.
LID_RENDER = """
    local lid_db = db_Open ['{mdb}', 'read'];
    local lid_ent, lid_lig, lid_k = 0;
    for lid_ent in db_Entries lid_db loop
        lid_k = lid_k + 1;
        lid_lig = cat cAtoms mol_Create first db_ReadFields [lid_db, lid_ent, 'mol'];
        WriteLigandInteractions [lid_lig, protein_atoms, twrite ['{prefix}{{}}.svg', lid_k]];
        oDestroy uniq aChain lid_lig;
    endloop
    db_Close lid_db;
"""


def interactions_svg_path(out_sdf: Path, pose_index: int) -> Path:
    """Where the docking run leaves pose `pose_index`'s (1-based) MOE Ligand
    Interactions diagram for the poses file `out_sdf`."""
    return out_sdf.with_name(f"{out_sdf.stem}_lid_{pose_index}.svg")


def browser_svg(svg: str) -> str:
    """MOE embeds an SVG font ("sf") that browsers ignore, so name real ones;
    and give the transparent drawing the white page its colors are meant
    for, since the app's panels are dark."""
    svg = svg.replace('font-family="sf"', 'font-family="Helvetica, Arial, sans-serif"')
    return re.sub(r"(<svg\b[^>]*>)", r'\1<rect width="100%" height="100%" fill="#ffffff"/>', svg, count=1)


# MOE's Dock panel defaults: placement scored with London dG, 30 placement
# poses carried into refinement, refined poses rescored with GBVI/WSA dG.
PLACEMENT_SCORING = "London dG"
PLACEMENT_POSES = 30
DEFAULT_RESCORING = "GBVI/WSA dG"


def build_svl_script(
    receptor: Path,
    ligand: Path,
    out_mdb: Path,
    out_sdf: Path,
    center: Tuple[float, float, float],
    size: Tuple[float, float, float],
    num_poses: int = 5,
    scoring: str = DEFAULT_RESCORING,
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
        PLACEMENT_POSES=max(PLACEMENT_POSES, num_poses),
        PLACEMENT_SCORING=PLACEMENT_SCORING,
        SCORING=scoring,
        PLACEMENT=placement,
        REFINEMENT=refinement,
        PH4_FILE=str(ph4_file.absolute()) if ph4_file else '',
        COVALENT_OPTIONS=covalent_options,
        LID_FUNCTION=LID_FUNCTION,
        LID_RENDER=LID_RENDER.format(
            mdb=out_mdb.absolute(),
            prefix=str(interactions_svg_path(out_sdf, 0))[:-len("0.svg")],
        ),
    )

    return script


QUICKPREP_TEMPLATE = """#svl
function _LIGX_Execute;

global function main []
    ReadPDB '{RECEPTOR_FILE}';

    // QuickPrep keeps receptor atoms within 8 A of the ligand free during
    // minimization and fixes the rest; naming the site ligand stops it from
    // guessing (it would otherwise pick by class and drop known cofactors).
    // With no ligand at all every atom is fixed and minimization is a no-op.
    local atoms = Atoms [];
    local site_lig = atoms | rName aResidue atoms == '{SITE_RESNAME}'
        and rUID aResidue atoms == {SITE_RESSEQ};
    if length site_lig then
        oSetCollection ['LIGX_LIGAND', site_lig];
    endif

    // MOE's saved QuickPrep panel settings are used, as in the GUI.
    _LIGX_Execute [verbose: 0];

    oSetCollection ['LIGX_LIGAND', []];
    WritePDB '{OUTPUT_PDB}';
    Close [];
endfunction
"""


def quickprep(
    receptor_pdb: Path,
    out_pdb: Path,
    site_ligand: Optional[Tuple[str, int]] = None,
    log_path: Optional[Path] = None,
    timeout: int = 1800,
) -> Path:
    """Run MOE QuickPrep (Structure Preparation, Protonate3D, distant solvent
    removal, tethered minimization) on a receptor and write it back as PDB,
    hydrogens and formal charges included.

    Args:
        receptor_pdb: Receptor PDB, co-crystal ligand(s) included
        out_pdb: Prepared receptor output
        site_ligand: (resname, resseq) of the HETATM ligand defining the site
        log_path: Optional log file path
        timeout: Max runtime in seconds
    """
    moebatch = config.find_moebatch()
    if not moebatch:
        raise MoeBatchNotFound("moebatch not found - install MOE")

    resname, resseq = site_ligand or ("", 0)
    script = QUICKPREP_TEMPLATE.format(
        RECEPTOR_FILE=str(receptor_pdb.absolute()),
        OUTPUT_PDB=str(out_pdb.absolute()),
        SITE_RESNAME=resname,
        SITE_RESSEQ=resseq,
    )
    script_path = out_pdb.with_suffix(".svl")
    script_path.write_text(script)

    cmd = [moebatch, "-load", str(script_path), "-exec", "main []"]
    try:
        proc = _run_moebatch(cmd, timeout, cwd=out_pdb.parent)
    except subprocess.TimeoutExpired:
        raise MoeDockingFailed(f"QuickPrep timed out after {timeout}s (script: {script_path})")

    if log_path:
        log_path.write_text(f"$ {' '.join(cmd)}\n\n{proc.stdout}\n{proc.stderr}")

    if proc.returncode != 0 or not out_pdb.exists():
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()
        reason = tail[-1] if tail else f"moebatch exited {proc.returncode}"
        raise MoeDockingFailed(f"QuickPrep failed: {reason} (script: {script_path})")

    script_path.unlink()
    return out_pdb


RENDER_INTERACTIONS_TEMPLATE = """#svl
function Open;
{LID_FUNCTION}
global function main []
    ReadPDB '{RECEPTOR_FILE}';
    local rec = Atoms [];
    Open '{POSE_FILE}';
    WriteLigandInteractions [diff [Atoms [], rec], rec, '{OUTPUT_SVG}'];
endfunction
"""


def render_interactions(receptor_pdb: Path, pose_sdf: Path, out_svg: Path, timeout: int = 120) -> Path:
    """Draw MOE's Ligand Interactions diagram for one pose (an SDF holding just
    that pose) against a PDB receptor, outside of a docking run - for poses
    docked before diagrams were drawn at docking time. Takes a license slot
    like any moebatch run."""
    moebatch = config.find_moebatch()
    if not moebatch:
        raise MoeBatchNotFound("moebatch not found - install MOE")
    script_path = out_svg.with_suffix(".svl")
    script_path.write_text(RENDER_INTERACTIONS_TEMPLATE.format(
        LID_FUNCTION=LID_FUNCTION,
        RECEPTOR_FILE=str(receptor_pdb.absolute()),
        POSE_FILE=str(pose_sdf.absolute()),
        OUTPUT_SVG=str(out_svg.absolute()),
    ))
    try:
        proc = _run_moebatch([moebatch, "-load", str(script_path), "-exec", "main []"], timeout, cwd=out_svg.parent)
    except subprocess.TimeoutExpired:
        raise MoeDockingFailed(f"Ligand Interactions timed out after {timeout}s")
    finally:
        script_path.unlink(missing_ok=True)
    if proc.returncode != 0 or not out_svg.exists():
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()
        raise MoeDockingFailed(tail[-1] if tail else f"moebatch exited {proc.returncode}")
    return out_svg


def dock(
    receptor: Path,
    ligand: Path,
    out_mdb: Path,
    out_sdf: Path,
    center: Tuple[float, float, float],
    size: Tuple[float, float, float],
    num_poses: int = 5,
    scoring: str = DEFAULT_RESCORING,
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
        scoring: Rescoring function applied after refinement (GBVI/WSA dG,
            London dG, Affinity dG, ASE); placement always uses London dG
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
    moebatch = config.find_moebatch()
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
        # Per ligand: ligands of one job share `receptor`, and each run
        # deletes its converted copy when done.
        receptor_moe = out_mdb.with_name(f"{out_mdb.stem}_receptor.moe")

        # Inline SVL for PDB to MOE conversion
        # Note: Uses ReadPDB/WriteMOE/Close - simple one-liner via -exec
        svl_cmd = f"ReadPDB '{receptor.absolute()}'; WriteMOE '{receptor_moe.absolute()}'; Close [];"

        try:
            result = _run_moebatch([moebatch, '-exec', svl_cmd], timeout=60, cwd=out_mdb.parent)

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

    proc = _run_moebatch(cmd, timeout, cwd=out_mdb.parent)

    if log_path:
        log_path.write_text(
            f"$ {' '.join(cmd)}\n\n{proc.stdout}\n{proc.stderr}"
        )

    # The SDF is written before the interaction diagrams are drawn, so a
    # diagram failure doesn't throw away the docking result; the endpoint
    # falls back to drawing the diagram some other way.
    if proc.returncode != 0 and not out_sdf.exists():
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
    parser.add_argument("--scoring", default=DEFAULT_RESCORING,
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
