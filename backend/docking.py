"""Thin wrapper around the smina command line."""

from __future__ import annotations

import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from rdkit import Chem, RDLogger

from . import config
from . import pharmacophore as pharm

RDLogger.DisableLog("rdApp.*")

SCORING_FUNCTIONS = ("vina", "vinardo", "dkoes_scoring", "ad4_scoring")
# MOE scoring: London dG, GBVI/WSA dG, Affinity dG, ASE, Alpha HB
MOE_SCORING_FUNCTIONS = ("London dG", "GBVI/WSA dG", "Affinity dG", "ASE", "Alpha HB")
ENGINES = ("smina", "gnina", "moe")


class SminaNotFound(RuntimeError):
    pass


class GninaNotFound(RuntimeError):
    pass


class MoeBatchNotFound(RuntimeError):
    pass


class SminaFailed(RuntimeError):
    pass


class MoeDockingFailed(RuntimeError):
    pass


# Raised for a nonzero exit from either engine; kept under the historical
# name since callers/tests already match on it.
DockingFailed = SminaFailed


@dataclass
class Pose:
    rank: int
    affinity: float | None
    rmsd_lb: float | None
    rmsd_ub: float | None
    # Populated only when dock()/parse_poses() is given constraints; otherwise
    # left at these defaults and the pose order is untouched.
    constraint_satisfied: bool | None = None
    constraint_detail: list[dict] = field(default_factory=list)
    # 0-based position this pose was read from in its source SDF, kept
    # through any re-sort so minimize()'s reorder pass (see
    # _rewrite_in_rank_order) can find its way back to the original block.
    orig_index: int = -1


def _find_engine(engine: str) -> str:
    if engine == "gnina":
        finder = config.find_gnina
        exc = GninaNotFound
    elif engine == "moe":
        finder = config.find_moebatch
        exc = MoeBatchNotFound
    else:
        finder = config.find_smina
        exc = SminaNotFound
    exe = finder()
    if not exe:
        raise exc(f"{engine} executable not found; run ./setup.sh")
    return exe


def engine_version(engine: str = "smina") -> str:
    exe = _find_engine(engine)
    out = subprocess.run(
        [exe, "--version"], capture_output=True, text=True, timeout=60, env=config.subprocess_env(exe)
    )
    return (out.stdout or out.stderr).strip().splitlines()[0] if (out.stdout or out.stderr) else engine


# Kept for callers that only ever cared about smina.
def smina_version() -> str:
    return engine_version("smina")


def build_command(
    receptor: Path,
    ligand: Path,
    out_sdf: Path,
    center: tuple[float, float, float],
    size: tuple[float, float, float],
    exhaustiveness: int = 8,
    num_modes: int = config.DEFAULT_NUM_MODES,
    energy_range: float = 4.0,
    seed: int = 42,
    scoring: str = "vina",
    cpu: int = 1,
    engine: str = "smina",
    use_gpu: bool = True,
) -> list[str]:
    if engine not in ENGINES:
        raise ValueError(f"unknown docking engine: {engine}")
    exe = _find_engine(engine)
    if scoring not in SCORING_FUNCTIONS:
        raise ValueError(f"unknown scoring function: {scoring}")
    cmd = [
        exe,
        "--receptor", str(receptor),
        "--ligand", str(ligand),
        "--out", str(out_sdf),
        "--center_x", f"{center[0]:.3f}",
        "--center_y", f"{center[1]:.3f}",
        "--center_z", f"{center[2]:.3f}",
        "--size_x", f"{size[0]:.3f}",
        "--size_y", f"{size[1]:.3f}",
        "--size_z", f"{size[2]:.3f}",
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--energy_range", str(energy_range),
        "--seed", str(seed),
        "--scoring", scoring,
        "--cpu", str(cpu),
        "--quiet",
    ]
    # smina has no GPU path at all; gnina uses one by default when a CUDA
    # device is present unless told not to.
    if engine == "gnina" and not use_gpu:
        cmd.append("--no_gpu")
    return cmd


def build_minimize_command(
    receptor: Path,
    ligand: Path,
    out_sdf: Path,
    center: tuple[float, float, float],
    size: tuple[float, float, float],
    scoring: str = "vina",
    cpu: int = 1,
    engine: str = "smina",
    use_gpu: bool = True,
) -> list[str]:
    """`--minimize` mode: no search, just a local energy minimization of
    whatever pose(s) are already in `ligand` against the receptor (each
    conformer in a multi-pose SDF is relaxed independently, same order and
    count out as in). Meant as a post-docking relax pass on top of the
    minimization vina/smina already applies per pose during the search
    itself, e.g. when a low exhaustiveness/energy_range left clashes behind.
    """
    if engine not in ENGINES:
        raise ValueError(f"unknown docking engine: {engine}")
    exe = _find_engine(engine)
    if scoring not in SCORING_FUNCTIONS:
        raise ValueError(f"unknown scoring function: {scoring}")
    cmd = [
        exe,
        "--receptor", str(receptor),
        "--ligand", str(ligand),
        "--out", str(out_sdf),
        "--center_x", f"{center[0]:.3f}",
        "--center_y", f"{center[1]:.3f}",
        "--center_z", f"{center[2]:.3f}",
        "--size_x", f"{size[0]:.3f}",
        "--size_y", f"{size[1]:.3f}",
        "--size_z", f"{size[2]:.3f}",
        "--scoring", scoring,
        "--cpu", str(cpu),
        "--minimize",
        "--quiet",
    ]
    if engine == "gnina" and not use_gpu:
        cmd.append("--no_gpu")
    return cmd


def parse_poses(sdf_path: Path, constraints: list[pharm.Constraint] | None = None) -> list[Pose]:
    """Read the scores smina writes as SDF properties on each pose.

    When `constraints` is given, each pose is also checked geometrically
    against every constraint (see backend.pharmacophore) and the list is then
    reordered: poses meeting every *required* constraint sort first, ties
    broken by affinity plus a penalty for how far the soft/failed ones miss
    by. Ranks are renumbered to match the new order.
    """
    poses: list[Pose] = []
    if not sdf_path.exists():
        return poses
    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=False)
    penalties: list[float] = []
    for i, mol in enumerate(supplier, start=1):
        if mol is None:
            continue

        def prop(name: str) -> float | None:
            if not mol.HasProp(name):
                return None
            try:
                return float(mol.GetProp(name))
            except ValueError:
                return None

        pose = Pose(
            rank=i,
            affinity=prop("minimizedAffinity"),
            rmsd_lb=prop("minimizedRMSD"),
            rmsd_ub=prop("RMSD"),
            orig_index=i - 1,
        )
        pen = 0.0
        if constraints:
            # The supplier skips sanitization (smina's SDF output can trip
            # RDKit's valence checks); the donor/acceptor SMARTS need ring
            # and valence info to match, so sanitize defensively per-pose and
            # fall back to "unknown" constraint state if that fails.
            try:
                Chem.SanitizeMol(mol)
                results = pharm.evaluate(mol, constraints)
                pose.constraint_detail = [r.as_dict() for r in results]
                pose.constraint_satisfied = pharm.all_required_met(results)
                pen = pharm.penalty(results, constraints)
            except Exception:
                pose.constraint_satisfied = None
        penalties.append(pen)
        poses.append(pose)

    if constraints:
        order = sorted(
            range(len(poses)),
            key=lambda idx: (
                not (poses[idx].constraint_satisfied in (True, None)),
                (poses[idx].affinity or 0.0) + penalties[idx],
            ),
        )
        poses = [poses[idx] for idx in order]
        for rank, pose in enumerate(poses, start=1):
            pose.rank = rank
    return poses


def _rewrite_in_rank_order(sdf_path: Path, poses: list[Pose]) -> None:
    """Physically reorders the SDF on disk to match `poses`' current order.

    Every other place a pose is looked up by position (the poses.sdf
    download, the pose bar, interactions2d) is purely positional against
    this file - a Pose list re-sorted only in memory (parse_poses' own
    constraint-first sort, or minimize()'s affinity re-sort below) would
    otherwise silently point "pose 1" at the wrong block.
    """
    blocks = sdf_path.read_text().split("$$$$\n")
    if blocks and blocks[-1].strip() == "":
        blocks = blocks[:-1]  # trailing empty chunk after the last delimiter
    sdf_path.write_text("".join(blocks[p.orig_index] + "$$$$\n" for p in poses))


def dock(
    receptor: Path,
    ligand_sdf: Path,
    out_sdf: Path,
    center: tuple[float, float, float],
    size: tuple[float, float, float],
    log_path: Path | None = None,
    constraints: list[pharm.Constraint] | None = None,
    engine: str = "smina",
    **kwargs,
) -> list[Pose]:
    """Dock one ligand and return its poses, best first (or best-constrained-first, see parse_poses)."""
    # Route to MOE if specified
    if engine == "moe":
        return dock_with_moe(
            receptor, ligand_sdf, out_sdf, center, size, log_path=log_path, **kwargs
        )

    # Smina/Gnina workflow
    cmd = build_command(receptor, ligand_sdf, out_sdf, center, size, engine=engine, **kwargs)
    proc = subprocess.run(
        cmd, capture_output=True, text=True, timeout=config.SMINA_TIMEOUT, env=config.subprocess_env(cmd[0])
    )
    if log_path is not None:
        log_path.write_text(
            "$ " + " ".join(cmd) + "\n\n" + (proc.stdout or "") + "\n" + (proc.stderr or "")
        )
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()
        raise SminaFailed(tail[-1] if tail else f"smina exited {proc.returncode}")
    poses = parse_poses(out_sdf, constraints=constraints)
    if not poses:
        raise SminaFailed("smina produced no poses (ligand may not fit the box)")
    return poses


def minimize(
    receptor: Path,
    poses_sdf: Path,
    out_sdf: Path,
    center: tuple[float, float, float],
    size: tuple[float, float, float],
    log_path: Path | None = None,
    constraints: list[pharm.Constraint] | None = None,
    **kwargs,
) -> list[Pose]:
    """Relax an already-docked poses.sdf in place (see build_minimize_command).

    Each pose is minimized independently, so the affinity order coming out
    isn't guaranteed to still be best-first the way dock()'s own output is -
    re-sorted here by affinity (constraint-first-then-affinity when
    `constraints` is given, same as parse_poses; plain affinity otherwise) so
    pose 1 stays "the best pose" for callers same as after dock().
    """
    cmd = build_minimize_command(receptor, poses_sdf, out_sdf, center, size, **kwargs)
    proc = subprocess.run(
        cmd, capture_output=True, text=True, timeout=config.SMINA_TIMEOUT, env=config.subprocess_env(cmd[0])
    )
    if log_path is not None:
        # Appended, not overwritten - this normally runs right after dock()
        # against the same log_path, and the docking search's own log is
        # worth keeping alongside the minimize pass that followed it.
        with log_path.open("a") as f:
            f.write("\n\n--- minimize ---\n\n$ " + " ".join(cmd) + "\n\n" + (proc.stdout or "") + "\n" + (proc.stderr or ""))
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()
        raise SminaFailed(tail[-1] if tail else f"smina exited {proc.returncode}")
    poses = parse_poses(out_sdf, constraints=constraints)
    if not poses:
        raise SminaFailed("smina produced no poses after minimization")
    if not constraints:
        poses.sort(key=lambda p: p.affinity if p.affinity is not None else 0.0)
        for rank, pose in enumerate(poses, start=1):
            pose.rank = rank
    _rewrite_in_rank_order(out_sdf, poses)
    return poses


def dock_with_moe(
    receptor: Path,
    ligand_sdf: Path,
    out_sdf: Path,
    center: tuple[float, float, float],
    size: tuple[float, float, float],
    log_path: Path | None = None,
    num_modes: int = config.DEFAULT_NUM_MODES,
    scoring: str = "London dG",
    placement: str = "Triangle Matcher",
    refinement: str = "Rigid Receptor",
    timeout: int = 3600,
    constraints: list[pharm.Constraint] | None = None,
    **kwargs,
) -> list[Pose]:
    """MOE docking workflow: calls moebatch SVL script that exports to SDF directly."""
    from . import moe_dock

    moebatch = config.find_moebatch()
    if not moebatch:
        raise MoeBatchNotFound("moebatch not found - install MOE")

    # MOE creates temp .mdb, SVL script exports to SDF and deletes MDB
    out_mdb = out_sdf.with_suffix(".mdb")

    # Build covalent reaction file path if enabled
    covalent_reaction_file = None
    if kwargs.get('covalent_enabled'):
        # Custom .rxn file takes precedence over built-in template
        if kwargs.get('custom_rxn_id'):
            custom_path = config.RUNS_DIR / "custom_reactions" / f"{kwargs['custom_rxn_id']}.rxn"
            if custom_path.exists():
                covalent_reaction_file = custom_path
        elif kwargs.get('covalent_reaction'):
            rxn_dir = config.find_moe_reactions()
            if rxn_dir:
                covalent_reaction_file = rxn_dir / f"{kwargs['covalent_reaction']}.rxn"

    moe_dock.dock(
        receptor=receptor,
        ligand=ligand_sdf,
        out_mdb=out_mdb,
        out_sdf=out_sdf,
        center=center,
        size=size,
        num_poses=num_modes,
        scoring=scoring,
        placement=placement,
        refinement=refinement,
        log_path=log_path,
        timeout=timeout,
        auto_detect_site=True,  # Use HETATM center with 10Å radius if found
        ph4_file=kwargs.get('ph4_file'),  # Allow direct .ph4 file upload
        covalent_reaction_file=covalent_reaction_file,
        covalent_residue=kwargs.get('covalent_residue'),
    )

    # Parse SDF with MOE scores (written as SD properties by SVL script)
    poses = _parse_moe_sdf(out_sdf, constraints=constraints)
    if not poses:
        raise MoeDockingFailed("MOE produced no poses")
    return poses


def _parse_moe_sdf(sdf_path: Path, constraints: list[pharm.Constraint] | None = None) -> list[Pose]:
    """Parse MOE-generated SDF with scores as SD properties."""
    poses: list[Pose] = []
    if not sdf_path.exists():
        return poses

    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=False)
    penalties: list[float] = []

    for i, mol in enumerate(supplier, start=1):
        if mol is None:
            continue

        def prop(name: str) -> float | None:
            if not mol.HasProp(name):
                return None
            try:
                return float(mol.GetProp(name))
            except ValueError:
                return None

        # MOE writes primary score as minimizedAffinity (matches smina format)
        affinity = prop("minimizedAffinity")
        rmsd = prop("minimizedRMSD")

        pose = Pose(
            rank=i,
            affinity=affinity,
            rmsd_lb=rmsd,
            rmsd_ub=rmsd,
            orig_index=i - 1,
        )

        pen = 0.0
        if constraints:
            try:
                Chem.SanitizeMol(mol)
                results = pharm.evaluate(mol, constraints)
                pose.constraint_detail = [r.as_dict() for r in results]
                pose.constraint_satisfied = pharm.all_required_met(results)
                pen = pharm.penalty(results, constraints)
            except Exception:
                pose.constraint_satisfied = None

        penalties.append(pen)
        poses.append(pose)

    # Re-rank by constraints if present
    if constraints:
        order = sorted(
            range(len(poses)),
            key=lambda idx: (
                not (poses[idx].constraint_satisfied in (True, None)),
                (poses[idx].affinity or 0.0) + penalties[idx],
            ),
        )
        poses = [poses[idx] for idx in order]
        for rank, pose in enumerate(poses, start=1):
            pose.rank = rank

    return poses
