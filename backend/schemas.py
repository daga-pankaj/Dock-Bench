"""Request models shared by the HTTP API (backend.main) and the Streamlit UI."""

from __future__ import annotations

from pydantic import BaseModel, Field

from . import config


class FetchRequest(BaseModel):
    pdb_id: str


class QuickPrepRequest(BaseModel):
    # het_ligands key (resname:chain:resseq) of the ligand defining the site;
    # QuickPrep keeps receptor atoms near it flexible during minimization.
    site_ligand_key: str | None = None


class LigandInput(BaseModel):
    name: str
    smiles: str | None = None
    molblock: str | None = None


class ParseRequest(BaseModel):
    text: str


class ConstraintIn(BaseModel):
    """One pharmacophoric constraint: an H-bond donor/acceptor, an aromatic
    ring (pi-stacking proxy), a metal-coordinating atom, any heavy atom
    ("positional"), or any atom matching a user SMARTS ("smarts", needs the
    `smarts` field) - checked for proximity to `radius` of `point` (a
    receptor-atom or hand-picked coordinate). Checked against each pose after
    docking, not during smina/gnina's own search - see
    backend/pharmacophore.py."""

    kind: str
    point: tuple[float, float, float]
    radius: float = Field(3.5, ge=0.5, le=15.0)
    required: bool = True
    weight: float = Field(1.0, ge=0.0, le=20.0)
    label: str = ""
    smarts: str | None = None  # only used (and required) when kind == "smarts"


class JobRequest(BaseModel):
    receptor_id: str
    center: tuple[float, float, float]
    size: tuple[float, float, float]
    ligands: list[LigandInput]
    exhaustiveness: int = Field(8, ge=1, le=64)
    num_modes: int = Field(config.DEFAULT_NUM_MODES, ge=1, le=20)
    energy_range: float = Field(4.0, ge=0.5, le=20.0)
    seed: int = 42
    scoring: str = "vina"
    engine: str = "smina"
    use_gpu: bool = True
    minimize_after_docking: bool = False
    keep_waters: bool = False
    remove_ligand: bool = True
    exclude_chains: list[str] = Field(default_factory=list)
    trim_distant_chains: bool = False
    chain_cutoff: float = Field(15.0, ge=5.0, le=60.0)
    constraints: list[ConstraintIn] = Field(default_factory=list)
    # Self-docking validation: the key (from het_ligands, resname:chain:resseq)
    # of a bound ligand to compare each docked pose against by RMSD (see
    # backend/rmsd.py). Read from this job's receptor before it's cleaned, so
    # it still works with remove_ligand=True stripping that same ligand out
    # of the receptor actually used for docking.
    reference_ligand_key: str | None = None
    # MOE-specific parameters
    moe_placement: str = "Triangle Matcher"
    moe_refinement: str = "Rigid Receptor"
    ph4_id: str | None = None  # MOE pharmacophore file ID from /api/pharmacophore/upload
    # Covalent docking (MOE-only)
    covalent_enabled: bool = False
    covalent_reaction: str | None = None  # e.g., "michael_acceptor_s" (built-in template)
    covalent_residue: str | None = None   # e.g., "A:CYS:189"
    custom_rxn_id: str | None = None      # Custom .rxn file ID from /api/reactions/upload
