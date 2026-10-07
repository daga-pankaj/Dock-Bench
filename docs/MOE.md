# MOE integration

How Docking Bench drives MOE, and what we learned getting the batch results to
match the MOE GUI. The user-facing feature list is in the top-level README.
This file is for working on the code.

MOE is optional. `config.find_moebatch()` looks for `moebatch` under
`/Applications/moe*/bin/`, then on `PATH`. Set `MOEBATCH_BIN` to override.
Without it the MOE engine, QuickPrep and the MOE diagrams are unavailable,
and everything else works.

## How a MOE docking run works

All MOE work runs headless through `moebatch -load <script>.svl -exec "main []"`.
The code is in `backend/moe_dock.py`.

1. RDKit builds the ligand 3D structure from SMILES, or the uploaded molblock
   is used as is (same as smina/gnina). **The ligand is not protonated or
   washed**, so write charged groups explicitly (e.g. `[NH3+]`). This changes
   both the poses and the interactions a lot (see DPP-4 / sitagliptin
   below).
2. The job's cleaned receptor PDB is converted to `.moe` with
   `ReadPDB`/`WriteMOE`. Each ligand gets its own copy
   (`<ligand>_poses_receptor.moe`) because ligands share the receptor and
   each run deletes its copy when it finishes.
3. `build_svl_script()` fills in `backend/moe_dock_template.svl`:
   - `Add_H` runs on the receptor before any atom indices are taken (see
     "Covalent docking" for why).
   - The site is every receptor atom inside the sphere around the docking box.
   - Options follow MOE's Dock panel defaults. Placement poses are scored
     with London dG and the top 30 (`PLACEMENT_POSES`) are carried into
     refinement. Refined poses are rescored with the user's function
     (default GBVI/WSA dG) and the best `num_poses` are kept.
   - A `.ph4` file, if one was uploaded, is read with `ph4_aType` and passed
     as `opt.ph4`.
   - `Dock` writes an `.mdb`, and `db_ExportSD` exports it to the poses SDF
     with MOE's own fields (`S`, `E_place`, `E_score1`, `E_refine`,
     `E_score2`, `E_conf`, `rmsd_refine`; see `docking.MOE_SCORE_FIELDS`).
   - Each pose's Ligand Interactions diagram is drawn next
     (`LID_FUNCTION`/`LID_RENDER`), while the receptor is still loaded.
4. `docking._parse_moe_sdf()` reads the SDF. `S` becomes the pose score, and
   every MOE field is kept per pose for the results table.

If the run fails after the SDF has been written (i.e. while drawing a
diagram), the docking result still counts. The diagram endpoint then falls
back to drawing the diagram on demand, or to ProLIF.

### Scores are only comparable with the same rescoring function

`S` is whatever rescoring function the job used. London dG and GBVI/WSA dG
values can't be compared with each other. Rescoring also affects pose
quality: on 3PBS covalent docking, London dG rescoring put the top pose
7.6 Å from the crystal pose, while GBVI/WSA dG got 1.1–1.4 Å. That's why the
default is GBVI/WSA dG, as in the GUI.

## License queue

The site has 6 MOE license tokens and each `moebatch` run takes 2, so with
unlimited parallelism the 4th ligand onward failed with
`Can't check out license ... Not enough tokens available`.

- Every `moebatch` call goes through `moe_dock._run_moebatch`. A server-wide
  semaphore caps concurrent runs at `MOE_MAX_CONCURRENT` (default 1, which
  leaves tokens for other users).
- MOE ligands run on their own thread pool, sized to the same cap
  (`JobManager._moe_pool`). Ligands still waiting show as "queued", and MOE
  jobs don't hold up smina/gnina workers.
- If the license server still has no tokens (other users hold them), the run
  retries every `MOE_LICENSE_RETRY` seconds for up to `MOE_LICENSE_WAIT`
  (default 3600 s), then gives up. No other kind of failure is retried.
- QuickPrep and on-demand diagrams use the same queue. They wait behind a
  running docking ligand.

## QuickPrep

Panel 1's **MOE QuickPrep** button calls `POST /api/receptor/{id}/quickprep`,
which calls `moe_dock.quickprep()`.

- The batch entry point is the global SVL function `_LIGX_Execute opt` in
  `$MOE/svl/calc.svl/ligx_u.svl`. It is the QuickPrep panel's OK action.
  Called with no options, it uses the user's saved QuickPrep settings
  (`LIGX_SETUP` defaults plus the RC file), as the GUI does. The steps are
  Structure Preparation, Protonate3D, deleting waters more than 4.5 Å from
  receptor and ligand, tethers, fixing atoms more than 8 Å from the ligand,
  and minimizing to RMS gradient 0.1.
- The raw upload goes in, header included, so Structure Preparation can use
  SEQRES to build short missing loops.
- **The flexible region is defined by the ligand.** We put the selected site
  ligand in the `LIGX_LIGAND` collection; otherwise MOE guesses by
  `aClassRLS` and skips known cofactors. With no ligand at all QuickPrep
  fixes every atom, so the minimization does nothing (2479 of 2480 atoms
  fixed in a test).
- `WritePDB` keeps chain IDs, residue numbers, hydrogens and formal charges
  (columns 79–80), and `ReadPDB` reads the charges back. The prepared PDB is
  stored as a new receptor with a `<id>.json` sidecar pointing at its
  source. `clean_receptor(keep_hydrogens=True)` keeps its hydrogens and
  charges at job time. RMSD references and LINK-record lookups still read
  the source upload, because the minimization moves the crystal ligand and
  MOE doesn't write LINK records.
- `prep.quickprep_report()` compares the upload and the prepared PDB:
  hydrogens added, residues built, waters removed, heavy-atom shifts, and
  protonation states that differ from pH 7 defaults (every HIS tautomer,
  ASH/GLH/LYN/ARN, CYM/TYM). Residues within 6 Å of the site ligand come
  first.
- Timing on DPP-4 1x70 (about 6k heavy atoms): 80–116 s on the raw upload
  with waters, about 40 s on a cleaned receptor.

## Ligand Interactions diagram

MOE's 2D diagram (the GUI's Ligand Interactions panel, `$MOE/svl/run/prolig2d.svl`)
can be drawn in batch:

- `DrawLigandInteractions [[lig], [rec], opt]` (wrapper in
  `calc.svl/prolig_u.svl`) returns a graphics object (grob). The legend
  comes from `call ['_DrawLigandInteractionsLegend', [[], []], 'prolig2d.svl']`.
  We stack the two in one grob and write it with
  `gr_fwrite [file, grob, 'image/svg+xml', [0,0,1,1]]`. That takes about 1 s
  and needs no GUI.
- A MOE docking run draws one diagram per pose
  (`poses/<ligand>_poses_lid_<n>.svg`) at no extra license cost. For
  smina/gnina poses, or older MOE jobs, `render_interactions()` draws one on
  demand. That is a short moebatch run in the license queue, and the result
  is cached next to the poses.
- `browser_svg()` replaces MOE's embedded SVG font ("sf"), which browsers
  ignore, with a system sans-serif, and adds a white background because the
  app panels are dark.
- `GET .../interactions2d.svg?style=auto|moe|prolif`. The `X-Diagram-Source`
  header says which renderer was used. If MOE fails, the endpoint falls back
  to ProLIF.

## Covalent docking

### Options must match the GUI's

The first version invented keys (`opt.covalent`, `opt.rxn_file`,
`opt.tether = 'A:CYS:189'`) and MOE ignored them. The GUI-generated
`run.sh` showed the real structure, which `build_svl_script()` now emits:

```svl
opt.dockmode = 'Covalent';
opt.placement = 'Covalent Docking';
opt.placement_opt = [rxnFile: '<file>.rxn', rxsite: 'asel', mode: 'Covalent', ...];
opt.refine_opt = [fixrec: 'Fix', tether: 10, ...];   // tether is a force constant
opt.rxsite = rxsite_atoms(1);                         // the nucleophile atom
```

### Choosing the reactive atom

- Nucleophile per residue: CYS SG, SER OG, THR OG1, TYR OH, LYS NZ
  (`moe_dock._get_nucleophile_atom_name`).
- The filter is `rUID aResidue protein_atoms == <resnum> and aName ... == '<atom>'`.
  `aResidue` returns residue *keys*, not numbers, so comparing it to 436
  never matches. `rUID` is the PDB residue number. The chain is ignored
  because MOE renames chains on PDB import (the chain tag becomes the file
  name).
- **Dock runs `Add_H` internally after `opt.rxsite` has become an atom
  index**, so the index shifts and points at the wrong atom. The template
  runs `Add_H` itself before taking indices. An empirical bisection against
  the GUI's `run.sh` found this; guessing at option changes didn't.
- A GUI rxsite is an index into the receptor sorted by `aNumber`. An earlier
  GUI run targeted Ser349 (index 4702) instead of Ser294 (3830), so check the
  GUI side too when results differ.
- In the UI, the **Pick** button next to "target residue" lets you click the
  nucleophile in the viewer. It fills in `chain:RES:num` and warns about
  non-nucleophilic residues or the wrong atom.

### Reactions

- There are 22 built-in templates in `$MOE/lib/reactions/covalent/*.rxn`,
  found by `config.list_covalent_reactions()`. Common ones:

  | Target    | Template                                                                   |
  |-----------|----------------------------------------------------------------------------|
  | Cys       | `michael_acceptor_s`, `epoxide_s`, `alpha_halocarbonyl_s`, `vinyl_sulfone_s` |
  | Ser/Thr   | `michael_acceptor_o`, `epoxide_o`, `carbonyl_o`                            |
  | Lys       | `sulfonyl_halide`, `carbamate`                                             |
  | Tyr       | `sulfonyl_halide_tyr`                                                      |
  | Cys–Cys   | `thiols`                                                                   |

- Custom reactions: draw one in Ketcher, export it as `.rxn`, and upload it
  (`/api/reactions/upload`, stored in `runs/custom_reactions/`).

### Check the ligand's reactant form first

When a covalent run differs from the GUI, the most common cause is ligand
chemistry, not the script. For 3PBS (AZR on Ser294), the PDB CCD SMILES is
the cut adduct, so the reacting carbon is an aldehyde. The ester-forming
reaction needs a carboxylic acid, so MOE acylated a side-chain carboxylate
instead.

- `covalent_check.reactant_warning()` warns when a ligand matches the
  reaction's ligand reactant zero times or more than once.
- `covalent_check.acid_form_from_link()`: when a bound ligand is LINKed to
  the protein as an acyl adduct, "self-dock" rebuilds it as the acid it
  reacts from.
- GUI-side trap: a `lig.moe` without stereo makes Dock enumerate every
  stereoisomer ("chiral 8"). An SDF input keeps the stereo.

## SVL pitfalls

- `rAtom` doesn't exist. Built-ins (`aName`, `aChain`, `rUID`, ...) are
  used without a `function` declaration. Only declare functions from loaded
  modules, e.g. `function Dock;`.
- `oGetCollection 'X'` returns the atoms themselves, so
  `oDestroy oGetCollection 'X'` deletes those atoms. To clear the
  collection, use `oSetCollection ['X', []]`.
- Functions in `$MOE/svl/run/*.svl` (e.g. `_DrawLigandInteractionsLegend`)
  aren't loaded in batch. Call them with `call [fn, arg, 'file.svl']`.
- The templates go through Python `str.format`, so literal SVL braces must
  be doubled (`{{}}`). Values substituted in aren't re-parsed.
- Dock writes scratch `dock_lig_*_bound_*.mdb` files into the working
  directory, so moebatch runs from the job's folder.

## Debugging approach

When batch results differ from the GUI:

1. Get a GUI baseline and keep the GUI-generated `run.sh`.
2. Rerun that `run.sh` in a copy and swap our inputs in one at a time
   (receptor, ligand, site, options) to find which one makes the
   difference.
3. Read MOE's own sources (`$MOE/svl/calc.svl/dock_u.svl`, `ligx_u.svl`,
   `run/prolig2d.svl`). Load an instrumented copy with `moebatch -load` to
   print internal state.

Confirm the root cause before changing the script. Guessed option changes
wasted several rounds early on.

## Example: DPP-4 / sitagliptin (1x70, ligand 715)

- With the SMILES as given (neutral amine), the MOE diagram shows almost
  nothing apart from one fluorine contact to Arg358. With `[NH3+]` it shows
  the Glu206 H-bond and the Glu205 ionic contact, as in the GUI.
- QuickPrep with ligand 715 as the site: His740 (catalytic triad) comes out
  HIP, Asp663 neutral (ASH). About 70 heavy atoms move more than 0.5 Å.

## Ideas not done yet

- Ligand protonation/wash step for MOE jobs (pH 7)
- Auto-detect nucleophilic residues near the site; suggest a reaction
  template from the ligand's warhead
- Multi-residue attachment (crosslinkers)
