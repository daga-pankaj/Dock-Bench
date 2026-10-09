---
marp: true
theme: default
paginate: true
---

# MOE Docking Integration
## Extending the Docking Bench Platform

**Multi-Engine Docking Platform**
Smina · Gnina · MOE

---

## Overview

**Goal:** Integrate MOE (Molecular Operating Environment) docking engine into existing smina/gnina platform

**Why MOE?**
- Commercial-grade docking with advanced placement algorithms
- Induced-fit refinement options
- Native pharmacophore constraint support
- Alternative scoring functions (London dG, GBVI/WSA dG, Affinity dG)

**Result:** Unified interface for 3 docking engines with engine-specific optimizations

---

## Architecture

```
Frontend (React-like vanilla JS)
    ↓
FastAPI Backend
    ├── smina/gnina wrapper (existing)
    └── MOE wrapper (new)
         ├── moe_dock.py       ← Python orchestration
         ├── moe_dock_template.svl ← SVL script template
         └── moebatch CLI       ← MOE batch processor
```

**Key decision:** Template-based SVL generation instead of direct Python-MOE API
- Cross-platform compatibility
- No MOE Python bindings required
- Full control over docking parameters

---

## MOE Batch Execution Flow

1. **PDB → MOE conversion** (inline SVL via `-exec`)
   ```bash
   moebatch -exec "ReadPDB 'file.pdb'; WriteMOE 'file.moe'; Close [];"
   ```

2. **Template instantiation** (`build_svl_script()`)
   - Replace placeholders: `{RECEPTOR_FILE}`, `{CENTER_X}`, etc.
   - Generate per-job SVL script

3. **Docking execution**
   ```bash
   moebatch -load script.svl -exec "main []"
   ```

4. **Export to SDF** (unified format for all engines)

---

## SVL Script Template

```svl
#svl
function Open; function Close; function Atoms; 
function Dock; function db_ExportSD; function ph4_aType;

global function main []
    local site_atoms, rec, lig, site, opt, ph4_data;
    
    Open '{RECEPTOR_FILE}';
    all_atoms = Atoms [];
    protein_atoms = all_atoms | aResidue all_atoms <> '';
    
    distances = sqrt add pow [sub [aPos protein_atoms, 
                [{CENTER_X}, {CENTER_Y}, {CENTER_Z}]], 2];
    site_atoms = protein_atoms | distances <= {SITE_RADIUS};
    
    rec = [type: 'MOE', atom_sel: protein_atoms];
    lig = [type: 'SD File', file: '{LIGAND_FILE}'];
    
    // Pharmacophore support
    ph4_data = [];
    if '{PH4_FILE}' <> '' then
        ph4_data = ph4_aType '{PH4_FILE}';
    endif;
    
    opt = [run_mode: 'run', placement: '{PLACEMENT}', 
           rescoring: '{SCORING}', refine: '{REFINEMENT}', 
           maxpose: {NUM_POSES}, ph4: ph4_data];
    
    Dock [rec, site, lig, '{OUTPUT_MDB}', opt];
    db_ExportSD ['{OUTPUT_MDB}', '{OUTPUT_SDF}'];
    Close [];
endfunction
```

---

## Auto-Detection Features

### Binding Site Auto-Detection
**Problem:** User shouldn't manually define binding site for co-crystallized structures

**Solution:** Extract HETATM ligands from PDB, use center + 10Å radius
```python
if auto_detect_site and receptor.suffix.lower() == '.pdb':
    het_ligands = prep.list_het_ligands(pdb_text, min_atoms=4)
    site_ligand = next((lig for lig in het_ligands if not lig.is_ion), None)
    if site_ligand:
        center = tuple(site_ligand.center)
        size = (20.0, 20.0, 20.0)  # 10Å radius = 20Å box
```

**Benefit:** One-click docking for PDB structures with bound ligands

---

## Binding Site Definition

**MOE requires atom selection, not just box coordinates**

Initial approach (failed):
```svl
site = [type: 'box', center: [x,y,z], size: [sx,sy,sz]]
// ❌ Poses scattered everywhere
```

Working approach:
```svl
distances = sqrt add pow [sub [aPos protein_atoms, [x,y,z]], 2];
site_atoms = protein_atoms | distances <= radius;
site = site_atoms;  // actual atom selection
// ✅ Poses confined to binding pocket
```

**Learning:** MOE site definition semantics differ from smina's geometric box

---

## Pharmacophore Constraints

### Smina/Gnina Approach
- Manual point-based constraints (H-bond donor/acceptor, positional)
- Post-docking filtering and reranking
- No search bias (constraints don't guide docking)

### MOE Approach
- Native .ph4 file support
- Upload pharmacophore query file
- Constraints guide placement algorithm
- Only available with "Pharmacophore" placement method

**Implementation:**
- Added `/api/pharmacophore/upload` endpoint
- Store .ph4 files in `runs/pharmacophores/{ph4_id}.ph4`
- Pass to SVL via `ph4_aType` function
- UI shows .ph4 upload only when MOE engine selected

---

## Parameter Mapping

| Concept | Smina/Gnina | MOE |
|---------|-------------|-----|
| **Search effort** | exhaustiveness (1-64) | placement method |
| **Refinement** | --minimize flag | refine: Rigid/Induced Fit/None |
| **Scoring** | vina, vinardo, ad4 | London dG, GBVI/WSA dG, Affinity dG |
| **Constraints** | post-filter points | native .ph4 pharmacophore |
| **Poses kept** | num_modes | maxpose |
| **Energy range** | energy_range | (implicit in rescoring) |

**Challenge:** Different paradigms require conditional UI/backend logic

---

## Technical Challenges Solved

### 1. SVL Syntax Errors
**Problem:** "Expressions only allowed in functions" with inline `-exec`

**Solution:** Use script file + `-load` then `-exec "main []"` pattern
- Forward declare functions at top
- Wrap all logic in `global function main []`

### 2. PDB Conversion Failures
**Problem:** Empty error messages from moebatch

**Solution:** Check both stderr AND stdout, add fallback message
```python
stderr = (e.stderr or e.stdout or "").strip()
raise MoeDockingFailed(f"moebatch error: {stderr or 'unknown error'}")
```

---

## Technical Challenges Solved (cont.)

### 3. Site Selection Expression
**Problem:** `'protein'` keyword not recognized in MOE

**Solution:** Manual residue-based selection
```svl
protein_atoms = all_atoms | aResidue all_atoms <> '';
```

### 4. Pharmacophore API Discovery
**Problem:** `ph4_Read` function doesn't exist

**Solution:** Use `ph4_aType` for atom-type pharmacophore queries
```svl
function ph4_aType;  // declare first
ph4_data = ph4_aType '{PH4_FILE}';  // load .ph4 file
```

### 5. Config Path Inconsistency
**Problem:** `config.RUNS` vs `config.RUNS_DIR` attribute errors

**Solution:** Standardized on `config.RUNS_DIR` throughout codebase

---

## Job Naming Convention

**Old system:** Random UUID hex (`8146a68fa69e`)

**New system:** `{engine}_{yyyymmddhhmmss}`
- `moe_20261001143530`
- `smina_20261001144215`
- `gnina_20261001150033`

**Benefits:**
- Instantly identify docking engine
- Chronological sorting
- Human-readable job directories
- Easier debugging and log review

---

## Frontend Integration

### Engine Selection
```javascript
$("engineSel").onchange = () => {
  const isMoe = $("engineSel").value === "moe";
  
  // Show/hide engine-specific params
  $("moeParamsRow").style.display = isMoe ? "grid" : "none";
  $("sminaGninaParamsRow").style.display = isMoe ? "none" : "flex";
  
  // Switch constraint UI
  $("sminaGninaConstraints").style.display = isMoe ? "none" : "block";
  $("moeConstraints").style.display = isMoe ? "block" : "none";
  
  updateScoringOptions(engine);
};
```

**Result:** Context-aware UI that shows only relevant options per engine

---

## MOE-Specific UI Elements

1. **Placement method dropdown**
   - Triangle Matcher (default)
   - Alpha PMI
   - Pharmacophore (requires .ph4 file)
   - Proxy Triangle
   - Alpha Triangle

2. **Refinement dropdown**
   - Rigid Receptor (default)
   - Induced Fit
   - Rigid Body
   - None

3. **Pharmacophore upload** (in constraint panel)
   - File validation: checks for `#moe:ph4que` or `#pharmacophore` header
   - Returns `ph4_id` for backend reference

---

## User Workflow

### Standard MOE Docking
1. Fetch/upload PDB receptor (auto-detects binding site from HETATM)
2. Upload ligand SMILES or SDF
3. Select **MOE** engine
4. Choose placement method (default: Triangle Matcher)
5. Choose refinement (default: Rigid Receptor)
6. Select scoring function
7. Run

### Pharmacophore-Constrained Docking
1-2. Same as above
3. Select **MOE** engine
4. **Upload .ph4 file** in pharmacophore panel
5. Set placement to **Pharmacophore**
6. Configure refinement/scoring
7. Run

---

<style scoped>section { font-size: 22px; } pre { font-size: 0.8em; }</style>

## What's New

Three features added since the initial MOE integration:

| Feature | What it does |
|---------|--------------|
| **Covalent docking** | Dock a warhead onto a chosen Cys/Ser/Thr/Tyr/Lys with MOE reaction templates |
| **2D protein–ligand interaction diagram** | MOE's Ligand Interactions plot per pose, with a ProLIF fallback |
| **MOE license queue** | Server-wide cap on `moebatch` runs; ligands wait as "queued" instead of failing |

Also: **MOE QuickPrep** for receptor preparation, per-pose MOE score tables, pose stepping, 2D ligand overlay on the 3D view

---

# Covalent Docking

---

<style scoped>section { font-size: 22px; } pre { font-size: 0.8em; }</style>

## Covalent Docking: Overview

**Inputs (MOE engine only):**
- Target residue as `chain:RES:num` — typed, or **Pick** an atom in the 3D viewer
- Reaction: one of 22 built-in templates (`$MOE/lib/reactions/covalent/*.rxn`) or an uploaded `.rxn` (draw in Ketcher, export)

**Nucleophile per residue:** CYS SG · SER OG · THR OG1 · TYR OH · LYS NZ

| Target | Common templates |
|--------|------------------|
| Cys | `michael_acceptor_s`, `epoxide_s`, `alpha_halocarbonyl_s`, `vinyl_sulfone_s` |
| Ser/Thr | `michael_acceptor_o`, `epoxide_o`, `carbonyl_o` |
| Lys | `sulfonyl_halide`, `carbamate` |
| Tyr | `sulfonyl_halide_tyr` |

**Results:** each pose's MOE scores plus the ligand group it is bonded through

---

<style scoped>section { font-size: 22px; } pre { font-size: 0.8em; }</style>

## Covalent Docking: Matching the GUI's Options

**First attempt:** made-up keys that MOE silently ignored
```svl
opt.covalent = 1; opt.rxn_file = '...'; opt.tether = 'A:CYS:189';   // ❌ no effect
```

**Fix:** copy the structure from the GUI-generated `run.sh`
```svl
opt.dockmode      = 'Covalent';
opt.placement     = 'Covalent Docking';
opt.placement_opt = [rxnFile: '<file>.rxn', rxsite: 'asel', mode: 'Covalent', ...];
opt.refine_opt    = [fixrec: 'Fix', tether: 10, ...];   // tether = force constant
opt.rxsite        = rxsite_atoms(1);                     // the nucleophile atom
```

**Selecting the reactive atom:**
```svl
rxsite_atoms = protein_atoms | rUID aResidue protein_atoms == 294
                             and aName protein_atoms == 'OG';
```
- `aResidue` returns residue *keys*; `rUID` gives the PDB residue number
- Chain is ignored: MOE renames chains on PDB import

---

## Covalent Docking: The Reactive-Site Bug

**Symptom:** `Reaction did not match protein` on unprepared receptors

**Root cause:** `Dock` runs `Add_H` internally **after** `opt.rxsite` is resolved to an atom index
→ added hydrogens shift the indices → rxsite points at an unrelated atom

**Fix:** add hydrogens up front, before taking any index
```svl
Open '{RECEPTOR_FILE}';
Add_H Atoms [];          // keeps opt.rxsite stable
all_atoms = Atoms [];
```

**How it was found:** reran the GUI's `run.sh`, swapped in our inputs one at a time, and printed state from an instrumented copy of `dock_u.svl`. Guessing at option changes didn't find it.

---

<style scoped>section { font-size: 22px; } pre { font-size: 0.8em; }</style>

## Covalent Docking: Ligand Chemistry Matters

**Case study: 3PBS (ligand AZR on Ser294)**
- The PDB CCD SMILES is the *bound adduct*, so the reacting carbon is an aldehyde
- The ester-forming reaction needs a carboxylic acid, so MOE acylated a side-chain carboxylate instead

**Safeguards added (`covalent_check.py`):**
- `reactant_warning()` warns when the ligand matches the reaction's reactant **0 or >1 times**
- `acid_form_from_link()` rebuilds a LINKed acyl adduct as its acid form for self-docking
- Detects which ligand group each pose is attached through

**Rescoring changes pose quality:**

| Rescoring | Top-pose RMSD to crystal (3PBS) |
|-----------|-------------------------------|
| London dG | 7.6 Å |
| GBVI/WSA dG (default) | 1.1–1.4 Å |

---

# 2D Protein–Ligand Interaction Diagram

---

<style scoped>section { font-size: 22px; } pre { font-size: 0.8em; }</style>

## Ligand Interactions Diagram

**MOE's 2D diagram (GUI: Ligand Interactions panel) drawn headless in `moebatch`**

```svl
grob   = DrawLigandInteractions [[lig], [rec], opt];
legend = call ['_DrawLigandInteractionsLegend', [[], []], 'prolig2d.svl'];
gr_fwrite [file, <grob + legend stacked>, 'image/svg+xml', [0,0,1,1]];
```
- About 1 s per diagram, no GUI needed
- `$MOE/svl/run/*.svl` functions aren't loaded in batch → call them with `call [...]`

**When diagrams are drawn:**
- **MOE jobs:** one per pose, inside the docking run, while the receptor is still loaded
  → `poses/<ligand>_poses_lid_<n>.svg`, **no extra license checkout**
- **smina/gnina or older MOE jobs:** drawn on demand (short queued `moebatch` run), then cached

---

<style scoped>section { font-size: 22px; } pre { font-size: 0.8em; }</style>

## Ligand Interactions: Serving & Fallback

**Endpoint:** `GET .../interactions2d.svg?style=auto|moe|prolif`
- `X-Diagram-Source` response header reports which renderer was used
- If MOE fails or isn't installed → falls back to the **ProLIF** diagram
- MOE/ProLIF switch in the diagram panel

**Browser fixes (`browser_svg()`):**
- Replaces MOE's embedded SVG font, which browsers ignore, with system sans-serif
- Adds a white background (app panels are dark)

**Example: DPP-4 / sitagliptin (1x70)**
- Neutral amine SMILES → almost no interactions shown
- `[NH3+]` → Glu206 H-bond and Glu205 ionic contact, matching the GUI
- Ligands aren't protonated automatically, so write charged groups explicitly

---

# MOE License Queue

---

## License Queue: The Problem

**Site license:** 6 MOE tokens; each `moebatch` run checks out 2

**Before:** ligands docked in parallel on the shared thread pool
```
ligand 1 ✅  ligand 2 ✅  ligand 3 ✅
ligand 4 ❌  Can't check out license ... Not enough tokens available
```

Also:
- Other users on the same license server compete for tokens
- QuickPrep and on-demand diagrams also need `moebatch`
- Long MOE jobs blocked smina/gnina workers

---

<style scoped>section { font-size: 22px; } pre { font-size: 0.8em; }</style>

## License Queue: Design

**Every `moebatch` call goes through `moe_dock._run_moebatch`:**

```
docking ligands ───┐
QuickPrep ─────────┼──► semaphore (MOE_MAX_CONCURRENT) ──► moebatch
on-demand diagrams ┘                                          │
                         license checkout error? ◄────────────┘
                           yes → sleep MOE_LICENSE_RETRY, retry (≤ MOE_LICENSE_WAIT)
                           no  → done (other errors fail immediately)
```

- **Dedicated MOE thread pool** (`JobManager._moe_pool`) sized to the cap; smina/gnina keep their own workers
- Waiting ligands show status **"queued"** in the UI
- Only license-checkout failures are retried (every 30 s, up to 3600 s by default)
- Each ligand gets its own receptor `.moe` copy, which fixes a shared-file race
- `moebatch` runs from the job folder so Dock's scratch `.mdb` files stay out of the repo root

---

## License Queue: Configuration

| Variable | Default | Meaning |
|----------|---------|---------|
| `MOE_MAX_CONCURRENT` | 1 | Max simultaneous `moebatch` runs (leaves tokens for other users) |
| `MOE_LICENSE_WAIT` | 3600 s | How long to keep retrying while tokens are taken |
| `MOE_LICENSE_RETRY` | 30 s | Delay between retries |
| `MOEBATCH_BIN` | auto | Path to `moebatch` if not under `/Applications/moe*/bin/` |

```bash
MOE_MAX_CONCURRENT=2 ./run.sh   # 2 runs × 2 tokens = 4 of 6 tokens
```

---

<style scoped>section { font-size: 22px; } pre { font-size: 0.8em; }</style>

## MOE QuickPrep

**Receptor panel button → `POST /api/receptor/{id}/quickprep`**

- Batch entry point: `_LIGX_Execute` (the QuickPrep panel's OK action), using your saved QuickPrep settings
- Structure Preparation → Protonate3D → remove distant waters → tethered minimization
- **Selected site ligand goes in the `LIGX_LIGAND` collection**: it defines the flexible region. With no ligand, every atom is fixed and minimization does nothing
- Prepared receptor is stored as a new receptor, shown in the viewer, and used for docking (its hydrogens are kept)

**Report:** hydrogens added, residues built, waters removed, atom shifts, non-default protonation states (site residues first)

**DPP-4 1x70:** His740 → HIP, Asp663 → ASH, ~70 heavy atoms move > 0.5 Å · ~40 s on a cleaned receptor

---

## Code Organization

```
backend/
├── moe_dock.py              # moebatch wrapper, ~550 lines
│   ├── _run_moebatch()      # License queue: semaphore + retry on token errors
│   ├── build_svl_script()   # Template → SVL (incl. covalent options, diagram render)
│   ├── dock()               # Docking orchestration
│   ├── quickprep()          # MOE QuickPrep (_LIGX_Execute)
│   └── render_interactions()# On-demand Ligand Interactions SVG
├── moe_dock_template.svl    # SVL script template, ~60 lines
├── covalent_check.py        # Reactant-form checks, attachment detection
├── docking.py               # Multi-engine dispatcher, MOE SDF/score parsing
├── jobs.py                  # Job store; separate MOE thread pool
├── prep.py                  # Receptor cleanup, QuickPrep report
└── main.py                  # FastAPI endpoints
    ├── /api/pharmacophore/upload, /api/reactions/upload
    ├── /api/receptor/{id}/quickprep
    └── .../poses/{n}/interactions2d.svg
frontend/
├── index.html               # MOE params, covalent controls, QuickPrep button
└── src/app.js               # Engine switching, residue Pick, diagram panel
```

---

## Scoring Functions

### MOE Native Scoring
- **London dG** — Fast empirical scoring (default)
- **GBVI/WSA dG** — Generalized Born/volume integral with weighted surface area
- **Affinity dG** — Calibrated binding affinity prediction
- **ASE** — Atomic solvation energy

### Cross-Engine Compatibility
MOE can also use smina scoring via backend wrapper:
- Vina
- Vinardo
- AutoDock4

**UI:** Scoring dropdown dynamically updates based on selected engine

---

## Performance Characteristics

| Engine | Placement Speed | Refinement | GPU Support | Best For |
|--------|----------------|------------|-------------|----------|
| **smina** | Fast | --minimize | No | High-throughput screening |
| **gnina** | Fast | CNN scoring | Yes | ML-based rescoring |
| **MOE** | Medium | Induced fit | No | Accurate pose prediction, flexible receptors |

**MOE advantages:**
- Induced-fit refinement (side-chain flexibility)
- Better handling of large/flexible ligands
- Pharmacophore-guided placement
- Commercial support & validation

---

## Testing & Validation

### Test Case: DPP-4 (PDB: 1X70)
- Receptor: DPP-4 protein (11,989 atoms)
- Native ligand: 715 (28 atoms)
- Binding site: Auto-detected from HETATM

### Validation Metrics
- RMSD to crystal pose (self-docking)
- Docking score distribution
- Pose clustering
- Interaction fingerprint preservation

### Cross-Engine Comparison
Run same ligand set with smina/gnina/MOE, compare:
- Top-ranked pose RMSD
- Score correlation
- Computational time
- Success rate (RMSD < 2.0Å)

---

## File Cleanup Strategy

**Generated per job:**
- `{ligand}_dock.svl` — Generated SVL script
- `{ligand}_poses.mdb` — MOE database (temp)
- `{ligand}_poses.sdf` — Final output (kept)

**Auto-cleanup on success:**
```python
if script_path.exists():
    script_path.unlink()
if out_mdb.exists():
    out_mdb.unlink()
if receptor_for_moe != receptor:  # converted PDB
    receptor_for_moe.unlink()
```

**Kept on failure:** Script file preserved for debugging
```
(script: /path/to/job/poses/ligand_dock.svl)
```

---

## Error Handling

### Graceful Degradation
```python
try:
    moebatch = find_moebatch()
except MoeBatchNotFound:
    # MOE engine disabled in UI, not an error
    engines_available['moe'] = False
```

### Informative Error Messages
```python
if proc.returncode != 0:
    tail = (proc.stderr or proc.stdout or "").strip().splitlines()
    raise MoeDockingFailed(
        f"{tail[-1]} (script: {script_path})"
    )
```

**User sees:** Exact MOE error + script path for debugging

---

## Configuration

### MOE Discovery
```python
def find_moebatch() -> Optional[str]:
    # 1. Environment variable
    if path := os.environ.get("MOE_BATCH"):
        return path if Path(path).exists() else None
    
    # 2. Common install locations (macOS)
    for base in ["/Applications", Path.home() / "Applications"]:
        for moe_dir in sorted(base.glob("moe*"), reverse=True):
            moebatch = moe_dir / "bin" / "moebatch"
            if moebatch.exists():
                return str(moebatch)
    
    return None
```

**Flexibility:** Works with any MOE version, no hardcoded paths

---

## Future Enhancements

### Short-Term
- [ ] MOE-specific scoring function selection (moe vs smina)
- [ ] Placement method presets (speed vs accuracy)
- [x] Parallel multi-ligand MOE jobs → license-aware queue
- [ ] .ph4 file preview/validation before upload

### Medium-Term
- [ ] Induced-fit receptor flexibility visualization
- [x] MOE-native interaction diagram integration
- [ ] Custom SVL script upload (power users)
- [ ] Benchmark suite (RMSD distribution across PDBbind)

### Long-Term
- [ ] MOE Cloud execution (license pooling)
- [ ] Fragment-based docking workflow
- [ ] Water placement analysis
- [x] Covalent docking support
- [ ] Ligand protonation/wash step for MOE jobs (pH 7)
- [ ] Auto-detect nucleophiles near the site; suggest a reaction from the warhead
- [ ] Multi-residue covalent attachment (crosslinkers)

---

## Lessons Learned

### Technical
1. **Template approach > API bindings** for cross-platform stability
2. **SVL execution model** requires function wrapping, not top-level code
3. **MOE site semantics** differ fundamentally from geometric boxes
4. **Error message archaeology** — check stdout when stderr empty
5. **Bisect against the GUI's `run.sh`**, don't guess option changes
6. **Check ligand chemistry first** when covalent results differ

### Architectural
1. **Engine abstraction layer** critical for multi-engine support
2. **Conditional UI rendering** keeps interface clean
3. **Job ID conventions** improve operational visibility
4. **Unified SDF output** enables cross-engine comparison

---

## Comparison: Integration Effort

| Task | Smina | Gnina | MOE |
|------|-------|-------|-----|
| Wrapper code | ~150 lines | ~200 lines | ~270 lines |
| Template/config | Command args | Command args + GPU | SVL template |
| Auto-detection | Box coords | Box coords | Atom selection |
| Constraints | Post-filter | Post-filter | Native .ph4 |
| UI changes | Minimal | GPU toggle | Placement/refinement |
| **Total effort** | 1 day | 2 days | 4 days |

**MOE complexity drivers:**
- SVL scripting learning curve
- Binding site atom selection semantics
- Pharmacophore API discovery
- Multi-paradigm parameter mapping

---

## Success Metrics

### Technical Achievement
✅ Three docking engines in unified platform
✅ Auto-detection works across all engines
✅ Engine-specific optimizations preserved
✅ Unified SDF output for cross-comparison

### User Experience
✅ Single interface for all engines
✅ Context-aware parameter visibility
✅ Pharmacophore support (MOE-native)
✅ Informative error messages
✅ Job naming shows engine at a glance

### Code Quality
✅ Template-based generation (maintainable)
✅ Graceful degradation (MOE optional)
✅ Comprehensive error handling
✅ Auto-cleanup on success

---

## Demo Setup

### Prerequisites
```bash
# MOE installation (commercial license required)
export MOE_BATCH=/Applications/moe2024.0604/bin/moebatch

# Start server
cd /Users/pdaga/GIT/Smina_docking_bench
python3 -m uvicorn backend.main:app --port 8000

# Open browser
open http://127.0.0.1:8000
```

### Quick Test
1. Fetch PDB: `1x70`
2. Load ligand SMILES: `CC(C)C[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)c1nc(N)nc2[nH]cnc12)C(=O)N[C@@H](CC(N)=O)C(F)(F)F`
3. Select engine: **MOE**
4. Placement: **Triangle Matcher**
5. Refinement: **Rigid Receptor**
6. Run

---

## Summary

### What We Built
- Full MOE integration into existing docking platform
- SVL template system for flexible scripting
- Auto-detection of binding sites from HETATM
- Native pharmacophore constraint support
- Systematic job naming by engine + timestamp
- Covalent docking with MOE reaction templates and reactant-form checks
- MOE Ligand Interactions diagrams per pose (ProLIF fallback)
- License-aware MOE queue and QuickPrep receptor preparation

### Impact
- **Users:** Multi-engine comparison in one interface
- **Platform:** Extensible architecture for future engines
- **Science:** Cross-validation of docking predictions

### Code Stats
- **Initial integration:** ~400 lines (backend + frontend)
- **`moe_dock.py` today:** ~550 lines Python, ~60 lines SVL template

---

## Questions?

**Repository:** `/Users/pdaga/GIT/Smina_docking_bench`

**Key Files:**
- `backend/moe_dock.py` — MOE wrapper
- `backend/moe_dock_template.svl` — SVL template
- `backend/covalent_check.py` — Covalent reactant checks
- `docs/MOE.md` — Full developer notes

**Contact:** Check README for issues/feedback

---

## Appendix: SVL Debugging Tips

### Run SVL script manually
```bash
/Applications/moe2024.0604/bin/moebatch \
  -load /path/to/script.svl \
  -exec "main []"
```

### Common errors
```
"Function 'X' does not exist"
→ Add function declaration at top

"Expressions only allowed in functions"
→ Wrap code in global function main []

"Site selection expression is invalid"
→ Use atom selection, not geometric box

"This placement method requires a pharmacophore"
→ Either upload .ph4 or change placement method

"Reaction did not match protein"
→ rxsite index shifted (Add_H before taking indices) or wrong residue

"Not enough tokens available"
→ Handled by the license queue; lower MOE_MAX_CONCURRENT
```

### Debug output
Check job log: `runs/jobs/{job_id}/logs/{ligand}.log`
