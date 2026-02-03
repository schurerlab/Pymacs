# 🧬 PyMACS
**PyMACS: A Python-Based Automation Suite for GROMACS Molecular Dynamics Setup, Simulation, and Analysis**

*A modular GROMACS automation, simulation, analysis, and figurebook-generation toolkit for molecular dynamics workflows built around CHARMM36 + CGenFF support.*

**Authors:** Joseph M. Schulz¹, Robert C. Reynolds², Stephan C. Schürer*¹˒³˒⁴  
¹ Department of Molecular and Cellular Pharmacology, University of Miami, Miami, FL, USA  
² O’Neal Comprehensive Cancer Center, University of Alabama at Birmingham, Birmingham, AL 35205, USA  
³ Sylvester Comprehensive Cancer Center, University of Miami Miller School of Medicine, Miami, FL, USA  
⁴ Frost Institute for Data Science & Computing, University of Miami, Miami, FL, USA  
*Corresponding author:* sschurer@miami.edu  

**Repo name:** `PyMACs`  
**Status:** Active development (manuscript + software in development)

---

## 📖 Overview

**PyMACS** is a Python-based toolkit for **automating molecular dynamics (MD) setup**, **running MD simulations**, **trajectory analysis**, and **publication-ready figure/report generation** using **GROMACS** and **CHARMM36**.

It is designed to support:
- **Protein-only systems**
- **Protein–ligand systems**
- **Protein–protein systems**
- **PROTAC / multi-component workflows** (optional utilities)

Core goals:
- Reproducibility and clear logging
- Scriptable, headless execution (works well in `tmux`)
- Scalable/batch-friendly analysis across many simulations
- Clean downstream visualization and reporting

---

## 📝 Manuscript (Tentative / In Development)

**Title:** *PyMACS: A Python-Based Automation Suite for GROMACS Molecular Dynamics Setup, Simulation, and Analysis*  
**Authors:** Joseph M. Schulz¹, Robert C. Reynolds², Stephan C. Schürer*¹˒³˒⁴  
*Corresponding author:* sschurer@miami.edu  

> Manuscript is currently in development. Software behavior, CLI/options, and defaults may evolve.

---

## 📦 What’s in this repo (Current directory layout)

This reflects the **current** repository structure (matches `tree` in this repo):

```text
PyMACs/
├── 1_AutomateGromacs.py
├── 2_AutomateGromacs.py
├── 3A_AutomateGromacs.py
├── 3B_NETWORX.py
├── 4PDF4MD.py
├── 4_MDfigs.txt
├── cgenff_charmm2gmx_py3_nx2.py
├── charmm36.ff/
├── charmm36_ljpme-jul2022.ff/
├── Example/
│   ├── A1D/
│   │   ├── A1D.cgenff.mol2
│   │   ├── A1D.err
│   │   ├── A1D.mol2
│   │   └── A1D.str
│   ├── CPD32_9G94.pdb
│   └── MD_ANALYSIS_FIGUREBOOK.pdf
├── em.mdp
├── ions.mdp
├── md.mdp
├── npt.mdp
├── nvt.mdp
├── environment_cgenff.yml
├── environment_mdanalysis.yml
├── recreate_envs.sh
└── README.md
```

> Your clone may include additional helper scripts/files as the project evolves.

---

## ⚙️ System requirements

- **OS:** Linux or WSL2 recommended  
  - macOS may be fine for **analysis-only**, but GROMACS install and GPU support vary.
- **GROMACS:** 2022+ (non-MPI binary expected)
  - Must be callable as: `gmx`
  - If you use `gmx_mpi`, you’ll need to adapt script calls
- **Python:** managed by conda environments
- **Conda/Mamba:** recommended for reproducibility
- **tmux:** optional but strongly recommended for long runs
- **GPU:** optional but recommended for production MD

---

## 🧪 Environment setup

Create the environments from the YAMLs in the repo root:

```bash
conda env create -f environment_cgenff.yml
conda env create -f environment_mdanalysis.yml
```

Or use the helper:

```bash
bash recreate_envs.sh
```

---

## 🧬 Force fields included (and why BOTH are shipped)

This repo ships with **two** CHARMM36 variants:

- `charmm36.ff` → **classic CHARMM36** (recommended default for many CGenFF-ligand workflows)
- `charmm36_ljpme-jul2022.ff` → **CHARMM36 LJ-PME** variant

### ✅ Important: PyMACS may automatically try BOTH force fields

PyMACS is designed to be practical for “real lab inputs”. In particular, for some systems (especially **ligand/CGenFF** runs), the pipeline may:
1. attempt setup using one force field variant (often **LJ-PME**), then
2. **fall back** to the other (often **classic CHARMM36**) if the first attempt fails due to known parameter/compatibility issues.

**Therefore:**
- **Repo root should contain BOTH force field folders** so the fallback logic is available.
- For per-system run folders, you have two options:
  - **Auto-fallback enabled (recommended for new users):** copy **both** force field folders into the run directory.
  - **Manual/controlled:** copy only the force field you intend to use (advanced), and disable/avoid fallback behavior in your workflow.

---

## 🧪 Ligand parameterization and CGenFF support

PyMACS supports CGenFF ligands, but **not all CGenFF output formats are interchangeable** across conversion styles and force field assumptions.

### ✅ What PyMACS is built around (canonical path)

PyMACS is designed around **CGenFF CHARMM-format outputs**:

- `<LIG>.cgenff.mol2`
- `<LIG>.str`

and a local conversion step via:

- `cgenff_charmm2gmx_py3_nx2.py` → producing GROMACS-compatible ligand include files (e.g., `.itp/.prm`)

### ✅ CGenFF version requirement (validated)

PyMACS workflows are currently validated against **CGenFF 4.6**.

- **Use CGenFF 4.6** for generating `.str` + `.cgenff.mol2` inputs.
- Newer CGenFF releases may produce outputs that are **not currently compatible** with the conversion + include conventions used here.

> If you hit unexpected `grompp`/include/atomtype issues and you’re on a newer CGenFF release, the first troubleshooting step is to regenerate ligand files using **CGenFF 4.6**.

### ⚠️ Common failure mode

Pre-converted “GROMACS bundles” (from servers or third-party converters) can differ in:
- include ordering
- atomtype declarations
- duplication/overlap with force field parameters
- assumptions about LJ-PME vs classic CHARMM36

PyMACS includes multiple modes/workarounds to accommodate common lab realities.

---

## ✅ The 3 CGenFF modes (PyMACS-supported)

### Mode 1 — Import *pre-converted* GROMACS ligand files (advanced / bring-your-own)

Use this when you already have GROMACS-formatted ligand files (e.g., `.itp` + `.prm`, plus supporting includes).

**Rules to avoid breakage:**
1. Make sure your run directory has a consistent force field context (classic vs LJ-PME) with your ligand parameters.
2. Keep include duplication under control:
   - exactly one ligand `.itp` include
   - any required ligand parameter include(s) (if used)
   - one `[ molecules ]` entry for the ligand

> This mode is powerful, but it assumes you understand the parameter topology ecosystem you’re bringing in.

---

### Mode 2 — CGenFF server outputs only (`.str + .cgenff.mol2`) (**recommended**)

Use this when you can obtain from the CGenFF server (or local CGenFF):

- `<LIG>.cgenff.mol2`
- `<LIG>.str`

Place them in your run directory before running PyMACS. PyMACS will convert locally using the included converter.

**Naming must match your PDB residue name exactly** (typically 3 letters):
- PDB resname: `A1D`
- files: `A1D.cgenff.mol2` and `A1D.str`

---

### Mode 3 — Fully automated ligand processing (local SILCSBio/CGenFF)

Use this when SILCSBio/CGenFF is installed locally and available in `PATH`. PyMACS will attempt to:
- extract ligand coordinates from the input PDB,
- prepare/convert ligand,
- generate CGenFF `.str + .cgenff.mol2`,
- convert to GROMACS includes,
- integrate into `topol.top`.

If using SILCSBio, a typical setup looks like:

```bash
export SILCSBIO_HOME=$HOME/silcsbio.2024.1
export PATH="$SILCSBIO_HOME/programs:$SILCSBIO_HOME:$PATH"
```

---

## ✅ Output behavior: overwrite-by-default (intentional)

PyMACS is designed for iterative workflows. **Re-running steps overwrites outputs by default** (e.g., regenerated `.tpr`, refreshed plots, updated PDF figurebooks).

- If you want to preserve a run state, copy/rename the run folder before re-running.
- Best practice is to keep runs in separate folders, e.g.:
  - `RUNS/SystemA_run01/`
  - `RUNS/SystemA_run02/`

---

## 🚀 Quick start templates (how to run things)

Below are **two short templates** you can copy/paste.

### Template A — Run a new system (end-to-end)

> Recommended: do not run in the repo root. Make a per-system folder.

```bash
mkdir -p RUNS/MySystem_01
cd RUNS/MySystem_01

# Copy core scripts
cp ../../1_AutomateGromacs.py .
cp ../../2_AutomateGromacs.py .
cp ../../3A_AutomateGromacs.py .
cp ../../3B_NETWORX.py .
cp ../../4PDF4MD.py .
cp ../../4_MDfigs.txt .
cp ../../cgenff_charmm2gmx_py3_nx2.py .

# Copy MDP templates
cp ../../em.mdp .
cp ../../ions.mdp .
cp ../../nvt.mdp .
cp ../../npt.mdp .
cp ../../md.mdp .

# Copy BOTH forcefields so auto-fallback works
cp -r ../../charmm36.ff .
cp -r ../../charmm36_ljpme-jul2022.ff .

# Copy your input structure
cp /path/to/your/input.pdb .
```

Then run:

```bash
# Step 1: setup (+ ligand integration if present)
conda activate cgenff
python 1_AutomateGromacs.py

# Step 2: equilibration + production
conda activate mdanalysis
python 2_AutomateGromacs.py

# Step 3: analysis
conda activate mdanalysis
python 3A_AutomateGromacs.py

# (optional) network analysis
python 3B_NETWORX.py

# Step 4: figurebook PDF
python 4PDF4MD.py
```

---

### Template B — Analysis-only (when MD is already done)

If you already have trajectories/topologies (e.g., a finished run folder) and only want plots + PDF:

```bash
conda activate mdanalysis
python 3A_AutomateGromacs.py
python 4PDF4MD.py
```

---

## 🧩 Pipeline stages (scripts)

### 🧱 Step 1 — Setup (`1_AutomateGromacs.py`)
Typical actions:
- reads the PDB and detects chains
- runs `pdb2gmx` for the protein topology
- detects non-water HETATM residues as ligand candidates
- integrates ligand parameters using one of the CGenFF modes
- solvates + ionizes
- writes helper maps like `atomIndex.txt` (workflow-dependent)

### ⚙️ Step 2 — Equilibration + Production (`2_AutomateGromacs.py`)
Typical actions:
- EM → NVT → NPT → Production MD
- GPU-ready execution (if available)
- thread tuning via environment and script prompts


## 🔁 Restarting / Resuming Production MD (Checkpoint-Safe Mode)

PyMACS supports **checkpoint-safe resumption** of production molecular dynamics runs using standard GROMACS restart mechanics. This is intended for runs interrupted by walltime limits, node failures, or manual termination.

> **Key idea:** the production length is stored in `md_0_1.tpr`.  
> On restart, **GROMACS uses the existing `.tpr` + `.cpt`**, and **does not read `md.mdp`**.

---

### ✅ Requirements (files that must exist)

To resume production, your run folder must contain:

- `md_0_1.tpr` — production input (created by `gmx grompp`)
- `md_0_1.cpt` — production checkpoint (written during `gmx mdrun`)

Optional but common:
- `md_0_1.xtc`, `md_0_1.edr`, `md_0_1.log` (existing outputs)

---

### 🚀 Resume production only (GPU)

Use this when you want to **skip EM/NVT/NPT** and simply continue production.

```bash
conda activate mdanalysis

python 2_AutomateGromacs.py   --mode ligand   --ligand UNK   --compute GPU   --gpu 0   --production_only   --resume   --headless   --ns 250
```

#### Do I need `--ns` for restarts?
- **Conceptually: NO.** A restart uses `md_0_1.tpr`, so the simulation continues with whatever length is already encoded in that `.tpr`.
- **Practically (current CLI behavior): YES for headless runs.** In the current script, `--ns` avoids an interactive prompt.  
  **To ensure PyMACS does not alter your `.tpr`**, set `--ns` to the *original* runtime (e.g., `250`) or **any value that is ≤ your existing `.tpr` length**.

> If you do not want PyMACS to even *attempt* any TPR extension logic, you can:
> - keep `--ns` equal to the original run length (recommended), OR
> - remove the extension call in `maybe_resume_production_only()` (dev option), OR
> - add a future flag such as `--no_extend` (recommended enhancement).

---

### 🧠 Resume production only (CPU)

```bash
conda activate mdanalysis

python 2_AutomateGromacs.py   --mode ligand   --ligand UNK   --compute CPU   --production_only   --resume   --headless   --ns 250
```

---

### 🧯 Force a clean restart (ignore checkpoint)

If you want PyMACS to **not resume**, even if a checkpoint exists:

```bash
conda activate mdanalysis

python 2_AutomateGromacs.py \
  --mode ligand \
  --ligand UNK \
  --ns 250 \
  --compute GPU \
  --gpu 0 \
  --headless \
  --force_restart

```

This will run the normal pipeline stages again (EM → NVT → NPT → Production).

---

### 🔎 Optional: verify planned production length in the current `.tpr`

If you want to confirm what your `md_0_1.tpr` is configured to run:

```bash
gmx dump -s md_0_1.tpr | grep -E "delta_t|nsteps"
```

Total production time (ns) is:

- **time (ns) = nsteps × delta_t(ps) / 1000**

---

### 📊 Step 3A — Analysis (`3A_AutomateGromacs.py`)
Typical actions (workflow-dependent):
- RMSD / RMSF
- contacts / pocket extraction
- additional structure/interaction analysis depending on your installed tools

### 🕸 Step 3B — Network analysis (`3B_NETWORX.py`)
Optional residue/contact network utilities.

### 🧾 Step 4 — PDF figurebook (`4PDF4MD.py`)
Compiles analysis outputs into a PDF report, ordered via `4_MDfigs.txt`.

---

## 🧭 Guided walkthrough template (interactive)

This section is intentionally **verbose** and mirrors the “what you actually do” experience, including typical prompts.  
Use it as a **copy/paste** scaffold for your own systems.

### 📁 Step 0 — Prepare a clean run directory (recommended)

```bash
mkdir -p RUNS/MySystem_01
cd RUNS/MySystem_01

# Copy core scripts
cp ../../1_AutomateGromacs.py .
cp ../../2_AutomateGromacs.py .
cp ../../3A_AutomateGromacs.py .
cp ../../3B_NETWORX.py .
cp ../../4PDF4MD.py .
cp ../../4_MDfigs.txt .
cp ../../cgenff_charmm2gmx_py3_nx2.py .

# Copy MDP templates
cp ../../em.mdp .
cp ../../ions.mdp .
cp ../../nvt.mdp .
cp ../../npt.mdp .
cp ../../md.mdp .

# Copy BOTH forcefields so auto-fallback works
cp -r ../../charmm36.ff .
cp -r ../../charmm36_ljpme-jul2022.ff .

# Copy your input structure(s)
cp /path/to/your/input.pdb .
```

---

### 🧱 Step 1 — System preparation (`1_AutomateGromacs.py`)

```bash
conda activate cgenff
python 1_AutomateGromacs.py
```

Typical prompt flow (example):

**Select input PDB**
```
Available PDB files:
1) input.pdb
Select a PDB file by number: 1
```

**Ligand detection (if present)**
```
Detected ligand candidates: A1D
Use A1D as ligand? [Y/n]: Y
```

**Chain naming / mapping**
```
Detected chains: A B
Enter descriptive name for chain A: ProteinA
Enter descriptive name for chain B: ProteinB
```

---

### ⚙️ Step 2 — EM → NVT → NPT → Production (`2_AutomateGromacs.py`)

```bash
conda activate mdanalysis
python 2_AutomateGromacs.py
```

Typical flow:
- energy minimization
- equilibration (NVT/NPT)
- production MD (duration often selected interactively)

---

> **Need to resume a crashed/paused production run?** See **“🔁 Restarting / Resuming Production MD (Checkpoint-Safe Mode)”** above.


### 📊 Step 3A — Analysis (`3A_AutomateGromacs.py`)

```bash
conda activate mdanalysis
python 3A_AutomateGromacs.py
```

Typical prompts:
```
Detected ligand candidate: A1D
Use A1D as ligand? [Y/n]: Y
Enter compound display name (or press ENTER to use resname): CPD32
Detected 24 CPU threads available
Enter number of threads to use [ENTER = auto]:
```

---

### 📄 Step 4 — PDF figurebook (`4PDF4MD.py`)

```bash
conda activate mdanalysis
python 4PDF4MD.py
```

Figure ordering is controlled via:
- `4_MDfigs.txt`

---

## 🧪 Reproduce the included Example figurebook (demo)

The repo includes:

- `Example/CPD32_9G94.pdb` (demo input)
- `Example/A1D/A1D.str` + `Example/A1D/A1D.cgenff.mol2` (demo ligand inputs)
- `Example/MD_ANALYSIS_FIGUREBOOK.pdf` (the expected figurebook style/output)

### Goal
Run the pipeline using the included Example inputs and generate a figurebook PDF in your run directory that matches the **structure and content** of the included example.

> Note: exact byte-for-byte identity can vary across OS/matplotlib versions, but the **figure order and content types** should match when using the same settings.

### Steps

```bash
mkdir -p RUNS/Example_CPD32
cd RUNS/Example_CPD32

# Copy scripts + templates
cp ../../1_AutomateGromacs.py .
cp ../../2_AutomateGromacs.py .
cp ../../3A_AutomateGromacs.py .
cp ../../3B_NETWORX.py .
cp ../../4PDF4MD.py .
cp ../../4_MDfigs.txt .
cp ../../cgenff_charmm2gmx_py3_nx2.py .
cp ../../em.mdp .
cp ../../ions.mdp .
cp ../../nvt.mdp .
cp ../../npt.mdp .
cp ../../md.mdp .

# Copy BOTH forcefields (required for auto-fallback demo)
cp -r ../../charmm36.ff .
cp -r ../../charmm36_ljpme-jul2022.ff .

# Copy Example inputs
cp ../../Example/CPD32_9G94.pdb .
mkdir -p A1D
cp ../../Example/A1D/A1D.str ./A1D/
cp ../../Example/A1D/A1D.cgenff.mol2 ./A1D/
```

Run the pipeline:

```bash
conda activate cgenff
python 1_AutomateGromacs.py

conda activate mdanalysis
python 2_AutomateGromacs.py
python 3A_AutomateGromacs.py
python 4PDF4MD.py
```

### What to compare
- Your generated figurebook (typically in the run folder) vs:
  - `../../Example/MD_ANALYSIS_FIGUREBOOK.pdf`

If your PDF is missing sections:
- confirm the analysis stage completed successfully
- confirm `4_MDfigs.txt` points to plots that exist in your run folder
- confirm you are using CGenFF **4.6**-generated ligand files (or the included example ones)

---

## 🧠 tmux quick tips (recommended)

| Action | Command |
|---|---|
| Start a named session | `tmux new -s MDRUN` |
| Detach | `Ctrl+b`, then `d` |
| List sessions | `tmux ls` |
| Re-attach | `tmux attach -t MDRUN` |
| Kill | `tmux kill-session -t MDRUN` |

---

## 🧯 Troubleshooting notes

Common causes of `grompp` failures:
- force field mismatch (classic vs LJ-PME)
- duplicated includes / atomtypes when importing pre-converted ligand files
- unexpected CGenFF output conventions (often from newer CGenFF versions)

If a run fails early:
- keep both force fields present and allow the fallback logic to run
- validate ligand input files are present and named correctly
- ensure you are using CGenFF **4.6** if you are generating new ligand `.str/.mol2`

---

## 📌 Best practices

- One system per folder; keep run folders clean.
- Let the pipeline overwrite outputs in a run folder (intended behavior), but **clone the folder** if you need an archive snapshot.
- For CGenFF ligands, prefer classic `charmm36.ff` unless you have validated LJ-PME end-to-end.
- Archive for reproducibility:
  - input PDB
  - ligand `.str/.cgenff.mol2` and generated `.itp/.prm` (or your imported ligand files)
  - script logs
  - final `.tpr` and `.cpt`

---

## 🧬 Citation

If you use PyMACS in academic work, please cite the manuscript (in preparation) and/or the software release:

> Schulz, J. M., Reynolds, R. C., & Schürer, S. C. **PyMACS: A Python-Based Automation Suite for GROMACS Molecular Dynamics Setup, Simulation, and Analysis.** (in preparation)

---

## 📬 Contact

- **Joseph M. Schulz:** University of Miami  
- **Stephan C. Schürer (corresponding author):** sschurer@miami.edu
