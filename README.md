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

This reflects the **current** repository structure (matches `tree`):

```text
PyMACs/
├── 00_GenerateMDP.sh
├── 00_RUNMDfull.sh
├── 1_AutomateGromacs.py
├── 2_AutomateGromacs.py
├── 2A_AutoGMXrestart.py
├── 3A_AutomateGromacs.py
├── 3B_NETWORX.py
├── 4PDF4MD.py
├── 4_Files.txt
├── 4_MDfigs.txt
├── README.md
├── cgenff_charmm2gmx_py3_nx2.py
├── charmm36.ff/
├── charmm36_ljpme-jul2022.ff/
├── MDPs/
│   ├── em.mdp
│   ├── ions.mdp
│   ├── md.mdp
│   ├── npt.mdp
│   └── nvt.mdp
├── Additional/
│   ├── OLD2_AutomateMD.py
│   └── PROTAC_AnalysisScripts.py
├── envs/
│   ├── ENV_MD.yml
│   ├── GMXanalysis.md
│   ├── environment_cgenff.yml
│   └── environment_mdanalysis.yml
├── gromacs_install.md
├── mdanalysis_environment.yml
├── mdanalysis2_environment.yml
├── recreate_envs.sh
└── step3Fix.pdf
```

---

## ⚙️ System requirements

- **OS:** Linux or WSL2 recommended (macOS may work for analysis-only; GROMACS install varies)
- **GROMACS:** 2022+ (non-MPI binary expected)
  - Must be callable as: `gmx`
  - If you use `gmx_mpi`, you’ll need to adapt scripts
- **Python:** managed by conda environments
- **Conda/Mamba:** recommended for reproducibility
- **tmux:** optional but strongly recommended for long runs
- **GPU:** optional but recommended for production MD

---

## 🧪 Environment setup

Create the environments from the YAMLs in `envs/`:

```bash
conda env create -f envs/environment_cgenff.yml
conda env create -f envs/environment_mdanalysis.yml
```

Or use the helper:

```bash
bash recreate_envs.sh
```

See also: `gromacs_install.md`.

---

## 🧬 Force fields included

This repo ships with **two** CHARMM36 variants:

- `charmm36.ff` → classic CHARMM36 (**recommended default for many CGenFF-ligand workflows**)
- `charmm36_ljpme-jul2022.ff` → CHARMM36 **LJ-PME** variant

### Important (force field folder hygiene)

If you keep **multiple** force field folders in the same **run directory**, auto-detection and include-order can become ambiguous.

**Best practice:**
- Run PyMACS in a **clean per-system project folder**
- Copy in **only the force field folder you intend to use** (and keep the original “golden” repo clone untouched)

---

## 🧪 Ligand parameterization and CGenFF support

PyMACS supports CGenFF ligands, but **not all CGenFF output formats are interchangeable** across force field variants and conversion styles.

### ✅ What PyMACS is built around (canonical path)

PyMACS is designed around **CGenFF CHARMM-format outputs**:

- `<LIG>.cgenff.mol2`
- `<LIG>.str`

and a local conversion step via:

- `cgenff_charmm2gmx_py3_nx2.py`  → producing GROMACS-compatible ligand files (e.g., `.itp/.prm` and includes)

### ⚠️ Common failure mode

Pre-converted “GROMACS bundles” from the CGenFF server or third-party converters can differ in:
- include ordering
- atomtype declarations
- duplication/overlap with force field parameters
- assumptions about LJ-PME vs classic CHARMM36

To support real-world lab usage, PyMACS includes **three** workarounds/modes.

---

## ✅ The 3 CGenFF workarounds (PyMACS-supported modes)

### Mode 1 — Import *pre-converted* GROMACS ligand files (advanced / bring-your-own)

Use this when you already have GROMACS-formatted ligand files (e.g., `.itp` + `.prm`, plus any supporting includes).

**Critical rule:**  
If you import already-converted GROMACS ligand files, you must avoid force field duplication and ambiguity.

**Do this in your run folder:**
1. **Keep only the intended force field folder** (recommended: `charmm36.ff` for CGenFF ligand systems)
2. **Remove/move any other force field folders**, especially if you downloaded a conversion bundle that already contains FF-like includes.
   - Example: if you are using classic CHARMM36, remove/move:
     - `charmm36_ljpme-jul2022.ff/`
3. Ensure your `topol.top` includes:
   - exactly one ligand `.itp` include
   - any required ligand parameter include(s) (if used)
   - one `[ molecules ]` entry for the ligand

> If you import ligand files and also keep multiple force fields or conversion folders around, you can get silent parameter conflicts or `grompp` errors.

---

### Mode 2 — CGenFF server outputs only (`.str + .cgenff.mol2`) (recommended / reproducible)

Use this when you can obtain from the CGenFF server:

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

## 🚀 Quick start (interactive)

### 1) Make launchers executable

```bash
chmod +x 00_RUNMDfull.sh 00_GenerateMDP.sh
```

### 2) Run the full pipeline

```bash
./00_RUNMDfull.sh
```

The launcher guides you through:
- selecting an input PDB
- ligand detection/confirmation (if applicable)
- chain naming / mapping
- simulation length
- long-run execution workflow (often with `tmux`)

---

## 🧩 Pipeline stages (scripts)

### 🧱 Step 1 — Setup (`1_AutomateGromacs.py`)
Typical actions:
- reads the PDB and detects chains
- runs `pdb2gmx` for the protein topology
- detects non-water HETATM residues as ligand candidates
- integrates ligand parameters using one of the **3 CGenFF modes**
- solvates + ionizes
- writes helper maps like `atomIndex.txt`

### ⚙️ Step 2 — Equilibration + Production (`2_AutomateGromacs.py`)
Typical actions:
- EM → NVT → NPT → Production MD
- GPU-ready execution (if available)
- thread tuning via environment and script prompts

### 🔁 Step 2A — Restart + Finalize (`2A_AutoGMXrestart.py`)
Typical actions:
- detects `.cpt` checkpoints and resumes
- exports cleaned/final trajectories (e.g., `Final_Trajectory.*`)

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

## 🧭 Guided walkthrough (interactive template)

This section is intentionally **verbose** and mirrors the “what you actually do” experience, including typical prompts.  
Use it as a **copy/paste** scaffold for your own systems.

### 📁 Step 0 — Prepare a clean run directory (recommended)

**Do not run in the repo root.**  
Make a per-system directory and copy in what you need. This prevents:
- accidental forcefield ambiguity
- parameter collisions
- issues if you choose a workflow that removes/moves folders during fallback/repair

Example:

```bash
mkdir -p RUNS/MySystem_01
cd RUNS/MySystem_01

# Copy scripts
cp ../../1_AutomateGromacs.py .
cp ../../2_AutomateGromacs.py .
cp ../../2A_AutoGMXrestart.py .
cp ../../3A_AutomateGromacs.py .
cp ../../3B_NETWORX.py .
cp ../../4PDF4MD.py .
cp ../../4_MDfigs.txt .
cp ../../4_Files.txt .
cp ../../00_RUNMDfull.sh .
cp ../../00_GenerateMDP.sh .
cp ../../cgenff_charmm2gmx_py3_nx2.py .

# Copy MDP templates
mkdir -p MDPs
cp -r ../../MDPs/* ./MDPs/

# Copy exactly ONE force field folder (recommended for CGenFF): classic CHARMM36
cp -r ../../charmm36.ff .

# Copy your input structure(s)
cp /path/to/your/input.pdb .
```

If you intend to use LJ-PME for a protein-only system, copy `charmm36_ljpme-jul2022.ff` instead.

---

### 🧱 Step 1 — System preparation (`1_AutomateGromacs.py`)

Activate the ligand/setup environment:

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

**CGenFF handling (choose your mode)**

- **Mode 2 (recommended):** ensure `A1D.cgenff.mol2` and `A1D.str` exist in the run folder beforehand
- **Mode 3 (automated):** ensure SILCSBio is installed and in PATH
- **Mode 1 (import):** ensure your ligand `.itp/.prm` are present and you removed extra FF folders/includes that could collide

Expected outputs (names may vary slightly by workflow and system):
- `topol.top`
- solvated/ionized coordinate files (`*.gro`)
- `index.ndx`
- a minimization-ready `.tpr` (or the files needed to build it)
- helper mapping files (e.g., `atomIndex.txt`)

---

### ⚙️ Step 2 — EM → NVT → NPT → Production (`2_AutomateGromacs.py`)

Activate the MD environment:

```bash
conda activate mdanalysis
python 2_AutomateGromacs.py
```

Typical flow:
- energy minimization
- equilibration (NVT/NPT)
- production MD (duration often selected interactively)

Outputs commonly include:
- `md_0_1.xtc`
- `md_0_1.tpr`
- `.edr`, `.log`, and checkpoint files (`.cpt`)

---

### 🔁 Step 2A — Restart/finalize (`2A_AutoGMXrestart.py`)

If your run was interrupted or you want “final cleaned exports”:

```bash
conda activate mdanalysis
python 2A_AutoGMXrestart.py
```

Typical exports:
- `Final_Trajectory.xtc`
- `Final_Trajectory.pdb`
- additional convenience files depending on your workflow

---

### 📊 Step 3A — Analysis (`3A_AutomateGromacs.py`)

Run the analysis:

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

Common outputs (workflow-dependent), often stored under:
- `Analysis_Results/`

---

### 🕸 Step 3B — Networks (`3B_NETWORX.py`) (optional)

```bash
conda activate mdanalysis
python 3B_NETWORX.py
```

---

### 📄 Step 4 — PDF figurebook (`4PDF4MD.py`)

Compile plots into a PDF figurebook/report:

```bash
conda activate mdanalysis
python 4PDF4MD.py
```

Figure ordering is controlled via:
- `4_MDfigs.txt`

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

- For known setup/parameter issues and fixes, see: `step3Fix.pdf`
- If you see `grompp` errors about missing bonded/dihedral types, re-check:
  - force field folder selection (`charmm36.ff` vs LJ-PME)
  - CGenFF mode used
  - duplicated includes / atomtypes when importing pre-converted ligand files

---

## 📌 Best practices

- One system per folder; keep run folders clean.
- For CGenFF ligands, prefer **classic `charmm36.ff`** unless you have validated LJ-PME end-to-end.
- Archive for reproducibility:
  - input PDB
  - ligand `.str/.cgenff.mol2` and generated `.itp/.prm` (or your imported ligand files)
  - `mdrun.log`/script logs
  - final `.tpr` and `.cpt`

---

## 🧬 Citation

If you use PyMACS in academic work, please cite the manuscript (in preparation) and/or the software release:

> Schulz, J. M., Reynolds, R. C., & Schürer, S. C. **PyMACS: A Python-Based Automation Suite for GROMACS Molecular Dynamics Setup, Simulation, and Analysis.** (in preparation)

---

## 📬 Contact

- **Joseph M. Schulz:** University of Miami  
- **Stephan C. Schürer (corresponding author):** sschurer@miami.edu
