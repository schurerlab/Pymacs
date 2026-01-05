# 🧬 **GROMACS Automation Suite (GMX-Automation)**

### End-to-End Molecular Dynamics Workflow for Proteins, Ligands, and PROTAC Systems

**Author:** Joseph-Michael Schulz (University of Miami)
**Version:** 3.0 — October 2025

---

## 📖 Overview

This repository provides a **fully automated molecular dynamics (MD)** pipeline built around **GROMACS 2023+**, **CHARMM36**, and **Conda-based environments**.
The workflow prepares, equilibrates, simulates, and analyzes complex biomolecular systems (proteins, protein–ligand, or PROTAC ternary systems) with a **single command**.

**Note:** This workflow assumes that the `gmx` binary (non-MPI version) is available within your active environment.
It is *not* intended for `gmx_mpi` or cluster MPI execution unless explicitly adapted.

---

## 🧩 Directory Structure

```
GMX-Automation/
│
├── 00_RUNMDfull.sh            # Master launcher: setup → simulation → logging (tmux compatible)
├── 1_AutomateGromacs.py       # Step 1: system setup (chains, topology, solvation, ions)
├── 2_AutomateGromacs.py       # Step 2: equilibration + production MD (GPU-ready)
├── 2A_AutoGMXrestart.py       # Step 2B: restart & postprocessing (resume MD + wrap trajectory + pocket)
├── 3A_AutomateGromacs.py      # Step 3: advanced RMSD/RMSF + contact/pocket analysis
│
├── MDPs/                      # Energy minimization, NVT, NPT, MD parameter templates
├── charmm36_ljpme-jul2022.ff/ # CHARMM36 forcefield directory
│
├── cgenff_charmm2gmx_py3_nx2.py  # Ligand topology conversion helper
├── PROTAC_AnalysisScripts.py     # (optional) Advanced PROTAC-specific utilities
│
├── environment_cgenff.yml     # Conda environment for ligand parameterization
├── environment_mdanalysis.yml # Conda environment for MD + analysis
│
└── GMXanalysis.md             # (optional) Notes or prior run summaries
```

---

## ⚙️ System Requirements

* **OS:** Linux or WSL2
* **GROMACS:** 2022+ (GPU build, non-MPI) — must be accessible as `gmx`
* **Python:** ≥3.9
* **Conda Environments:**

  * `cgenff` → ligand parameterization (Open Babel, SILCSBio, CGenFF)
  * `mdanalysis` → simulation & postprocessing (GROMACS, MDAnalysis, MDTraj, NumPy, etc.)
* **tmux:** optional but recommended for background runs
  Install: `sudo apt install tmux`

---

## 🚀 Quick Start (Full Pipeline)

### 1️⃣ Make the launcher executable:

```bash
chmod +x 00_RUNMDfull.sh
```

### 2️⃣ Run interactively:

```bash
./00_RUNMDfull.sh
```

You will be prompted to:

* Choose your input PDB file
* Select simulation type (Ligand–Protein, Protein–Protein, PROTAC, or Protein-only)
* Confirm or specify ligand (3-letter code)
* Assign chain names
* Define simulation time (e.g., 50 ns)

The launcher automatically:

* Activates required Conda environments (`cgenff` → setup, `mdanalysis` → MD run)
* Executes both setup (`1_AutomateGromacs.py`) and simulation (`2_AutomateGromacs.py`)
* Logs all output to `mdrun.log`
* Runs within a persistent **tmux session** (`MDRUN` by default)

---

## 🧠 Headless / Batch Mode Usage

You can skip all prompts by predefining environment variables:

```bash
SESSION=MDRUN \
PROJECT_DIR=$PWD \
LIG_CODE=GDP \
SIM_TYPE=1 \
MD_NS=100 \
./00_RUNMDfull.sh
```

Or export them beforehand:

```bash
export SESSION=Rap1_MD
export PROJECT_DIR=$PWD
export LIG_CODE=GDP
export SIM_TYPE=1
./00_RUNMDfull.sh
```

---

## 🧩 Step-by-Step Breakdown

### 🧱 Step 1 — System Setup (`1_AutomateGromacs.py`)

Configures your molecular system: chain mapping, topology generation, solvation, and ionization.

**Features:**

* Automatic ligand detection (non-water HETATM records)
* CHARMM36 forcefield selection (auto-detected)
* CGenFF-based ligand parameterization (SILCSBio)
* Solvation and neutralization (0.15 M NaCl)
* Writes processed topology and coordinate files
* Generates `atomIndex.txt` mapping chain names and indices

**Example (headless):**

```bash
python 1_AutomateGromacs.py \
  --pdb Model_2.pdb \
  --chain-map "A:Rap1B,B:Rap1GAP" \
  --ligand GDP
```

**Output Example:**

```
Model_2/
├── topol.top
├── complex.gro
├── solv_ions.gro
├── GDP.itp / GDP.prm
├── index.ndx
└── em.tpr
```

---

### ⚙️ Step 2 — Equilibration & Production (`2_AutomateGromacs.py`)

Performs the **Energy Minimization → NVT → NPT → MD** workflow in sequence.

**Modes:**

* `protein` — single-protein simulation
* `ligand` — protein–ligand complex
* `protac` — ternary (E3–linker–target) complex

**Example (100 ns GPU run):**

```bash
python 2_AutomateGromacs.py --mode ligand --ligand GDP --ns 100 --gpu 0
```

**Outputs:**

```
em.gro, nvt.gro, npt.gro
md_0_1.xtc, md_0_1.tpr, md_0_1.log, md_0_1.cpt
```

**Key Features:**

* Automatic GPU detection (`nvidia-smi`)
* Adaptive multithreading (`OMP_NUM_THREADS`)
* Dynamic update of tc-grps (protein/ligand)
* Automatic `nsteps` scaling based on nanoseconds
* Ligand restraint integration into topology

---

### 🔁 Step 2B — Restart & Finalization (`2A_AutoGMXrestart.py`)

Handles continuation from checkpoints and trajectory cleanup.

**Example:**

```bash
python 2A_AutoGMXrestart.py --ligand GDP --gpu 0
```

**Actions:**

* Detects `md_0_1.cpt` checkpoint and resumes simulation
* Wraps and recenters the trajectory (`Final_Trajectory.xtc`)
* Extracts binding pocket (`binding_pocket_only.pdb/.xtc`)
* Generates visual-ready `Final_Trajectory.pdb`

**Outputs:**

```
Final_Trajectory.xtc / Final_Trajectory.pdb
binding_pocket_only.xtc / binding_pocket_only.pdb
```

---

### 📊 Step 3 — Trajectory Analysis (`3A_AutomateGromacs.py`)

Performs comprehensive RMSD/RMSF analysis and interaction network mapping.

**Example:**

```bash
python 3A_AutomateGromacs.py \
  --topo Final_Trajectory.pdb \
  --traj Final_Trajectory.xtc \
  --ligand GDP \
  --threads 16
```

**Outputs:**

```
Analysis_Results/
├── RMSD_RMSF_Plots/
├── residue_contact_frequency.png
├── GDP_interaction_network.png
└── binding_pocket_only.pdb/.xtc
```

**Features:**

* Parallelized frame alignment (multi-core)
* GPU-accelerated distance computations (CuPy)
* Ligand–protein contact persistence heatmaps
* RMSD, RMSF, and residue-wise interaction plots

---

## 🧪 Environment Setup (Reference)

Two reproducible environments are provided:

| Environment  | Purpose                 | YAML File                    |
| ------------ | ----------------------- | ---------------------------- |
| `cgenff`     | Ligand parameterization | `environment_cgenff.yml`     |
| `mdanalysis` | Simulation & analysis   | `environment_mdanalysis.yml` |

Recreate them on another system with:

```bash
conda env create -f environment_cgenff.yml
conda env create -f environment_mdanalysis.yml
```

---

## 🧭 tmux Session Management

| Action            | Command                      |
| ----------------- | ---------------------------- |
| View sessions     | `tmux ls`                    |
| Attach to session | `tmux attach -t MDRUN`       |
| Detach safely     | `Ctrl + B`, then `D`         |
| Kill session      | `tmux kill-session -t MDRUN` |

All logs and progress are written to `mdrun.log` in real time.

---

## 📁 Example Output Layout

```
MyProtein_MD/
├── mdrun.log
├── .mdrun_payload.sh
├── protein_processed.gro
├── topol.top
├── em.gro → nvt.gro → npt.gro → md_0_1.xtc
├── md_0_1.cpt (checkpoint)
├── Final_Trajectory.xtc / .pdb
├── binding_pocket_only.pdb / .xtc
└── Analysis_Results/
```

---

## 🧠 Best Practices

* Ensure `gmx` (non-MPI) is accessible from your environment path.
* Always include `charmm36_ljpme-jul2022.ff` in your project directory.
* Run inside `tmux` to prevent job interruption.
* Use `2A_AutoGMXrestart.py` for long simulations to enable checkpoint recovery.
* Prefer GPU-enabled environments for production runs.
* Annotate all run details (ligand, time, GPU, date) in `mdrun.log` for reproducibility.

---

## 🧬 Citation

If you use this pipeline in your work, please cite:

> Schulz, J.-M. *Automated GROMACS MD Workflow for Multi-Protein and PROTAC Systems* (University of Miami, 2025).
