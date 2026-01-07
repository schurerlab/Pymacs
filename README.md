# 🧬 PyMACS  
**Python Molecular Analysis & Clustering Suite**

*A modular GROMACS automation, analysis, and figure-generation toolkit for molecular dynamics workflows*

**Author:** Joseph-Michael Schulz  
**Affiliation:** University of Miami  
**Status:** Active development (official repository)

---

## 📖 Overview

**PyMACS** is a Python-based toolkit for **automating molecular dynamics (MD) simulations**, **trajectory analysis**, and **publication-ready figure generation** using **GROMACS** and the **CHARMM36 force field**.

The project is designed to support **protein**, **protein–ligand**, and **PROTAC / ternary complex** workflows with a strong emphasis on:

- Reproducibility  
- Scriptable, headless execution  
- Scalable analysis across many simulations  
- Clean downstream visualization and reporting  

PyMACS is intended for **computational structural biology**, **drug discovery**, and **targeted protein degradation (TPD)** research.

---

## 🧩 Repository Structure

```

Pymacs/
│
├── 1_AutomateGromacs.py        # System preparation & topology generation
├── 2_AutomateGromacs.py        # Energy minimization, equilibration, production MD
├── 3A_AutomateGromacs.py       # Trajectory analysis (RMSD, RMSF, contacts)
├── 3B_NETWORX.py               # Interaction network & contact graph analysis
│
├── 4PDF4MD.py                  # Automated PDF report / figure compilation
├── 4_MDfigs.txt                # Figure configuration / figure list
│
├── charmm36_ljpme-jul2022.ff/  # CHARMM36 force field (LJ-PME)
│
├── cgenff_charmm2gmx_py3_nx2.py # CGenFF → GROMACS ligand conversion helper
│
├── em.mdp                      # Energy minimization parameters
├── nvt.mdp                     # NVT equilibration parameters
├── npt.mdp                     # NPT equilibration parameters
├── md.mdp                      # Production MD parameters
├── ions.mdp                    # Ion generation parameters
│
├── environment_cgenff.yml      # Conda env: ligand parameterization
├── environment_mdanalysis.yml  # Conda env: MD & analysis
├── recreate_envs.sh            # Convenience script for env recreation
│
├── Example/                    # Reproducible example system (HRP-2 PWWP + CPD32)
│   ├── CPD32_9G94.pdb           # Crystal structure 9G94
│   └── A1D/                     # CGenFF ligand parameters (compound 32)
│       ├── A1D.cgenff.mol2
│       ├── A1D.mol2
│       ├── A1D.str
│       └── A1D.err
│
└── README.md
```

---
## ⚙️ System Requirements

- **OS:** Linux or WSL2  
- **Python:** ≥ 3.9  
- **GROMACS:** 2022+ (non-MPI build, available as `gmx`)  
- **GPU:** Optional but recommended for production MD  
- **Conda / Mamba:** Required  

> ⚠️ This repository assumes **`gmx` (non-MPI)** execution.  
> MPI / cluster usage requires adaptation.

---

## 🧪 Ligand Parameterization Requirements (CGenFF / SILCSBio)

Ligand parameterization in PyMACS relies on **CGenFF** and supports **two equivalent workflows**:

---

### ✅ Option 1 — Local SILCSBio Installation (Recommended)

For fully automated ligand processing, **SILCSBio** must be installed locally and available in your `PATH`.  
SILCSBio is **free for academic users** and can be obtained from:

👉 https://cgenff.com/

After installation, ensure the following environment variables are set:

```bash
export SILCSBIO_HOME=$HOME/silcsbio.2024.1
export PATH="$SILCSBIO_HOME/programs:$SILCSBIO_HOME:$PATH"
```

This configuration enables PyMACS to:

- Automatically extract ligands from the input PDB  
- Generate `.mol2` files  
- Run CGenFF locally  
- Convert CHARMM stream files into GROMACS-compatible formats  

This is the **preferred workflow** for headless execution and high-throughput studies.

---

### 🔁 Option 2 — CGenFF Online Server (No Local Installation Required)

If SILCSBio cannot be installed locally, PyMACS can still be used by leveraging the **CGenFF online server**.

In this case, users must manually:

1. Upload the ligand structure to the CGenFF server  
2. Download the resulting files:  
   - `<LIG>.cgenff.mol2`  
   - `<LIG>.str`  
3. Place both files in the working directory **before running PyMACS**

#### ⚠️ Critical naming and formatting requirements

To ensure correct ligand detection and topology generation:

- The ligand **residue name** in the PDB (three-letter code, e.g. `A1D`) **must match exactly**  
- The MOL2 file **must be named**:

```
A1D.cgenff.mol2
```

- The CHARMM stream file **must be named**:

```
A1D.str
```

- The uploaded ligand structure **must contain hydrogens**  
- Atom naming and bonding must be chemically consistent

Failure to meet these requirements may cause ligand parsing or force-field generation to fail.

Once the correctly named `.cgenff.mol2` and `.str` files are present, PyMACS will automatically detect and use them, bypassing the need for a local SILCSBio installation.

---

### 🧠 Summary

| Workflow | SILCSBio Required | Automation Level |
|--------|------------------|------------------|
| Local SILCSBio | ✅ Yes | Fully automated |
| Online CGenFF | ❌ No | Semi-manual (ligand only) |

Both approaches generate **identical CGenFF parameters** and are fully compatible with downstream PyMACS workflows.

---
## 🧪 Conda Environments

Two reproducible environments are provided:

| Environment | Purpose | File |
|------------|-------|------|
| `cgenff` | Ligand parameterization (CGenFF, Open Babel) | `environment_cgenff.yml` |
| `mdanalysis` | MD execution & analysis | `environment_mdanalysis.yml` |

Create both with:

```bash
conda env create -f environment_cgenff.yml
conda env create -f environment_mdanalysis.yml
```

Or recreate automatically:

```bash
bash recreate_envs.sh
```

---

## 🧱 Step 1 — System Setup  
**`1_AutomateGromacs.py`**

Prepares a complete GROMACS-ready system from an input PDB.

### Capabilities

- Chain parsing and mapping  
- CHARMM36 topology generation  
- Automatic ligand detection (HETATM, non-water)  
- CGenFF-based ligand parameterization  
- Solvation and ionization (physiological salt)  
- Index and atom mapping generation

### Example

```bash
conda activate cgenff
python 1_AutomateGromacs.py 
```

---

## ⚙️ Step 2 — Equilibration & Production MD  
**`2_AutomateGromacs.py`**

Runs the full MD protocol:

**Energy Minimization → NVT → NPT → Production MD**

### Features

- GPU-aware execution  
- Automatic `nsteps` scaling from nanoseconds  
- Protein / ligand temperature coupling groups  
- Checkpoint-safe execution

### Example

```bash
conda activate mdanalysis
python 2_AutomateGromacs.py 
```

---

## 📊 Step 3 — Trajectory Analysis  
**`3A_AutomateGromacs.py`**

Performs structural and dynamical analysis of MD trajectories.

### Outputs

- RMSD and RMSF profiles  
- Residue-wise contact frequencies  
- Ligand–protein interaction persistence

```bash
python 3A_AutomateGromacs.py 
```

---

## 🕸 Step 3B — Interaction Networks  
**`3B_NETWORX.py`**

Generates interaction networks and contact graphs for protein–ligand or ternary systems.

- Residue–residue and residue–ligand networks  
- Graph-based interaction persistence  
- Visual-ready network outputs

---

## 📄 Automated Figures & Reports  
**`4PDF4MD.py`**

Assembles analysis outputs into a **single, publication-ready PDF**.

- Figure ordering controlled via `4_MDfigs.txt`  
- Batch-safe execution across many systems  
- Suitable for supplements and internal reports

```bash
python 4PDF4MD.py 
```
---

## 🧠 Best Practices

- Ensure `gmx` is available in your active environment  
- Keep `charmm36_ljpme-jul2022.ff` in the project root  
- Separate setup (`cgenff`) and MD (`mdanalysis`) environments  
- Retain logs and parameter files for reproducibility

---


## 🧭 Guided Example Walkthrough (Interactive Template)

This section provides a complete, reproducible walkthrough of the PyMACS pipeline using the included example system (HRP-2 PWWP domain bound to compound 32; `CPD32_9G94.pdb`, ligand code `A1D`). The steps below mirror the **actual interactive prompts** encountered when running each script, allowing users to follow along exactly or adapt the workflow for their own systems.

---

### 📁 Step 0 — Prepare the Working Directory

Create a clean working directory and copy all required files into a single location. All PyMACS scripts, parameter files, and inputs are expected to reside in the same directory during execution.

```bash
mkdir MD_Example
cd MD_Example

cp ../Pymacs/*.py .
cp ../Pymacs/*.mdp .
cp -r ../Pymacs/charmm36_ljpme-jul2022.ff .
cp ../Pymacs/Example/CPD32_9G94.pdb .
cp -r ../Pymacs/Example/A1D .
```

This directory now contains:

* The crystal structure (`CPD32_9G94.pdb`)
* Pre-generated CGenFF ligand parameters (`A1D.cgenff.mol2`, `A1D.str`, etc.)
* All PyMACS scripts and MDP templates
* The CHARMM36/LJ-PME force field directory

---

### 🧱 Step 1 — System Preparation & Parameterization (`1_AutomateGromacs.py`)

Run the preparation script in interactive mode:

```bash
python 1_AutomateGromacs.py
```

Typical prompts and recommended responses for the example system are shown below.

**Select input structure**

```
Available PDB files:
1) CPD32_9G94.pdb

Select a PDB file by number: 1
```

**Ligand detection and confirmation**

```
Detected ligand candidates: A1D
Use A1D as ligand? [Y/n]: Y
```

**Chain naming**

```
Detected chains: A
Enter descriptive name for chain A: LEDGF
```

**pH and structure fixing**

```
Enter pH for protonation [default 7.4]:
Fix missing residues and atoms? [Y/n]: Y
```

**CGenFF handling**

* If SILCSBio/CGenFF is available in the PATH, the ligand is parameterized automatically.
* If not, the script will detect the pre-generated `A1D.cgenff.mol2` and `A1D.str` files and proceed without reparameterization.

At completion, the directory contains a fully solvated and ionized system, a master topology (`topol.top`), an index file (`index.ndx`), and an energy-minimization-ready input (`em.tpr`).

---

### ⚙️ Step 2 — Energy Minimization, Equilibration & Production (`2_AutomateGromacs.py`)

Run the simulation driver:

```bash
python 2_AutomateGromacs.py
```

The script executes the following stages sequentially:

* Steepest-descent energy minimization
* NVT equilibration with positional restraints
* NPT equilibration with pressure coupling and gradual restraint release
* Production MD for the user-defined simulation length

Thread count and GPU usage are auto-detected, with interactive confirmation when applicable. Checkpoint files are written automatically, allowing interrupted simulations to be restarted safely.

Successful completion produces the final trajectory (`md_0_1.xtc`), topology (`md_0_1.tpr`), and associated energy and log files.

---

### 📊 Step 3 — Trajectory Analysis & Interaction Mapping (`3A_AutomateGromacs.py`)

Run the analysis script:

```bash
python 3A_AutomateGromacs.py
```

Key interactive prompts include:

**Ligand selection**

```
Detected ligand candidate: A1D
Use A1D as ligand? [Y/n]: Y
```

**Compound display name**

```
Enter compound display name (or press ENTER to use resname): CPD32
```

**Thread configuration**

```
Detected 24 CPU threads available
Enter number of threads to use [ENTER = auto]:
```

The script automatically:

* Recenters the trajectory and extracts a ligand-focused subset
* Identifies the dynamic binding pocket (5 Å cutoff)
* Computes global and per-chain RMSD/RMSF
* Quantifies ligand stability and radius of gyration
* Generates secondary-structure timelines (DSSP)
* Builds residue-level contact maps and interaction frequency plots

All results are written to the `Analysis_Results/` directory.

---

### 🧬 Step 4 — Interaction Networks & Reporting (`3B_NETWORX.py`, `4PDF4MD.py`)

Interaction network diagrams are generated automatically at the end of Step 3 or can be rerun independently. Optional prompts control layout geometry and edge visibility.

Finally, assemble all plots into a single report:

```bash
python 4PDF4MD.py
```

This produces a multi-panel, publication-ready PDF summarizing structural stability, ligand dynamics, interaction persistence, and binding-site organization.

---

### ✅ Outcome

Following this walkthrough yields a fully reproducible MD simulation and analysis of the HRP-2 PWWP–CPD32 complex, serving as a validated reference for adapting PyMACS to new protein–ligand or multi-component systems.

---
## 📘 Reference Output: Example Analysis Figurebook

To further support reproducibility and validation, the `Example/` directory includes a **pre-generated analysis report** produced by running the full PyMACS pipeline on the example system:

```
Example/
├── A1D/
├── CPD32_9G94.pdb
└── MD_ANALYSIS_FIGUREBOOK.pdf
```

**`MD_ANALYSIS_FIGUREBOOK.pdf`** is a consolidated, publication-style figurebook generated by `4PDF4MD.py` following trajectory analysis and interaction mapping. It contains representative outputs from each major analysis stage, including:

* Protein and ligand RMSD and RMSF profiles
* Radius of gyration time series
* Secondary-structure evolution (DSSP)
* Ligand–protein contact maps and interaction frequencies
* Interaction network diagrams and binding-site summaries

This document serves as a **ground-truth reference** for users running the example workflow locally. When executing the guided walkthrough described above, users should obtain qualitatively similar plots and trends (allowing for minor numerical differences due to platform, threading, or stochastic effects). Deviations beyond expected variation may indicate differences in simulation parameters, environment configuration, or execution order.

Including this reference output allows users to:

* Verify that their PyMACS installation is functioning correctly
* Confirm correct ligand detection, topology generation, and analysis behavior
* Compare newly generated results against a known, validated run
* Rapidly understand the scope and structure of PyMACS analysis outputs without first running long simulations

---

### 🔍 Recommended Use

* **First-time users**: Review `MD_ANALYSIS_FIGUREBOOK.pdf` before running simulations to understand expected outputs.
* **Validation runs**: Compare locally generated figures against the reference after completing the walkthrough.
* **Method development**: Use the figurebook as a visual baseline when extending or modifying analysis scripts.


## 🧬 Citation

If you use PyMACS in academic work, please cite:

> Schulz, J.-M. *PyMACS: A Modular GROMACS Automation and Analysis Suite* (University of Miami, 2026).

---

## 📬 Contact

For questions, issues, or contributions, please open a GitHub Issue or contact the jmschulz@med.miami.edu directly.




