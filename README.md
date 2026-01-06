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

## 🧬 Citation

If you use PyMACS in academic work, please cite:

> Schulz, J.-M. *PyMACS: A Modular GROMACS Automation and Analysis Suite* (University of Miami, 2026).

---

## 📬 Contact

For questions, issues, or contributions, please open a GitHub Issue or contact the jmschulz@med.miami.edu directly.




