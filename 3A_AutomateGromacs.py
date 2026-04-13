#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=======================================================================
💠 Advanced MD Trajectory Analyzer — RMSD/RMSF + Pocket Interaction Maps
=======================================================================

Author:  Joseph-Michael Schulz (PhD Candidate, UM Biochemistry & Molecular Biology)
Version: 2.2 — October 2025
Dependencies: MDAnalysis, MDTraj, CuPy (optional), NumPy, Pandas, Matplotlib, Seaborn, NetworkX, tqdm

📘 Overview
-----------
This script performs *comprehensive structural trajectory analysis* of molecular
dynamics simulations, integrating **MDTraj**, **MDAnalysis**, and **CuPy** for GPU acceleration.

It automatically detects the ligand, extracts the binding pocket, and computes
high-resolution RMSD/RMSF profiles, residue contact maps, and protein–ligand
interaction networks — all in a reproducible, step-wise workflow with full multithreading
support and detailed CLI logging.

💡 Key Features
---------------
1️⃣ **Automatic Ligand Detection**
    • Detects the largest non-protein, non-water residue in the system (or use --ligand)
    • Optional fallback via atom indices (--lig-atom-start / --lig-atom-end)

2️⃣ **OpenMP / Thread Optimization**
    • Interactive detection of CPU threads with auto-scaling to 66 % of total cores
    • Environment variable propagation for MKL / OpenBLAS / NumExpr
    • Compatible with headless or bash-pipeline execution

3️⃣ **GPU Acceleration (CuPy-native)**
    • Vectorized ligand–protein distance and contact calculations on GPU
    • Automatic fallback to NumPy on CPU when GPU unavailable

4️⃣ **Dynamic Binding Pocket Refinement**
    • Chain-aware residue selection within 5 Å of the ligand
    • Writes reduced pocket-only trajectory and topology for downstream analysis

5️⃣ **Comprehensive RMSD/RMSF Analysis**
    • Global backbone RMSD
    • Ligand RMSD and per-atom RMSF
    • Per-protein RMSD/RMSF from atomIndex.txt
    • Bound-complex and full-complex RMSD time series

6️⃣ **Residue Contact & Interaction Analysis**
    • Per-frame ligand–protein contact mapping (with tqdm progress bars)
    • Heatmap of contact persistence over simulation time
    • Contact frequency histogram (residue-wise)
    • Classification of H-bond, hydrophobic, and ionic interactions
    • Automatically colored network graph (percent of frames per residue)

7️⃣ **Robust CLI & Logging**
    • Modular `phase()` / `tick()` / `done()` functions for live status feedback
    • Auto-timestamps, elapsed phase timing, and parallel progress visualization
    • Clean failure handling for missing files, ligand misnaming, or topology issues

📊 Output Summary
----------------
All results are written to the output directory (`--outdir`, default: `Analysis_Results/`):

│── RMSD / RMSF Data
│   ├── rmsd_over_time.csv / .png
│   ├── {LIG}_rmsf.csv / _rmsf_plot.png
│   ├── RMSD_{CHAIN}.csv / .png  (per chain)
│   ├── RMSF_{CHAIN}.csv / .png
│   ├── {LIG}_BoundComplex_RMSD.csv / .png
│   ├── GlobalComplex_RMSD.csv / .png
│
│── Pocket Contact & Interaction Analysis
│   ├── {LIG}_contact_map_residue_time.png
│   ├── residue_contact_frequency.png
│   ├── interaction_stacked_normalized.png
│   └── {LIG}_interaction_network_percent_in_nodes_pastel.png
│
│── Pocket & Subset Trajectories
│   ├── Final_Trajectory.xtc / .pdb (Protein|Ligand subset)
│   ├── binding_pocket_only.xtc / .pdb
│
└── Log Output
    └── Detailed phase-by-phase progress with timestamps and thread summary

⚙️ Example Usage
----------------
Interactive (auto-detect threads):
    $ python analyze_md.py --topo md_0_1.gro --traj md_0_1.xtc -l GDP

Headless (predefined threading):
    $ export OMP_NUM_THREADS=16
    $ python analyze_md.py --threads 16 --topo md_0_1.gro --traj md_0_1.xtc

🧩 Workflow Summary
-------------------
STEP A  — Recenter & Subset (Protein|Ligand)
STEP B  — Dynamic Pocket Refinement (5 Å cutoff)
STEP C  — Identify Binding Chains (via atomIndex.txt)
STEP D1 — Global Protein RMSD
STEP D2 — Ligand RMSD / RMSF
STEP D3 — Per-Protein RMSD/RMSF
STEP D4 — Bound Complex RMSD
STEP D5 — Full Complex RMSD
STEP E1 — Pocket Contact Analysis (MDAnalysis)
STEP E2 — Interaction Classification (H-bond / Hydrophobic / Ionic)
STEP E3 — Network Graph Generation

✅ Final Deliverable
-------------------
A complete post-processing suite producing quantitative and visual analyses
of ligand binding, residue dynamics, and interaction persistence — ideal for
supplementary figures, QC validation, and comparative complex evaluation.

=======================================================================
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import gc, MDAnalysis as mda
import mdtraj as md
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import networkx as nx
import multiprocessing
from tqdm import tqdm
import time
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor
import mdtraj as md
import matplotlib.pyplot as plt
from matplotlib.patches import Patch




def phase(msg):
    print("\n" + "="*70)
    print(f"🧩 {msg.upper()} — {datetime.now().strftime('%H:%M:%S')}")
    print("="*70)
    time.sleep(0.3)

def tick(msg):
    print(f"   ⏳ {msg}...", end="", flush=True)

def done(msg="Done"):
    print(f" ✅ {msg}")

# ------------------------- CLI -------------------------
parser = argparse.ArgumentParser(
    description="Generate RMSD/RMSF and pocket interaction plots (robust ligand detection)."
)

parser.add_argument(
    "-l", "--ligand", type=str, default=None,
    help="Ligand resname (e.g., PTC, GDP). If omitted, auto-detects largest non-protein residue."
)

parser.add_argument(
    "--lig-atom-start", type=int, default=None,
    help="Fallback: starting atom index (global, 0-based) for ligand in MDTraj topology."
)

parser.add_argument(
    "--lig-atom-end", type=int, default=None,
    help="Fallback: ending atom index (inclusive) for ligand in MDTraj topology."
)

# 👇 Combined alias: --topo OR --gro (both map to args.topo)
parser.add_argument(
    "--topo", "--gro", default="Final_Trajectory.pdb",
    help="Topology for MDTraj (e.g., .gro or .pdb). (--gro is an alias for convenience.)"
)

parser.add_argument(
    "--traj", default="Final_Trajectory.xtc",
    help="Trajectory for MDTraj (e.g., .xtc)."
)

parser.add_argument(
    "--pocket_pdb", default="binding_pocket_only.pdb",
    help="Pocket-only PDB for MDAnalysis."
)

parser.add_argument(
    "--pocket_xtc", default="binding_pocket_only.xtc",
    help="Pocket-only XTC for MDAnalysis."
)

parser.add_argument(
    "--outdir", default="Analysis_Results",
    help="Output directory for results and plots."
)

parser.add_argument(
    "--contact_cutoff", type=float, default=4.0,
    help="Contact cutoff distance in Å (default 4.0)."
)

parser.add_argument(
    "--min_contact_frac", type=float, default=0.10,
    help="Minimum fraction of frames a residue must contact ligand to appear in plots (default 0.10 = 10%)."
)

parser.add_argument(
    "--threads", type=int, default=None,
    help="Number of CPU threads to use (default: auto-detect 66%% of available cores). "
         "If empty in interactive mode, it will auto-fill."
)

parser.add_argument(
    "--compound-name",
    type=str,
    default=None,
    help="Readable compound name to use in output (e.g., Paprotrain or CPD32_Inhibitor)"
)

parser.add_argument(
    "--numbering-template",
    type=str,
    default=None,
    help="Optional PDB file to use as residue-numbering template."
)

# ------------------ NETWORX-compatible flags ------------------
parser.add_argument(
    "--minfrac", type=float, default=None,
    help="Minimum fraction of frames a residue must contact ligand (NETWORX threshold, default 0.10 = 10%)."
)

parser.add_argument(
    "--net-out", type=str, default=None,
    help="Optional output base name for NETWORX figure (e.g., A1D_network)."
)

parser.add_argument(
    "--ellipse-rx", type=float, default=None,
    help="Override ellipse horizontal radius (NETWORX)."
)

parser.add_argument(
    "--ellipse-ry", type=float, default=None,
    help="Override ellipse vertical radius (NETWORX)."
)

parser.add_argument(
    "--no-edges", action="store_true",
    help="If set, NETWORX expanded panel will hide all edges (node-only mode)."
)


args = parser.parse_args()
NUMBERING_TEMPLATE = args.numbering_template
MIN_CONTACT_FRAC = args.min_contact_frac



def list_pdb_files():
    return sorted([
        f for f in os.listdir(".")
        if f.lower().endswith(".pdb") and os.path.isfile(f)
    ])


def select_numbering_template_interactive():
    pdbs = list_pdb_files()

    if not pdbs:
        print("⚠️ No PDB files found in current directory.")
        return None

    print("\n🧬 Residue numbering template selection")
    print("Do you want to use a PDB file as a residue-numbering template?")
    resp = ask("Use PDB template for residue numbering? [y/N]: ", str, default="n")

    if not resp.lower().startswith("y"):
        return None

    print("\n📂 Available PDB files:")
    for i, f in enumerate(pdbs, 1):
        print(f"   {i}) {f}")

    choice = ask(
        "Select a PDB by number (or press ENTER to cancel): ",
        str,
        default=""
    )

    if not choice:
        print("⏭️ No template selected — skipping residue remapping.")
        return None

    if choice.isdigit() and 1 <= int(choice) <= len(pdbs):
        selected = pdbs[int(choice) - 1]
        print(f"✅ Using {selected} as residue numbering template.")
        return selected

    print("⚠️ Invalid selection — skipping residue remapping.")
    return None





# ------------------------------------------------------------
# INTERACTIVE PROMPTS FOR NETWORX ARGUMENTS (if missing)
# ------------------------------------------------------------
import string

def ask(prompt, cast=str, default=None):
    """
    Ask user for input interactively with validation.
    - Strips non-printable characters
    - Retries on invalid input
    - Falls back to default on ENTER or EOF
    """
    while True:
        try:
            raw = input(prompt)

            # Remove non-printable / control characters
            v = "".join(ch for ch in raw if ch in string.printable).strip()

            if v == "":
                return default

            try:
                return cast(v)
            except ValueError:
                print(
                    f"❌ Invalid input: '{raw}'.\n"
                    f"   Please enter a valid {cast.__name__} "
                    f"(digits only, e.g. 0.2)"
                )

        except EOFError:
            return default

# ------------------------------------------------------------
# ALWAYS FORCE NETWORX OUTPUT INTO Analysis_Results/NETWORX
# ------------------------------------------------------------
lig_code = args.ligand or os.path.basename(args.topo).split('.')[0]

NET_OUTDIR = os.path.join("Analysis_Results", "NETWORX")
os.makedirs(NET_OUTDIR, exist_ok=True)

# The base output stem used for all NETWORX panel images (A, B, C)
args.net_out = os.path.join(NET_OUTDIR, lig_code)


# If no --minfrac provided explicitly → ask
if args.minfrac is None:
    args.minfrac = ask(
        "📊 Minimum fraction of frames for residue inclusion [0.0–1.0] (default 0.10): ",
        float,
        default=0.10
    )

# If ellipse radii were not provided → ask interactively
if args.ellipse_rx is None:
    args.ellipse_rx = ask(
        "⚪ NETWORX ellipse horizontal radius (rx) (press ENTER to use auto layout): ",
        float,
        default=None
    )

if args.ellipse_ry is None:
    args.ellipse_ry = ask(
        "⚪ NETWORX ellipse vertical radius (ry) (press ENTER to use auto layout): ",
        float,
        default=None
    )

# If --no-edges was not passed → ask yes/no
# Default behavior: KEEP edges.
# Only hide edges if:
#  1) user passed --no-edges
#  2) or they explicitly say "y" to hiding

if not args.no_edges:
    v = ask("Hide edges in NETWORX graph? (y/n, default n): ", str, default="n")
    if v.lower().startswith("y"):
        args.no_edges = True
    else:
        args.no_edges = False







PLOT_TYPES = ["H-bond", "Hydrophobic", "Ionic"]
# ------------------------------------------------------------
# Define interaction types and colors for plotting
# ------------------------------------------------------------



# Default color map (safe + readable)
DEFAULT_COLORS = {
    "H-bond": "#1f77b4",        # blue
    "Hydrophobic": "#ff7f0e",   # orange
    "Ionic": "#2ca02c",         # green
}

# Assign colors, fallback to gray if unknown type
COLORS = {
    t: DEFAULT_COLORS.get(t, "#7f7f7f")
    for t in PLOT_TYPES
}

print("🎨 Plot interaction types:", PLOT_TYPES)
print("🎨 Assigned colors:", COLORS)


TOPO_FILE   = args.topo
TRAJ_FILE   = args.traj
POCKET_PDB  = args.pocket_pdb
POCKET_XTC  = args.pocket_xtc
OUTPUT_DIR  = args.outdir
CONTACT_CUTOFF = float(args.contact_cutoff)
MIN_CONTACT_FRAC = float(args.min_contact_frac)
MAX_WORKERS = max(1, os.cpu_count() // 2)
# ============================================================
# MIAMI COLOR PALETTE (LOCKED)
# ============================================================
MIAMI_ORANGE = "#F47321"   # Protein
MIAMI_GREEN  = "#005030"   # Ligand
MIAMI_BLUE   = "#0033A0"   # Water

# ============================================================


os.makedirs(OUTPUT_DIR, exist_ok=True)
# Create subfolder for distance/contact histograms
DIST_DIR = os.path.join(OUTPUT_DIR, "RESIDUE_Contact_Distance")
os.makedirs(DIST_DIR, exist_ok=True)

PIE_DIR = os.path.join(OUTPUT_DIR, "RESIDUE_Pies")
os.makedirs(PIE_DIR, exist_ok=True)

# ============================================================
# 🧠 OpenMP Thread Management
# ============================================================
if args.threads is None:
    try:
        # Interactive mode: show available cores
        total_cores = multiprocessing.cpu_count()
        print(f"\n🧩 Detected {total_cores} total CPU threads available.")
        user_input = input("⚙️ Enter number of threads to use [press ENTER for auto (66%)]: ").strip()
        if user_input:
            threads = int(user_input)
        else:
            threads = max(1, round(total_cores * 0.66))
    except EOFError:
        # Non-interactive / headless run (e.g., bash pipeline)
        threads = max(1, round(multiprocessing.cpu_count() * 0.66))
else:
    threads = args.threads

# Use either environment or user-provided override
threads_env = os.environ.get("OMP_NUM_THREADS")
if threads_env and threads_env.isdigit():
    print(f"🧠 OMP_NUM_THREADS already set in environment = {threads_env}")
    threads = int(threads_env)
else:
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["MKL_NUM_THREADS"] = str(threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(threads)
    os.environ["NUMEXPR_NUM_THREADS"] = str(threads)
    os.environ["OMP_PROC_BIND"] = "spread"
    os.environ["OMP_PLACES"] = "cores"

print(f"⚙️ Using {threads} CPU threads for RMSD/RMSF calculations.")

# --- Optional GPU acceleration (CuPy) ---
try:
    import cupy as cp
    GPU_AVAILABLE = True
    print("⚡ CuPy GPU acceleration enabled for distance calculations.")
except ImportError:
    import numpy as cp  # fallback to NumPy API
    GPU_AVAILABLE = False
    print("⚠️ CuPy not found — using CPU (NumPy) for distance calculations.")


# ============================================================
# 🧩 STARTUP: Ligand Handling, Recentered Subset, and Pocket Extraction
# ============================================================
import subprocess, os
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np

# ------------------------- Helper functions -------------------------
def file_exists(path):
    return os.path.isfile(path) and os.path.getsize(path) > 0

def run_cmd(cmd, cwd="."):
    """Run a shell command and stream output live."""
    process = subprocess.Popen(cmd, shell=True, cwd=cwd,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               text=True, universal_newlines=True)
    for line in process.stdout:
        print(line, end="")
    process.wait()
    if process.returncode != 0:
        print(f"⚠️ Command failed (continuing): {cmd}")

# ============================================================
# 🔍 SMART LIGAND DETECTION BLOCK (Auto + User Confirmation)
# ============================================================

def detect_ligand_candidates():
    """Scan directory and topology to find possible ligand candidates."""
    candidates = set()

    # 1) Any file like X_ini.pdb → extract residue code X
    for f in os.listdir("."):
        if f.endswith("_ini.pdb") and len(f.split("_")[0]) <= 4:
            candidates.add(f.split("_")[0].upper())

    # 2) Scan topology for residues NOT standard amino acids or water
    non_standard_res = set()
    try:
        topo = mda.Universe(TOPO_FILE)
        for res in topo.residues:
            if res.resname not in {
                "ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS",
                "ILE","LEU","LYS","MET","PHE","PRO","SER","THR",
                "TRP","TYR","VAL",
                "HOH","WAT","SOL"
            }:
                non_standard_res.add(res.resname.strip().upper())
    except Exception:
        pass

    candidates.update(non_standard_res)
    return sorted(list(candidates))


# ------------------------- Ligand Input -------------------------
if args.ligand:
    ligand_code = args.ligand.strip().upper()
    print(f"🧬 Using ligand from CLI: {ligand_code}")

else:
    candidates = detect_ligand_candidates()

    if len(candidates) == 1:
        auto = candidates[0]
        print(f"\n🔍 Detected ligand candidate: {auto}")
        resp = input(f"Use {auto} as ligand? [Y/n]: ").strip().lower()
        if resp in ("", "y", "yes"):
            ligand_code = auto
        else:
            ligand_code = input("Enter ligand resname manually: ").strip().upper()

    elif len(candidates) > 1:
        print("\n🔍 Multiple ligand candidates detected:")
        for i, c in enumerate(candidates, 1):
            print(f"   {i}) {c}")

        choice = input("Select a ligand or press ENTER to type manually: ").strip()
        if choice.isdigit() and 1 <= int(choice) <= len(candidates):
            ligand_code = candidates[int(choice)-1]
        else:
            ligand_code = input("Enter ligand resname manually: ").strip().upper()

    else:
        # No detectable ligand → manual fallback
        ligand_code = input("❓ Enter ligand resname (e.g., PTC, GDP): ").strip().upper()

    if not ligand_code:
        print("❌ ERROR: No ligand provided and auto-detection failed.")
        exit(1)

print(f"🧬 Final ligand selected: {ligand_code}")

# ===========================================================
# FOLLOW-UP QUESTION: USE THIS NAME OR ENTER A LONG NAME?
# ===========================================================
print("\n📌 You may now choose a full compound name for labeling graphs/output.")
print("   • Press ENTER to use the ligand code above (recommended only for short names).")
print("   • Or type a more descriptive compound name (e.g., Paprotrain, CPD32_Inhibitor).")

compound_input = input("🏷️  Enter compound display name (or press ENTER to use resname): ").strip()

if compound_input:
    compound_name = compound_input
    print(f"🏷️ Using custom compound name for all plots + summaries: {compound_name}")
else:
    compound_name = ligand_code
    print(f"🏷️ Using ligand resname for all plots + summaries: {compound_name}")

# Make sure both variables exist globally for all later steps
print(f"\n🧬 RESNAME for topology selections: {ligand_code}")
print(f"🏷️ Compound name for figures/output: {compound_name}\n")


# ============================================================
# 🧬 OPTIONAL RESIDUE NUMBERING TEMPLATE SELECTION
# ============================================================

if NUMBERING_TEMPLATE:
    if not os.path.exists(NUMBERING_TEMPLATE):
        print(f"⚠️ Provided numbering template not found: {NUMBERING_TEMPLATE}")
        NUMBERING_TEMPLATE = None
    else:
        print(f"🧬 Using residue numbering template from CLI: {NUMBERING_TEMPLATE}")
else:
    NUMBERING_TEMPLATE = select_numbering_template_interactive()




# ------------------------- Ensure index group -------------------------
def ensure_index_has_protein_ligand_group(ligand_code):
    index_file = "index.ndx"
    group_name = f"Protein_{ligand_code}"

    if not os.path.exists(index_file):
        raise SystemExit("❌ index.ndx not found. Run MD setup script first.")

    with open(index_file) as f:
        if any(group_name.lower() in line.lower() for line in f if "[" in line):
            print(f"✅ Found existing group [{group_name}] in index.ndx.")
            return

    print(f"⚠️ Group [{group_name}] not found — creating it now...")
    proc = subprocess.run("gmx make_ndx -f md_0_1.tpr -o temp.ndx <<< 'q'",
                          shell=True, capture_output=True, text=True)
    lig_num = None
    for line in proc.stdout.splitlines():
        if ligand_code and ligand_code in line:
            lig_num = line.strip().split()[0]
            break
    lig_num = lig_num or "13"
    run_cmd(f"echo '1 | {lig_num}\\nname 21 {group_name}\\nq' | gmx make_ndx -f md_0_1.tpr -o index.ndx")
    print(f"✅ Created [{group_name}] (1 | {lig_num}) in index.ndx.")

ensure_index_has_protein_ligand_group(ligand_code)





import re

def residue_number(res):
    """Extract numeric residue index from labels like ARG329"""
    match = re.search(r"(\d+)", res)
    return int(match.group(1)) if match else -1



from collections import defaultdict

def build_residue_number_map(template_pdb):
    """
    Build a simple ordered list of residue numbers
    in the order they appear in the template PDB.
    """
    resseqs = []
    last_seen = None

    # Only consider ATOM records (protein atoms). This prevents HETATM
    # entries (ligands, solvent, ions) in the template from shifting
    # the residue ordering used for remapping.
    with open(template_pdb) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue

            # columns 23-26 are residue sequence number in PDB (1-indexed field width 4)
            try:
                resseq = int(line[22:26])
            except ValueError:
                # skip malformed lines
                continue

            if resseq != last_seen:
                resseqs.append(resseq)
                last_seen = resseq

    return resseqs



def apply_template_numbering(pdb_in, pdb_out, resseqs):
    """
    Apply residue numbers sequentially to pdb_in
    using the order defined in resseqs.
    """
    res_idx = -1
    last_resseq = None

    with open(pdb_in) as fin, open(pdb_out, "w") as fout:
        for line in fin:
            # Leave non-coordinate lines intact
            if not line.startswith(("ATOM", "HETATM")):
                fout.write(line)
                continue

            # Do NOT consume template residue numbers for HETATM records.
            # HETATM (ligand/ion) entries should be preserved as-is so they
            # don't shift the protein residue mapping.
            if line.startswith("HETATM"):
                fout.write(line)
                continue

            # Only ATOM records are remapped (protein atoms)
            try:
                cur_resseq = int(line[22:26])
            except ValueError:
                fout.write(line)
                continue

            if cur_resseq != last_resseq:
                res_idx += 1
                if res_idx >= len(resseqs):
                    raise RuntimeError(
                        f"❌ Template has fewer residues ({len(resseqs)}) "
                        f"than trajectory ({res_idx+1})"
                    )
                new_resseq = resseqs[res_idx]
                last_resseq = cur_resseq

            fout.write(
                line[:22] + f"{new_resseq:4d}" + line[26:]
            )





# ============================================================
# STEP A+B — PREPROCESS FINAL TRAJECTORY + POCKET EXTRACTION
# (with correct ATOM → HETATM ligand rewriting)
# ============================================================

phase("STEP A+B — Generate Final_Trajectory.* and binding_pocket_only.*")

if ligand_code is None:
    raise SystemExit("❌ Ligand code not provided and auto-detection not permitted.")

lig = ligand_code.upper()
group = f"Protein_{lig}"

required = ["index.ndx", "md_0_1.tpr", "md_0_1.xtc"]
for f in required:
    if not os.path.exists(f):
        raise SystemExit(f"❌ Missing required file: {f}")


# ------------------------------------------------------------
# Helper — rewrite ligand to HETATM
# ------------------------------------------------------------
def convert_ligand_to_hetatm(pdb_in, pdb_out, ligand_resname):
    ligand_resname = ligand_resname.upper()
    with open(pdb_in, "r") as fin, open(pdb_out, "w") as fout:
        for line in fin:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                resn = line[17:20].strip().upper()
                if resn == ligand_resname:
                    fout.write("HETATM" + line[6:])
                else:
                    fout.write(line)
            else:
                fout.write(line)
    print(f"🔧 Rewrote ligand {ligand_resname} as HETATM → {pdb_out}")


# ============================================================
# STEP A — Final_Trajectory.xtc / Final_Trajectory.pdb
# ============================================================
if not (file_exists("Final_Trajectory.xtc") and file_exists("Final_Trajectory.pdb")):

    print("🎞️ Creating Final_Trajectory.xtc/pdb ...")

    run_cmd(
        f"printf '{group}\\n{group}\\n' | "
        f"gmx trjconv -s md_0_1.tpr -f md_0_1.xtc "
        f"-o Final_Trajectory.xtc -n index.ndx -pbc mol -center"
    )

    run_cmd(
        f"printf '{group}\\n' | "
        f"gmx trjconv -s md_0_1.tpr -f Final_Trajectory.xtc "
        f"-o Final_Trajectory_RAW.pdb -n index.ndx -dump 0"
    )

    convert_ligand_to_hetatm(
        "Final_Trajectory_RAW.pdb",
        "Final_Trajectory_TMP.pdb",
        lig
    )

    if NUMBERING_TEMPLATE and os.path.exists(NUMBERING_TEMPLATE):
        resmap = build_residue_number_map(NUMBERING_TEMPLATE)
        apply_template_numbering(
            "Final_Trajectory_TMP.pdb",
            "Final_Trajectory.pdb",
            resmap
        )
    else:
        os.rename("Final_Trajectory_TMP.pdb", "Final_Trajectory.pdb")

    os.remove("Final_Trajectory_RAW.pdb")

    print("✅ Final_Trajectory.xtc/pdb written.\n")
    print("🧬 First 5 residues after renumbering:")
    with open("Final_Trajectory.pdb") as f:
        seen = set()
        for line in f:
            if line.startswith("ATOM"):
                key = line[22:26]
                if key not in seen:
                    print("   ", line[17:26])
                    seen.add(key)
                if len(seen) == 5:
                    break


else:
    print("⏭️ Final_Trajectory.* already exists — skipping.")



# ============================================================
# STEP B — binding_pocket_only.xtc / .pdb
# ============================================================
if not (file_exists("binding_pocket_only.xtc") and file_exists("binding_pocket_only.pdb")):

    print("🔬 Extracting binding pocket using MDTraj...")

    traj = md.load("Final_Trajectory.xtc", top="Final_Trajectory.pdb")

    # IMPORTANT: quote resname (needed for digit-leading ligands like 4PT)
    lig_sel = "resname == " + repr(lig)
    lig_atoms = traj.topology.select(lig_sel)

    prot_atoms = traj.topology.select("protein")


    if len(lig_atoms) == 0:
        raise SystemExit(f"❌ Ligand {lig} not found in Final_Trajectory.pdb.")

    neighbors = md.compute_neighbors(
        traj,
        cutoff=0.5,   # 0.5 nm = 5 Å
        query_indices=lig_atoms,
        haystack_indices=prot_atoms
    )

    if len(neighbors) == 0:
        raise SystemExit("❌ No protein atoms within 5 Å of ligand.")

    pocket_atoms = np.unique(np.concatenate(neighbors))

    pocket_res = sorted(
        {traj.topology.atom(idx).residue.index for idx in pocket_atoms}
    )

    print(f"🎯 Binding pocket contains {len(pocket_res)} residues.")

    sel = traj.topology.select(
        "resid " + " ".join(map(str, pocket_res)) + " or " + ("resname == " + repr(lig))
    )


    pocket = traj.atom_slice(sel)

    pocket.save("binding_pocket_only.xtc")
    pocket[0].save("binding_pocket_only_RAW.pdb")

    # rewrite ligand → HETATM
    convert_ligand_to_hetatm("binding_pocket_only_RAW.pdb",
                             "binding_pocket_only.pdb", lig)
    os.remove("binding_pocket_only_RAW.pdb")

    print("✅ binding_pocket_only.xtc/pdb written.\n")

else:
    print("⏭️ binding_pocket_only.* already exists — skipping.")

phase("STEP A+B complete")






# ============================================================
# STEP C — Determine Which Chains Interact with the Ligand
# ============================================================
phase("STEP C — Identify Binding Chains for Downstream RMSD")

print("\n🔍 Identifying which protein chains appear in the final binding pocket...")

# We now rely ONLY on binding_pocket_only.pdb from the new STEP B
u_pocket = mda.Universe("binding_pocket_only.pdb")

# Extract ligand
lig_atomgroup = u_pocket.select_atoms(f"resname {ligand_code}")
if lig_atomgroup.n_atoms == 0:
    # Debug: show what resnames MDAnalysis thinks exist
    print("DEBUG pocket resnames:", sorted(set(u_pocket.residues.resnames)))
    raise SystemExit(f"❌ Could not find ligand '{ligand_code}' in binding_pocket_only.pdb")


# All protein residues in the pocket
protein_residues = u_pocket.select_atoms("protein").residues

# Load atomIndex.txt to define chain ranges
chain_map = {}
if os.path.exists("atomIndex.txt"):
    with open("atomIndex.txt") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            name, rng = line.split()
            start, end = [int(x) - 1 for x in rng.split("-")]
            chain_map[name] = (start, end)
else:
    print("⚠️ atomIndex.txt not found — treating entire protein as one chain")
    chain_map["Protein"] = (protein_residues.atoms[0].index,
                            protein_residues.atoms[-1].index)

# Determine active (binding) chains
active_chains = {}

for cname, (start, end) in chain_map.items():
    # Does any residue in the pocket fall within this chain range?
    overlapping = [
        r for r in protein_residues
        if start <= r.atoms[0].index <= end
    ]
    if overlapping:
        active_chains[cname] = (start, end)
        print(f"✅ {cname} participates in ligand binding ({len(overlapping)} pocket residues).")
    else:
        print(f"⚪ {cname} does not interact with ligand.")

if not active_chains:
    print("⚠️ No protein chains interact with ligand!")
else:
    print("📎 Active binding chains:", ", ".join(active_chains.keys()))





# ============================================================
# ✅ Ready for downstream RMSD/RMSF and contact analysis
# ============================================================


# ============================================================
# ✅ Ready for downstream RMSD/RMSF and contact analysis
# ============================================================
TOPO_FILE = "Final_Trajectory.pdb"
TRAJ_FILE = "Final_Trajectory.xtc"



# ------------------------- Helpers -------------------------
def savefig(path):
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()


def infer_element(atom_name: str):
    c = atom_name.strip()[:1]
    return c if c in "ONCHSP" else None

# ============================================================
# Utilities for Checkpointing and RMSD Compatibility
# ============================================================

def checkpoint(path):
    with open(path, "w") as f:
        f.write("done")

def done_checkpoint(path):
    return os.path.exists(path)

def safe_rmsd(traj, ref, frame, atom_indices):
    """Try modern MDTraj rst syntax; fall back when unsupported."""
    try:
        return md.rmsd(traj, ref, frame, atom_indices=atom_indices, parallel=True)
    except TypeError:
        return md.rmsd(traj, ref, frame, atom_indices=atom_indices)

def align_chunk(args):
        traj, bb, idxs = args
        sub = traj.slice(idxs)
        sub.superpose(traj, 0, atom_indices=bb)
        return sub

# ============================================================
# GLOBAL TRAJECTORY LOAD + ALIGNMENT + LIGAND DETECTION
# (runs ONCE, reused by D1/D2/D3/D4/D5)
# ============================================================

_traj_cache = {
    "traj": None,
    "bb": None,
    "ligand_atoms": None,
    "lig_resname": None,
}

def get_aligned_traj_and_ligand():
    global ligand_code

    if _traj_cache["traj"] is not None:
        return (_traj_cache["traj"],
                _traj_cache["bb"],
                _traj_cache["ligand_atoms"],
                _traj_cache["lig_resname"])

    print("\n🌍 Loading trajectory & performing backbone alignment…")
    traj = md.load(TRAJ_FILE, top=TOPO_FILE)
    bb = traj.topology.select("backbone")

    # -------------- FIXED LIGAND LOGIC --------------
    if ligand_code is None:
        raise RuntimeError("❌ No ligand resname provided and auto-detection disabled.")

    lig_resname = ligand_code.upper()

    atoms = traj.topology.select("resname == " + repr(lig_resname))
    if len(atoms) == 0:
        raise RuntimeError(
            f"❌ Ligand '{lig_resname}' not found in trajectory.\n"
            "Check Final_Trajectory.pdb for residue naming."
        )

    ligand_atoms = atoms

    print(f"🔍 Using user-provided ligand: {lig_resname} ({len(ligand_atoms)} atoms)")
    # -------------------------------------------------

    _traj_cache.update({
        "traj": traj,
        "bb": bb,
        "ligand_atoms": ligand_atoms,
        "lig_resname": lig_resname,
    })

    return traj, bb, ligand_atoms, lig_resname


# ============================================================
# STEP D1 — GLOBAL PROTEIN RMSD
# ============================================================

# ------------------ STEP D1 replacement (Protein RMSD) ------------------
if not done_checkpoint("chk_D1.txt"):
    phase("STEP D1 — Global Protein RMSD")

    traj, bb, ligand_atoms, lig_resname = get_aligned_traj_and_ligand()

    print(f"📈 Computing protein backbone RMSD ({traj.n_frames} frames)…")
    prot_rmsd = safe_rmsd(traj, traj, 0, atom_indices=bb) * 10.0

    df = pd.DataFrame({"Time_ps": traj.time, "Protein_RMSD": prot_rmsd})
    df.to_csv(os.path.join(OUTPUT_DIR, "Protein_RMSD.csv"), index=False)

    plt.figure(figsize=(9,5))
    plt.plot(traj.time/1000.0, prot_rmsd, lw=1.8, color=MIAMI_ORANGE)   # use Miami orange for protein
    plt.xlabel("Time (ns)")
    plt.ylabel("RMSD (Å)")
    plt.title("Backbone Residue RMSD")
    savefig(os.path.join(OUTPUT_DIR, "Protein_RMSD.png"))

    checkpoint("chk_D1.txt")
    print("✅ STEP D1 complete.")
else:
    print("⏭️ STEP D1 skipped (checkpoint found).")


# ============================================================
# STEP D2 — LIGAND RMSD & RMSF
# ============================================================

# ------------------ STEP D2 replacement (Ligand RMSD & RMSF) ------------------
if not done_checkpoint("chk_D2.txt"):
    phase("STEP D2 — Ligand RMSD & RMSF")

    traj, bb, ligand_atoms, lig_resname = get_aligned_traj_and_ligand()

    print("📈 Ligand RMSD…")
    lig_rmsd = safe_rmsd(traj, traj, 0, atom_indices=ligand_atoms) * 10.0

    print("📈 Ligand RMSF…")
    mean_pos = traj.xyz[:, ligand_atoms, :].mean(axis=0)
    rmsf = np.sqrt(
        ((traj.xyz[:, ligand_atoms, :] - mean_pos)**2).sum(axis=2).mean(axis=0)
    ) * 10.0

    pd.DataFrame({"Time_ps": traj.time, "Ligand_RMSD": lig_rmsd}).to_csv(
        os.path.join(OUTPUT_DIR, f"{compound_name}_Ligand_RMSD.csv"), index=False
    )
    pd.DataFrame({"Atom": ligand_atoms, "RMSF": rmsf}).to_csv(
        os.path.join(OUTPUT_DIR, f"{compound_name}_Ligand_RMSF.csv"), index=False
    )

    plt.figure(figsize=(9,5))
    plt.plot(traj.time/1000.0, lig_rmsd, color=MIAMI_GREEN, lw=1.6)   # ligand in Miami green
    plt.title(f"{compound_name} Ligand RMSD")
    plt.xlabel("Time (ns)")
    plt.ylabel("RMSD (Å)")
    savefig(os.path.join(OUTPUT_DIR, f"{lig_resname}_Ligand_RMSD.png"))

    plt.figure(figsize=(8,4))
    plt.plot(rmsf, color=MIAMI_GREEN, lw=1.6)   # ligand RMSF in Miami green
    plt.title(f"{compound_name} Ligand RMSF")
    plt.xlabel("Ligand Atom Index")
    plt.ylabel("RMSF (Å)")
    savefig(os.path.join(OUTPUT_DIR, f"{lig_resname}_Ligand_RMSF.png"))

    checkpoint("chk_D2.txt")
    print("✅ STEP D2 complete.")
else:
    print("⏭️ STEP D2 skipped (checkpoint found).")


# ============================================================
# STEP D3 — PER-CHAIN RMSD / RMSF (ATOM + RESIDUE + PERSISTENCE)
# ============================================================

if not done_checkpoint("chk_D3.txt"):
    phase("STEP D3 — Per-Chain RMSD / RMSF")

    traj, bb, ligand_atoms, lig_resname = get_aligned_traj_and_ligand()
    n_frames = traj.n_frames

    # --------------------------------------------------------
    # Load chain map
    # --------------------------------------------------------
    chain_map = {}
    if os.path.exists("atomIndex.txt"):
        with open("atomIndex.txt") as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                name, r = line.split()
                s, e = [int(x) - 1 for x in r.split("-")]
                chain_map[name] = (s, e)

    # --------------------------------------------------------
    # Loop over chains
    # --------------------------------------------------------
    for chain, (start, end) in chain_map.items():
        sel = np.arange(start, end + 1)
        if sel.max() >= traj.n_atoms:
            print(f"⚠️ {chain}: invalid atom range {start}-{end}")
            continue

        # =========================
        # RMSD (Å)
        # =========================
        rmsd = safe_rmsd(traj, traj, 0, atom_indices=sel) * 10

        pd.DataFrame({
            "Time_ps": traj.time,
            "RMSD_Å": rmsd
        }).to_csv(
            os.path.join(OUTPUT_DIR, f"RMSD_{chain}.csv"),
            index=False
        )

        plt.figure(figsize=(9, 5))
        plt.plot(traj.time / 1000.0, rmsd, color=MIAMI_ORANGE, lw=1.8)   # protein chain in Miami orange
        plt.title(f"{chain} RMSD")
        plt.xlabel("Time (ns)")
        plt.ylabel("RMSD (Å)")
        savefig(os.path.join(OUTPUT_DIR, f"RMSD_{chain}.png"))
        plt.close()

        # =========================
        # ATOM-LEVEL RMSF (Å)
        # =========================
        mean_pos = traj.xyz[:, sel, :].mean(axis=0)
        rmsf_atom = np.sqrt(
            ((traj.xyz[:, sel, :] - mean_pos) ** 2).sum(axis=2).mean(axis=0)
        ) * 10

        pd.DataFrame({
            "Atom_Index": sel,
            "RMSF_Å": rmsf_atom
        }).to_csv(
            os.path.join(OUTPUT_DIR, f"RMSF_{chain}.csv"),
            index=False
        )

        plt.figure(figsize=(9, 4))
        plt.plot(rmsf_atom, color=MIAMI_ORANGE, lw=1.6)   # protein RMSF in Miami orange
        plt.title(f"{chain} RMSF (per atom)")
        plt.xlabel("Atom Index")
        plt.ylabel("RMSF (Å)")
        savefig(os.path.join(OUTPUT_DIR, f"RMSF_{chain}.png"))
        plt.close()

        # =========================
        # RESIDUE-LEVEL RMSF (Cα)
        # =========================
        ca_indices = [
            atom.index for atom in traj.topology.atoms
            if atom.name == "CA" and start <= atom.index <= end
        ]

        if not ca_indices:
            print(f"⚠️ {chain}: no CA atoms found")
            continue

        res_ids = [
            traj.topology.atom(i).residue.resSeq
            for i in ca_indices
        ]

        mean_ca = traj.xyz[:, ca_indices, :].mean(axis=0)
        rmsf_res = np.sqrt(
            ((traj.xyz[:, ca_indices, :] - mean_ca) ** 2)
            .sum(axis=2)
            .mean(axis=0)
        ) * 10

        pd.DataFrame({
            "Residue": res_ids,
            "RMSF_Å": rmsf_res
        }).to_csv(
            os.path.join(OUTPUT_DIR, f"RMSF_{chain}_per_residue.csv"),
            index=False
        )

        # =========================
        # CONTACT PERSISTENCE
        # =========================
        INTERACTION_CUTOFF = args.contact_cutoff  # Å
        contact_counts = {r: 0 for r in res_ids}

        for frame_i in range(n_frames):
            lig_xyz = traj.xyz[frame_i, ligand_atoms, :]
            ca_xyz = traj.xyz[frame_i, ca_indices, :]

            d = np.linalg.norm(
                lig_xyz[:, None, :] - ca_xyz[None, :, :],
                axis=2
            )

            min_d = d.min(axis=0) * 10  # Å

            for r, dist in zip(res_ids, min_d):
                if dist < INTERACTION_CUTOFF:
                    contact_counts[r] += 1

        contact_frac = {
            r: contact_counts[r] / n_frames
            for r in res_ids
        }

        persistent_residues = {
            r for r, frac in contact_frac.items()
            if frac >= MIN_CONTACT_FRAC
        }

        # =========================
        # BUILD RMSF LINE WITH TRUE GAPS
        # =========================
        res_plot = []
        rmsf_plot = []

        for i, r in enumerate(res_ids):
            if i > 0 and r - res_ids[i - 1] > 1:
                res_plot.append(np.nan)
                rmsf_plot.append(np.nan)
            res_plot.append(r)
            rmsf_plot.append(rmsf_res[i])

        res_plot = np.array(res_plot)
        rmsf_plot = np.array(rmsf_plot)

        from matplotlib.lines import Line2D

        # =========================
        # RESIDUE RMSF PLOT
        # =========================
        plt.figure(figsize=(10, 4))

        # RMSF curve (broken at gaps)
        plt.plot(res_plot, rmsf_plot, color=MIAMI_ORANGE, lw=1.8)  # residue RMSF in orange

        # Persistent contact markers (hashed vertical lines)
        for r, rmsf_val in zip(res_ids, rmsf_res):
            if r in persistent_residues:
                plt.vlines(
                    x=r,
                    ymin=0,
                    ymax=rmsf_val,
                    color=MIAMI_GREEN,   # ligand contact marker color
                    linestyles="dashed",
                    alpha=0.6,
                    lw=1.2
                )

        plt.title(f"{chain} RMSF per Residue")
        plt.xlabel("Residue Number")
        plt.ylabel("RMSF (Å)")

        # Legend (contact-only)
        legend_elements = [
            Line2D(
                [0], [0],
                color=MIAMI_GREEN,
                lw=1.5,
                linestyle="dashed",
                label=f"Ligand contact ≥ {MIN_CONTACT_FRAC:.0%} frames"
            )
        ]

        plt.legend(
            handles=legend_elements,
            loc="upper right",
            frameon=False
        )

        plt.tight_layout()
        savefig(
            os.path.join(
                OUTPUT_DIR,
                f"RMSF_{chain}_per_residue_contacts.png"
            )
        )
        plt.close()

    checkpoint("chk_D3.txt")
    print("✅ STEP D3 complete.")

else:
    print("⏭️ STEP D3 skipped (checkpoint found).")



# ============================================================
# STEP D4 — BOUND COMPLEX RMSD (Dual-Axis Protein + Ligand)
# ============================================================

# ------------------ STEP D4 replacement (Bound Complex RMSD) ------------------
if not done_checkpoint("chk_D4.txt"):
    phase("STEP D4 — Bound Complex RMSD (Dual-Axis)")

    traj, bb, ligament_atoms, lig_resname = get_aligned_traj_and_ligand()

    # Protein RMSD
    prot_rmsd = safe_rmsd(traj, traj, 0, atom_indices=bb) * 10.0  # Å

    # Ligand RMSD
    lig_rmsd = safe_rmsd(traj, traj, 0, atom_indices=ligand_atoms) * 10.0  # Å

    time_ns = traj.time / 1000.0

    df = pd.DataFrame({
        "Time_ns": time_ns,
        "Protein_RMSD": prot_rmsd,
        "Ligand_RMSD": lig_rmsd
    })

    df.to_csv(os.path.join(OUTPUT_DIR,
                           f"{lig_resname}_Complex_RMSD_Overlay.csv"),
                           index=False)

    # Dual-Axis Plot — protein = orange, ligand = green
    fig, ax1 = plt.subplots(figsize=(10, 6))

    ax1.plot(time_ns, prot_rmsd, color=MIAMI_ORANGE, label="Protein RMSD")
    ax1.set_ylabel("Protein RMSD (Å)", color=MIAMI_ORANGE)
    ax1.tick_params(axis="y", labelcolor=MIAMI_ORANGE)

    ax2 = ax1.twinx()
    ax2.plot(time_ns, lig_rmsd, color=MIAMI_GREEN, label="Ligand RMSD")
    ax2.set_ylabel("Ligand RMSD (Å)", color=MIAMI_GREEN)
    ax2.tick_params(axis="y", labelcolor=MIAMI_GREEN)

    ax1.set_title(f"RMSD of {compound_name} bound with {chain}")
    ax1.set_xlabel("Time (ns)")

    fig.tight_layout()
    savefig(os.path.join(OUTPUT_DIR,
                         f"{lig_resname}_Complex_RMSD_Overlay.png"))

    checkpoint("chk_D4.txt")
    print("✅ STEP D4 complete.")
else:
    print("⏭️ STEP D4 skipped (checkpoint found).")






# ============================================================
# STEP D4B — Radius of Gyration (Protein / Ligand / Complex)
# ============================================================

if not done_checkpoint("chk_D4B.txt"):
    phase("STEP D4B — Radius of Gyration (Protein, Ligand, Complex)")

    traj, bb, ligand_atoms, lig_resname = get_aligned_traj_and_ligand()

    print("📐 Computing radius of gyration for protein, ligand, and full complex…")

    # ------------------------------
    # Slice trajectories (API-safe)
    # ------------------------------
    traj_protein = traj.atom_slice(bb)
    traj_ligand  = traj.atom_slice(ligand_atoms)
    traj_complex = traj  # full system

    # ------------------------------
    # Compute Rg (nm → Å)
    # ------------------------------
    rg_protein = md.compute_rg(traj_protein) * 10.0
    rg_ligand  = md.compute_rg(traj_ligand)  * 10.0
    rg_complex = md.compute_rg(traj_complex) * 10.0

    time_ns = traj.time / 1000.0

    # ------------------------------
    # Save CSV
    # ------------------------------
    df = pd.DataFrame({
        "Time_ns": time_ns,
        "Rg_Protein_Å": rg_protein,
        "Rg_Ligand_Å": rg_ligand,
        "Rg_Complex_Å": rg_complex
    })

    df.to_csv(
        os.path.join(OUTPUT_DIR, "Radius_of_Gyration_Protein_Ligand_Complex.csv"),
        index=False
    )

    # ------------------------------
    # Overlay plot
    # ------------------------------
    plt.figure(figsize=(10,6))

    plt.plot(time_ns, rg_protein, lw=2.0, label="Protein Rg", color="#1E88E5")
    plt.plot(time_ns, rg_ligand,  lw=2.0, label="Ligand Rg",  color="#D32F2F")
    plt.plot(time_ns, rg_complex, lw=2.0, label="Complex Rg", color="#2E7D32")

    plt.xlabel("Time (ns)")
    plt.ylabel("Radius of Gyration (Å)")
    plt.title("Radius of Gyration — Protein, Ligand, and Complex")
    plt.legend(frameon=False)
    plt.grid(alpha=0.25)

    savefig(os.path.join(
        OUTPUT_DIR,
        "Radius_of_Gyration_Overlay.png"
    ))

    checkpoint("chk_D4B.txt")
    print("✅ STEP D4B complete.")

else:
    print("⏭️ STEP D4B skipped (checkpoint found).")







# ============================================================
# STEP D5 — FULL SYSTEM RMSD
# ============================================================

# ------------------ STEP D5 replacement (Full System RMSD) ------------------
if not done_checkpoint("chk_D5.txt"):
    phase("STEP D5 — Full Complex RMSD")

    traj, bb, ligand_atoms, lig_resname = get_aligned_traj_and_ligand()

    full_rmsd = safe_rmsd(traj, traj, 0, atom_indices=np.arange(traj.n_atoms))*10
    pd.DataFrame({"Time_ps": traj.time, "RMSD_Å": full_rmsd}).to_csv(
        os.path.join(OUTPUT_DIR, "FullComplex_RMSD.csv"), index=False
    )

    plt.figure(figsize=(9,5))
    plt.plot(traj.time/1000., full_rmsd, color=MIAMI_ORANGE, lw=1.6)  # full-system (protein-majority) in orange
    plt.title("Backbone RMSD of ALL Chains in Complex")
    plt.ylabel("RMSD (Å)")
    plt.xlabel("Time (ns)")

    savefig(os.path.join(OUTPUT_DIR, "FullComplex_RMSD.png"))

    checkpoint("chk_D5.txt")
    print("✅ STEP D5 complete.")
else:
    print("⏭️ STEP D5 skipped (checkpoint found).")


# ============================================================
# STEP D6 — Secondary Structure Evolution (DSSP)
# ============================================================

if not done_checkpoint("chk_D6.txt"):
    phase("STEP D6 — Secondary Structure Evolution (DSSP)")

    traj, bb, ligand_atoms, lig_resname = get_aligned_traj_and_ligand()

    print("📐 Computing DSSP for all frames...")
    dssp = md.compute_dssp(traj)   # shape = (n_frames, n_residues)

    # ============================================================
    # 1) SSE Encoding (Helix / Sheet / Coil)
    # ============================================================

    def encode_sse(c):
        """Convert DSSP letter to numeric SSE class."""
        if c in {"H", "G", "I"}:   # Helical family
            return 2               # Helix
        if c in {"E", "B"}:        # Beta family
            return 1               # Sheet
        return 0                   # Everything else = coil

    int_dssp = np.vectorize(encode_sse)(dssp)

    # Save raw and encoded SSE
    pd.DataFrame(dssp).to_csv(
        os.path.join(OUTPUT_DIR, "DSSP_raw.csv"), index=False
    )
    pd.DataFrame(int_dssp).to_csv(
        os.path.join(OUTPUT_DIR, "DSSP_encoded_numeric.csv"), index=False
    )

    # ============================================================
    # 2) HEATMAP — SSE Timeline (Residue × Time)
    # ============================================================

    from matplotlib.colors import ListedColormap, BoundaryNorm

    print("🎨 Plotting DSSP heatmap...")

    # Coil = yellow, Sheet = blue, Helix = red
    cmap = ListedColormap(["#FFD700", "#1E88E5", "#D32F2F"])
    bounds = [-0.5, 0.5, 1.5, 2.5]          # 0,1,2 centered
    norm = BoundaryNorm(bounds, cmap.N)

    fig, ax = plt.subplots(figsize=(16, 6))

    im = ax.imshow(
        int_dssp.T,
        aspect="auto",
        cmap=cmap,
        norm=norm,
        origin="upper"  # keep residue 0 at the top
    )

    ax.set_xlabel("Frame")
    ax.set_ylabel("Residue Index")
    ax.set_title("Secondary Structure Evolution (Coil / Sheet / Helix)")

    # Colorbar with fixed ticks/labels
    cbar = fig.colorbar(im, ax=ax, boundaries=bounds, ticks=[0, 1, 2])
    cbar.ax.set_yticklabels(["Coil", "Sheet", "Helix"])
    cbar.set_label("SSE", rotation=270, labelpad=15)

    fig.tight_layout()
    savefig(os.path.join(OUTPUT_DIR, "DSSP_Heatmap.png"))

    # ============================================================
    # 3) FRACTIONS — Helix / Sheet / Coil Over Time
    # ============================================================

    print("📊 Computing SSE fractions per frame...")

    helix_codes = {"H", "G", "I"}
    sheet_codes = {"E", "B"}
    coil_codes  = {"S", "T", "C", " "}

    helix_frac = []
    sheet_frac = []
    coil_frac  = []

    for row in dssp:
        n = len(row)
        helix_frac.append(sum(c in helix_codes for c in row) / n)
        sheet_frac.append(sum(c in sheet_codes for c in row) / n)
        coil_frac.append(sum(c in coil_codes  for c in row) / n)

    df_sse = pd.DataFrame({
        "Time_ns": traj.time / 1000.0,
        "Helix": helix_frac,
        "Sheet": sheet_frac,
        "Coil":  coil_frac,
    })

    df_sse.to_csv(os.path.join(OUTPUT_DIR, "SSE_Fractions.csv"), index=False)

    print("🎨 Plotting SSE fraction timeline...")
    plt.figure(figsize=(12, 6))
    plt.plot(df_sse["Time_ns"], df_sse["Helix"], label="Helix", lw=2, color="#D32F2F")
    plt.plot(df_sse["Time_ns"], df_sse["Sheet"], label="Sheet", lw=2, color="#1E88E5")
    plt.plot(df_sse["Time_ns"], df_sse["Coil"],  label="Coil",  lw=2, color="#FFD700")

    plt.xlabel("Time (ns)")
    plt.ylabel("Fraction of Protein")
    plt.title("Secondary Structure Fractions Over Time")
    plt.legend()

    savefig(os.path.join(OUTPUT_DIR, "SSE_Fractions.png"))

    checkpoint("chk_D6.txt")
    print("✅ STEP D6 complete — DSSP analysis finished.\n")

else:
    print("⏭️ STEP D6 skipped (checkpoint found).")



# ============================================================
# STEP D6B — Per-Chain Secondary Structure Evolution (DSSP)
# ============================================================

import string

def chain_index_to_letter(idx):
    letters = string.ascii_uppercase
    if idx < len(letters):
        return letters[idx]
    else:
        return f"{letters[idx // 26 - 1]}{letters[idx % 26]}"



if not done_checkpoint("chk_D6B.txt"):
    phase("STEP D6B — Per-Chain Secondary Structure Evolution")

    traj, bb, ligand_atoms, lig_resname = get_aligned_traj_and_ligand()

    print("📐 Loading DSSP (reusing D6 results if available)...")
    dssp = md.compute_dssp(traj)

    # Map each residue to its chain index
    residue_chains = np.array([res.chain.index for res in traj.topology.residues])
    unique_chains = np.unique(residue_chains)

    from matplotlib.colors import ListedColormap, BoundaryNorm

    # Shared encoder
    def encode_sse(c):
        if c in {"H","G","I"}: return 2   # Helix
        if c in {"E","B"}:     return 1   # Sheet
        return 0                          # Coil

    # Shared color mapping
    cmap = ListedColormap(["#FFD700", "#1E88E5", "#D32F2F"])
    bounds = [-0.5, 0.5, 1.5, 2.5]
    norm = BoundaryNorm(bounds, cmap.N)

    for chain_id in unique_chains:
        chain_label = chain_index_to_letter(chain_id)
        print(f"\n🔍 Processing chain {chain_label} ...")


        mask = residue_chains == chain_id
        dssp_chain = dssp[:, mask]          # shape: (frames, residues_in_chain)

        int_sse = np.vectorize(encode_sse)(dssp_chain)

        # Save numeric matrix
        pd.DataFrame(int_sse).to_csv(
            os.path.join(OUTPUT_DIR, f"DSSP_chain_{chain_label}.csv"),
            index=False
        )


        # ---------------------------
        # Heatmap per chain
        # ---------------------------
        fig, ax = plt.subplots(figsize=(16, 6))
        im = ax.imshow(
            int_sse.T,
            aspect="auto",
            cmap=cmap,
            norm=norm,
            origin="upper"
        )

        ax.set_ylabel("Residue Index (Chain-local)")
        ax.set_xlabel("Frame")
        ax.set_title(f"Chain {chain_label} — Secondary Structure Evolution")


        cbar = fig.colorbar(im, ax=ax, boundaries=bounds, ticks=[0, 1, 2])
        cbar.ax.set_yticklabels(["Coil", "Sheet", "Helix"])
        cbar.set_label("SSE", rotation=270, labelpad=15)

        fig.tight_layout()
        savefig(os.path.join(
            OUTPUT_DIR, f"DSSP_chain_{chain_label}_heatmap.png"
        ))


        # ---------------------------
        # Chain-level fractions
        # ---------------------------
        helix_codes = {"H","G","I"}
        sheet_codes = {"E","B"}
        coil_codes  = {"S","T","C"," "}

        helix_frac = []
        sheet_frac = []
        coil_frac = []

        for row in dssp_chain:
            n = len(row)
            helix_frac.append(sum(c in helix_codes for c in row) / n)
            sheet_frac.append(sum(c in sheet_codes for c in row) / n)
            coil_frac.append(sum(c in coil_codes  for c in row) / n)

        df = pd.DataFrame({
            "Time_ns": traj.time/1000,
            "Helix": helix_frac,
            "Sheet": sheet_frac,
            "Coil": coil_frac
        })

        df.to_csv(
            os.path.join(OUTPUT_DIR, f"SSE_chain_{chain_label}_fractions.csv"),
            index=False
        )


        plt.figure(figsize=(10,5))
        plt.plot(df["Time_ns"], df["Helix"], label="Helix", color="#D32F2F")
        plt.plot(df["Time_ns"], df["Sheet"], label="Sheet", color="#1E88E5")
        plt.plot(df["Time_ns"], df["Coil"],  label="Coil",  color="#FFD700")
        plt.title(f"Chain {chain_label} — Secondary Structure Fractions Over Time")
        plt.xlabel("Time (ns)")
        plt.ylabel("Fraction")
        plt.legend()

        savefig(os.path.join(
            OUTPUT_DIR, f"SSE_chain_{chain_label}_fractions.png"
        ))


    checkpoint("chk_D6B.txt")
    print("✅ STEP D6B complete.\n")

else:
    print("⏭️ STEP D6B skipped (checkpoint found).")










# ============================================================
# STEP D6C — Residue-Level SSE Persistence Ranking
# ============================================================

if not done_checkpoint("chk_D6C.txt"):
    phase("STEP D6C — Residue-Level SSE Persistence Ranking")

    traj, bb, ligand_atoms, lig_resname = get_aligned_traj_and_ligand()
    dssp = md.compute_dssp(traj)

    helix_codes = {"H","G","I"}
    sheet_codes = {"E","B"}
    coil_codes  = {"S","T","C"," "}

    helix_persist = []
    sheet_persist = []
    coil_persist = []

    for col in range(dssp.shape[1]):
        s = dssp[:, col]
        helix_persist.append(sum(c in helix_codes for c in s) / len(s))
        sheet_persist.append(sum(c in sheet_codes for c in s) / len(s))
        coil_persist.append(sum(c in coil_codes for c in s) / len(s))

    df = pd.DataFrame({
        "Residue_Index": np.arange(dssp.shape[1]),
        "Helix_Persistence": helix_persist,
        "Sheet_Persistence": sheet_persist,
        "Coil_Persistence": coil_persist
    })

    df.to_csv(os.path.join(OUTPUT_DIR, "Residue_SSE_Persistence.csv"), index=False)

    # Plot ranking (Helix)
    plt.figure(figsize=(14,6))
    plt.bar(df["Residue_Index"], df["Helix_Persistence"], color="red")
    plt.title("Residue Helix Persistence")
    plt.xlabel("Residue Index")
    plt.ylabel("Fraction of Time in Helix")
    savefig(os.path.join(OUTPUT_DIR, "Residue_Helix_Persistence.png"))

    # Plot ranking (Sheet)
    plt.figure(figsize=(14,6))
    plt.bar(df["Residue_Index"], df["Sheet_Persistence"], color="blue")
    plt.title(" Residue Sheet Persistence")
    plt.xlabel("Residue Index")
    plt.ylabel("Fraction of Time in Sheet")
    savefig(os.path.join(OUTPUT_DIR, "Residue_Sheet_Persistence.png"))

    # Flexibility = 1 - max(structured)
    df["Flexibility"] = 1 - df[[
        "Helix_Persistence",
        "Sheet_Persistence"
    ]].max(axis=1)

    plt.figure(figsize=(14,6))
    plt.bar(df["Residue_Index"], df["Flexibility"], color="green")
    plt.title(" Residue Flexibility Index (1 = Highly Flexible)")
    plt.xlabel("Residue Index")
    plt.ylabel("Flexibility")
    savefig(os.path.join(OUTPUT_DIR, "Residue_Flexibility.png"))

    checkpoint("chk_D6C.txt")
    print("✅ STEP D6C complete.\n")

else:
    print("⏭️ STEP D6C skipped (checkpoint found).")







# ============================================================
# STEP E — FULL BINDING-POCKET ANALYSIS WITH CSV EXPORTS
# ============================================================
_, _, _, lig_resname = get_aligned_traj_and_ligand()
phase("STEP E — Binding Pocket Analysis")

u = mda.Universe("binding_pocket_only.pdb", "binding_pocket_only.xtc")
pocket_protein = u.select_atoms("protein")
pocket_ligand  = u.select_atoms(f"resname {lig_resname}")
total_frames = len(u.trajectory)

CONTACT_CUTOFF = float(args.contact_cutoff)
MIN_CONTACT_FRAC = float(args.min_contact_frac)

# Storage tables
interaction_data = []          # (Time_ps, Residue)
framewise_interaction_rows = []  # (Time_ns, Residue, Type)

# ============================================================
# STEP E1 — CONTACT DETECTION (RESIDUE LEVEL)
# ============================================================

phase("STEP E1 — Pocket Contact Detection")

interaction_data = []

try:
    if pocket_ligand.n_atoms == 0:
        print(f"⚠️ No ligand atoms found for resname {lig_resname}. Skipping STEP E1/E1A/E1B.")
    elif pocket_protein.n_atoms == 0:
        print("⚠️ No protein atoms found in binding pocket selection. Skipping STEP E1/E1A/E1B.")
    else:
        for ts in tqdm(u.trajectory, desc="🔍 Detecting contacts"):
            time_ps = ts.time

            # Compute protein–ligand distances
            dist_mat = distances.distance_array(
                pocket_ligand.positions,
                pocket_protein.positions
            )

            # For each protein residue, find closest ligand atom
            for res in pocket_protein.residues:
                atom_indices = res.atoms.indices

                # Guard against empty residues
                if len(atom_indices) == 0:
                    continue

                dists = dist_mat[:, atom_indices]  # shape = (lig_atoms, res_atoms)

                # Guard against malformed empty distance slices
                if dists.size == 0:
                    continue

                ligand_atom_idx, protein_atom_idx = np.unravel_index(
                    np.argmin(dists),
                    dists.shape
                )

                min_dist = dists[ligand_atom_idx, protein_atom_idx]

                if min_dist < CONTACT_CUTOFF:
                    ligand_atom = pocket_ligand.atoms[ligand_atom_idx]
                    protein_atom = res.atoms[protein_atom_idx]

                    interaction_data.append((
                        time_ps,
                        f"{res.resname}{res.resid}",
                        ligand_atom.name,
                        ligand_atom.index,
                        protein_atom.name,
                        protein_atom.index,
                        float(min_dist)
                    ))

except Exception as e:
    print(f"⚠️ STEP E1 contact detection failed — continuing pipeline.")
    print(f"   ↳ {type(e).__name__}: {e}")


# ============================================================
# DATAFRAME CONSTRUCTION
# ============================================================

if len(interaction_data) == 0:
    print("⚠️ No ligand–protein contacts detected under current cutoff.")
    interaction_df = pd.DataFrame(columns=[
        "Time_ps", "Residue",
        "LigAtomName", "LigAtomIndex",
        "ProtAtomName", "ProtAtomIndex",
        "Distance", "Time_ns", "ResidNum"
    ])
else:
    interaction_df = pd.DataFrame(
        interaction_data,
        columns=[
            "Time_ps", "Residue",
            "LigAtomName", "LigAtomIndex",
            "ProtAtomName", "ProtAtomIndex",
            "Distance"
        ]
    )

    interaction_df["Time_ns"] = interaction_df["Time_ps"] / 1000.0

    # Extract numeric residue number for sorting
    interaction_df["ResidNum"] = (
        interaction_df["Residue"]
        .str.extract(r"(\d+)")
        .astype(int)
    )

# Save full framewise contacts whether empty or not
interaction_df.to_csv(
    os.path.join(OUTPUT_DIR, "AllContacts_Framewise.csv"),
    index=False
)

print(f"📊 interaction_df rows: {len(interaction_df)}")
print(f"📊 unique contacting residues: {interaction_df['Residue'].nunique() if not interaction_df.empty else 0}")


# ============================================================
# FILTER BY CONTACT PERSISTENCE
# ============================================================

if interaction_df.empty:
    keep_res = []
    filtered_df = interaction_df.copy()
else:
    counts = interaction_df["Residue"].value_counts()
    keep_res = counts[counts >= total_frames * MIN_CONTACT_FRAC].index
    filtered_df = interaction_df[interaction_df["Residue"].isin(keep_res)].copy()

filtered_df.to_csv(
    os.path.join(OUTPUT_DIR, "FilteredContacts_Framewise.csv"),
    index=False
)

print(f"📊 residues passing persistence filter ({MIN_CONTACT_FRAC:.1%}): {len(keep_res)}")
print(f"📊 filtered_df rows: {len(filtered_df)}")


# ============================================================
# STEP E1A — CONTACT HEATMAP (NUMERICALLY SORTED)
# ============================================================

phase("STEP E1A — Contact Heatmap")

try:
    if filtered_df.empty:
        print("⚠️ No persistent contacts found — skipping STEP E1A.")
    else:
        pivot = filtered_df.pivot_table(
            index=["ResidNum", "Residue"],
            columns="Time_ns",
            aggfunc="size",
            fill_value=0
        )

        if pivot.empty or pivot.shape[0] == 0 or pivot.shape[1] == 0:
            print("⚠️ Empty pivot table — skipping STEP E1A.")
        else:
            # Sort residues numerically
            pivot = pivot.sort_index(level="ResidNum")

            # Drop numeric index level (keep clean labels)
            pivot.index = pivot.index.get_level_values("Residue")

            plt.figure(figsize=(16, 9))
            sns.heatmap(pivot, cmap="coolwarm")
            plt.title(f"{compound_name}–{chain} Contact Map (Binding Pocket Only)")
            plt.xlabel("Time (ns)")
            plt.ylabel("Residue")

            savefig(
                os.path.join(
                    OUTPUT_DIR,
                    f"{lig_resname}_contact_map_residue_time.png"
                )
            )

            # Optional: save pivot used for plotting
            pivot.to_csv(
                os.path.join(
                    OUTPUT_DIR,
                    f"{lig_resname}_contact_map_residue_time.csv"
                )
            )

            print("✅ STEP E1A complete.")

except Exception as e:
    print("⚠️ STEP E1A failed — continuing pipeline.")
    print(f"   ↳ {type(e).__name__}: {e}")


# ============================================================
# STEP E1B — CONTACT FREQUENCY BARPLOT (MATCHING ORDER)
# ============================================================

phase("STEP E1B — Contact Frequency Barplot")

try:
    if filtered_df.empty:
        print("⚠️ No persistent contacts found — skipping STEP E1B.")
    else:
        freq = (
            filtered_df
            .drop_duplicates(["Residue", "ResidNum", "Time_ns"])
            .groupby(["ResidNum", "Residue"])
            .size()
            .sort_index(level="ResidNum")
        )

        if freq.empty:
            print("⚠️ Frequency table is empty — skipping STEP E1B.")
        else:
            # Clean index for plotting
            freq.index = freq.index.get_level_values("Residue")

            # Convert counts → percent of total frames
            freq_pct = (freq.astype(float) / float(total_frames)) * 100.0

            # Save numeric table
            freq_pct.to_csv(
                os.path.join(
                    OUTPUT_DIR,
                    f"{lig_resname}_contact_frequency_filtered.csv"
                ),
                header=["PercentFrames"]
            )

            plt.figure(figsize=(8, 12))
            ax = freq_pct.plot(kind="barh", color=MIAMI_ORANGE)
            ax.set_title(f"{compound_name}-{chain} Contact Frequency")
            ax.set_xlabel("Percent of frames (%)")
            ax.set_xlim(0, 100)

            # Annotate bars
            for p in ax.patches:
                width = p.get_width()
                if np.isfinite(width):
                    ax.text(
                        width + 0.8,
                        p.get_y() + p.get_height() / 2,
                        f"{width:.1f}%",
                        va="center",
                        fontsize=8
                    )

            savefig(
                os.path.join(
                    OUTPUT_DIR,
                    f"{lig_resname}_contact_frequency_filtered.png"
                )
            )

            print("✅ STEP E1B complete.")

except Exception as e:
    print("⚠️ STEP E1B failed — continuing pipeline.")
    print(f"   ↳ {type(e).__name__}: {e}")

# ============================================================
# STEP E2 — INTERACTION TYPE CLASSIFICATION (SAFE SINGLE THREAD)
# ============================================================

phase("STEP E2 — Classifying Interaction Types")

def classify_frame(ts):
    time_ns = ts.time / 1000.0

    dist_mat = distances.distance_array(
        pocket_ligand.positions, pocket_protein.positions
    )
    rows = []

    for res in pocket_protein.residues:
        rid = f"{res.resname}{res.resid}"
        atom_indices = res.atoms.indices

        # compute per–atom nearest distances   
        dists = dist_mat[:, atom_indices]
        lig_i, prot_j = np.unravel_index(np.argmin(dists), dists.shape)
        min_dist = dists[lig_i, prot_j]

        lig_atom = pocket_ligand.atoms[lig_i]
        prot_atom = res.atoms[prot_j]

        itype = None

        # H-bond detection
        if infer_element(prot_atom.name) in ("O","N") and infer_element(lig_atom.name) in ("O","N"):
            if min_dist < 3.5:
                itype = "H-bond"

        # Hydrophobic
        if itype is None and res.resname in {"PHE","LEU","ILE","VAL","MET","TRP","TYR","ALA"} and min_dist < 4.5:
            itype = "Hydrophobic"

        # Ionic
        if itype is None and res.resname in {"ASP","GLU","ARG","LYS","HIS"} and min_dist < 6.0:
            itype = "Ionic"

        if itype:
            rows.append((
                time_ns, rid, itype,
                lig_atom.name, lig_atom.index,
                prot_atom.name, prot_atom.index,
                float(min_dist)
            ))

    return rows



# Execute (safe, interruptible)
for ts in tqdm(u.trajectory, desc="🔬 Classifying interactions"):
    framewise_interaction_rows.extend(classify_frame(ts))


# Convert & save
int_df = pd.DataFrame(
    framewise_interaction_rows,
    columns=[
        "Time_ns", "Residue", "Type",
        "LigAtomName", "LigAtomIndex",
        "ProtAtomName", "ProtAtomIndex",
        "Distance"
    ]
)

int_df.to_csv(os.path.join(OUTPUT_DIR, "InteractionTypes_Framewise.csv"), index=False)

summary = int_df.groupby(["Residue","Type"]).size().unstack(fill_value=0)
summary["Total"] = summary.sum(axis=1)
summary.to_csv(os.path.join(OUTPUT_DIR, "InteractionTypes_Summary.csv"))

stacked_frac = summary.drop(columns=["Total"], errors="ignore").div(total_frames).fillna(0)

import re

def residue_sort_key(res):
    # Extract trailing residue number (e.g. ARG329 → 329)
    m = re.search(r"(\d+)$", res)
    return int(m.group(1)) if m else 0

# Sort stacked_frac by residue number
stacked_frac = stacked_frac.loc[
    sorted(stacked_frac.index, key=residue_sort_key)
]

# ============================================================
# STEP E2A — Grouped bar plot (comparison view)
# ============================================================

fig, ax = plt.subplots(figsize=(16, 6), dpi=150)

stacked_frac.plot(
    kind="bar",
    ax=ax,
    stacked=False,
    width=0.75,
    color=[COLORS[c] for c in PLOT_TYPES],
    edgecolor="black",
    linewidth=0.4,
)

ax.set_title(
    f"Normalized {chain}–{compound_name} Interaction Frequencies",
    fontsize=16,
    pad=20
)
ax.set_ylabel("Fraction of Frames", fontsize=12)
ax.set_xlabel("Residue", fontsize=12)

# Comparison emphasis
ax.yaxis.grid(True, linestyle="--", alpha=0.4)
ax.set_axisbelow(True)

ax.tick_params(axis="x", rotation=90, labelsize=8)
ax.tick_params(axis="y", labelsize=10)

# Legend ABOVE axes, right side of title (figure-level)
fig.legend(
    handles=ax.get_legend_handles_labels()[0],
    labels=ax.get_legend_handles_labels()[1],
    title="Interaction Type",
    fontsize=11,
    title_fontsize=12,
    ncol=len(PLOT_TYPES),
    loc="upper right",
    bbox_to_anchor=(0.98, 1.01),
    frameon=False
)

# Remove axis legend
ax.get_legend().remove()

ax.margins(x=0.01)
plt.tight_layout(rect=[0, 0, 1, 0.88])

savefig(os.path.join(OUTPUT_DIR, "interaction_grouped_normalized.png"))
plt.close()


# ============================================================
# STEP E2B — Stacked bar plot (composition view)
# ============================================================

fig, ax = plt.subplots(figsize=(16, 6), dpi=150)

stacked_frac.plot(
    kind="bar",
    ax=ax,
    stacked=True,
    width=0.75,
    color=[COLORS[c] for c in PLOT_TYPES],
    edgecolor="white",
    linewidth=0.6,
)

ax.set_title(
    f"Normalized {chain}–{compound_name} Interaction Composition",
    fontsize=16,
    pad=20
)
ax.set_ylabel("Fraction of Frames", fontsize=12)
ax.set_xlabel("Residue", fontsize=12)

# Composition emphasis
ax.yaxis.grid(False)

ax.tick_params(axis="x", rotation=90, labelsize=8)
ax.tick_params(axis="y", labelsize=10)

# Legend ABOVE axes, right side of title (figure-level)
fig.legend(
    handles=ax.get_legend_handles_labels()[0],
    labels=ax.get_legend_handles_labels()[1],
    title="Interaction Type",
    fontsize=11,
    title_fontsize=12,
    ncol=len(PLOT_TYPES),
    loc="upper right",
    bbox_to_anchor=(0.98, 1.01),
    frameon=False
)

# Remove axis legend
ax.get_legend().remove()

ax.margins(x=0.01)
plt.tight_layout(rect=[0, 0, 1, 0.88])

savefig(os.path.join(OUTPUT_DIR, "interaction_stacked_normalized.png"))
plt.close()



# ============================================================
# STEP E2.5 — Call NETWORX (NEWNETWORX_noLP.py)
# ============================================================
phase("STEP E3 —Running Advanced NETWORX Renderer")

NETWORX_SCRIPT = "3B_NETWORX.py"   # or 3B_NETWORX.py if preferred

if not os.path.exists(NETWORX_SCRIPT):
    print(f"❌ NETWORX script {NETWORX_SCRIPT} not found — skipping NETWORX rendering.")
else:
    output_name = args.net_out or f"{lig_resname}"
    minfrac = args.minfrac           # same threshold as analyzer
    ellipse_rx = args.ellipse_rx
    ellipse_ry = args.ellipse_ry
    
    # Build CLI call string
    cmd = f"python {NETWORX_SCRIPT} -l {lig_resname} -f {minfrac} -o {output_name}"
    
    if ellipse_rx:
        cmd += f" --ellipse-rx {ellipse_rx}"
    if ellipse_ry:
        cmd += f" --ellipse-ry {ellipse_ry}"
    if args.no_edges:
        cmd += " --no-edges"

    print("➡️ Running NETWORX with command:")
    print("   ", cmd)

    run_cmd(cmd)
    print("✅ NETWORX rendering complete.")


# ============================================================
# STEP E3 — NETWORK GRAPH
# ============================================================

phase("STEP E3 — Network Graph")

filtered_summary = summary[summary["Total"] >= total_frames * MIN_CONTACT_FRAC]

G = nx.Graph()
center_label = lig_resname
G.add_node(center_label, size=3000)

for res, row in filtered_summary.iterrows():
    pct = 100 * row["Total"] / total_frames
    node_label = f"{res}\n{pct:.1f}%"

    G.add_node(res, size=2600, label=node_label)
    G.add_edge(center_label, res, weight=row["Total"])

pos = nx.spring_layout(G, seed=42)



# ============================================================
# STEP E4 — INTERACTION TIMELINE HEATMAP (Priority Encoded)
# ============================================================

phase("STEP E4 — Interaction Timeline Heatmap")

# -------------------------------
# Priority encoding
# -------------------------------
PRIORITY_MAP = {
    "H-bond": 3,
    "Ionic": 2,
    "Hydrophobic": 1
}

# Encode interaction strength numerically
int_df["InteractionCode"] = int_df["Type"].map(PRIORITY_MAP).fillna(0)

# ------------------------------------------------------------
# OPTIONAL (HIGHLY RECOMMENDED): bin time to reduce sparsity
# ------------------------------------------------------------
TIME_BIN_NS = 2.0  # adjust if needed
int_df["Time_bin"] = (
    (int_df["Time_ns"] / TIME_BIN_NS)
    .round(0)
    * TIME_BIN_NS
)

# ------------------------------------------------------------
# Pivot using MAX priority per residue per time bin
# ------------------------------------------------------------
timeline_matrix = (
    int_df
    .pivot_table(
        index="Residue",
        columns="Time_bin",
        values="InteractionCode",
        aggfunc="max",
        fill_value=0
    )
    .sort_index(
        key=lambda x: x.map(residue_number),
        ascending=False   # 🔥 highest residue number on top
    )
)

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
plt.figure(figsize=(20, 10))

sns.heatmap(
    timeline_matrix,
    cmap="viridis",
    cbar_kws={
        "label": "Interaction Strength (3=H-bond, 2=Ionic, 1=Hydrophobic)"
    },
    linewidths=0.05
)

plt.title(f"{compound_name} Interaction Timeline Heatmap")
plt.xlabel("Time (ns)")
plt.ylabel("Residue")
plt.tight_layout()

savefig(os.path.join(OUTPUT_DIR, f"{lig_resname}_interaction_timeline.png"))

print("✅ STEP E4 complete — interaction timeline heatmap created.")


# ============================================================
# STEP E4B — INTERACTION EVENT TIMELINE (Figure 24)
# ============================================================

phase("STEP E4B — Interaction Event Timeline (Figure 25)")

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Sort data
df = int_df.sort_values(["Residue", "Time_ns"])

# Residue ordering
residues = sorted(
    df["Residue"].unique(),
    key=residue_number,
    reverse=False   # 🔥 highest number at top
)

res_to_y = {res: i for i, res in enumerate(residues)}

fig, ax = plt.subplots(figsize=(22, 12))
BAR_HEIGHT = 0.8

# Color map
TYPE_COLORS = {
    "H-bond": "#FDE725",
    "Ionic": "#5EC962",
    "Hydrophobic": "#31688E"
}

# Gap threshold (ns) to define event break
GAP_NS = 1.5  # ns; maximum allowed gap between frames for a continuous event

for residue, g in df.groupby("Residue"):
    y = res_to_y[residue]
    g = g.sort_values("Time_ns")

    current_type = None
    start_time = None
    prev_time = None

    for _, row in g.iterrows():
        t = row["Time_ns"]
        itype = row["Type"]

        # Start new event
        if start_time is None:
            start_time = t
            current_type = itype

        # Break event if gap or type change
        elif (t - prev_time > GAP_NS) or (itype != current_type):
            ax.broken_barh(
                [(start_time, prev_time - start_time)],
                (y - BAR_HEIGHT / 2, BAR_HEIGHT),
                facecolors=TYPE_COLORS.get(current_type, "gray")
            )
            start_time = t
            current_type = itype

        prev_time = t

    # Final segment
    if start_time is not None:
        ax.broken_barh(
            [(start_time, prev_time - start_time)],
            (y - BAR_HEIGHT / 2, BAR_HEIGHT),
            facecolors=TYPE_COLORS.get(current_type, "gray")
        )

# Axis formatting
ax.set_yticks(range(len(residues)))
ax.set_yticklabels(residues)
ax.set_xlabel("Time (ns)")
ax.set_ylabel("Residue")
ax.set_title(f"{compound_name} Interaction Event Timeline")

ax.set_ylim(-1, len(residues))
ax.set_xlim(df["Time_ns"].min(), df["Time_ns"].max())

# Legend
legend_elements = [
    Patch(facecolor=c, label=t) for t, c in TYPE_COLORS.items()
]
ax.legend(
    handles=legend_elements,
    title="Interaction Type",
    loc="upper right"
)

plt.tight_layout()
savefig(os.path.join(OUTPUT_DIR, f"{lig_resname}_interaction_event_timeline.png"))

print("✅ STEP E4B complete — Figure 25 interaction event timeline created.")



# ============================================================
# STEP E5 — Per-Residue Interaction Composition Pie Charts
# ============================================================

phase("STEP E5 — Residue Pie Charts")
print("\n📊 Generating pie charts for per-residue interaction composition...")

top_residues = summary.sort_values("Total", ascending=False).head(20).index
# top_residues = summary.sort_values("Total", ascending=False).index

for res in top_residues:
    print(f"   • Creating pie chart for {res} ...")

    # Extract full row
    vals = summary.loc[res][["H-bond", "Hydrophobic", "Ionic"]]

    # Remove zero-valued categories
    vals = vals[vals > 0]

    # If a residue somehow has no counted interactions (very rare)
    if len(vals) == 0:
        print(f"     ⚠️  Skipping {res} — no nonzero interactions.")
        continue

    # Make pie chart
    plt.figure(figsize=(4,4))
    plt.pie(
        vals,
        labels=vals.index,
        autopct="%1.1f%%",
        startangle=90,
        pctdistance=0.75
    )
    plt.title(f"{res} Interaction Composition")
    plt.tight_layout()

    savefig(os.path.join(PIE_DIR, f"pie_{res}.png"))

print("✅ STEP E5 complete — Pie charts saved.\n")


# ============================================================
# STEP E6 — Interaction Distance Distributions
# ============================================================

phase("STEP E6 — Interaction Distance Distributions")
print("\n📈 Generating contact distance histograms per residue...")

for res in filtered_df["Residue"].unique():
    d = filtered_df.loc[filtered_df["Residue"] == res, "Distance"]

    plt.figure(figsize=(5,4))
    sns.histplot(d, bins=30, kde=True)

    plt.xlabel("Distance (Å)")
    plt.ylabel("Count")
    plt.title(f"{res} Contact Distance Distribution")

    savefig(os.path.join(DIST_DIR, f"dist_{res}.png"))
    plt.close()

print("✅ STEP E6 complete — Distance histograms saved.\n")



# ============================================================
# STEP E7 — Residue–Ligand Contact Persistence Distribution
# ============================================================

phase("STEP E7 — Interaction Degree Distribution")
print("\n📊 Generating residue–ligand contact persistence distribution...")
print("\n🔍 int_df columns:")
print(int_df.columns.tolist())

# --- Any-contact persistence ---
contact_persistence = (
    int_df
    .groupby("Residue")["Time_ns"]
    .nunique()
)

contact_persistence = contact_persistence[
    contact_persistence >= total_frames * MIN_CONTACT_FRAC
]

# --- Hydrogen-bond persistence ---
hbond_df = int_df[int_df["Type"].str.lower().isin(["hbond", "hydrogen_bond"])]


hbond_persistence = (
    hbond_df
    .groupby("Residue")["Time_ns"]
    .nunique()
)

# Apply same persistence threshold
hbond_persistence = hbond_persistence[
    hbond_persistence >= total_frames * MIN_CONTACT_FRAC
]

plt.figure(figsize=(6, 4))

sns.histplot(
    contact_persistence,
    bins=15,
    color="steelblue",
    alpha=0.6,
    label="Any ligand contact"
)

sns.histplot(
    hbond_persistence,
    bins=15,
    color="darkorange",
    alpha=0.6,
    label="Hydrogen bond contact"
)

plt.xlabel("Number of Simulation Frames with Contact")
plt.ylabel("Number of Residues")
plt.title(f"{compound_name} Residue–Ligand Contact Persistence")

plt.legend(frameon=False)
savefig(os.path.join(DIST_DIR, "interaction_degree_distribution.png"))
plt.close()

print("✅ STEP E7 complete — Contact persistence distribution saved.\n")



# ============================================================
# STEP E7B — Composite Distance Histogram
# ============================================================

phase("STEP E7B — Composite Distance Histogram")
print("\n📊 Generating composite histogram of ALL interaction distances...")

plt.figure(figsize=(10, 7))

for itype, color in zip(
    ["H-bond", "Hydrophobic", "Ionic"],
    ["#4DB6AC", "#64B5F6", "#FFD54F"]
):
    d = int_df.loc[int_df["Type"] == itype, "Distance"]
    sns.histplot(d, bins=40, kde=True, label=itype, color=color, alpha=0.6)

plt.xlabel("Interaction Distance (Å)")
plt.ylabel("Count")
plt.title(f"{compound_name} – Composite Interaction Distance Distributions")
plt.legend()

savefig(os.path.join(OUTPUT_DIR, "Composite_Distance_Distributions.png"))
plt.close()

print("✅ STEP E7B complete — Composite distance histogram saved.\n")


# ============================================================
# STEP E8 — Distance Violin Plot by Interaction Type
# ============================================================

phase("STEP E8 — Violin Plot by Interaction Type")
print("\n🎻 Generating violin plot for interaction type distance distributions...")

plt.figure(figsize=(9,6))
sns.violinplot(
    data=int_df[int_df["Type"].isin(["H-bond","Hydrophobic","Ionic"])],
    x="Type",
    y="Distance",
    palette="Set2"
)

plt.title(f"{compound_name} — Distance Distribution by Interaction Type")
plt.ylabel("Distance (Å)")

savefig(os.path.join(OUTPUT_DIR, "interaction_type_violin.png"))
plt.close()

print("✅ STEP E8 complete — Violin plot saved.\n")

# ============================================================
# STEP E9 — Interaction Heatmap (Residue × Type)
# ============================================================

phase("STEP E9 — Interaction Heatmap")
print("\n🔥 Generating interaction frequency heatmap...")

PLOT_TYPES = ["H-bond", "Hydrophobic", "Ionic"]

# Normalize counts → fraction of frames
heatmap_df = summary[PLOT_TYPES].div(total_frames)

# Order residues by overall interaction frequency (NOT a Total column)
order = heatmap_df.sum(axis=1).sort_values(ascending=False).index
heatmap_df = heatmap_df.loc[order]

plt.figure(figsize=(10, 8))
sns.heatmap(
    heatmap_df,
    cmap="magma",
    vmin=0,
    vmax=1,
    cbar_kws={
        "label": "Interaction Frequency\n(fraction of simulation frames)"
    }
)

plt.title(f"{compound_name} Interaction Type vs Residue Heatmap")
plt.xlabel("Interaction Type")
plt.ylabel("Residue")

savefig(os.path.join(OUTPUT_DIR, "interaction_heatmap_normalized.png"))
plt.close()

print("✅ STEP E9 complete — Normalized interaction heatmap saved.\n")






# ============================================================
# STEP E10 — Residue Binding Importance Ranking
# ============================================================

phase("STEP E10 — Binding Importance Ranking")
print("\n🏆 Computing and plotting binding importance ranking...")

importance = (
    3*summary["H-bond"] +
    2*summary["Ionic"] +
    1*summary["Hydrophobic"]
).sort_values(ascending=False)

plt.figure(figsize=(10,6))
importance.plot(kind="bar", color="slateblue")
plt.title("Residue Binding Importance Score")
plt.ylabel("Weighted Interaction Score")
plt.xlabel("Residue")

savefig(os.path.join(OUTPUT_DIR, "binding_importance_ranking.png"))

print("✅ STEP E10 complete — Binding importance chart saved.\n")




# ============================================================
# STEP E11 — PERSISTENCE VS DISTANCE SCATTER (FIXED)
# ============================================================

phase("STEP E11 — Persistence vs Distance Scatter")

# --- Compute per-residue mean interaction distance ---
mean_dist_series = int_df.groupby("Residue")["Distance"].mean()

# --- Compute persistence fraction for the same residues ---
persistence_series = summary["Total"] / total_frames

# --- Align indices so x and y match ---
common_residues = mean_dist_series.index.intersection(persistence_series.index)

mean_dist = mean_dist_series.loc[common_residues]
frac = persistence_series.loc[common_residues]

# --- Make sure lengths match ---
print("Residues used:", len(common_residues))

plt.figure(figsize=(10, 6))
plt.scatter(mean_dist, frac, s=120, edgecolor="k")

for r in common_residues:
    plt.text(mean_dist[r], frac[r], r, fontsize=8)

plt.xlabel("Mean Minimum Distance (Å)")
plt.ylabel("Persistence Fraction")
plt.title(f"{compound_name} — Persistence vs Distance")

savefig(os.path.join(OUTPUT_DIR, f"{lig_resname}_persistence_vs_distance.png"))

print("✅ STEP E11 complete.")





# ============================================================
# STEP D6D — FULL SSE ↔ Ligand Contact Dynamics Analysis
# Includes:
#  - Stacked SSE + Contacts plots
#  - Pearson & Spearman correlations
#  - Rolling correlations aligned to trajectory
#  - Change-point detection using rupture-Pelt
#  - Segment statistics export
# ============================================================

if not done_checkpoint("chk_D6D.txt"):
    phase("STEP D6D — FULL SSE–Ligand Contact Coupling Analysis")
    # Load aligned trajectory and ligand info
    traj, bb, ligand_atoms, lig_resname = get_aligned_traj_and_ligand()


    import ruptures as rpt
    from scipy.stats import pearsonr, spearmanr

    # ============================================================
    # 1) LOAD DSSP (from STEP D6)
    # ============================================================
    dssp = md.compute_dssp(traj)

    helix_codes = {"H", "G", "I"}
    sheet_codes = {"E", "B"}
    coil_codes  = {"S", "T", "C", " "}

    helix_frac, sheet_frac, coil_frac = [], [], []

    for row in dssp:
        total = len(row)
        helix_frac.append(sum(c in helix_codes for c in row) / total)
        sheet_frac.append(sum(c in sheet_codes for c in row) / total)
        coil_frac.append(sum(c in coil_codes  for c in row) / total)

    time_ns = traj.time / 1000.0

    # ============================================================
    # 2) LOAD CONTACTS FROM STEP E2
    # ============================================================
    if "int_df" not in globals():
        raise RuntimeError("❌ int_df is missing. Run STEP E2 before STEP D6D.")

    contact_counts = (
        int_df.groupby("Time_ns")["Residue"]
        .count()
        .reindex(time_ns, fill_value=0)
    )

    # Smooth contacts (helps interpretability)
    smooth_window = 50
    contacts_smooth = contact_counts.rolling(smooth_window, center=True).mean()

    # ============================================================
    # 3) CORRELATION ANALYSIS
    # ============================================================
    pear_corr, _ = pearsonr(contacts_smooth.fillna(0), sheet_frac)
    spear_corr, _ = spearmanr(contacts_smooth.fillna(0), sheet_frac)

    print("\n📈 SSE ↔ Contact Correlations")
    print(f"   Pearson  (contacts vs sheet): {pear_corr:.3f}")
    print(f"   Spearman (contacts vs sheet): {spear_corr:.3f}")

    # ============================================================
    # 4) FIXED, PROPER ROLLING CORRELATION
    # ============================================================
    contacts_series = pd.Series(contacts_smooth.values, index=time_ns)
    sheet_series    = pd.Series(sheet_frac, index=time_ns)

    roll_window_frames = 500  # ~1–2 ns depending on dt

    rolling_corr = (
        contacts_series
        .rolling(roll_window_frames, center=True)
        .corr(sheet_series)
    )

    # Align to time axis strictly
    rolling_corr = rolling_corr.reindex(time_ns)

    # Fix NaN edges
    rolling_corr = rolling_corr.fillna(method="bfill").fillna(method="ffill")

    # ============================================================
    # 5) CHANGE-POINT DETECTION
    # ============================================================
    print("\n🔍 Detecting change-points in ligand contacts...")

    algo = rpt.Pelt(model="rbf").fit(contacts_smooth.fillna(0).values)
    change_points = algo.predict(pen=5)

    # Remove the final endpoint (ruptures always adds last index)
    if change_points and change_points[-1] == len(time_ns):
        change_points = change_points[:-1]

    print(f"   Change-points detected at frames: {change_points}")

    # ============================================================
    # 6) CHANGE-POINT SUMMARY CSV
    # ============================================================
    cp_out = os.path.join(OUTPUT_DIR, f"{lig_resname}_ChangePoint_Summary.csv")
    with open(cp_out, "w") as f:
        f.write("Segment,StartFrame,EndFrame,MeanContacts,MeanHelix,MeanSheet,MeanCoil\n")

        start = 0
        for i, cp in enumerate(change_points + [len(time_ns)]):
            seg_contacts = contacts_smooth.iloc[start:cp].mean()
            seg_helix = np.mean(helix_frac[start:cp])
            seg_sheet = np.mean(sheet_frac[start:cp])
            seg_coil  = np.mean(coil_frac[start:cp])

            f.write(f"{i},{start},{cp},{seg_contacts:.3f},{seg_helix:.3f},"
                    f"{seg_sheet:.3f},{seg_coil:.3f}\n")

            start = cp

    print(f"📄 Change-point summary saved: {cp_out}")

    # ============================================================
    # 7) PLOTTING — SSE TOP, CONTACTS BOTTOM
    # ============================================================

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(16, 12),
        sharex=True, gridspec_kw={"height_ratios": [2, 1]}
    )

    # ------------------ TOP PLOT: SSE FRACTIONS ------------------
    ax1.plot(time_ns, helix_frac, lw=2, color="blue", label="Helix")
    ax1.plot(time_ns, sheet_frac, lw=2, color="orange", label="Sheet")
    ax1.plot(time_ns, coil_frac,  lw=2, color="green", label="Coil")

    ax1.set_ylabel("Secondary Structure Fraction")
    ax1.set_title(
        f"{compound_name} — SSE vs Ligand Contacts\n"
        f"Pearson={pear_corr:.3f}   Spearman={spear_corr:.3f}"
    )
    ax1.legend(loc="upper right", frameon=False)

    # Change-point markers
    for cp in change_points:
        ax1.axvspan(time_ns[cp] - 0.2, time_ns[cp] + 0.2,
                    color="gray", alpha=0.2)

    # ------------------ BOTTOM PLOT: CONTACTS ------------------
    ax2.plot(time_ns, contacts_smooth, color="red", lw=1.8,
             label="Ligand Contacts (smoothed)")
    ax2.fill_between(time_ns, contacts_smooth, color="red", alpha=0.2)

    ax2.plot(time_ns, contact_counts, color="black", alpha=0.2,
             label="Raw Contacts")

    ax2.set_ylabel("Ligand Contacts / Frame")
    ax2.set_xlabel("Time (ns)")
    ax2.legend(loc="upper right", frameon=False)

    # Same change-points
    for cp in change_points:
        ax2.axvspan(time_ns[cp] - 0.2, time_ns[cp] + 0.2,
                    color="gray", alpha=0.2)

    # Rolling correlation overlay
    ax2_2 = ax2.twinx()
    ax2_2.plot(time_ns, rolling_corr, color="purple", lw=1.5,
               label="Rolling Corr")
    ax2_2.set_ylabel("Rolling Corr (Sheet vs Contacts)", color="purple")
    ax2_2.tick_params(axis="y", labelcolor="purple")

    fig.tight_layout()

    out_plot = os.path.join(OUTPUT_DIR, f"{lig_resname}_SSE_vs_Contacts_FULL.png")
    savefig(out_plot)

    checkpoint("chk_D6D.txt")
    print(f"\n✅ STEP D6D complete — full analysis saved to:\n   {out_plot}\n")

else:
    print("⏭️ STEP D6D skipped (checkpoint found).")








del u
gc.collect()
print("\n✅ Pocket analysis complete — all memory released.")
print("✅ All analyses complete. Outputs saved to:", os.path.abspath(OUTPUT_DIR))

print("\n🧾 Pipeline Summary:")
print("✅ Step A: Final_Trajectory.xtc/pdb — OK")
print("✅ Step B: binding_pocket_only.xtc/pdb — OK")
print("✅ Step C: RMSD/RMSF — OK")
print("✅ Step D: Contact + Network Analysis — OK")
print(f"📂 Outputs saved to {os.path.abspath(OUTPUT_DIR)}")
