#!/usr/bin/env python3
"""
====================================================================
🏗️  Automated GROMACS MD Pipeline — Step 2: Equilibration + Production
====================================================================

Author:  Joseph-Michael Schulz (University of Miami)
Created: 2025-10-08
Version: 2.3 (GPU-Adaptive Build)
Dependencies:
    • GROMACS 2022+ (compiled with CUDA support)
    • Python 3.9+
    • MDAnalysis (optional, for binding pocket extraction)
    • NumPy
    • NVIDIA CUDA drivers installed and visible via `nvidia-smi`

--------------------------------------------------------------------
📖  OVERVIEW
--------------------------------------------------------------------
This script automates **equilibration and production molecular dynamics (MD)**
runs for protein, ligand, and PROTAC systems prepared in Step 1 of the workflow.

It dynamically configures:
• GPU and CPU threading
• MDP files (`tc-grps`, `nsteps`)
• Topology inclusion for ligands
• Ligand restraint generation
• Energy minimization → NVT → NPT → Production stages
• Optional post-run pocket extraction via MDAnalysis

The script supports three main simulation modes:
1️⃣  `protein` — single-protein systems (no ligand)
2️⃣  `ligand`  — small-molecule or peptide bound to a protein
3️⃣  `protac`  — ternary complex systems (E3 ligase + PROTAC + target)

--------------------------------------------------------------------
⚙️  KEY FEATURES
--------------------------------------------------------------------
• **Automatic GPU detection** (via `nvidia-smi`)
• **Thread scaling** via detected CPU cores
• **Dynamic tc-grps** updates to include ligand groups
• **Automatic nsteps scaling** (500,000 steps per ns)
• **Automatic topology patching** for ligand restraints
• **Energy minimization, NVT, NPT, and production MD** in series
• **Pocket builder** — extracts binding-site trajectories (optional)

--------------------------------------------------------------------
📦  REQUIRED INPUT FILES
--------------------------------------------------------------------
Expected to exist in the working directory (output of Step 1):
    topol.top
    em.mdp, nvt.mdp, npt.mdp, md.mdp
    protein_processed.gro / complex.gro
    ions.mdp
    forcefield directory (e.g., charmm36-mar2023.ff)
Optional (for ligand/protac):
    LIG.cgenff.mol2, LIG.str, LIG.gro, posre_LIG.itp, index_LIG.ndx

--------------------------------------------------------------------
🧠  USAGE
--------------------------------------------------------------------
Run directly from your system or cluster terminal after Step 1 setup:

    python 2_AutomateMD.py [options]

Examples:

▶ Protein-only simulation (default 50 ns):
    $ python 2_AutomateMD.py --mode protein

▶ Ligand-bound system (custom ligand and runtime):
    $ python 2_AutomateMD.py --mode ligand --ligand PTC --ns 25

▶ PROTAC ternary system (automatically includes atomIndex.txt):
    $ python 2_AutomateMD.py --mode protac --ligand PTX --ns 100

▶ Specify GPU manually (default = auto-detect):
    $ python 2_AutomateMD.py --mode ligand --ligand HDG --gpu 1

--------------------------------------------------------------------
🧩  OPTIONAL: BUILD POCKET TRAJECTORY
--------------------------------------------------------------------
After MD completion, extract binding-site residues interacting
with the ligand using:

    >>> from 2_AutomateMD import build_binding_pocket
    >>> build_binding_pocket("./", ligand_code="PTC", cutoff_ang=5.0)

Outputs:
    binding_pocket_only.pdb
    binding_pocket_only.xtc

--------------------------------------------------------------------
🧰  OUTPUTS
--------------------------------------------------------------------
Energy minimization:
    em.gro, em.log, em.edr
NVT/NPT equilibration:
    nvt.gro, npt.gro, nvt.cpt, npt.cpt, nvt.log, npt.log
Production MD:
    md_0_1.tpr, md_0_1.trr, md_0_1.xtc, md_0_1.log, md_0_1.edr
Optional (Ligand systems):
    posre_LIG.itp, index_LIG.ndx
Optional (Pocket extraction):
    binding_pocket_only.pdb, binding_pocket_only.xtc

--------------------------------------------------------------------
🔬  NOTES
--------------------------------------------------------------------
• Uses single precision and GPU-accelerated nonbonded + PME computations.
• Adjusts tc-grps automatically for two-temperature coupling groups:
    `Protein_LIGAND` and `Water_and_ions`
• Designed for reproducible multi-GPU use on HPC systems (e.g., Pegasus).
• All parameters editable via the MDP templates in the same directory.
• Will exit safely on any command failure with logged STDERR.

--------------------------------------------------------------------
🧭  FUTURE EXTENSIONS
--------------------------------------------------------------------
• Add adaptive sampling (REMD, REST2)
• Integrate with slurm/LSF submission scheduler
• Implement real-time energy drift monitoring
• Export JSON summaries for DockQ/PLIP analysis pipelines

====================================================================
"""


import os
import subprocess
import re
import multiprocessing
import argparse
import time
import multiprocessing, os, subprocess, time
import shutil
import os

def parse_args():
    p = argparse.ArgumentParser(
        description="Automated GROMACS Step 2: MD equilibration + production run"
    )
    p.add_argument("--mode", choices=["ligand", "protein", "protac"], required=False)
    p.add_argument("--ligand", type=str, help="3-letter ligand code (required for ligand or protac modes)")
    p.add_argument("--ns", type=float, default=None, help="Production length in nanoseconds (default 50)")
    p.add_argument("--gpu", type=int, help="GPU ID (optional; auto-detected if omitted)")
    p.add_argument("--compute", choices=["CPU", "GPU"], default="CPU", help="Computation mode")
    p.add_argument("--threads", type=int, default=None, help="Number of CPU threads to use (default = auto)")
    p.add_argument("--pinoffset", type=int, default=None, help="Starting CPU core index for thread pinning")
    p.add_argument("--headless", action="store_true",
                   help="Non-interactive mode; skip all prompts and use provided flags")
    p.add_argument("--production_only", action="store_true",
               help="Skip EM/NVT/NPT and run production only (resume if checkpoint exists).")
    p.add_argument("--resume", action="store_true",
                help="If md checkpoint exists, resume production from it (non-interactive).")
    p.add_argument("--force_restart", action="store_true",
                help="Do NOT resume even if checkpoint exists (start pipeline normally).")


    args = p.parse_args()

    # Interactive fallback only when not headless
    if not args.headless and not args.mode:
        print("\nPlease select the simulation mode:")
        print("  1. ligand")
        print("  2. protein")
        print("  3. protac")
        choice = input("Enter choice (1–3): ").strip()
        args.mode = {"1": "ligand", "2": "protein", "3": "protac"}.get(choice, "protein")
        print(f"🧠 Selected mode: {args.mode}")

    return args


def system_has_virtual_sites(tpr_path: str) -> bool:
    """
    Detect whether the system contains virtual sites by inspecting `gmx dump`.
    """
    try:
        out = subprocess.check_output(
            f"gmx dump -s {tpr_path} | grep -i virtual",
            shell=True,
            text=True
        )
        return bool(out.strip())
    except subprocess.CalledProcessError:
        return False

GPU_PROFILES = [
    # 🚀 SAFE for ALL systems (virtual sites OK)
    {
        "name": "gpu_safe",
        "flags": "-nb gpu -pme gpu -bonded gpu -update cpu"
    },

    # 🚀 FULL GPU (ONLY if no virtual sites)
    {
        "name": "gpu_full",
        "flags": "-nb gpu -pme gpu -bonded gpu -update gpu"
    },
]



def run_mdrun_with_fallback(
    base_cmd,
    cwd,
    pin_offset,
    threads,
    gpu_id,
    use_gpu=True,
    tpr_check=None
):
    """
    Try best GPU profile first; if virtual sites exist, force CPU update.
    Fall back to CPU-only if all GPU attempts fail.
    """

    has_vs = False
    if tpr_check:
        tpr_path = os.path.join(cwd, tpr_check)
        if os.path.exists(tpr_path):
            has_vs = system_has_virtual_sites(tpr_path)
            if has_vs:
                print("⚠️ Virtual sites detected → forcing -update cpu (no GPU update).")

    # Profile order:
    #   - virtual sites: only safe
    #   - no virtual sites: full first, then safe
    if has_vs:
        profiles = [p for p in GPU_PROFILES if "-update cpu" in p["flags"]]
    else:
        profiles = []
        # prefer full if available
        for p in GPU_PROFILES:
            if "-update gpu" in p["flags"]:
                profiles.append(p)
        for p in GPU_PROFILES:
            if "-update cpu" in p["flags"]:
                profiles.append(p)

    # --------------------
    # GPU attempts
    # --------------------
    if use_gpu and gpu_id is not None and gpu_id >= 0:
        for profile in profiles:
            cmd = (
                f"{base_cmd} "
                f"-pin on -pinoffset {pin_offset} "
                f"-ntmpi 1 -ntomp {threads} "
                f"-gpu_id 0 {profile['flags']}"
            )

            print(f"\n🚀 Trying {profile['name']} → {cmd}")
            rc = run_command_check_rc(cmd, cwd=cwd)

            if rc == 0:
                print(f"✅ Success using {profile['name']}")
                return

            print(f"⚠️ {profile['name']} failed, trying next fallback...")

    # --------------------
    # CPU fallback
    # --------------------
    print("\n🧯 Falling back to CPU-only execution")
    cpu_cmd = (
        f"{base_cmd} "
        f"-pin on -pinoffset {pin_offset} "
        f"-ntmpi 1 -ntomp {threads}"
    )

    rc = run_command_check_rc(cpu_cmd, cwd=cwd)
    if rc != 0:
        print("❌ CPU fallback also failed. Aborting.")
        exit(rc)

    print("✅ CPU-only execution successful")



def run_command(command, cwd=None, input_text=None):
    """Runs a shell command in the given directory and exits on failure."""
    print(f"\n🔹 Running command: {command}")
    result = subprocess.run(command, shell=True, cwd=cwd, text=True, capture_output=True, input=input_text)
    if result.stdout:
        print(f"📜 STDOUT:\n{result.stdout}")
    if result.returncode != 0:
        print(f"❌ Error running command: {command}")
        print(f"🔺 STDERR: {result.stderr}")
        exit(1)


import subprocess
import os


import shutil
import re
import subprocess
import os

def read_tpr_nsteps_dt(tpr_path: str):
    """
    Pull dt (ps) and nsteps from a TPR by streaming `gmx dump`,
    stopping as soon as we find both.
    """
    cmd = f"gmx dump -s {tpr_path}"
    proc = subprocess.Popen(cmd, shell=True, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)

    dt_ps = None
    nsteps = None

    for line in proc.stdout:
        if dt_ps is None and "delta_t" in line:
            m = re.search(r"delta_t\s*=\s*([0-9.]+)", line)
            if m:
                dt_ps = float(m.group(1))

        if nsteps is None and "nsteps" in line:
            m = re.search(r"nsteps\s*=\s*(\d+)", line)
            if m:
                nsteps = int(m.group(1))

        if dt_ps is not None and nsteps is not None:
            proc.kill()
            break

    proc.wait()
    return nsteps, dt_ps


def tpr_total_ns(tpr_path: str):
    nsteps, dt_ps = read_tpr_nsteps_dt(tpr_path)
    if nsteps is None or dt_ps is None:
        return None
    # time = nsteps * dt (ps) -> convert ps to ns
    return (nsteps * dt_ps) / 1000.0


def extend_tpr_to_target_ns(directory: str, tpr_name: str, target_ns: float):
 
    tpr_path = os.path.join(directory, tpr_name)
    current_ns = tpr_total_ns(tpr_path)

    if current_ns is None:
        print("⚠️ Could not read current TPR length; skipping extension.")
        return

    if target_ns <= current_ns + 1e-6:
        print(f"ℹ️ TPR already set to ~{current_ns:.3f} ns (>= {target_ns} ns). No extension needed.")
        return

    extend_ns = target_ns - current_ns
    extend_ps = extend_ns * 1000.0

    tmp_tpr = tpr_path + ".tmp"
    print(f"🧩 Extending TPR: {current_ns:.3f} ns → {target_ns:.3f} ns (extend {extend_ns:.3f} ns)")

    run_command_cpu(
        f"gmx convert-tpr -s {tpr_name} -extend {extend_ps:.3f} -o {os.path.basename(tmp_tpr)}",
        cwd=directory
    )
    shutil.move(tmp_tpr, tpr_path)
    print("✅ TPR extension complete.")



def run_command_check_rc(command, cwd=None, input_text=None):
    """
    Special version of run_command:
    - DOES NOT exit on failure
    - Returns the command's exit code
    - Captures stdout/stderr for debugging
    """
    print(f"\n🔸 Running (check_rc): {command}")
    result = subprocess.run(
        command, shell=True, cwd=cwd, text=True,
        capture_output=True, input=input_text
    )
    if result.stdout:
        print(result.stdout)
    if result.stderr:
        print(result.stderr)
    return result.returncode  # <-- IMPORTANT





def maybe_resume_production_only(directory, args, gpu_id, pin_offset, threads, use_gpu, target_ns, md_deffnm="md_0_1"):
    md_tpr = os.path.join(directory, f"{md_deffnm}.tpr")
    md_cpt = os.path.join(directory, f"{md_deffnm}.cpt")
    md_prev_cpt = os.path.join(directory, f"{md_deffnm}_prev.cpt")

    # Prefer md_0_1.cpt, fallback to md_0_1_prev.cpt if needed
    cpt_to_use = md_cpt if os.path.exists(md_cpt) else (md_prev_cpt if os.path.exists(md_prev_cpt) else None)

    if args.force_restart:
        return False

    if not (os.path.exists(md_tpr) and cpt_to_use):
        return False

    # Decide whether to resume
    if args.production_only or args.resume:
        do_resume = True
    elif args.headless:
        # sensible default for headless runs: resume if possible
        do_resume = True
    else:
        ans = input(
            f"\n♻️ Found production checkpoint ({os.path.basename(cpt_to_use)}).\n"
            f"Resume production ONLY (skip EM/NVT/NPT)? [Y/n]: "
        ).strip().lower()
        do_resume = (ans in ("", "y", "yes"))

    if not do_resume:
        return False

    print("\n🚀 Resuming PRODUCTION only...")
    # If user requested a longer run, extend the TPR so it doesn’t stop early
    if target_ns is not None:
        extend_tpr_to_target_ns(directory, f"{md_deffnm}.tpr", float(target_ns))

    base_cmd = f"gmx mdrun -v -deffnm {md_deffnm} -cpi {os.path.basename(cpt_to_use)} -append"

    # Reuse your existing GPU→GPU-lite→CPU fallback runner
    run_mdrun_with_fallback(
        base_cmd=base_cmd,
        cwd=directory,
        pin_offset=pin_offset,
        threads=threads,
        gpu_id=gpu_id,
        use_gpu=use_gpu,
        tpr_check=f"{md_deffnm}.tpr"   # ✅ ADD THIS
    )


    print("✅ Production resume complete.\n")
    return True



def run_command_cpu(command, cwd=None, input_text=None):
    """Run a shell command (CPU version). Exits on failure."""
    print(f"\n🔹 [CPU] Running command: {command}")
    result = subprocess.run(command, shell=True, cwd=cwd, text=True, capture_output=True, input=input_text)
    if result.stdout:
        print(result.stdout)
    if result.returncode != 0:
        print(f"❌ Error running command: {command}")
        print(f"🔺 STDERR: {result.stderr}")
        exit(result.returncode)


def run_command_gpu(command, cwd=None, input_text=None, gpu_id=0, threads=None):
    """Run a GPU-accelerated GROMACS command without overwriting the existing CUDA isolation."""
    threads = threads or os.cpu_count()

    # Do NOT reset CUDA_VISIBLE_DEVICES here — it is already set globally in environment setup
    os.environ["OMP_NUM_THREADS"] = str(threads)

    print(f"\n🚀 [GPU] Running on visible GPU 0 (physical GPU mapped via CUDA_VISIBLE_DEVICES={os.environ.get('CUDA_VISIBLE_DEVICES')})")
    print(f"🔹 Command: {command}")

    process = subprocess.Popen(
        command,
        shell=True,
        cwd=cwd,
        text=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        bufsize=1,
        universal_newlines=True
    )
    if input_text:
        process.stdin.write(input_text)
        process.stdin.close()
    for line in process.stdout:
        print(line, end="")
    process.wait()

    if process.returncode != 0:
        print(f"❌ GPU command failed: {command}")
        exit(process.returncode)



def run_command_auto(command, cwd=None, input_text=None, use_gpu=True, gpu_id=0):
    if use_gpu:
        run_command_gpu(command, cwd=cwd, input_text=input_text, gpu_id=gpu_id)
    else:
        run_command_cpu(command, cwd=cwd, input_text=input_text)



def modify_topology_for_ligand(topol_path, ligand_code):
    """Ensures ligand position restraints are correctly added to topol.top (if ligand present)."""
    with open(topol_path, "r") as f:
        lines = f.readlines()
    if ligand_code:
        for i, line in enumerate(lines):
            if f'#include "{ligand_code.lower()}.itp"' in line:
                insert_index = i + 1
                break
        else:
            print(f"⚠️ `#include \"{ligand_code.lower()}.itp\"` not found in topol.top. Skipping modification.")
            return
        posres_block = [
            "; Ligand position restraints\n",
            "#ifdef POSRES\n",
            f'#include "posre_{ligand_code.lower()}.itp"\n',
            "#endif\n"
        ]
        if not any(f'#include "posre_{ligand_code.lower()}.itp"' in line for line in lines):
            lines[insert_index:insert_index] = posres_block
        with open(topol_path, "w") as f:
            f.writelines(lines)
        print(f"✅ Updated {topol_path} with ligand position restraints.")
    else:
        print("ℹ️ No ligand present; no topology modifications needed.")


def build_binding_pocket(directory: str, ligand_code: str, traj="Final_Trajectory.xtc",
                        topo_tpr="md_0_1.tpr", cutoff_ang=5.0,
                        out_pdb="binding_pocket_only.pdb", out_xtc="binding_pocket_only.xtc"):
    """
    Create pocket-only topology+trajectory containing protein residues within `cutoff_ang` Å of the ligand
    across the whole trajectory.
    """
    try:
        import MDAnalysis as mda
        from MDAnalysis.analysis import distances
        import numpy as np

        traj_path = os.path.join(directory, traj)
        tpr_path  = os.path.join(directory, topo_tpr)

        if not os.path.exists(traj_path):
            print(f"❌ Pocket build: trajectory not found: {traj_path}")
            return
        if not ligand_code:
            print("ℹ️ Pocket build skipped (no ligand code provided).")
            return

        print(f"🔍 Pocket build: loading {topo_tpr} + {traj} ...")
        try:
            u = mda.Universe(tpr_path, traj_path)
        except Exception:
            print("⚠️ TPR not supported, using GRO fallback...")
            u = mda.Universe(os.path.join(directory, "md_0_1.gro"), traj_path)

        # ====== main logic ======True
        protein = u.select_atoms("protein")
        ligand  = u.select_atoms(f"resname {ligand_code}")

        if ligand.n_atoms == 0:
            print(f"❌ Pocket build: no atoms for ligand resname '{ligand_code}'. Skipping.")
            return

        print(f"📏 Finding protein residues within {cutoff_ang:.1f} Å of ligand {ligand_code} over all frames...")
        close_resids = set()
        for ts in u.trajectory:
            D = distances.distance_array(ligand.positions, protein.positions)  # Å
            for res in protein.residues:
                if np.any(D[:, res.atoms.indices] < cutoff_ang):
                    close_resids.add(res.resid)

        if not close_resids:
            print("⚠️ Pocket build: no contacting residues found; writing ligand-only pocket.")
            selection = f"resname {ligand_code}"
        else:
            res_sel = " or ".join([f"resid {r}" for r in sorted(close_resids)])
            selection = f"resname {ligand_code} or ({res_sel})"

        pocket_atoms = u.select_atoms(selection)
        print(f"🧬 Pocket atoms: {pocket_atoms.n_atoms}")

        # Write reduced trajectory
        out_xtc_path = os.path.join(directory, out_xtc)
        out_pdb_path = os.path.join(directory, out_pdb)

        print(f"💾 Writing pocket trajectory → {out_xtc_path}")
        with mda.Writer(out_xtc_path, pocket_atoms.n_atoms) as w:
            for ts in u.trajectory:
                w.write(pocket_atoms)

        # Write a single-frame PDB topology (first frame)
        u.trajectory[0]
        print(f"💾 Writing pocket topology → {out_pdb_path}")
        pocket_atoms.write(out_pdb_path)

        print("✅ Pocket build complete.")

    except ImportError as e:
        print(f"❌ Pocket build: missing dependency → {e}. Install MDAnalysis in this environment.")


def print_numa_summary():
    print("\n🧩 NUMA/GPU core mapping summary:")
    try:
        gpu_list = subprocess.check_output(
            "nvidia-smi --query-gpu=name --format=csv,noheader",
            shell=True, text=True
        ).strip().splitlines()
    except Exception:
        gpu_list = []
    for i, name in enumerate(gpu_list):
        cores = get_numa_cores_for_gpu(i)
        print(f"  GPU {i} ({name.strip()}): cores {cores[0]}–{cores[-1]} ({len(cores)} cores)")



def set_gpu_env(gpu_id, threads):
    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    os.environ["OMP_NUM_THREADS"] = str(threads)

def clear_gpu_env():
    for k in [
        "CUDA_VISIBLE_DEVICES",
        "GMX_FORCE_GPU_ID",
        "GMX_GPU_DD_COMMS",
        "GMX_GPU_PME_DECOMPOSITION"
    ]:
        os.environ.pop(k, None)




def get_numa_cores_for_gpu(gpu_id: int):
    """
    Return a list of CPU cores local to a GPU (NUMA-aware if possible).

    Behavisor:
    • Uses `lscpu -e=CPU,NODE` to map physical cores to NUMA nodes.
    • If NUMA nodes ≥ number of GPUs → assign GPU N → Node N.
    • If fewer NUMA nodes than GPUs → split all cores evenly.
    • If NUMA info unavailable → fallback to even split across GPUs.
    """

    import subprocess, multiprocessing

    total_cores = multiprocessing.cpu_count()

    # --- Detect GPU count ---
    try:
        gpu_list = subprocess.check_output(
            "nvidia-smi --query-gpu=name --format=csv,noheader",
            shell=True, text=True
        ).strip().splitlines()
        gpu_count = len([g for g in gpu_list if g])
    except Exception:
        gpu_count = 1

    # --- Attempt NUMA-aware mapping ---
    try:
        output = subprocess.check_output("lscpu -e=CPU,NODE", shell=True, text=True)
        node_map = {}
        for line in output.splitlines()[1:]:
            parts = line.split()
            if len(parts) >= 2 and parts[1] != "-":
                cpu, node = int(parts[0]), int(parts[1])
                node_map.setdefault(node, []).append(cpu)

        # Sort nodes for deterministic order
        node_ids = sorted(node_map.keys())

        # If NUMA nodes ≥ GPUs, map 1:1
        if gpu_id < len(node_ids):
            mapped_cores = node_map[node_ids[gpu_id]]
            if mapped_cores:
                return mapped_cores

        # Otherwise fall through to even-split fallback
    except Exception:
        pass  # will fallback below

    # --- Fallback: even split across GPUs ---
    per_gpu = max(1, total_cores // max(1, gpu_count))
    start = gpu_id * per_gpu
    end = min(start + per_gpu, total_cores)
    return list(range(start, end))


import os, subprocess, psutil

def set_env(gpu_id:int, threads:int, offset:int):
    """Prepare environment vars for proper GPU/CPU pinning."""
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    env["OMP_NUM_THREADS"] = str(threads)
    core_range = f"{offset}-{offset + threads - 1}"
    env["OMP_PLACES"] = f"{{{core_range}}}"
    env["OMP_PROC_BIND"] = "close"
    env["GOMP_CPU_AFFINITY"] = core_range
    env["GMX_GPU_DD_COMMS"] = "CUDA"
    env["GMX_GPU_PME_DECOMPOSITION"] = "CUDA"
    return env




def detect_ligand_from_cgenff():
    """
    Detect ligand name from *.cgenff.mol2 file.
    Returns ligand code (uppercase) or None.
    """
    mol2_files = [
        f for f in os.listdir(".")
        if f.lower().endswith(".cgenff.mol2")
    ]

    if not mol2_files:
        return None

    if len(mol2_files) > 1:
        print("⚠️ Multiple CGenFF MOL2 files detected:")
        for i, f in enumerate(mol2_files):
            print(f"  [{i}] {f}")
        try:
            idx = int(input("Select ligand MOL2 index [default 0]: ").strip() or "0")
            fname = mol2_files[idx]
        except Exception:
            fname = mol2_files[0]
    else:
        fname = mol2_files[0]

    ligand = fname.split(".cgenff.mol2")[0]
    return ligand.upper()


def index_has_group(index_path, group_name):
    """Return True if index.ndx contains a group with the exact given name."""
    if not os.path.exists(index_path):
        return False

    with open(index_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("[") and line.endswith("]"):
                name = line[1:-1].strip()
                if name == group_name:
                    return True
    return False


def standardize_mdp_groups(mdp_path, coupling_group_str):
    """
    Force tc-grps, comm-mode, and comm-grps to match the expected group layout.
    This prevents stale values from previous ligand/protein runs.
    """
    if not os.path.exists(mdp_path):
        print(f"❌ ERROR: Missing MDP file: {mdp_path}")
        exit(1)

    with open(mdp_path, "r") as f:
        lines = f.readlines()

    filtered = []
    found_tc = False
    found_comm_mode = False
    found_comm_grps = False

    for line in lines:
        stripped = line.strip()

        if stripped.startswith("tc-grps"):
            filtered.append(f"tc-grps                 = {coupling_group_str}\n")
            found_tc = True
        elif stripped.startswith("comm-mode"):
            filtered.append("comm-mode               = Linear\n")
            found_comm_mode = True
        elif stripped.startswith("comm-grps"):
            filtered.append(f"comm-grps               = {coupling_group_str}\n")
            found_comm_grps = True
        else:
            filtered.append(line)

    if not found_tc:
        filtered.append(f"tc-grps                 = {coupling_group_str}\n")
    if not found_comm_mode:
        filtered.append("comm-mode               = Linear\n")
    if not found_comm_grps:
        filtered.append(f"comm-grps               = {coupling_group_str}\n")

    with open(mdp_path, "w") as f:
        f.writelines(filtered)


def setup_md(directory, gpu_id, ligand_code, simulation_type, simulation_time_ns, use_gpu, args):
    """
    Runs full MD simulation setup on a specified GPU, with dynamic MDP configuration.
    Includes automatic tc-grps update, nsteps scaling, ligand topology handling,
    and indexed group management.
    """

    # ============================================================
    # ##0##  INITIAL ENVIRONMENT + GPU SELECTION + NUMA-AWARE THREADING
    # ============================================================
    print("\n🧩 [STEP 0] Configuring environment and CPU/GPU threading ...")

    available_cores = multiprocessing.cpu_count()
    print(f"🧠 Detected {available_cores} total CPU cores on system.")

    # --- Detect GPUs ---
    try:
        gpu_list = subprocess.check_output(
            "nvidia-smi --query-gpu=name --format=csv,noheader",
            shell=True, text=True
        ).strip().splitlines()
        gpu_list = [g.strip() for g in gpu_list if g.strip()]
        gpu_count = len(gpu_list)
    except Exception:
        gpu_list, gpu_count = [], 0

    # --- Interactive GPU selection ---
    headless_mode = getattr(args, "headless", False)
    if not headless_mode and gpu_count > 0:
        print(f"\n🧠 Detected {gpu_count} GPU(s):")
        for i, name in enumerate(gpu_list):
            print(f"  [{i}] {name}")
        gpu_choice = input(f"\nSelect GPU number (0–{gpu_count-1}) [default 0]: ").strip()
        gpu_id = int(gpu_choice) if gpu_choice.isdigit() and int(gpu_choice) < gpu_count else 0
        print(f"🎯 Selected GPU {gpu_id}: {gpu_list[gpu_id]}")
    elif gpu_count > 0:
        gpu_id = args.gpu if getattr(args, "gpu", None) is not None else 0
        print(f"🤖 Headless mode: using GPU {gpu_id} ({gpu_list[gpu_id]})")
    else:
        gpu_id = -1
        print("⚠️ No GPUs detected — running in CPU mode.")

    # --- NUMA-aware core detection ---
    try:
        cores = get_numa_cores_for_gpu(gpu_id) if gpu_id >= 0 else list(range(available_cores))
    except Exception:
        cores = list(range(available_cores))

    if not cores:
        per_gpu = max(1, available_cores // max(1, gpu_count or 1))
        start = gpu_id * per_gpu
        end = min(start + per_gpu, available_cores)
        cores = list(range(start, end))

    suggested_threads = len(cores)
    suggested_offset = cores[0]
    print(f"📊 Suggested NUMA mapping for GPU {gpu_id}: cores {suggested_offset}–{cores[-1]} ({suggested_threads} cores)")

    # --- Interactive override ---
    if not headless_mode:
        accept = input(f"Use this mapping for GPU {gpu_id}? [Y/n]: ").strip().lower()
        if accept in ["n", "no"]:
            try:
                threads_in = input(f"Enter number of threads [default {suggested_threads}]: ").strip()
                threads = int(threads_in) if threads_in else suggested_threads
                offset_in = input(f"Enter pin offset [default {suggested_offset}]: ").strip()
                pin_offset = int(offset_in) if offset_in else suggested_offset
            except Exception:
                print("⚠️ Invalid input, reverting to suggested values.")
                threads, pin_offset = suggested_threads, suggested_offset
        else:
            threads, pin_offset = suggested_threads, suggested_offset
    else:
        threads, pin_offset = suggested_threads, suggested_offset
        print(f"🤖 Headless mode: auto-accepted NUMA mapping → threads={threads}, pin offset={pin_offset}")

    # --- Apply CLI overrides ---
    if getattr(args, "threads", None):
        threads = int(args.threads)
    if getattr(args, "pinoffset", None):
        pin_offset = int(args.pinoffset)

    available_threads = max(1, threads)

    # --- Environment setup (GPU isolation like 2A_AutoGMXrestart.py) ---
    if gpu_id is not None and gpu_id >= 0:
        # Map the chosen physical GPU (e.g., 1 or 2) into a single visible device (index 0)
        os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
        os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
        os.environ["GMX_FORCE_GPU_ID"] = "0"  # GROMACS now always sees it as device 0
        print(f"🎯 Forcing isolation: physical GPU {gpu_id} → visible GPU 0")
    else:
        # No GPU explicitly chosen: expose all GPUs
        try:
            visible = subprocess.check_output(
                "nvidia-smi --query-gpu=index --format=csv,noheader",
                shell=True, text=True
            ).strip().replace("\n", ",")
            os.environ["CUDA_VISIBLE_DEVICES"] = visible
            print(f"💻 No GPU ID specified; exposing all GPUs: {visible}")
        except Exception:
            os.environ["CUDA_VISIBLE_DEVICES"] = "0"
            print("⚠️ Could not query nvidia-smi; defaulting to GPU 0")

    # Set thread binding
    os.environ["OMP_NUM_THREADS"] = str(available_threads)



    print(f"\n🧠 Environment configured:")
    print(f"   • Detected total cores: {available_cores}")
    if gpu_id >= 0:
        print(f"   • GPU {gpu_id}: {gpu_list[gpu_id]} (cores {suggested_offset}–{cores[-1]})")
    else:
        print("   • CPU-only mode active.")
    print(f"   • Using {available_threads} threads, pin offset {pin_offset}")
    print(f"   • CUDA_VISIBLE_DEVICES={gpu_id if gpu_id >= 0 else 'None'}")
    print(f"   • OMP_NUM_THREADS={available_threads}")




    # after your Step 0 environment config (once pin_offset and available_threads are known)
    if maybe_resume_production_only(
        directory=directory,
        args=args,
        gpu_id=gpu_id,
        pin_offset=pin_offset,
        threads=available_threads,
        use_gpu=use_gpu,
        target_ns=simulation_time_ns,
        md_deffnm="md_0_1"
    ):
        return  # <- IMPORTANT: skip steps 1–7 entirely


    # ============================================================
    # 🧩 1. Update MDP Files (tc-grps + comm-grps + nsteps)
    # ============================================================
    print("\n🧩 [STEP 1] Updating MDP files (tc-grps, comm-grps, nsteps) ...")

    mdp_files = ["md.mdp", "nvt.mdp", "npt.mdp"]

    # Decide the exact thermostat / COM groups based on system type
    if ligand_code:
        coupling_group_str = f"Protein_{ligand_code} Water_and_ions"
    else:
        coupling_group_str = "Protein Water_and_ions"

    for mdp in mdp_files:
        mdp_path = os.path.join(directory, mdp)
        if not os.path.exists(mdp_path):
            print(f"❌ ERROR: Required file {mdp} not found.")
            exit(1)

        with open(mdp_path) as f:
            lines = f.readlines()

        filtered = []
        found_tc = False
        found_comm_mode = False
        found_comm_grps = False

        for line in lines:
            stripped = line.strip()

            if stripped.startswith("tc-grps"):
                filtered.append(f"tc-grps                 = {coupling_group_str}\n")
                found_tc = True
            elif stripped.startswith("comm-mode"):
                filtered.append("comm-mode               = Linear\n")
                found_comm_mode = True
            elif stripped.startswith("comm-grps"):
                filtered.append(f"comm-grps               = {coupling_group_str}\n")
                found_comm_grps = True
            else:
                filtered.append(line)

        # Ensure required entries exist even if missing from template
        if not found_tc:
            filtered.append(f"tc-grps                 = {coupling_group_str}\n")
        if not found_comm_mode:
            filtered.append("comm-mode               = Linear\n")
        if not found_comm_grps:
            filtered.append(f"comm-grps               = {coupling_group_str}\n")

        with open(mdp_path, "w") as f:
            f.writelines(filtered)

    print(f"🛠️ Standardized tc-grps/comm-grps to: {coupling_group_str}")

    # Update nsteps in md.mdp (dt=0.002 ps → 500,000 steps/ns)
    md_path = os.path.join(directory, "md.mdp")
    if not os.path.exists(md_path):
        print("❌ ERROR: md.mdp not found.")
        exit(1)

    nsteps = int(round(simulation_time_ns * 500000))
    with open(md_path) as f:
        lines = [
            f"nsteps                  = {nsteps}\n" if l.strip().startswith("nsteps") else l
            for l in f
        ]
    with open(md_path, "w") as f:
        f.writelines(lines)

    print(f"🧮 Set nsteps = {nsteps} (~{simulation_time_ns} ns).")
    print("✅ Finished MDP standardization.\n")

    # ============================================================
    # 🧱 2. Topology Preparation
    # ============================================================
    print("🧱 [STEP 2] Preparing topology ...")
    topol_path = os.path.join(directory, "topol.top")
    if ligand_code:
        print(f"🔧 Modifying topology to include ligand {ligand_code} ...")
        modify_topology_for_ligand(topol_path, ligand_code)
        print(f"✅ Ligand '{ligand_code}' integrated into topology.")
    else:
        print("ℹ️ No ligand present; skipping topology modification.")

    # ============================================================
    # ⚡ 3. Energy Minimization
    # ============================================================
    print("\n⚡ [STEP 3] Starting energy minimization ...")

    run_mdrun_with_fallback(
        base_cmd="gmx mdrun -v -deffnm em",
        cwd=directory,
        pin_offset=pin_offset,
        threads=available_threads,
        gpu_id=gpu_id,
        use_gpu=use_gpu
    )

    print("✅ Energy minimization complete.\n")


    


    # ============================================================
    # 🧬 4. Ligand Restraints (generate if missing)
    # ============================================================
    if ligand_code:
        posre_file = f"posre_{ligand_code.lower()}.itp"
        ligand_gro = f"{ligand_code.lower()}.gro"
        index_file = f"index_{ligand_code.lower()}.ndx"

        if not os.path.exists(os.path.join(directory, posre_file)):
            print(f"🧪 Generating position restraints for ligand: {ligand_code}")

            if os.path.exists(os.path.join(directory, ligand_gro)):

                # ------------------------------
                # TRY 1 →  H* wildcard``
                # TRY 2 →  explicit H
                # TRY 3 →  keep full ligand
                # ------------------------------

                cmd1 = (
                    f"printf \"0 & ! a H*\\nq\\n\" | "
                    f"gmx make_ndx -f {ligand_gro} -o {index_file}"
                )

                cmd2 = (
                    f"printf \"0 & ! a H\\nq\\n\" | "
                    f"gmx make_ndx -f {ligand_gro} -o {index_file}"
                )

                cmd3 = (
                    f"printf \"0\\nq\\n\" | "
                    f"gmx make_ndx -f {ligand_gro} -o {index_file}"
                )

                print("🔹 Attempt 1: Removing hydrogens using wildcard H* ...")
                rc1 = run_command_check_rc(cmd1, cwd=directory)

                if rc1 != 0:
                    print("⚠️ Attempt 1 failed — trying explicit hydrogen selector...")

                    rc2 = run_command_check_rc(cmd2, cwd=directory)

                    if rc2 != 0:
                        print("⚠️ Attempt 2 failed — using full ligand group...")

                        rc3 = run_command_check_rc(cmd3, cwd=directory)

                        if rc3 != 0:
                            print("❌ All index generation attempts failed.")
                            index_file = None


                # If we successfully generated an index → create restraints
                if index_file and os.path.exists(os.path.join(directory, index_file)):
                    print("🔧 Creating ligand position restraints (genrester)...")
                    run_command(
                        f"echo '3' | gmx genrestr -f {ligand_gro} -n {index_file} "
                        f"-o {posre_file} -fc 1000 1000 1000",
                        cwd=directory
                    )
                    print(f"✅ Created {posre_file}")

            else:
                print(f"⚠️ Ligand GRO file not found ({ligand_gro}). Skipping posre generation.")
        else:
            print(f"✅ Using existing ligand position restraints: {posre_file}")






    # ============================================================
    # 🧩 5. Build Master Index File (Protein / Ligand / PROTAC)
    # ============================================================
    print("\n🧩 [STEP 5] Building clean index.ndx ...")

    # ------------------------------------------------------------
    # 5A — Clean old index files
    # ------------------------------------------------------------
    cleanup_files = ["index.ndx", "index_default.ndx"]
    if ligand_code:
        cleanup_files.append(f"index_{ligand_code.lower()}.ndx")

    for f in cleanup_files:
        fp = os.path.join(directory, f)
        if os.path.exists(fp):
            os.remove(fp)

    # ------------------------------------------------------------
    # 5B — Ensure ligand restraint file exists if ligand is present
    # ------------------------------------------------------------
    if ligand_code:
        lig = ligand_code.upper()
        posre = os.path.join(directory, f"posre_{lig.lower()}.itp")
        lig_gro = os.path.join(directory, f"{lig.lower()}.gro")

        if not os.path.exists(posre):
            print(f"⚠️ posre_{lig.lower()}.itp missing — generating fallback restraints...")

            if not os.path.exists(lig_gro):
                print(f"❌ FATAL: {lig_gro} not found — cannot create restraints.")
                exit(1)

            print(f"🔍 Scanning {lig_gro} for heavy atoms...")
            heavy_atoms = []

            with open(lig_gro, "r") as g:
                for line in g:
                    if len(line) < 20:
                        continue
                    resname = line[5:10].strip()
                    if resname != lig:
                        continue
                    atom_name = line[10:15].strip()
                    atom_id = int(line[15:20])
                    if not atom_name.startswith("H"):
                        heavy_atoms.append(atom_id)

            print(f"🔎 Found {len(heavy_atoms)} heavy atoms")

            if heavy_atoms:
                with open(posre, "w") as out:
                    out.write("[ position_restraints ]\n")
                    out.write("; atom  type    fx    fy    fz\n")
                    for a in heavy_atoms:
                        out.write(f"{a:5d}   1   1000   1000   1000\n")
                print(f"✅ Wrote fallback restraint file: {posre}")
            else:
                print("⚠️ No heavy atoms found — posre will be empty.")
        else:
            print(f"✔ Using existing ligand restraints: {posre}")

    # ------------------------------------------------------------
    # 5C — MODE-SPECIFIC INDEX BUILD
    # ------------------------------------------------------------
    simtype = simulation_type.upper().strip()

    if (not ligand_code) and (simtype != "PROTAC"):
        # ========================================================
        # PROTEIN-ONLY MODE
        # ========================================================
        print("🧬 Protein-only system → building default index.ndx")

        run_command(
            "gmx make_ndx -f em.gro -o index.ndx",
            cwd=directory,
            input_text="q\n"
        )

        expected_group = "Protein"
        expected_water_group = "Water_and_ions"
        print("✅ index.ndx built (protein only)")

    elif simtype == "PROTAC":
        # ========================================================
        # PROTAC MODE
        # ========================================================
        print("🧬 PROTAC system → building custom multi-group index.ndx")

        if not ligand_code:
            print("❌ FATAL: PROTAC mode requires a ligand_code.")
            exit(1)

        atom_file = os.path.join(directory, "atomIndex.txt")
        if not os.path.exists(atom_file):
            print("❌ FATAL: atomIndex.txt missing for PROTAC mode.")
            exit(1)

        run_command(
            "gmx make_ndx -f em.gro -o index.ndx",
            cwd=directory,
            input_text=(
                "1 | 13\n"
                f"name 21 Protein_{ligand_code.upper()}\n"
                + "".join(
                    f"{line.strip()}\nname {22+i} "
                    f"{['Ligase','Target_A','Target_B','Target_C','Target_D'][i] if i < 5 else f'Target_{i}'}\n"
                    for i, line in enumerate(open(atom_file))
                    if line.strip()
                )
                + "q\n"
            )
        )

        expected_group = f"Protein_{ligand_code.upper()}"
        expected_water_group = "Water_and_ions"
        print(f"✅ PROTAC index.ndx built with merged group: {expected_group}")

    elif ligand_code:
        # ========================================================
        # LIGAND MODE
        # KEEP OLD LOGIC THAT PREVIOUSLY WORKED
        # ========================================================
        expected_group = f"Protein_{ligand_code.upper()}"
        expected_water_group = "Water_and_ions"

        print(f"💊 Ligand-bound system → building merged index group {expected_group}")

        ndx_input = (
            "1 | 13\n"
            f"name 21 {expected_group}\n"
            "q\n"
        )

        run_command(
            "gmx make_ndx -f em.gro -o index.ndx",
            cwd=directory,
            input_text=ndx_input
        )

        print(f"✅ index.ndx built cleanly: {expected_group}")
        print("   ✔ Old ligand merge preserved")
        print("   ✔ Protein-only logic kept separate")
        print("   ✔ tc-grps can safely use this group")

    else:
        print("❌ FATAL: Could not determine simulation mode for Step 5.")
        print(f"   ligand_code={ligand_code}")
        print(f"   simulation_type={simulation_type}")
        exit(1)

    # ------------------------------------------------------------
    # 5D — Final validation
    # ------------------------------------------------------------
    index_path = os.path.join(directory, "index.ndx")
    if not os.path.exists(index_path):
        print("❌ FATAL: index.ndx was not created.")
        exit(1)

    missing_groups = []

    if not index_has_group(index_path, expected_group):
        missing_groups.append(expected_group)

    if not index_has_group(index_path, expected_water_group):
        missing_groups.append(expected_water_group)

    if missing_groups:
        print(f"❌ FATAL: Missing expected index group(s): {', '.join(missing_groups)}")
        print("📋 Groups present in index.ndx:")
        with open(index_path, "r") as f:
            for line in f:
                s = line.strip()
                if s.startswith("[") and s.endswith("]"):
                    print(f"   {s}")
        exit(1)

    print(f"✅ Verified Step 5 index groups exist: {expected_group}, {expected_water_group}\n")



    # ============================================================
    # 🌡️ 6. NVT Equilibration
    # ============================================================
    print("🌡️ [STEP 6] Starting NVT equilibration ...")

    # ------------------------------------------------------------
    # 6A — Determine the exact coupling group expected for this system
    # ------------------------------------------------------------
    if ligand_code:
        expected_group = f"Protein_{ligand_code.upper()}"
        coupling_group_str = f"{expected_group} Water_and_ions"
    else:
        expected_group = "Protein"
        coupling_group_str = "Protein Water_and_ions"

    print(f"🧭 Expected NVT coupling/index group: {expected_group}")

    # ------------------------------------------------------------
    # 6B — Force nvt.mdp to match the actual simulation mode
    # ------------------------------------------------------------
    nvt_mdp_path = os.path.join(directory, "nvt.mdp")
    standardize_mdp_groups(nvt_mdp_path, coupling_group_str)
    print(f"✅ Standardized nvt.mdp groups to: {coupling_group_str}")

    # ------------------------------------------------------------
    # 6C — Validate that index.ndx contains the groups nvt.mdp expects
    # ------------------------------------------------------------
    index_path = os.path.join(directory, "index.ndx")
    if not os.path.exists(index_path):
        print("❌ FATAL: index.ndx not found before NVT grompp.")
        exit(1)

    if not index_has_group(index_path, expected_group):
        print(f"❌ FATAL: Expected index group '{expected_group}' not found in index.ndx")
        print("   The NVT MDP file and index file are out of sync.")
        print("\n📋 Groups found in index.ndx:")
        with open(index_path, "r") as f:
            for line in f:
                s = line.strip()
                if s.startswith("[") and s.endswith("]"):
                    print(f"   {s}")
        exit(1)

    if not index_has_group(index_path, "Water_and_ions"):
        print("❌ FATAL: Expected index group 'Water_and_ions' not found in index.ndx")
        print("   The default solvent/ion group is missing from the generated index.")
        print("\n📋 Groups found in index.ndx:")
        with open(index_path, "r") as f:
            for line in f:
                s = line.strip()
                if s.startswith("[") and s.endswith("]"):
                    print(f"   {s}")
        exit(1)

    print(f"✅ Verified index groups exist: {expected_group}, Water_and_ions")

    # ------------------------------------------------------------
    # 6D — Build NVT TPR
    # ------------------------------------------------------------
    run_command_cpu(
        "gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 2",
        cwd=directory
    )

    # ------------------------------------------------------------
    # 6E — Run NVT with GPU→fallback handling
    # ------------------------------------------------------------
    run_mdrun_with_fallback(
        base_cmd="gmx mdrun -v -deffnm nvt",
        cwd=directory,
        pin_offset=pin_offset,
        threads=available_threads,
        gpu_id=gpu_id,
        use_gpu=use_gpu
    )

    print("✅ NVT equilibration complete.\n")


    # ============================================================
    # 💧 7. NPT Equilibration
    # ============================================================
    
    print("💧 [STEP 7] Starting NPT equilibration ...")

    run_command_cpu(
        "gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr -maxwarn 2",
        cwd=directory
    )

    # ✅ GPU→GPU-lite→CPU fallback
    run_mdrun_with_fallback(
        base_cmd="gmx mdrun -v -deffnm npt",
        cwd=directory,
        pin_offset=pin_offset,
        threads=available_threads,
        gpu_id=gpu_id,
        use_gpu=use_gpu
    )

    print("✅ NPT equilibration complete.\n")


    # ============================================================
    # 🧠 8. Production MD (resume from checkpoint if present)
    # ============================================================
    print("🧠 [STEP 8] Starting production MD ...")

    md_deffnm = "md_0_1"
    md_tpr = os.path.join(directory, f"{md_deffnm}.tpr")
    md_cpt = os.path.join(directory, f"{md_deffnm}.cpt")

    # build / rebuild TPR if needed
    if not (os.path.exists(md_tpr) and os.path.exists(md_cpt)):
        run_command_cpu(
            "gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx "
            f"-o {md_deffnm}.tpr -maxwarn 2",
            cwd=directory
        )

    # base mdrun (resume if checkpoint exists)
    base_cmd = f"gmx mdrun -v -deffnm {md_deffnm}"
    if os.path.exists(md_cpt) and os.path.exists(md_tpr):
        print(f"♻️ Found checkpoint: {md_deffnm}.cpt → resuming with -cpi/-append")
        base_cmd += f" -cpi {md_deffnm}.cpt -append"

    # GPU→fallback runner (virtual-site aware)
    run_mdrun_with_fallback(
        base_cmd=base_cmd,
        cwd=directory,
        pin_offset=pin_offset,
        threads=available_threads,
        gpu_id=gpu_id,
        use_gpu=use_gpu,
        tpr_check=f"{md_deffnm}.tpr"
    )

    print("✅ Production MD complete.\n")







if __name__ == "__main__":
    args = parse_args()
    print_numa_summary()

    # GPU selection
    try:
        gpu_list = subprocess.check_output(
            "nvidia-smi --query-gpu=name --format=csv,noheader", shell=True
        ).decode().strip().split('\n')
        gpu_list = [g for g in gpu_list if g]
        gpu_count = len(gpu_list)
    except Exception:
        gpu_list = []
        gpu_count = 0

    if args.compute.upper() == "GPU":
        if args.gpu is not None:
            gpu_id = args.gpu
        elif args.headless:
            # Headless GPU but no gpu id given -> default 0 if available, else fail to CPU
            gpu_id = 0 if gpu_count > 0 else -1
        else:
            if gpu_count == 0:
                print("❌ No GPU detected. Ensure NVIDIA drivers and CUDA are properly installed.")
                exit(1)
            elif gpu_count == 1:
                print("\n🎯 Detected 1 GPU. Automatically selecting GPU 0.")
                gpu_id = 0
            else:
                print(f"\n🧠 Detected {gpu_count} GPUs:")
                for idx, name in enumerate(gpu_list):
                    print(f"  [{idx}] {name}")
                gpu_choice = input(f"\nEnter GPU number (0-{gpu_count - 1}): ").strip()
                if not gpu_choice.isdigit() or int(gpu_choice) not in range(gpu_count):
                    print("❌ Invalid GPU selection.")
                    exit(1)
                gpu_id = int(gpu_choice)
    else:
        gpu_id = -1  # CPU mode


    simulation_directory = os.getcwd()
    print(f"\n📂 Running MD setup in: {simulation_directory} on GPU {gpu_id}")

    # Simulation type
    sim_map = {"1": "Ligand:Protein", "2": "Protein:Protein/Peptide", "3": "PROTAC", "4": "Just Protein"}
    if args.mode:
        mode_map = {"ligand": "Ligand:Protein", "protein": "Just Protein", "protac": "PROTAC"}
        simulation_type = mode_map.get(args.mode.lower(), args.mode)

    else:
        print("\nPlease select the simulation type:")
        print("  1. Ligand:Protein")
        print("  2. Protein:Protein/Peptide")
        print("  3. PROTAC")
        print("  4. Just Protein")
        sim_choice = input("Enter choice (1-4): ").strip()
        if sim_choice not in ["1", "2", "3", "4"]:
            print("❌ Invalid selection. Exiting.")
            exit(1)
        simulation_type = sim_map[sim_choice]

    # ─── Ligand code (auto-detect + confirm) ──────────────────
    ligand_code = None

    if simulation_type.lower() in ["ligand:protein", "ligand", "protac"]:

        # CLI always wins
        if args.ligand:
            ligand_code = args.ligand.upper()
            print(f"🧬 Using ligand from CLI: {ligand_code}")

        else:
            detected = detect_ligand_from_cgenff()

            if detected:
                ans = input(
                    f"🔍 Detected ligand '{detected}' from CGenFF files. "
                    "Use this ligand? [Y/n]: "
                ).strip().lower()

                if ans in ("", "y", "yes"):
                    ligand_code = detected
                else:
                    ligand_code = input(
                        "Enter the 3-letter ligand code: "
                    ).strip().upper()

            else:
                ligand_code = input(
                    "Enter the 3-letter ligand code: "
                ).strip().upper()

        if not ligand_code:
            print("❌ Ligand code is required for this simulation type.")
            exit(1)

    print(f"🧬 Final ligand selection: {ligand_code}")


    # Production length
    if args.ns is not None:
        simulation_time_ns = float(args.ns)
    else:
        length_input = input("How long should the production MD run (in nanoseconds)? [Press Enter for 50 ns]: ").strip()
        if length_input == "":
            simulation_time_ns = 50.0
        else:
            try:
                simulation_time_ns = float(length_input)
            except ValueError:
                print("ℹ️ Invalid input. Defaulting to 50 ns.")
                simulation_time_ns = 50.0
    if simulation_time_ns <= 0:
        print("ℹ️ Non-positive simulation time provided. Defaulting to 50 ns.")
        simulation_time_ns = 50.0
    # Determine compute mode
    use_gpu = (args.compute.upper() == "GPU")
    print(f"⚙️ Compute mode selected: {'GPU-accelerated' if use_gpu else 'CPU-only'}")


    setup_md(simulation_directory, gpu_id, ligand_code, simulation_type, simulation_time_ns, use_gpu, args)






