"""
====================================================================
⚙️  1_AutomateGromacs.py — Automated GROMACS MD Setup & Run Pipeline
====================================================================

This script automates the setup and execution of molecular dynamics (MD)
simulations involving one or more proteins and an optional small-molecule
ligand. It supports both **interactive** (guided prompts) and **headless**
(non-interactive, fully automated) operation.

--------------------------------------------------------------
🚀 QUICKSTART SUMMARY
--------------------------------------------------------------
# Protein–ligand system (headless):
python 1_AutomateGromacs.py \
    --pdb Model_2.pdb \
    --chain-map "A:Rap1B,B:Rap1GAP" \
    --ligand GDP

# Protein–protein system (no ligand):
python 1_AutomateGromacs.py \
    --pdb Model_2.pdb \
    --chain-names "Rap1B,Rap1GAP" \
    --mode protein

--------------------------------------------------------------
🧭 INTERACTIVE MODE
--------------------------------------------------------------
Run with no arguments to launch guided mode:
    python 1_AutomateGromacs.py

You will be prompted to:
  • Select a PDB file
  • Indicate whether the system includes a ligand
  • Enter the ligand’s 3-letter code (if applicable)
  • Assign custom names to protein chains

--------------------------------------------------------------
🤖 HEADLESS MODE (Recommended for pipelines)
--------------------------------------------------------------
Skip all prompts by providing arguments:

  --pdb <FILENAME>             Select input PDB file directly
  --chain-map "A:Name,B:Name"  Explicit chain naming (recommended)
  --chain-names "Name1,Name2"  Shorthand naming (ordered by chain order)
  --ligand GDP                 Ligand 3-letter code (optional)
  --mode [ligand|protein]      Skip system-type prompt entirely

🧠 NOTE:
If `--ligand` is provided, the script will automatically
switch to *ligand mode* — specifying `--mode ligand` is optional.

Only one of `--chain-map` or `--chain-names` is required.

--------------------------------------------------------------
🧩 EXAMPLES
--------------------------------------------------------------
# 1️⃣ Full headless run (Rap1B : Rap1GAP + GDP)
python 1_AutomateGromacs.py \
    --pdb Model_2.pdb \
    --chain-map "A:Rap1B,B:Rap1GAP" \
    --ligand GDP

# 2️⃣ Protein-only run (Rap1B : Rap1GAP)
python 1_AutomateGromacs.py \
    --pdb Model_2.pdb \
    --chain-names "Rap1B,Rap1GAP" \
    --mode protein

# 3️⃣ Minimal ligand-only quick test
python 1_AutomateGromacs.py \
    --pdb test.pdb \
    --ligand ATP

--------------------------------------------------------------
🧠 NOTES
--------------------------------------------------------------
- Detected chain IDs are cached in `.chainmap.json` for reuse.
- Ligand topologies are automatically generated via SILCSBio CGenFF.
- Ion concentration defaults to 0.15 M (override with):
      export GMX_SALT_M=0.10
- Compatible with CHARMM36 forcefields and GROMACS ≥2023.
- All setup steps are logged to `mdrun.log` in the working directory.
- A runtime banner automatically summarizes the chosen mode,
  e.g. “🚀 Running headless: Rap1B–Rap1GAP + GDP [ligand mode]”.

====================================================================


--------------------------------------------------------------
📦 OUTPUT SUMMARY
--------------------------------------------------------------

After a successful run, the following files and folders will be generated
in your working directory (`Model_X/` or equivalent):

🧬 STRUCTURE & TOPOLOGY FILES
─────────────────────────────
protein.pdb                 → Ligand-stripped protein (from original PDB)
protein_processed.gro       → Processed protein structure (via pdb2gmx)
topol.top                   → Master topology including forcefield and ligand entries
atomIndex.txt               → Chain index mapping (e.g., Rap1B 1–2773, Rap1GAP 2774–9550)
.chainmap.json              → Cached chain-to-name map for future reuse

💊 LIGAND FILES (only if --ligand provided)
──────────────────────────────────────────
<lig>.cgenff.mol2           → Ligand structure for CGenFF parameterization
<lig>.str                   → CGenFF output parameter stream
<lig>.itp / <lig>.prm       → Ligand topology and parameters (auto-generated)
<lig>_ini.pdb               → Ligand reference coordinates (for pose grafting)
<lig>_pose_match.pdb        → Pose-aligned ligand in correct coordinates
<lig>.gro                   → Final GROMACS-format ligand coordinates

⚙️ SYSTEM ASSEMBLY FILES
────────────────────────
complex.gro                 → Combined protein + ligand structure
newbox.gro                  → Simulation box (1.0 nm buffer, cubic)
solv.gro                    → Solvated box with water molecules
ions.tpr                    → Pre-ionization input structure
solv_ions.gro               → Final solvated, neutralized system with 0.15 M NaCl

🧠 DIAGNOSTIC & LOG FILES
──────────────────────────
mdrun.log                   → Master run log (appends key runtime info)
grompp.mdp / ions.mdp       → Input parameter files used for setup
*.itp / *.prm               → Topology includes generated by pdb2gmx and CGenFF
_water.ndx (if used)        → Backup water-only index for genion
index.ndx                   → System index groups used by GROMACS
index_default.ndx           → Default group definitions before merging
chainmap.json               → Cached chain-to-name mapping for reruns

✅ FINAL PREPARED INPUT
────────────────────────
em.tpr                      → Ready-to-run energy minimization input file

This marks the completion of the **Setup Phase**.
The next script (e.g., `2_AutomateMD.py` or `2_AutomateGromacsRun.py`)
continues with equilibration (NVT, NPT) and production MD.

====================================================================
"""
# === END OF 1_AutomateGromacs.py HEADER ===


import os
import subprocess
# Define number of threads for GROMACS execution
import multiprocessing
import shutil  # For file operations

import re


import argparse


parser = argparse.ArgumentParser(description="Automate MD setup with optional chain naming and ligand support.")

parser.add_argument("--mode", choices=["ligand", "protein"], help="Skip interactive mode selection (headless mode).")
parser.add_argument("--chain-names", type=str, help="Comma-separated custom chain names (e.g. Rap1B,Rap1GAP)")
parser.add_argument("--chain-map", type=str, help="Explicit chainID:name pairs (e.g. A:Rap1B,B:Rap1GAP)")
parser.add_argument("--ligand", "-l", type=str, help="Ligand 3-letter code (e.g., PTC). If omitted, runs in protein-only mode.")
parser.add_argument("--pdb", type=str, help="Input PDB filename (e.g., Model_2.pdb)")

args = parser.parse_args()
# ============================================================
# 🚀 Runtime Summary Banner
# ============================================================
if args.mode or args.ligand or args.chain_names or args.chain_map:
    chain_display = ""
    if args.chain_map:
        # Convert "A:Rap1B,B:Rap1GAP" → "Rap1B–Rap1GAP"
        chain_display = "–".join([p.split(":")[1] for p in args.chain_map.split(",")])
    elif args.chain_names:
        chain_display = "–".join(args.chain_names.split(","))
    else:
        chain_display = "UnnamedChains"

    if args.mode == "ligand" or args.ligand:
        mode_label = "ligand mode"
        ligand_label = args.ligand or "UnknownLigand"
        banner = f"🚀 Running headless: {chain_display} + {ligand_label} [{mode_label}]"
    elif args.mode == "protein":
        mode_label = "protein-only mode"
        banner = f"🚀 Running headless: {chain_display} [{mode_label}]"
    else:
        banner = f"🚀 Running headless: {chain_display} [auto mode]"

    print(f"\n{banner}\n")
    # Optional: also log it to mdrun.log for traceability
    with open(os.path.join(os.getcwd(), "mdrun.log"), "a") as log:
        log.write(f"{banner}\n")


def comment_out_faulty_dihedral(log_text, directory):
    """
    Detects and comments out a faulty dihedral entry in the specified .itp file based on grompp error.
    Returns True if successful, False otherwise.
    """
    # Match error message like:
    # ERROR 1 [file topol_Protein_chain_A.itp, line 15833]
    match = re.search(r"ERROR \d+ \[file (\S+), line (\d+)\]", log_text)
    if not match:
        print("❌ No faulty dihedral line found in error log.")
        return False

    itp_file = match.group(1)
    line_num = int(match.group(2))

    itp_path = os.path.join(directory, itp_file)
    if not os.path.exists(itp_path):
        print(f"❌ Could not find .itp file: {itp_path}")
        return False

    print(f"🩹 Commenting out line {line_num} in {itp_file}...")
    with open(itp_path, "r") as file:
        lines = file.readlines()
    if line_num - 1 < len(lines):
        lines[line_num - 1] = "; " + lines[line_num - 1]  # Comment out
    else:
        print("❌ Line number out of range in file.")
        return False
    with open(itp_path, "w") as file:
        file.writelines(lines)
    print(f"✅ Successfully commented out problematic dihedral in {itp_file}")
    return True

def fix_missing_atoms(input_pdb, output_pdb, ph=7.4):
    """Uses PDBFixer to fill missing atoms and add hydrogens to the PDB file."""
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile

        fixer = PDBFixer(filename=input_pdb)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(pH=ph)
        PDBFile.writeFile(fixer.topology, fixer.positions, open(output_pdb, 'w'))
        print(f"✅ Fixed missing atoms and saved cleaned file to: {output_pdb}")
        return True
    except Exception as e:
        print(f"❌ Failed to fix atoms in {input_pdb}: {str(e)}")
        return False

def modify_topology_file(topol_path, ligand_code):
    """Updates the topology file by inserting required lines at the correct positions."""
    # Modified: ensure ligand parameters are included dynamically if ligand is present
    with open(topol_path, "r") as f:
        lines = f.readlines()
    # If ligand is present, include its parameter and topology files
    if ligand_code:
        # Find insertion point after forcefield includes
        insertion_index = None
        for i, line in enumerate(lines):
            if "; Include forcefield parameters" in line:
                insertion_index = i + 2
                break
        if insertion_index is not None:
            prm_line = f'#include "{ligand_code.lower()}.prm"\n'
            itp_line = f'; Include ligand topology\n#include "{ligand_code.lower()}.itp"\n'
            if prm_line not in lines[insertion_index:]:
                lines.insert(insertion_index, prm_line)
                lines.insert(insertion_index + 1, itp_line)
        # Append ligand molecule entry at end of file if not present
        ligand_entry = f"{ligand_code:<20}1\n"
        if not any(line.strip().startswith(ligand_code) for line in lines):
            lines.append(ligand_entry)
        print(f"✅ Added ligand includes and entry for {ligand_code} in topol.top")
    else:
        print("ℹ️ No ligand present; skipping ligand-specific topology modifications.")
    # Write changes back to topology file
    with open(topol_path, "w") as f:
        f.writelines(lines)

def merge_gro_files(protein_gro, ligand_gro, complex_gro):
    """Merges protein_processed.gro and ligand .gro into complex.gro while preserving the last row."""
    if not os.path.exists(protein_gro) or not os.path.exists(ligand_gro):
        print(f"❌ ERROR: Missing input files: {protein_gro} or {ligand_gro}")
        return
    # Read input .gro files
    with open(protein_gro, "r") as f:
        protein_lines = f.readlines()
    with open(ligand_gro, "r") as f:
        ligand_lines = f.readlines()
    correct_box = protein_lines[-1].strip()  # Box from protein
    ligand_body = ligand_lines[2:-1]         # Ligand atoms (skip header/count/box)
    total_atoms = int(protein_lines[1].strip()) + int(ligand_lines[1].strip())
    merged_lines = []
    merged_lines.append(protein_lines[0])                # Title
    merged_lines.append(f"{total_atoms}\n")              # Updated atom count
    merged_lines.extend(protein_lines[2:-1])             # Protein atoms
    merged_lines.extend(ligand_body)                     # Ligand atoms
    merged_lines.append(correct_box + "\n")              # Box dimensions
    with open(complex_gro, "w") as f:
        f.writelines(merged_lines)
    # Verify last line
    with open(complex_gro, "r") as f:
        final_lines = f.readlines()
        if final_lines[-1].strip() != correct_box:
            print("❌ ERROR: Final row of complex.gro does not match protein_processed.gro!")
            exit(1)
    print(f"✅ complex.gro merged successfully with {total_atoms} atoms.")

def run_command(command, cwd=None, input_text=None):
    """Runs a shell command in the given directory with optional input and exits on failure."""
    print(f"\n🔹 Running command: {command}")
    print(f"📂 In directory: {cwd if cwd else os.getcwd()}")
    try:
        result = subprocess.run(command, shell=True, cwd=cwd, text=True, capture_output=True, input=input_text)
        print(f"📜 STDOUT:\n{result.stdout}")
        if result.returncode != 0:
            print(f"❌ Error running command: {command}")
            print(f"🔺 STDERR: {result.stderr}")
            print("⚠️ Exiting script due to failure.")
            exit(1)
    except Exception as e:
        print(f"❌ Exception occurred while running command: {command}\n{str(e)}")
        print("⚠️ Exiting script due to failure.")
        exit(1)
    return True

def detect_forcefield(directory):
    """Detects if a CHARMM forcefield directory exists in the working directory."""
    for item in os.listdir(directory):
        if os.path.isdir(os.path.join(directory, item)) and item.startswith("charmm36"):
            forcefield_path = os.path.join(directory, item)
            print(f"✅ CHARMM forcefield detected at {forcefield_path}. Selecting automatically.")
            return forcefield_path
    print("❗ CHARMM forcefield not found in the current directory. User must select manually.")
    return None

def select_pdb_file(directory):
    """Prompts the user to select a PDB file from available options."""
    pdb_files = [f for f in os.listdir(directory) if f.endswith(".pdb")]
    if not pdb_files:
        print("❌ No PDB files found in the current directory.")
        exit(1)
    print("\n📁 Available PDB files:")
    for idx, file in enumerate(pdb_files, 1):
        print(f"{idx}. {file}")
    while True:
        selection = input("\n🔍 Select a PDB file by number (e.g., 1): ").strip()
        if selection.isdigit() and 1 <= int(selection) <= len(pdb_files):
            return pdb_files[int(selection) - 1]
        else:
            print("❗ Invalid selection. Please enter a valid number.")

def get_ligand_code():
    """Prompts user to enter ligand 3-letter code or use default."""
    choice = input("\n💊 Enter ligand 3-letter code or press ENTER to use default [PTC]: ").strip()
    return choice.upper() if choice else "PTC"

def automate_pdb2gmx(directory):
    """Automates GROMACS pdb2gmx selection for force field and termini settings."""
    omp_threads = os.environ.get("OMP_NUM_THREADS") or str(multiprocessing.cpu_count())
    os.environ["OMP_NUM_THREADS"] = omp_threads
    print(f"🧠 Using OMP_NUM_THREADS = {omp_threads}")
    protein_pdb = os.path.join(directory, "protein.pdb")
    print("\n🔍 Checking for protein.pdb in the working directory...")
    if not os.path.exists(protein_pdb):
        print(f"❌ Skipping {directory}: protein.pdb not found.")
        return False
    # Detect force field
    forcefield_path = detect_forcefield(directory)
    if forcefield_path is None:
        print("\n❌ CHARMM forcefield not found. Manual selection required.")
        return False
    print(f"✅ Using CHARMM force field from: {forcefield_path}")
    forcefield_option = "1"  # CHARMM
    chain_count = count_chains(protein_pdb)
    print(f"🔍 Detected {chain_count} chains in the protein.")
    if chain_count == 0:
        print("❌ No chains detected. Check protein.pdb formatting.")
        return False
    input_selections = [forcefield_option]
    for i in range(chain_count):
        input_selections.append("1" if i == 0 else "0")  # Start terminus: first chain NH3+, others PRO-NH2+
        input_selections.append("0")                    # End terminus: COO-
    input_text = "\n".join(input_selections) + "\n"
    print("\n✅ Auto-selecting force field and termini for pdb2gmx.")
    return run_command(f"gmx pdb2gmx -f protein.pdb -o protein_processed.gro -p topol.top -ignh -water spc -ter", cwd=directory, input_text=input_text)

def count_chains(protein_pdb):
    """Counts the number of unique chains in a PDB file."""
    chains = set()
    with open(protein_pdb, "r") as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain_id = line[21]  # Chain ID column
                chains.add(chain_id)
    return len(chains)

import os, json

def generate_atom_index_file(directory, chain_names_arg=None, chain_map_arg=None):
    """
    Generates atomIndex.txt with chain names.
    Priority:
      1. If chain_map_arg or chain_names_arg provided (from CLI or env), use that.
      2. If saved .chainmap.json exists, reuse mapping.
      3. Else prompt interactively.
    """
    protein_pdb = os.path.join(directory, "protein.pdb")
    protein_gro = os.path.join(directory, "protein_processed.gro")
    output_file = os.path.join(directory, "atomIndex.txt")
    cache_file = os.path.join(directory, ".chainmap.json")

    if not os.path.exists(protein_pdb) or not os.path.exists(protein_gro):
        print(f"❌ Missing required files (protein.pdb or protein_processed.gro) in {directory}")
        return False

    # --- Detect chain IDs
    chain_ids = []
    with open(protein_pdb) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                cid = line[21].strip() or "A"
                if cid not in chain_ids:
                    chain_ids.append(cid)

    if not chain_ids:
        print("⚠️ No chain IDs found in PDB. Defaulting to single-chain [A].")
        chain_ids = ["A"]

    print(f"🔍 Detected {len(chain_ids)} chain(s): {', '.join(chain_ids)}")

    # --- Priority 1: CLI overrides
    chain_map = {}
    if chain_map_arg:
        for pair in chain_map_arg.split(","):
            if ":" in pair:
                cid, name = pair.split(":", 1)
                chain_map[cid.strip()] = name.strip()
    elif chain_names_arg:
        names = [n.strip() for n in chain_names_arg.split(",")]
        for i, cid in enumerate(chain_ids):
            if i < len(names):
                chain_map[cid] = names[i]

    # --- Priority 2: existing .chainmap.json cache
    elif os.path.exists(cache_file):
        with open(cache_file) as f:
            saved = json.load(f)
        for cid in chain_ids:
            if cid in saved:
                chain_map[cid] = saved[cid]
        print(f"💾 Loaded existing chain map: {chain_map}")

    # --- Priority 3: interactive prompt
    # If provided via CLI, skip prompts entirely
    if args.chain_map:
        try:
            for pair in args.chain_map.split(","):
                cid, name = [p.strip() for p in pair.split(":")]
                chain_map[cid] = name
            print(f"✅ Chain map provided via CLI: {chain_map}")
        except Exception as e:
            print(f"⚠️ Failed to parse chain map '{args.chain_map}': {e}")
    elif args.chain_names:
        names = [n.strip() for n in args.chain_names.split(",")]
        for i, cid in enumerate(chain_ids):
            if i < len(names):
                chain_map[cid] = names[i]
        print(f"✅ Chain names provided via CLI: {chain_map}")
    else:
        print("💬 No chain mapping provided; entering interactive chain naming mode.")
        chain_map = {}

        for cid in chain_ids:
            while True:
                name = input(f"🧩 Enter descriptive name for chain {cid} (e.g. Rap1B): ").strip()
                if name:
                    chain_map[cid] = name
                    break
                else:
                    print("⚠️ Chain name cannot be empty. Please try again.")

        print("\n✅ You entered the following chain map:")
        for cid, name in chain_map.items():
            print(f"   {cid} → {name}")
        confirm = input("Confirm these names? [Y/n]: ").strip().lower()
        if confirm == "n":
            print("🔁 Restarting chain naming...")
            return generate_atom_index_file(directory)  # rerun prompt

    # else:
    #     print("💬 No chain mapping provided; running fallback auto-naming.")
    #     chain_map = {cid: f"Protein{i+1}" for i, cid in enumerate(chain_ids)}


    # Fallback default labels
    for i, cid in enumerate(chain_ids):
        if cid not in chain_map:
            chain_map[cid] = f"Protein{i+1}"

    print("✅ Final chain map:")
    for cid, name in chain_map.items():
        print(f"   {cid} → {name}")

    # --- Save cache for future reuse
    with open(cache_file, "w") as f:
        json.dump(chain_map, f, indent=2)

    # --- Parse GRO atom index ranges
    with open(protein_gro) as gro:
        lines = gro.readlines()

    atom_ranges = []
    atom_counter = 1
    last_resid = None
    current_start = 1
    chain_index = 0

    for i, line in enumerate(lines[2:-1]):
        try:
            resid = int(line[:5].strip())
        except ValueError:
            continue
        atom_counter = i + 1
        if last_resid and resid < last_resid:
            atom_ranges.append((chain_ids[chain_index], current_start, atom_counter))
            current_start = atom_counter + 1
            chain_index = min(chain_index + 1, len(chain_ids) - 1)
        last_resid = resid
    atom_ranges.append((chain_ids[chain_index], current_start, atom_counter + 1))

    # --- Write atomIndex.txt
    with open(output_file, "w") as f:
        for cid, start, end in atom_ranges:
            label = chain_map.get(cid, f"Protein{cid}")
            f.write(f"{label} {start}-{end}\n")

    print(f"✅ Generated atomIndex.txt in {directory}")
    return True


# ======== NEW HELPERS FOR REALTIME CGenFF PREP ========

def ensure_silcsbio_env():
    """Verify SILCSBio is on this machine and basic files exist."""
    home = os.environ.get("SILCSBIO_HOME", "")
    if not home or not os.path.isdir(home):
        print("❌ SILCSBIO_HOME is not set or directory missing. Make sure your ~/.bashrc SILCSBio setup is loaded.")
        return None, None, None
    cgenff_exec = os.path.join(home, "cgenff", "cgenff")
    prm36      = os.path.join(home, "cgenff", "par_all36_cgenff.prm")
    if not os.path.isfile(cgenff_exec):
        print("❌ Missing SILCSBio cgenff executable:", cgenff_exec); return None, None, None
    if not os.path.isfile(prm36):
        print("❌ Missing par_all36_cgenff.prm:", prm36); return None, None, None
    return cgenff_exec, prm36, home

def extract_ligand_from_pdb(pdb_path, ligand_code, out_pdb):
    """Extracts only the specified 3-letter ligand (HETATM) into out_pdb."""
    ligand_code = ligand_code.upper()
    wrote = False
    with open(pdb_path, "r") as fin, open(out_pdb, "w") as fout:
        for line in fin:
            if line.startswith("HETATM") and line[17:20].strip().upper() == ligand_code:
                fout.write(line); wrote = True
    if not wrote:
        print(f"❌ Could not find ligand {ligand_code} in {pdb_path}")
        return False
    print(f"✅ Extracted ligand {ligand_code} to {out_pdb}")
    return True

def ligand_has_hydrogens(pdb_file):
    """Detect if ligand PDB has hydrogens."""
    with open(pdb_file, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                aname = line[12:16].strip()
                elem  = line[76:78].strip() if len(line) >= 78 else ""
                if aname.upper().startswith("H") or elem.upper() == "H":
                    return True
    return False



def autodetect_ligands(pdb_path):
    """
    Returns a sorted list of unique ligand 3-letter codes in the PDB
    (excluding water, ions, and amino acids).
    """
    # ignore common non-standard atoms
    ignore = {
        "HOH", "WAT", "SOL",  # water
        "NA", "K", "CL", "CA", "MG", "ZN", "FE", "MN", "CU", "CO", "BR", "IOD",
        "SO4", "PO4", "PI", "ACE"  # ions/buffers
    }

    # amino acids, nucleotides, common cofactors
    ignore.update({
        "ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE",
        "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
        "DA","DT","DG","DC","A","T","G","C","U"
    })

    ligs = set()

    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("HETATM"):
                resname = line[17:20].strip()
                if resname.upper() not in ignore:
                    ligs.add(resname.upper())

    return sorted(ligs)




def graft_coords_onto_ini(ini_pdb, pose_pdb, out_pdb):
    """
    Preserve the atom ORDER from ini_pdb (matches .itp) but replace XYZ with coords
    from pose_pdb (extracted from complex). We match records by atom name within the only residue.
    """
    def parse_pdb_atoms(path):
        atoms = []
        with open(path) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    name = line[12:16].strip()
                    resn = line[17:20].strip()
                    resi = line[22:26].strip()
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                    atoms.append({"raw": line, "name": name, "resn": resn, "resi": resi, "x": x, "y": y, "z": z})
        return atoms

    ini_atoms  = parse_pdb_atoms(ini_pdb)
    pose_atoms = parse_pdb_atoms(pose_pdb)

    # Build simple name->list mapping for pose (to handle duplicates like H11/H12 etc.)
    from collections import defaultdict, deque
    by_name = defaultdict(deque)
    for a in pose_atoms:
        by_name[a["name"]].append(a)

    out_lines = []
    with open(ini_pdb) as f:
        ini_lines = f.readlines()

    for line in ini_lines:
        if line.startswith(("ATOM", "HETATM")):
            name = line[12:16].strip()
            if not by_name[name]:
                # fallback: if name missing, keep original ini coords (warn)
                # You could also relax to element-based matching here.
                out_lines.append(line)
                continue
            src = by_name[name].popleft()
            # write XYZ into fixed-column PDB fields (30-54)
            new = (line[:30] + f"{src['x']:8.3f}{src['y']:8.3f}{src['z']:8.3f}" + line[54:])
            out_lines.append(new)
        else:
            out_lines.append(line)

    with open(out_pdb, "w") as f:
        f.writelines(out_lines)
    return True






def build_cgenff_inputs_realtime(directory, ligand_code, source_complex_pdb):
    """
    1) Extract ligand from complex PDB -> {code}_raw.pdb
    2) Add hydrogens if needed and normalize via Open Babel -> {code}_h.pdb
    3) Convert to MOL2 with gen3d -> {code}.cgenff.mol2
    4) Set MOL2 molecule name to ligand_code
    5) Run SILCSBio cgenff (stdout -> {code}.str)
    Returns (mol2_file, str_file) or (None, None) on failure.
    """
    cgenff_exec, prm36, silcs_home = ensure_silcsbio_env()
    if not cgenff_exec:
        return None, None

    base    = ligand_code.lower()
    raw_pdb = os.path.join(directory, f"{base}_raw.pdb")
    h_pdb   = os.path.join(directory, f"{base}_h.pdb")
    mol2    = os.path.join(directory, f"{ligand_code}.cgenff.mol2")  # keep your expected name
    strout  = os.path.join(directory, f"{ligand_code}.str")

    # 1) Extract ligand from the complex PDB
    if not extract_ligand_from_pdb(os.path.join(directory, source_complex_pdb), ligand_code, raw_pdb):
        return None, None

    # 2) Add H if missing (or normalize anyway)
    if not ligand_has_hydrogens(raw_pdb):
        ok = run_command(f"obabel {os.path.basename(raw_pdb)} -O {os.path.basename(h_pdb)} -h", cwd=directory)
    else:
        ok = run_command(f"obabel {os.path.basename(raw_pdb)} -O {os.path.basename(h_pdb)}", cwd=directory)
    if not ok:
        return None, None

    # 3) Convert to MOL2 (gen 3D if needed)
    ok = run_command(f"obabel {os.path.basename(h_pdb)} -O {os.path.basename(mol2)}", cwd=directory)
    if not ok:
        return None, None

    # 4) Force the MOL2 molecule name to be the desired ligand_code
    #    (cgenff uses this name as the RESI name in the .str)
    try:
        mol2_lines = []
        with open(mol2, "r") as f:
            lines = f.readlines()
        i = 0
        while i < len(lines):
            mol2_lines.append(lines[i])
            if lines[i].strip().upper() == "@<TRIPOS>MOLECULE" and i + 1 < len(lines):
                # replace the very next line with ligand_code
                lines[i+1] = f"{ligand_code}\n"
                mol2_lines.append(lines[i+1])
                mol2_lines.extend(lines[i+2:])
                break
            i += 1
        else:
            # If no MOLECULE section found, just prepend one
            mol2_lines = ["@<TRIPOS>MOLECULE\n", f"{ligand_code}\n"] + lines
        with open(mol2, "w") as f:
            f.writelines(mol2_lines)
    except Exception as e:
        print(f"❌ Failed to set MOL2 molecule name: {e}")
        return None, None

    # 5) Run cgenff with ONLY the mol2 file and redirect stdout to .str
    cgenff_cmd = f"{cgenff_exec} {os.path.basename(mol2)} > {os.path.basename(strout)}"

    print("\n🧪 Running CGenFF parameterization...")
    print(f"⚙️ Command: {cgenff_cmd}")

    # Execute CGenFF manually — do NOT exit() yet
    result = subprocess.run(
        cgenff_cmd, shell=True, cwd=directory,
        capture_output=True, text=True
    )

    # SHOW ERRORS (stdout normally empty due to redirection)
    print(f"📜 STDERR:\n{result.stderr}")

    # FAILED?
    if result.returncode != 0:
        print("\n❌ **CGENFF PARAMETERIZATION FAILED**")
        print("   CGenFF could not assign parameters for this ligand.")
        print("   Most common causes:")
        print("   • unusual bonding or ring systems")
        print("   • missing hydrogens or incorrect valence")
        print("   • OpenBabel generated invalid geometry")
        print("   • unsupported atom types in SILCSBio CGenFF")

        if os.path.exists(strout):
            print("\n📄 Partial .str file produced:")
            print(open(strout).read())

        print("⚠️ Returning (None, None) to trigger fallback.\n")
        return None, None

    # CHECK EMPTY .str
    if not os.path.exists(strout) or os.path.getsize(strout) == 0:
        print("\n❌ CGENFF produced an EMPTY .str file — no parameters were generated.")
        print("⚠️ This means the ligand is unsupported or malformed.")
        return None, None

    print(f"✅ CGENFF successfully generated {os.path.basename(strout)}")
    print(f"✅ Generated {os.path.basename(mol2)} and {os.path.basename(strout)}")

    return mol2, strout



def parse_make_ndx_for_water(stdout_text):
    """
    Robustly parse `gmx make_ndx` listing to find the group number for SOL/Water/WAT.
    Returns the group number as a string if found, else None.
    Lines look like:
      15 SOL                 : 55782 atoms
    or
      14 Water               : 55782 atoms
    """
    candidates_exact = {"SOL", "Water", "WAT"}
    # regex: capture group number and group name up to the colon
    pat = re.compile(r'^\s*(\d+)\s+([A-Za-z0-9_\-\+ ]+?)\s*:', re.MULTILINE)
    matches = pat.findall(stdout_text or "")
    if not matches:
        return None

    # 1) exact priority: SOL > Water > WAT
    for idx, name in matches:
        name_clean = name.strip()
        if name_clean == "SOL":
            return idx
    for idx, name in matches:
        if name.strip() == "Water":
            return idx
    for idx, name in matches:
        if name.strip() == "WAT":
            return idx

    # 2) contains-based fallback
    for idx, name in matches:
        name_u = name.strip().upper()
        if any(tag in name_u for tag in ("SOL", "WATER", "WAT", "TIP3", "TIP4", "TIP5", "HOH")):
            return idx

    # 3) last resort: pick the largest water-like group by atom count (if parsable)
    #   (parse "<idx> <name> : <count> atoms")
    pat_count = re.compile(r'^\s*(\d+)\s+([A-Za-z0-9_\-\+ ]+?)\s*:\s*(\d+)\s+atoms', re.MULTILINE)
    best = None
    best_count = -1
    for idx, name, cnt in pat_count.findall(stdout_text or ""):
        name_u = name.strip().upper()
        if any(tag in name_u for tag in ("SOL", "WATER", "WAT", "TIP3", "TIP4", "TIP5", "HOH")):
            c = int(cnt)
            if c > best_count:
                best_count = c
                best = idx
    return best


def build_water_index_with_gmx_select(directory):
    """
    Fallback: build a dedicated water-only index file using gmx select.
    Returns (ndx_path, group_number) or (None, None) on failure.
    We use ions.tpr (exists after the ions grompp) and solv.gro.
    """
    ndx_path = os.path.join(directory, "_water.ndx")
    # select resnames commonly used for water:
    selection = 'resname SOL or resname WAT or resname HOH or resname TIP3 or resname TIP4 or resname TIP5'
    ok = run_command(
        f'gmx select -s ions.tpr -f solv.gro -select "{selection}" -on {os.path.basename(ndx_path)}',
        cwd=directory
    )
    if not ok:
        return None, None

    # gmx select typically writes one group named "Selection" -> group 0
    # (If multiple groups were made, the first is still a safe default)
    try:
        with open(ndx_path, "r") as f:
            content = f.read()
        # sanity check it actually contains a group
        if "[ " in content and "]" in content:
            return ndx_path, "0"
    except Exception:
        pass
    return None, None


def _mean(vals):
    return sum(vals)/len(vals) if vals else None

def centroid_from_gro(path, first=2, last=-1):
    xs = []; ys = []; zs = []
    try:
        with open(path) as f:
            lines = f.readlines()[first:last]
    except FileNotFoundError:
        return (None, None, None)
    for ln in lines:
        try:
            xs.append(float(ln[20:28])); ys.append(float(ln[28:36])); zs.append(float(ln[36:44]))
        except:
            pass
    return (_mean(xs), _mean(ys), _mean(zs))


# ======== END NEW HELPERS ========



def process_directory(directory, pdb_filename, ligand_code):
    """Processes a directory for MD setup using a selected PDB file and ligand code."""
    print(f"\n🔄 Processing directory: {directory}")
    pdb_path = os.path.join(directory, pdb_filename)
    print(f"\n🔍 Checking for {pdb_filename} in the working directory...")
    if not os.path.exists(pdb_path):
        print(f"⚠️ Skipping {directory}: {pdb_filename} not found.")
        return False

    # Step 1: Remove ligand or copy original PDB based on presence of ligand
    if ligand_code:
        print(f"🛠️ Removing ligand ({ligand_code}) to create protein.pdb...")
        if not run_command(f"grep -v '{ligand_code}' {pdb_filename} > protein.pdb", cwd=directory):
            return False
    else:
        print("🛠️ No ligand provided. Using the original PDB as protein.pdb")
        try:
            shutil.copy(pdb_path, os.path.join(directory, "protein.pdb"))
            print("✅ Copied original PDB to protein.pdb")
        except Exception as e:
            print(f"❌ Failed to copy {pdb_filename} to protein.pdb: {e}")
            return False

    # Step 1.5: Clean up missing atoms and sidechains
    protein_pdb_path = os.path.join(directory, "protein.pdb")
    fixed_pdb_path = os.path.join(directory, "protein_fixed.pdb")
    print(f"🧼 Fixing missing atoms and adding hydrogens with PDBFixer...")
    if not fix_missing_atoms(protein_pdb_path, fixed_pdb_path):
        return False
    os.replace(fixed_pdb_path, protein_pdb_path)

    # Step 2: Generate protein topology
    if not automate_pdb2gmx(directory):
        return False

    # Step 3: Generate atom index file
    if not generate_atom_index_file(base_directory,
                                chain_names_arg=args.chain_names,
                                chain_map_arg=args.chain_map):
        print("⚠️ Failed to generate atom index file.")

        return False
    
    # Step 4: Generate ligand topology (if applicable)
    if ligand_code:
        print("\n🛠️ Generating ligand topology...")
        forcefield_path = detect_forcefield(directory)
        if forcefield_path is None:
            return False

        try:
            # --- PRIMARY: Full SILCSBio + pose-preserving route ---
            mol2_file = os.path.join(directory, f"{ligand_code}.cgenff.mol2")
            str_file  = os.path.join(directory, f"{ligand_code}.str")

            # 4a) Ensure we HAVE the CGenFF inputs (.mol2 + .str)
            if not (os.path.exists(mol2_file) and os.path.exists(str_file)):
                print("ℹ️ ParamChem files not found — generating with local CGenFF from the input PDB...")
                built_mol2, built_str = build_cgenff_inputs_realtime(directory, ligand_code, pdb_filename)
                if not (built_mol2 and built_str):
                    raise RuntimeError("Local CGenFF generation failed — will trigger fallback.")
                mol2_file, str_file = built_mol2, built_str

            # 4b) Convert (.mol2 + .str) → .itp/.prm/_ini.pdb
            if not run_command(
                f"python cgenff_charmm2gmx_py3_nx2.py {ligand_code} "
                f"{os.path.basename(mol2_file)} {os.path.basename(str_file)} {forcefield_path}",
                cwd=directory
            ):
                raise RuntimeError("CGenFF to GROMACS conversion failed.")

            # 4c) Pose-matched coordinates (keep ligand where it is in original PDB)
            ini_pdb_file     = f"{ligand_code.lower()}_ini.pdb"
            pose_source_pdb  = f"{ligand_code.lower()}_h.pdb"          # created by build_cgenff_inputs_realtime
            pose_pdb         = f"{ligand_code.lower()}_pose_match.pdb" # grafted, but may be missing hydrogens
            ligand_gro_file  = f"{ligand_code.lower()}.gro"

            if not graft_coords_onto_ini(ini_pdb_file, pose_source_pdb, pose_pdb):
                raise RuntimeError("Pose grafting failed.")

            # ------------------------------------------------------------
            # 🔧 NEW PATCH: Ensure final pose has hydrogens before editconf
            # ------------------------------------------------------------
            pose_pdb_path = os.path.join(directory, pose_pdb)
            pose_h_pdb    = os.path.join(directory, f"{ligand_code.lower()}_pose_h.pdb")

            if not ligand_has_hydrogens(pose_pdb_path):
                print("🔧 Adding hydrogens to pose-matched ligand...")
                if not run_command(
                    f"obabel {os.path.basename(pose_pdb)} -O {os.path.basename(pose_h_pdb)} -h",
                    cwd=directory
                ):
                    raise RuntimeError("Failed to add hydrogens to pose-matched ligand.")
                pose_final = pose_h_pdb
            else:
                pose_final = pose_pdb_path

            # 4d) Generate final ligand .gro with hydrogens guaranteed
            if not run_command(
                f"gmx editconf -f {os.path.basename(pose_final)} -o {ligand_gro_file}",
                cwd=directory
            ):
                raise RuntimeError("Failed to generate ligand GRO file.")

            # Sanity centroid check
            prot_ctr = centroid_from_gro(os.path.join(directory, "protein_processed.gro"))
            lig_ctr  = centroid_from_gro(os.path.join(directory, ligand_gro_file))
            print(f"🔎 Centroids (nm) — protein: {prot_ctr}, ligand: {lig_ctr}")

        except Exception as e:
            # --- FALLBACK: Simple non-SILCSBio route ---
            print(f"⚠️ CGenFF/SILCSBio route failed ({e}). Falling back to basic ligand conversion...")

            mol2_file = f"{ligand_code}.cgenff.mol2"
            str_file  = f"{ligand_code}.str"
            ini_pdb_file = f"{ligand_code.lower()}_ini.pdb"
            ligand_gro_file = f"{ligand_code.lower()}.gro"

            if not run_command(
                f"python cgenff_charmm2gmx_py3_nx2.py {ligand_code} "
                f"{mol2_file} {str_file} {forcefield_path}",
                cwd=directory
            ):
                return False

            # Fallback does not preserve pose — but still needs hydrogens
            # editconf from _ini.pdb directly
            if not run_command(f"gmx editconf -f {ini_pdb_file} -o {ligand_gro_file}", cwd=directory):
                return False

            print(f"✅ Fallback ligand topology generated successfully for {ligand_code}.")

    else:
        print("ℹ️ No ligand present, skipping ligand topology generation.")



    # Step 5: Merge protein and ligand structures
    if ligand_code:
        print("\n🔄 Merging protein and ligand structures...")
        protein_gro = os.path.join(directory, "protein_processed.gro")
        ligand_gro = os.path.join(directory, ligand_gro_file)
        complex_gro = os.path.join(directory, "complex.gro")
        merge_gro_files(protein_gro, ligand_gro, complex_gro)
    else:
        print("ℹ️ Skipping merge step (protein-only system). Using protein structure as complex.")
        protein_gro = os.path.join(directory, "protein_processed.gro")
        complex_gro = os.path.join(directory, "complex.gro")
        shutil.copy(protein_gro, complex_gro)
        print("✅ Single protein structure saved as complex.gro")

    # Step 6: Modify topology
    if ligand_code:
        print("\n🛠️ Modifying topology file (topol.top)...")
        topol_path = os.path.join(directory, "topol.top")
        modify_topology_file(topol_path, ligand_code)
    else:
        print("ℹ️ Skipping topology file modifications for ligand.")

    # Step 7–11: GROMACS setup steps
    print("\n⚙️ Running GROMACS setup steps...")
    if not run_command("gmx editconf -f complex.gro -o newbox.gro -bt cubic -d 1.0", cwd=directory):
        return False
    if not run_command("gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro", cwd=directory):
        return False

    # Step 9: Prepare ions
    print("⚙️ Step 9: Running grompp for ion preparation...")
    grompp_result = subprocess.run(
        "gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2",
        shell=True, cwd=directory, text=True, capture_output=True
    )
    if grompp_result.returncode != 0:
        print("⚠️ Initial grompp for ions failed. Attempting to patch .itp file...")
        print("🔍 Analyzing grompp STDERR:\n", grompp_result.stderr)
        if comment_out_faulty_dihedral(grompp_result.stderr, directory):
            print("🔁 Retrying grompp after patching...")
            if not run_command("gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2", cwd=directory):
                return False
        else:
            print("❌ Could not automatically fix the grompp issue.")
            return False


    # Determine solvent (water) group index for genion automatically
    make_ndx_result = subprocess.run(
        "gmx make_ndx -f solv.gro",
        shell=True, cwd=directory, text=True, input="q\n", capture_output=True
    )

    sol_group = None
    ndx_arg = ""  # by default, we won't pass a custom index file

    if make_ndx_result.returncode == 0:
        # Robust parse of the printed group list
        sol_group = parse_make_ndx_for_water(make_ndx_result.stdout)

    if sol_group is None:
        print("❗ Could not confidently detect SOL group from make_ndx listing. Trying gmx select fallback...")
        # Build a dedicated water-only index and use that for genion
        ndx_path, group0 = build_water_index_with_gmx_select(directory)
        if ndx_path and group0 is not None:
            ndx_arg = f"-n {os.path.basename(ndx_path)}"
            sol_group = group0

    if sol_group is None:
        # As a final fallback, show user the printed groups and ask
        print("❌ Could not determine SOL group automatically.")
        print("📋 Here's what make_ndx returned:\n", make_ndx_result.stdout)
        sol_group = input("🔢 Please enter the group number corresponding to SOL/Water manually: ").strip()

    # make_ndx_result = subprocess.run("gmx make_ndx -f solv.gro", shell=True, cwd=directory, text=True, input="q\n", capture_output=True)
    # sol_group = None
    # if make_ndx_result.returncode == 0:
    #     for line in make_ndx_result.stdout.splitlines():
    #         parts = line.split()
    #         if parts and ("SOL" in parts[-1] or "Water" in parts[-1] or "WAT" in parts[-1]):
    #             sol_group = parts[0]
    #             break
    # if sol_group is None:
    #     print("❌ Could not determine SOL group from make_ndx.")
    #     print("📋 Here's what make_ndx returned:")
    #     print(make_ndx_result.stdout)
    #     sol_group = input("🔢 Please enter the group number corresponding to SOL manually: ").strip()

    # Use the detected group + optional custom index
    # if not run_command(
    #     f"gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral {ndx_arg}",
    #     cwd=directory,
    #     input_text=f"{sol_group}\n"
    # ):
    #     return False
    salt_m = os.environ.get("GMX_SALT_M", "0.15")  # default 0.15 M, override with env var
    if not run_command(
        f"gmx genion -s ions.tpr -o solv_ions.gro -p topol.top "
        f"-pname NA -nname CL -neutral -conc {salt_m} -seed 2025 {ndx_arg}",
        cwd=directory,
        input_text=f"{sol_group}\n"
    ):
        return False
    # salt_m = args.salt
    # if not run_command(
    #     f"gmx genion -s ions.tpr -o solv_ions.gro -p topol.top "
    #     f"-pname NA -nname CL -neutral -conc {salt_m} -seed 2025 {ndx_arg}",
    #     cwd=directory,
    #     input_text=f\"{sol_group}\\n\"
    # ):
    #     return False



    # if not run_command("gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral", cwd=directory, input_text=f"{sol_group}\n"):
    #     return False

    if not run_command("gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 2", cwd=directory):
        return False

    print("\n✅ All setup steps completed successfully.")
    return True




if __name__ == "__main__":
    if len(os.sys.argv) == 1:
        print("ℹ️ Running interactively — you can also provide:")
        print("   --chain-names 'Rap1B,Rap1GAP' or --chain-map 'A:Rap1B,B:Rap1GAP'\n")

    base_directory = os.getcwd()
    print("\n📂 Base Directory:", base_directory)
    if args.pdb:
        pdb_filename = args.pdb
        print(f"✅ Using specified PDB file: {pdb_filename}")
    else:
        pdb_filename = select_pdb_file(base_directory)
    # Ask user whether a ligand is present in the system
    
    # --- AUTO-DETECT LIGANDS FROM THE PDB ---
    detected_ligs = autodetect_ligands(os.path.join(base_directory, pdb_filename))

    ligand_code = None  # default

    if args.ligand:
        ligand_code = args.ligand.upper()
        print(f"\n✅ Ligand mode detected from CLI: {ligand_code}")

    elif args.mode == "protein":
        ligand_code = None
        print("\n✅ Protein-only mode forced by CLI. Skipping ligand detection.")

    else:
        # If script is interactive OR mode not pre-selected
        if detected_ligs:
            print("\n🔍 Auto-detected ligand(s) in the PDB:")
            for lig in detected_ligs:
                print(f"   • {lig}")

            if len(detected_ligs) == 1:
                choice = input(f"\nUse detected ligand '{detected_ligs[0]}'? [Y/n]: ").strip().lower()
                if choice in ("", "y", "yes"):
                    ligand_code = detected_ligs[0]
                else:
                    ligand_code = input("Enter 3-letter ligand code manually: ").strip().upper()
            else:
                print("\nMultiple ligands detected.")
                print("Enter one of these codes, or type a different ligand manually:")
                ligand_code = input(f"Choose ligand [{', '.join(detected_ligs)}]: ").strip().upper()

                if ligand_code == "":
                    ligand_code = detected_ligs[0]  # default first
        else:
            print("\nℹ️ No ligand detected in the PDB.")
            use_lig = input("Run as protein-only? [Y/n]: ").strip().lower()

            if use_lig in ("n", "no"):
                ligand_code = input("Enter 3-letter ligand code manually: ").strip().upper()
            else:
                ligand_code = None

    # Final report
    if ligand_code:
        print(f"\n✅ Using ligand: {ligand_code}")
    else:
        print("\n✅ No ligand selected for this simulation.")



    # Display chosen options
    print(f"\n✅ Selected PDB file: {pdb_filename}")
    if ligand_code:
        print(f"✅ Using ligand code: {ligand_code}")
    else:
        print("✅ No ligand selected for this simulation.")
    if not process_directory(base_directory, pdb_filename, ligand_code):
        print("⚠️ Error detected. Stopping execution.")
        exit(1)
