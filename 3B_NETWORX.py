#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
NETWORX v3 — Journal-Grade Ligand–Protein Interaction Diagram
--------------------------------------------------------------
Implements:
 • Two-panel figure: (a) stacked bar summary, (b) oval network
 • Oval layout for residues
 • Curved edges from ligand atoms to residues
 • Edge labels at 75% distance (fractions)
 • Bond order rendering (1/2/3 parallel lines)
 • Full halogen / sulfur / phosphorus coloring
 • Hydrogen bond edges: thick dashed maroon
 • Hydrophobic edges: soft dashed light-green
 • Water bridges: dashed blue
 • Ionic edges omitted in network
 • Final figure saved at high DPI

Assumptions:
 • Run from MD system root (where Analysis_Results/ lives)
 • Required files:
     Analysis_Results/InteractionTypes_Summary.csv
     Analysis_Results/FilteredContacts_Framewise.csv
 • Ligand MOL2:
     <ligand>.cgenff.mol2
"""

import os
import math
import argparse
import matplotlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from matplotlib.patches import Patch, FancyArrowPatch
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgba

from rdkit import Chem
from rdkit.Chem import rdDepictor
import re
import math

HALOGEN_ELEMENTS = {"F", "Cl", "Br", "I"}

def is_halogen(atom_name, atom_elements):
    """
    Return True if atom_name corresponds to a halogen (F, Cl, Br, I),
    allowing suffix digits (e.g. Br12, Cl3).
    """
    if atom_name not in atom_elements:
        return False

    elem = atom_elements[atom_name]  # e.g. "Br", "Cl"
    return elem in HALOGEN_ELEMENTS


def find_nearest_halogen(atom_positions, halogen_atoms, center):
    """
    If multiple halogens exist, pick the one *closest to center*.
    Otherwise return the single halogen available.
    """
    if not halogen_atoms:
        return None

    if len(halogen_atoms) == 1:
        return halogen_atoms[0]

    cx, cy = center
    best = None
    best_d = 9999

    for h in halogen_atoms:
        x, y = atom_positions[h]
        d = math.hypot(x - cx, y - cy)
        if d < best_d:
            best = h
            best_d = d

    return best




def auto_place_residues_guided(residues, residue_to_ligatom, atom_positions, center):
    """
    Place residues on an oval boundary, guided by the ligand atom they interact with.
    Ensures:
      - Residues sit closer to the ligand region they interact with
      - No residue overlaps another
      - No residue is allowed inside the molecular core
    """

    cx, cy = center

    # -----------------------------------------------------
    # 1. Determine ellipse radii based on ligand dimensions
    # -----------------------------------------------------
    xs = [p[0] for p in atom_positions.values()]
    ys = [p[1] for p in atom_positions.values()]

    ligand_w = (max(xs) - min(xs))
    ligand_h = (max(ys) - min(ys))

    # Base radii guided by ligand dimensions
    rx = ligand_w * 1.6   # horizontal radius
    ry = ligand_h * 1.8   # vertical radius

    # Enforce sensible minimums so very compact ligand layouts don't
    # place residues on top of the ligand. These constants keep behaviour
    # local to the network builders and won't affect other graphs.
    rx = max(rx, ELLIPSE_RX * 0.9)
    ry = max(ry, ELLIPSE_RY * 0.9)

    # -----------------------------------------------------
    # 2. Compute target angle for each residue
    # -----------------------------------------------------
    angles = {}

    for res in residues:

        lig_atoms = residue_to_ligatom.get(res, [])
        if isinstance(lig_atoms, str):
            lig_atoms = [lig_atoms]

        # Default angle: evenly distributed if no known atom
        if not lig_atoms:
            idx = residues.index(res)
            angles[res] = 2*np.pi * idx / len(residues)
            continue

        # Compute average angle from center to its interacting ligand atoms
        angs = []
        for la in lig_atoms:
            if la not in atom_positions:
                continue
            xL, yL = atom_positions[la]
            ang = np.arctan2(yL - cy, xL - cx)
            angs.append(ang)

        if angs:
            angles[res] = np.mean(angs)
        else:
            idx = residues.index(res)
            angles[res] = 2*np.pi * idx / len(residues)

    # -----------------------------------------------------
    # 3. Sort residues by angle for stable layout
    # -----------------------------------------------------
    sorted_res = sorted(residues, key=lambda r: angles[r])

    # -----------------------------------------------------
    # 4. Initial placement on ellipse
    # -----------------------------------------------------
    res_positions = {}

    for res in sorted_res:
        a = angles[res]
        x = cx + rx * np.cos(a)
        y = cy + ry * np.sin(a)
        res_positions[res] = [x, y]

    # -----------------------------------------------------
    # 5. Non-overlap correction by nudging angles + pairwise checks
    # -----------------------------------------------------
    min_dist = 2.0 * RESIDUE_RADIUS

    # Attempt multiple expand/relax cycles. In each cycle we relax angles
    # to reduce residue-residue overlaps, then check for any remaining
    # residue-residue or residue-ligand overlaps. If overlaps remain,
    # expand the ellipse radii slightly and retry.
    max_expand_passes = 10
    for expand_pass in range(max_expand_passes):
        # Relaxation: multiple short passes to nudge crowded residues apart
        for _pass in range(6):
            # check all pairs (not just adjacent) and nudge angles
            for i1 in range(len(sorted_res)):
                for i2 in range(i1 + 1, len(sorted_res)):
                    r1 = sorted_res[i1]
                    r2 = sorted_res[i2]

                    x1, y1 = res_positions[r1]
                    x2, y2 = res_positions[r2]

                    d = math.hypot(x1 - x2, y1 - y2)
                    if d < min_dist and d > 1e-6:
                        # amount to nudge scaled by proximity
                        overlap_frac = (min_dist - d) / min_dist
                        delta = 0.04 + 0.18 * overlap_frac

                        # separate angles: r1 clockwise, r2 counter-clockwise
                        angles[r1] -= delta
                        angles[r2] += delta

                        # recompute their positions on ellipse
                        a1 = angles[r1]
                        a2 = angles[r2]
                        res_positions[r1] = [cx + rx * np.cos(a1), cy + ry * np.sin(a1)]
                        res_positions[r2] = [cx + rx * np.cos(a2), cy + ry * np.sin(a2)]

        # After relaxation, check for any remaining overlaps
        any_overlap = False

        # 1) residue-residue pairwise
        for i1 in range(len(sorted_res)):
            for i2 in range(i1 + 1, len(sorted_res)):
                r1 = sorted_res[i1]
                r2 = sorted_res[i2]
                x1, y1 = res_positions[r1]
                x2, y2 = res_positions[r2]
                if math.hypot(x1 - x2, y1 - y2) < min_dist:
                    any_overlap = True
                    break
            if any_overlap:
                break

        # 2) residue-ligand overlap
        if not any_overlap:
            for (ax, ay) in atom_positions.values():
                for (rx0, ry0) in res_positions.values():
                    if math.hypot(ax - rx0, ay - ry0) < (RESIDUE_RADIUS + LIG_ATOM_RADIUS + 0.35):
                        any_overlap = True
                        break
                if any_overlap:
                    break

        if not any_overlap:
            return res_positions

        # gentle expansion and recompute positions
        rx *= 1.12
        ry *= 1.12
        for res in sorted_res:
            a = angles[res]
            res_positions[res] = [cx + rx * np.cos(a), cy + ry * np.sin(a)]

    # fallback: return best-effort placement
    return res_positions





def fix_edges_to_center_safe(residue_to_ligatom, atom_positions, atom_elements, center):
    """
    Only redirect ligand atoms that *remain unresolved* after
    hydrogen-parent mapping AND LP mapping.

    This prevents valid atoms like H5, H4, H6 from being replaced
    incorrectly by halogens. Also attempt similarity mapping for cases
    like 'N3' when possible.
    """

    def resolve_similar_atom(name, atom_positions, atom_elements, center=None):
        a = normalize_atom(name)
        if not a:
            return None
        if a in atom_positions:
            return a

        m = re.match(r"([A-Za-z]+)(\d*)", a)
        if not m:
            return None
        prefix, num = m.group(1), m.group(2)

        candidates = [nm for nm in atom_positions.keys()
                      if atom_elements.get(nm, "").upper() == prefix.upper() or nm.upper().startswith(prefix.upper())]
        if not candidates:
            candidates = [nm for nm in atom_positions.keys() if nm[0].upper() == prefix[0].upper()]

        if not candidates:
            return None
        if len(candidates) == 1:
            return candidates[0]

        if num:
            try:
                t = int(num)
                best, bestd = None, float("inf")
                for c in candidates:
                    dm = re.search(r"(\d+)$", c)
                    if dm:
                        v = int(dm.group(1))
                        d = abs(v - t)
                    else:
                        d = 10**6
                    if d < bestd:
                        best, bestd = c, d
                if best:
                    return best
            except Exception:
                pass

        if center is not None:
            cx, cy = center
            best, bestd = None, float("inf")
            for c in candidates:
                x, y = atom_positions[c]
                d = math.hypot(x - cx, y - cy)
                if d < bestd:
                    best, bestd = c, d
            return best

        return candidates[0]


    # Collect halogens
    halogens = [nm for nm, el in atom_elements.items() if el in {"F", "Cl", "Br", "I"}]

    if not halogens:
        print("⚠ No halogens available for fallback — attempting similarity mapping for unresolved atoms.")
        updated = residue_to_ligatom.copy()
        for res, ligatom in residue_to_ligatom.items():
            if not ligatom:
                continue
            ligatom_norm = normalize_atom(ligatom)
            if ligatom_norm in atom_positions:
                continue
            if ligatom_norm.startswith("H"):
                print(f"⚠ Unresolved hydrogen '{ligatom}' for {res} — keeping unresolved for now.")
                continue
            cand = resolve_similar_atom(ligatom_norm, atom_positions, atom_elements, center=center)
            if cand and cand in atom_positions:
                print(f"🔧 Mapping unresolved '{ligatom}' for {res} → similar atom '{cand}'")
                updated[res] = cand
            else:
                print(f"⚠ Could not map unresolved '{ligatom}' for {res}")
        return updated

    # fallback halogen path unchanged
    cx, cy = center
    halogen = min(
        halogens,
        key=lambda h: math.hypot(atom_positions[h][0] - cx, atom_positions[h][1] - cy)
    )

    updated = residue_to_ligatom.copy()
    for res, ligatom in residue_to_ligatom.items():
        if not ligatom:
            continue
        ligatom_norm = normalize_atom(ligatom)
        if ligatom_norm in atom_positions:
            continue
        if ligatom_norm.startswith("H"):
            print(f"⚠ Unresolved hydrogen '{ligatom}' for {res} — keeping center fallback.")
            continue
        print(f"🔧 Redirecting *UNRESOLVED* ligand atom '{ligatom}' for {res} → halogen {halogen}")
        updated[res] = halogen

    return updated


# --- small utility needed by several loaders (must be defined early) ---
def normalize_atom(a):
    """Clean ligand atom names (remove '.', spaces, etc.). Returns empty string for non-strings."""
    if not isinstance(a, str):
        return ""
    return "".join(c for c in a.strip() if c.isalnum())

def residue_number(resname: str) -> int:
    """Extract integer residue number from a label like 'LYS73' or 'ALA 123B'."""
    s = str(resname)
    digits = "".join(ch for ch in s if ch.isdigit())
    try:
        return int(digits) if digits else 0
    except Exception:
        return 0

def _resolve_h_to_heavy(atom, parent_map, atom_positions, atom_elements, center, max_heavy_dist=2.5):
    """
    Map a possibly-hydrogen ligand atom name to the corresponding heavy atom.
    Strategy (in order):
      1) If atom is already a heavy atom in atom_positions -> return it
      2) If parent_map has direct mapping (H -> heavy) -> use it
      3) If atom starts with 'H', try prefix matching against parent_map keys
      4) If hydrogen has 2D coords, pick nearest heavy within max_heavy_dist
      5) Fallback: return heavy atom nearest to ligand center
    Returns heavy_atom_name or None.
    """
    a = normalize_atom(atom)
    if not a:
        return None

    # 1) already heavy
    if a in atom_positions and not a.startswith("H"):
        return a

    # 2) direct parent_map mapping
    if a in parent_map:
        parent = parent_map[a]
        if parent in atom_positions:
            return parent

    # 3) prefix match for hydrogens (H3 -> H31/H32)
    if a.startswith("H"):
        root = a.rstrip("0123456789")
        for h_mol2, heavy in parent_map.items():
            if h_mol2.startswith(root) and heavy in atom_positions:
                return heavy

    # 4) if H coord exists, pick nearest heavy within cutoff
    if a in atom_positions:
        hx, hy = atom_positions[a]
        best = None
        best_d = float("inf")
        for nm, (x, y) in atom_positions.items():
            # skip hydrogens
            if atom_elements.get(nm, "").upper() == "H" or nm.startswith("H"):
                continue
            d = math.hypot(x - hx, y - hy)
            if d < best_d:
                best = nm
                best_d = d
        if best and best_d <= max_heavy_dist:
            return best

    # 5) fallback: pick heavy atom nearest ligand center
    cx, cy = center
    best = None
    best_d = float("inf")
    for nm, (x, y) in atom_positions.items():
        if atom_elements.get(nm, "").upper() == "H" or nm.startswith("H"):
            continue
        d = math.hypot(x - cx, y - cy)
        if d < best_d:
            best = nm
            best_d = d
    return best
# --------------------------------------------------------------------
# Global paths
# --------------------------------------------------------------------
ANALYSIS_DIR = "Analysis_Results"

# --------------------------------------------------------------------
# Color definitions
# --------------------------------------------------------------------

# Atom colors (Schrödinger-like)
ATOM_COLORS = {
    "C":  "#555555",   # dark gray
    "H":  "#D0D0D0",   # light gray
    "O":  "#E41A1C",   # red
    "N":  "#377EB8",   # blue
    "S":  "#FFD700",   # bright yellow
    "P":  "#FF69B4",   # hot pink

    "F":  "#7CFC00",   # light green
    "Cl": "#32CD32",   # lime green
    "Br": "#228B22",   # forest green
    "I":  "#800080",   # iodine purple
}

# Network edge colors
EDGE_HBOND = "#800000"    # maroon, dashed, thick
EDGE_HYDRO = "#90EE90"    # light green, dashed
EDGE_WATER = "#1E90FF"    # blue, dashed
    # FULL CONTACTS (for Panel C) — USE THE CORRECT FILE
# Bar plot colors
BAR_HBOND  = "#D62728"
BAR_HYDRO  = "#2CA02C"
BAR_IONIC  = "#FF7F0E"
BAR_WATER  = "#1F77B4"

# Radii for layout
RESIDUE_RADIUS = 0.55    # radius of residue node circle
LIG_ATOM_RADIUS = 0.18   # radius of ligand atom circle (increased for visual clearance)

# Oval parameters (Option B)
# Collapsed network only
ELLIPSE_RX = 5.6
ELLIPSE_RY = 3.4

# Ligand 2D layout scaling factor (increase to spread ligand atom nodes)
LIGAND_LAYOUT_SCALE = 1.25

RESIDUE_CATEGORIES = {
    "hydrophobic": {"ALA","VAL","LEU","ILE","MET","PHE","TRP","TYR","PRO"},
    "positive":    {"LYS","ARG","HIS"},       # HIS treated as positive
    "negative":    {"ASP","GLU"},
    "polar":       {"SER","THR","ASN","GLN","CYS"},
}

RESIDUE_COLORS = {
    "hydrophobic": "#C9F7C5",  # soft green
    "positive":    "#C5F1FF",  # soft cyan
    "negative":    "#FFCBC5",  # soft salmon
    "polar":       "#E8E8E8",  # soft gray
}

# --------------------------------------------------------------------
# Helper functions
# --------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Generate journal-ready PROTAC interaction diagrams",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    p.add_argument(
        "-l", "--ligand",
        type=str,
        help="Ligand residue name (e.g., A1D)"
    )

    p.add_argument(
        "-f", "--filter",
        dest="minfrac",
        type=float,
        default=None,
        help="Minimum contact fraction cutoff"
    )

    p.add_argument(
        "-o", "--output",
        type=str,
        help="Output base filename (without extension)"
    )


    # NEW — allow node-only panel
    p.add_argument(
        "--no-edges",
        action="store_true",
        help="Render network diagram without edges (node-only mode)."
    )


    return p.parse_args()


def ask_if_missing(args):
    print("\n🧬 NETWORX — Two-Panel Journal Mode\n")

    if not args.ligand:
        args.ligand = input("Ligand resname (e.g., A1D): ").strip()

    if args.minfrac is None:
        v = input("Minimum contact fraction [default 0.10]: ").strip()
        if v:
            v = v.replace(",", ".")
            args.minfrac = float(v)
        else:
            args.minfrac = 0.10

    # Keep output unset so we do not create a per-run subdirectory.
    if not args.output:
        args.output = None

    return args


# --------------------------------------------------------------------
# Ligand layout via MOL2 (2D drawing)
# --------------------------------------------------------------------
def get_ligand_layout(ligand_resname: str):
    """
    Read <ligand_resname>.cgenff.mol2 and extract a 2D layout WITHOUT
    modifying or renaming atoms. This preserves CGenFF atom names (H1, H2, H3…),
    ensures correct hydrogen → heavy mapping, and avoids RDKit sanitization
    that breaks name consistency.

    Returns:
        atom_positions : {atom_name: (x, y)}
        atom_elements  : {atom_name: element_symbol}
        bonds          : [(atom1, atom2, order)]
        center         : (cx, cy)
        parent_map     : {Hname: heavy_atom_name}
    """

    mol2_file = f"{ligand_resname}.cgenff.mol2"
    if not os.path.exists(mol2_file):
        raise FileNotFoundError(f"MOL2 file not found: {mol2_file}")

    # ---------------------------------------------------------------
    # 1) Load MOL2 EXACTLY as-is.
    # ---------------------------------------------------------------
    # sanitize=False prevents RDKit from renaming atoms or changing valence
    mol = Chem.MolFromMol2File(mol2_file, sanitize=False, removeHs=False)
    mol_arom = restore_aromaticity(mol)

    if mol is None:
        raise ValueError(f"RDKit failed to load {mol2_file}")

    # ---------------------------------------------------------------
    # 2) Fix any missing atomic numbers (rare but possible)
    # ---------------------------------------------------------------
    pt = Chem.GetPeriodicTable()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            sym = atom.GetSymbol()
            atom.SetAtomicNum(pt.GetAtomicNumber(sym))

    # ---------------------------------------------------------------
    # 3) DO NOT SANITIZE — THIS BREAKS ATOM NAMES
    # ---------------------------------------------------------------
    # No SanitizeMol, No SetAromaticity, No Kekulize
    # These cause H13/H14 renaming, destroy mapping.

    # ---------------------------------------------------------------
    # 4) Generate 2D coordinates (safe & does NOT rename anything)
    # ---------------------------------------------------------------
    rdDepictor.Compute2DCoords(mol)
    conf = mol.GetConformer()

    # ---------------------------------------------------------------
    # 5) Extract atom positions exactly as named in MOL2
    # ---------------------------------------------------------------
    atom_positions = {}
    atom_elements = {}
    atom_names = []
    parent_map = {}

    for i, atom in enumerate(mol.GetAtoms()):
        pos = conf.GetAtomPosition(i)

        # The REAL, ORIGINAL atom name from the MOL2 file
        info = atom.GetMonomerInfo()
        raw_name = info.GetName().strip() if info else f"{atom.GetSymbol()}{i+1}"
        atom_name = normalize_atom(raw_name)

        atom_positions[atom_name] = (pos.x, pos.y)
        atom_elements[atom_name] = atom.GetSymbol()
        atom_names.append(atom_name)

    # ---------------------------------------------------------------
    # 6) Bonds + parent map for hydrogens
    # ---------------------------------------------------------------
    bonds = []

    for b in mol.GetBonds():
        i = b.GetBeginAtomIdx()
        j = b.GetEndAtomIdx()

        a1 = atom_names[i]
        a2 = atom_names[j]

        bt = b.GetBondType()
        if bt == Chem.BondType.SINGLE:
            order = 1
        elif bt == Chem.BondType.DOUBLE:
            order = 2
        elif bt == Chem.BondType.TRIPLE:
            order = 3
        elif bt == Chem.BondType.AROMATIC:
            # Alternate single/double bond pattern for drawing ONLY.
            # Does NOT modify molecule, atom indices, or mapping.
            if (i + j) % 2 == 0:
                order = 2
            else:
                order = 1


        bonds.append((a1, a2, order))

        # Build hydrogen → heavy-atom parent map
        sym1 = mol.GetAtomWithIdx(i).GetSymbol()
        sym2 = mol.GetAtomWithIdx(j).GetSymbol()

        if sym1 == "H" and sym2 != "H":
            parent_map[a1] = a2
        elif sym2 == "H" and sym1 != "H":
            parent_map[a2] = a1
    

    

    # ---------------------------------------------------------------
    # 7) Ligand center (for fallback)
    # ---------------------------------------------------------------
    xs = [p[0] for p in atom_positions.values()]
    ys = [p[1] for p in atom_positions.values()]
    center = (float(np.mean(xs)), float(np.mean(ys)))

    # ---------------------------------------------------------------
    # 8) Optionally scale ligand atom layout relative to center to spread nodes
    # ---------------------------------------------------------------
    if LIGAND_LAYOUT_SCALE != 1.0:
        cx, cy = center
        for nm, (x, y) in list(atom_positions.items()):
            dx = x - cx
            dy = y - cy
            atom_positions[nm] = (cx + dx * LIGAND_LAYOUT_SCALE, cy + dy * LIGAND_LAYOUT_SCALE)

    # ---------------------------------------------------------------
    # Debug
    # ---------------------------------------------------------------
    print("\n🔬 Ligand atom names from MOL2 (PRESERVED):")
    print(sorted(atom_positions.keys()))

    return atom_positions, atom_elements, bonds, center, parent_map, mol_arom, atom_names


def resolve_lp_mapping(ligand, residue_to_ligatom):
    """
    Replace LP atoms (LP1, LP2...) with nearest heavy atom based on ligand_init.pdb geometry.

    Requirements:
        ligand_init.pdb exists in working directory (e.g. A1D_ini.pdb)
    """
    import math

    pdb_file = f"{ligand}_ini.pdb"
    if not os.path.exists(pdb_file):
        print(f"[LP-MAP] No {pdb_file} found — skipping LP mapping.")
        return residue_to_ligatom

    # Load coordinates from PDB
    atom_coords = {}  # atom_name → (x,y,z)
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                nm = line[12:16].strip()   # atom name
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atom_coords[nm] = (x, y, z)

    # Identify heavy atoms (exclude hydrogens + LP atoms)
    heavy_atoms = [nm for nm in atom_coords if not nm.startswith("H") and not nm.startswith("LP")]

    def dist(a, b):
        return math.dist(atom_coords[a], atom_coords[b])

    updated_map = residue_to_ligatom.copy()

    for res, at in residue_to_ligatom.items():
        if not at or not at.startswith("LP"):
            continue

        if at not in atom_coords:
            continue  #

        # Find nearest heavy atom
        lp_coord = atom_coords[at]
        best = None
        best_d = 9999

        for hvy in heavy_atoms:
            d = dist(at, hvy)
            if d < best_d:
                best_d = d
                best = hvy

        if best:
            print(f"[LP-MAP] {res}: {at} → {best}  (distance {best_d:.2f} Å)")
            updated_map[res] = best

    return updated_map

def shorten_to_circle(x_center, y_center, x_out, y_out, r):
    """
    Shorten line from outer point (x_out, y_out) toward center (x_center, y_center)
    so that it ends at radius r from center.
    """
    dx = x_out - x_center
    dy = y_out - y_center
    d = math.hypot(dx, dy)
    if d == 0:
        return x_center, y_center
    scale = r / d
    return x_center + dx * scale, y_center + dy * scale


def _atom_node_color(elem: str):
    """Return the display color for a ligand atom element (node color).

    Carbon -> black, S -> yellow, P -> pink, halogens -> Miami green,
    O -> red, N -> blue, otherwise fall back to ATOM_COLORS mapping.
    """
    if not elem:
        return ATOM_COLORS.get("C", "#000000")
    if elem == "C":
        return "#000000"
    if elem == "S":
        return "#FFD700"
    if elem == "P":
        return "#FF69B4"
    if elem in ("F", "Cl", "Br", "I"):
        return "#007849"
    if elem == "O":
        return "#D72638"
    if elem == "N":
        return "#1E90FF"
    return ATOM_COLORS.get(elem, "#888")


def _draw_bicolor_arc(ax, x1, y1, x2, y2, rad, c1, c2, lw=2.0, linestyle='--', zorder=2, nseg=48):
    """Draw a curved edge from (x1,y1) to (x2,y2) with a color gradient from c1->c2.

    The curve is approximated as a quadratic Bezier using a single control
    point computed from `rad` (mimicking `arc3,rad=`). The curve is sampled
    and a LineCollection is used with per-segment colors interpolated between
    the two endpoint colors.
    """
    dx = x2 - x1
    dy = y2 - y1
    midx = 0.5 * (x1 + x2)
    midy = 0.5 * (y1 + y2)
    dist = math.hypot(dx, dy)
    if dist == 0:
        cp = (midx, midy)
    else:
        perp_x = -dy / dist
        perp_y = dx / dist
        cp = (midx + perp_x * rad * dist, midy + perp_y * rad * dist)

    ts = np.linspace(0.0, 1.0, nseg)
    pts = []
    for t in ts:
        x = (1 - t) ** 2 * x1 + 2 * (1 - t) * t * cp[0] + t ** 2 * x2
        y = (1 - t) ** 2 * y1 + 2 * (1 - t) * t * cp[1] + t ** 2 * y2
        pts.append((x, y))

    segments = [[pts[i], pts[i + 1]] for i in range(len(pts) - 1)]

    rgba1 = np.array(to_rgba(c1))
    rgba2 = np.array(to_rgba(c2))
    frac_ts = np.linspace(0.0, 1.0, len(segments))
    colors = [tuple(rgba1 * (1 - t) + rgba2 * t) for t in frac_ts]

    lc = LineCollection(segments, colors=colors, linewidths=lw, linestyles=linestyle, zorder=zorder, capstyle='round')
    ax.add_collection(lc)


def label_point_frac(x1, y1, x2, y2, f=0.75):
    """Return point at fraction f along line from (x1,y1) → (x2,y2)."""
    return (x1 + f * (x2 - x1), y1 + f * (y2 - y1))


def classify_residue(res_name: str):
    # Extract 3-letter code
    code = "".join([c for c in res_name if c.isalpha()]).upper()

    for cat, members in RESIDUE_CATEGORIES.items():
        if code in members:
            return cat
    return "polar"   # fallback


def auto_place_residues_no_overlap(residue_list, atom_positions,
                                   start_rx=4.2, start_ry=2.3,
                                   scale_factor=1.12):
    """
    Expands ellipse radii until no residue node circle overlaps ligand atoms.
    Returns:
        res_positions, final_rx, final_ry
    """
    rx, ry = start_rx, start_ry

    while True:
        # First: try layout
        n = len(residue_list)
        angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
        res_positions = {}
        for res, t in zip(residue_list, angles):
            x = rx * math.cos(t)
            y = ry * math.sin(t)
            res_positions[res] = (x, y)

        # Check for overlaps
        overlap = False
        for (ax, ay) in atom_positions.values():
            for (rx0, ry0) in res_positions.values():
                dist = math.hypot(ax - rx0, ay - ry0)
                if dist < (RESIDUE_RADIUS + LIG_ATOM_RADIUS + 0.2):
                    overlap = True
                    break
            if overlap:
                break

        if not overlap:
            return res_positions, rx, ry

        # Expand ellipse
        rx *= scale_factor
        ry *= scale_factor



def extract_mol2_bond_orders(mol2_file):
    """
    Parse the MOL2 file directly and extract:
        atom_name → atom_name → bond_order (1,2,3)
    Aromatic bonds ('ar') are treated as alternating 1/2 bond pattern.
    """

    bonds = []
    in_bond_section = False

    with open(mol2_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("@<TRIPOS>BOND"):
                in_bond_section = True
                continue
            if line.startswith("@<TRIPOS>") and in_bond_section:
                break

            if in_bond_section and line:
                parts = line.split()
                if len(parts) < 4:
                    continue

                _, a1, a2, order_raw = parts[:4]

                # convert bond order
                if order_raw == "1":
                    order = 1
                elif order_raw == "2":
                    order = 2
                elif order_raw == "3":
                    order = 3
                elif order_raw in ("ar", "4"):
                    # assign aromatic alternating pattern later
                    order = "AROMATIC"
                else:
                    order = 1

                bonds.append((a1, a2, order))

    return bonds



def make_residue_legend(ax):
    """
    Unified legend for Panel C — identical style to Panel B.
    Includes:
      • Edge types (H-bond, Hydrophobic, Water bridge, Ionic)
      • Residue category colors (Hydrophobic, Positive, Negative, Polar)
    """

    edge_handles = [
        Line2D([0],[0],color=EDGE_HBOND, linestyle="--", linewidth=2.8, label="H-bond"),
        Line2D([0],[0],color=EDGE_HYDRO, linestyle="--", linewidth=1.6, label="Hydrophobic"),
        Line2D([0],[0],color=EDGE_WATER, linestyle="--", linewidth=2.0, label="Water bridge"),
        Line2D([0],[0],color="#C39BD3", linestyle="--", linewidth=2.2, label="Ionic"),
    ]

    residue_handles = [
        Patch(facecolor=RESIDUE_COLORS["hydrophobic"], edgecolor="black", label="Hydrophobic"),
        Patch(facecolor=RESIDUE_COLORS["positive"],   edgecolor="black", label="Positive"),
        Patch(facecolor=RESIDUE_COLORS["negative"],   edgecolor="black", label="Negative"),
        Patch(facecolor=RESIDUE_COLORS["polar"],      edgecolor="black", label="Polar"),
    ]

    handles = edge_handles + residue_handles

    ax.legend(
        handles=handles,
        fontsize=9,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.02),
        ncol=4
    )



def find_active_residues(df_full, minfrac):
    """
    Return a sorted list of residues that have ≥1 interaction
    with Fraction >= minfrac.
    """
    active = (
        df_full[df_full["Fraction"] >= minfrac]["Residue"]
        .unique()
        .tolist()
    )
    active = sorted(active, key=residue_number)
    return active




def safe_kekulize_and_get_rings(mol, atom_names):
    """
    Try to kekulize the molecule and return:
        - valid_rings: list of rings (list of atom names)
        - bad_atoms: atoms RDKit said were not kekulizable

    Only rings fully composed of OK atoms will be returned.
    """
    m2 = Chem.Mol(mol)  # copy
    bad_atoms = set()

    try:
        Chem.SanitizeMol(
            m2,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_KEKULIZE |
                        Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
        )
        print("✓ Kekulization succeeded — ring detection active.")
    except Exception as e:
        msg = str(e)
        print("⚠ Kekulization failed:", msg)

        # Extract atoms RDKit complains about
        # Format usually: "Can't kekulize mol. Unkekulized atoms: 6 7 8 9 10"
        m = re.search(r"atoms:\s*([\d\s]+)", msg)
        if m:
            idxs = m.group(1).split()
            bad_atoms = {atom_names[int(i)].strip() for i in idxs}
            print("⚠ Problematic atoms:", bad_atoms)

        # Continue but with no valid rings
        return [], bad_atoms

    # If kekulized: get rings
    ri = m2.GetRingInfo().AtomRings()

    valid_rings = []
    for ring in ri:
        # Convert indices → atom names
        ring_atoms = [atom_names[i] for i in ring]

        # Exclude rings containing bad atoms
        if any(a in bad_atoms for a in ring_atoms):
            continue

        valid_rings.append(ring_atoms)

    return valid_rings, bad_atoms



def fix_kekulization_errors(mol, bonds, atom_names):
    """
    Post-processing correction that:

    • Attempts a safe kekulization on a COPY of the mol
    • Restores proper bond orders (1/2/3) for aromatic / conjugated rings
    • Leaves atom names, indices, coordinates untouched
    • Updates the existing 'bonds' list IN PLACE

    Returns: corrected_bonds (new list of (atom1, atom2, order))
    """

    # Always work on a copy — do NOT touch original
    m2 = Chem.Mol(mol)

    # --- Attempt sanitization to assign proper bond orders ---
    try:
        Chem.SanitizeMol(
            m2,
            sanitizeOps = Chem.SanitizeFlags.SANITIZE_KEKULIZE |
                          Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
        )
        print("✓ Kekulization successful.")
    except Exception as e:
        print("⚠ Kekulization failed softly:", e)
        # Fallback: only aromaticity perception
        try:
            Chem.SanitizeMol(m2, sanitizeOps=Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
            print("✓ Aromaticity restored without kekulization.")
        except Exception as e2:
            print("⚠ Aromaticity restoration failed:", e2)
            return bonds   # keep original

    # --- Extract corrected orders ---
    corrected = []
    for b in m2.GetBonds():
        i = b.GetBeginAtomIdx()
        j = b.GetEndAtomIdx()

        a1 = atom_names[i]
        a2 = atom_names[j]

        bt = b.GetBondType()
        if bt == Chem.BondType.SINGLE:
            order = 1
        elif bt == Chem.BondType.DOUBLE:
            order = 2
        elif bt == Chem.BondType.TRIPLE:
            order = 3
        else:
            # AROMATIC now has explicit Kekulé form → either 1 or 2
            order = 1

        corrected.append((a1, a2, order))

    print("✓ Bond orders updated using kekulized molecule.")
    return corrected





def draw_curved_edge(ax, x1, y1, x2, y2,
                     color="#000", lw=1.8, style="solid",
                     curvature=0.25, zorder=1):
    """
    Draws a quadratic Bezier curved edge between two points.
    Identical behavior to Panel B.
    """

    # Midpoint
    mx, my = (x1 + x2) / 2.0, (y1 + y2) / 2.0

    # Perpendicular offset for curvature
    dx, dy = x2 - x1, y2 - y1
    nx, ny = -dy, dx
    norm = math.sqrt(nx*nx + ny*ny) or 1.0
    nx /= norm; ny /= norm

    cx = mx + curvature * nx
    cy = my + curvature * ny

    path = matplotlib.path.Path([
        (x1, y1),
        (cx, cy),
        (x2, y2)
    ], [matplotlib.path.Path.MOVETO,
        matplotlib.path.Path.CURVE3,
        matplotlib.path.Path.CURVE3])

    patch = matplotlib.patches.PathPatch(
        path,
        edgecolor=color,
        linewidth=lw,
        linestyle=style,
        fill=False,
        zorder=zorder
    )

    ax.add_patch(patch)



from rdkit.Chem.rdchem import BondType





def restore_aromaticity(mol):
    """
    Safely restore aromaticity perception *without* touching atom names,
    coordinates, or bond orders. Sanitization steps that would modify
    the molecule (e.g., kekulization, valence checking) are skipped.
    """

    params = Chem.SanitizeFlags.SANITIZE_SETAROMATICITY | \
             Chem.SanitizeFlags.SANITIZE_KEKULIZE

    m2 = Chem.Mol(mol)  # shallow copy

    try:
        Chem.SanitizeMol(m2, sanitizeOps=params)
    except:
        # KEKULIZE fails for some CGenFF heteroaromatics — run only aromaticity
        try:
            Chem.SanitizeMol(
                m2,
                sanitizeOps=Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
            )
        except Exception as e:
            print("⚠ Aromaticity restoration failed:", e)

    return m2



def parse_mol2_bonds(mol2_path):
    """
    Reads @<TRIPOS>BOND entries from a MOL2 file and returns:
        [(atom_name1, atom_name2, bond_order), ...]
    where bond_order is integer 1/2/3 or "AR" for aromatic.
    """

    bonds = []
    atoms = {}
    reading_atoms = False
    reading_bonds = False

    with open(mol2_path, "r") as f:
        for line in f:
            line = line.strip()

            if line.startswith("@<TRIPOS>ATOM"):
                reading_atoms = True
                reading_bonds = False
                continue

            if line.startswith("@<TRIPOS>BOND"):
                reading_atoms = False
                reading_bonds = True
                continue

            if line.startswith("@<TRIPOS>"):
                reading_atoms = False
                reading_bonds = False
                continue

            # Parse atom block
            if reading_atoms and line:
                parts = line.split()
                # MOL2 atoms: index name x y z type subst subst_id charge
                idx = int(parts[0])
                name = parts[1]
                atoms[idx] = name

            # Parse bond block
            if reading_bonds and line:
                parts = line.split()
                if len(parts) < 4:
                    continue

                _, a1_idx, a2_idx, order_raw = parts[:4]
                a1_idx = int(a1_idx)
                a2_idx = int(a2_idx)

                atom1 = atoms[a1_idx]
                atom2 = atoms[a2_idx]

                # Assign bond order
                if order_raw == "1":
                    order = 1
                elif order_raw == "2":
                    order = 2
                elif order_raw == "3":
                    order = 3
                elif order_raw in ("ar", "4"):
                    order = "AR"
                else:
                    order = 1

                bonds.append((atom1, atom2, order))

    return bonds




def normalize_bond_endpoint(a, atom_names):
    """
    Convert RDKit index -> atom name, or return atom name directly.
    """
    # Case 1: integer index
    if isinstance(a, int):
        return atom_names[a]

    # Case 2: string atom name (MOL2)
    if isinstance(a, str):
        return a

    print(f"[WARN] Unexpected bond endpoint type: {a}")
    return None




def find_aromatic_rings(mol):
    """
    Return rings where all atoms are aromatic.
    """
    aromatic_rings = []
    ri = mol.GetRingInfo().AtomRings()

    for ring in ri:
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            aromatic_rings.append(ring)

    return aromatic_rings



def ring_center(positions, atom_indices):
    xs = [positions[i][0] for i in atom_indices]
    ys = [positions[i][1] for i in atom_indices]
    return (np.mean(xs), np.mean(ys))


def draw_bond(ax, x1, y1, x2, y2, order, elem1="C", elem2="C", a1_name=None, a2_name=None):
    """
    Draw a single, double, or triple bond with correct parallel offsets.
    Quiet, robust: handle None names gracefully and only draw valid geometry.
    """
    dx = x2 - x1
    dy = y2 - y1
    L = math.hypot(dx, dy)
    if L == 0:
        return

    # Unit / perpendicular
    ux = dx / L
    uy = dy / L
    px = -uy
    py = ux

    # Default styling
    base_offset = 0.18
    # Draw ligand bonds below ligand atoms: use low bond_zorder
    bond_zorder = 1.0
    lw = 3.6

    # determine endpoint colors from element types
    c1 = _atom_node_color(elem1)
    c2 = _atom_node_color(elem2)

    def plot_split_line(xa, ya, xb, yb, col_a, col_b, linewidth=lw, z=bond_zorder):
        # split the line exactly at midpoint and draw two colored halves
        mx = 0.5 * (xa + xb)
        my = 0.5 * (ya + yb)
        ax.plot([xa, mx], [ya, my], color=col_a, linewidth=linewidth, zorder=z)
        ax.plot([mx, xb], [my, yb], color=col_b, linewidth=linewidth, zorder=z)

    # For triple bonds: increase spacing and make outer lines more prominent
    if order == 3:
        offset = 0.34
        outer_lw = 3.4
        mid_lw = 3.6
        # middle line
        plot_split_line(x1, y1, x2, y2, c1, c2, linewidth=mid_lw, z=bond_zorder)
        # outer parallel lines
        plot_split_line(x1 + px*offset, y1 + py*offset, x2 + px*offset, y2 + py*offset, c1, c2, linewidth=outer_lw, z=bond_zorder)
        plot_split_line(x1 - px*offset, y1 - py*offset, x2 - px*offset, y2 - py*offset, c1, c2, linewidth=outer_lw, z=bond_zorder)
        return

    # Double bond
    if order == 2:
        offset = base_offset
        plot_split_line(x1 + px*offset, y1 + py*offset, x2 + px*offset, y2 + py*offset, c1, c2, linewidth=lw, z=bond_zorder)
        plot_split_line(x1 - px*offset, y1 - py*offset, x2 - px*offset, y2 - py*offset, c1, c2, linewidth=lw, z=bond_zorder)
        return

    # Single bond
    if order == 1:
        plot_split_line(x1, y1, x2, y2, c1, c2, linewidth=lw, z=bond_zorder)
        return

    # fallback
    plot_split_line(x1, y1, x2, y2, c1, c2, linewidth=lw, z=bond_zorder)



# --------------------------------------------------------------------
# Panel A — stacked bar plot
# --------------------------------------------------------------------
def draw_contact_bar_panel(ax, df, total_frames, minfrac):
    """
    Draw stacked bar plot of interaction fractions (Panel a).
    Uses InteractionTypes_Summary.csv.
    """
    df_filt = df[df["Fraction"] >= minfrac].copy()
    if df_filt.empty:
        raise ValueError(f"No residues exceed minfrac={minfrac} for bar plot")

    df_filt = df_filt.sort_values("Residue", key=lambda s: s.map(residue_number))

    # Identify interaction columns present
    type_cols = []
    if "H-bond" in df_filt.columns:
        type_cols.append("H-bond")
    if "Hydrophobic" in df_filt.columns:
        type_cols.append("Hydrophobic")
    if "Ionic" in df_filt.columns:
        type_cols.append("Ionic")
    if "Water bridge" in df_filt.columns:
        type_cols.append("Water bridge")

    if not type_cols:
        raise ValueError("No known interaction type columns found for bar plot")

    frac_data = {}
    for col in type_cols:
        frac_data[col] = df_filt[col].fillna(0) / float(total_frames)

    residues = df_filt["Residue"].tolist()
    x = np.arange(len(residues))
    bottom = np.zeros(len(residues))

    color_map = {
        "H-bond":      (BAR_HBOND, "H-bonds"),
        "Hydrophobic": (BAR_HYDRO, "Hydrophobic"),
        "Ionic":       (BAR_IONIC, "Ionic"),
        "Water bridge": (BAR_WATER, "Water bridges"),
    }

    for col in type_cols:
        frac_vals = frac_data[col].values
        c, label = color_map.get(col, ("#999999", col))
        ax.bar(x, frac_vals, bottom=bottom,
               color=c, edgecolor="black", linewidth=0.6, label=label)
        bottom += frac_vals

    ax.set_xticks(x)
    ax.set_xticklabels(residues, rotation=60, ha="right", fontsize=7.5)
    ax.set_ylabel("Interaction fraction", fontsize=9)
    ax.set_title("Protein–Ligand Contacts", fontsize=12, weight="bold", pad=6)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Legend on the right, small
    handles = [
        Patch(facecolor=BAR_HBOND, edgecolor="black", label="H-bonds"),
        Patch(facecolor=BAR_HYDRO, edgecolor="black", label="Hydrophobic"),
        Patch(facecolor=BAR_IONIC, edgecolor="black", label="Ionic"),
        Patch(facecolor=BAR_WATER, edgecolor="black", label="Water bridges"),
    ]

    # Compact legend placed slightly left of center above the title so it doesn't overlap.
    # Move x a small amount (0.44) to shift the legend a few pixels left in axes coords.
    ax.legend(
        handles=handles,
        fontsize=7.5,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.84, 1.12),   # small left shift from center
        ncol=2,
        handletextpad=0.2,
        columnspacing=0.6,
        borderaxespad=0.2
    )

# --------------------------------------------------------------------
# Panel B — network diagram
# --------------------------------------------------------------------
def draw_network_panel(ax, ligand, df_filt, residue_to_ligatom):
    """
    Schrödinger-style ligand–protein interaction network.
    Multi-edge + hydrogen-safe + LP-mapping + aromatic rings preserved.

    Features:
      ✔ Multi-atom edges per residue (NO COLLAPSE)
      ✔ Hydrogen → heavy-atom remapping (prefix logic)
      ✔ LP1/LP2 → heavy-atom fallback + halogen fallback
      ✔ No center-collapsed edges unless truly unresolved
      ✔ Multiple arcs per residue with alternating curvature
      ✔ Correct % labels for EACH interaction
    """

    # ---------------------------------------------------------
    # 1. Load ligand layout (atom positions, bonds, parents)
    # ---------------------------------------------------------
    atom_positions, atom_elements, bonds, center, parent_map, mol, atom_names = \
        get_ligand_layout(ligand)
    rings, bad_atoms = safe_kekulize_and_get_rings(mol, atom_names)

    print("Valid rings:", rings)
    print("Skipped (bad) atoms:", bad_atoms)

    # ---------------------------------------------------------
    # Optional safety: if ligand atom layout is crowded (many atoms
    # overlapping visually), gently scale ligand atom positions away
    # from the ligand center so atoms/bonds and residue edges don't
    # occlude each other. This preserves original MOL2 coords but
    # improves figure readability when RDKit produced compact layouts.
    # ---------------------------------------------------------
    try:
        # consider only non-hydrogen atoms for crowding metric
        non_h_coords = [pos for nm, pos in atom_positions.items() if not nm.startswith("H")]
        if len(non_h_coords) > 1:
            min_d = float("inf")
            for i in range(len(non_h_coords)):
                xi, yi = non_h_coords[i]
                for j in range(i+1, len(non_h_coords)):
                    xj, yj = non_h_coords[j]
                    d = math.hypot(xi - xj, yi - yj)
                    if d < min_d:
                        min_d = d
            # if minimum inter-atom distance is small, scale outward
            if min_d < 0.9:  # threshold (Angstrom-like units from RDKit 2D coords)
                cx, cy = center
                # compute a conservative scale: either module default or enough
                # to push min_d to ~1.1
                safe_scale = max(LIGAND_LAYOUT_SCALE, 1.1 / max(min_d, 1e-6))
                for nm, (x, y) in list(atom_positions.items()):
                    dx = x - cx
                    dy = y - cy
                    atom_positions[nm] = (cx + dx * safe_scale, cy + dy * safe_scale)
    except Exception:
        # non-fatal: if something goes wrong here, keep original positions
        pass


    # ---------------------------------------------------------
    # 2. Fix LP atoms → nearest heavy atoms
    # ---------------------------------------------------------
    residue_to_ligatom = resolve_lp_mapping(ligand, residue_to_ligatom)

    # ---------------------------------------------------------
    # 3. Fallback unresolved atoms → halogen closest to center
    # ---------------------------------------------------------
    residue_to_ligatom = fix_edges_to_center_safe(
        residue_to_ligatom, atom_positions, atom_elements, center
    )

    # Prepare canvas
    ax.set_aspect("equal")
    ax.axis("off")

    # Draw ligand bonds earlier/later no longer required — zorder handles stacking.
    # ---------------------------------------------------------
    # 4. Draw ligand bonds (1/2/3-line rendering)
    # ---------------------------------------------------------
    for a1_idx, a2_idx, order in bonds:

        a1 = normalize_bond_endpoint(a1_idx, atom_names)
        a2 = normalize_bond_endpoint(a2_idx, atom_names)

        if a1 is None or a2 is None:
            continue

        # Skip hydrogens
        if a1.startswith("H") or a2.startswith("H"):
            continue

        x1, y1 = atom_positions[a1]
        x2, y2 = atom_positions[a2]

        elem1 = atom_elements.get(a1, "C")
        elem2 = atom_elements.get(a2, "C")

        draw_bond(ax, x1, y1, x2, y2, order, elem1, elem2)



    # ---------------------------------------------------------
    # 5. Aromatic rings
    # ---------------------------------------------------------
    # try:
    #     aromatic_rings = find_aromatic_rings(mol)
    #     for ring in aromatic_rings:
    #         ring_names = [atom_names[i] for i in ring]
    #         coords = [atom_positions[n] for n in ring_names]

    #         cx = np.mean([c[0] for c in coords])
    #         cy = np.mean([c[1] for c in coords])

    #         distances = [
    #             math.hypot(coords[i][0] - coords[(i+1) % len(coords)][0],
    #                        coords[i][1] - coords[(i+1) % len(coords)][1])
    #             for i in range(len(coords))
    #         ]
    #         r = np.mean(distances) * 0.38

    #         ax.add_patch(
    #             plt.Circle((cx, cy), r, fill=False, linewidth=3,
    #                        color="#555", zorder=1.5)
    #         )
    # except Exception as e:
    #     print("Aromatic ring drawing failed:", e)

    # ---------------------------------------------------------
    # 6. Draw ligand atom spheres (increase zorder so atoms are on top)
    # ---------------------------------------------------------
    for nm, (x, y) in atom_positions.items():

        # 🔥 HIDE ALL HYDROGEN ATOM NODES
        if nm.startswith("H"):
            continue

        elem = atom_elements.get(nm, "C")
        # Node color: carbon -> black, otherwise element-specific highlights
        if elem == "C":
            node_col = "#000000"
        elif elem == "S":
            node_col = "#FFD700"  # sulfur -> yellow
        elif elem == "P":
            node_col = "#FF69B4"  # phosphorus -> pink
        elif elem in ("F", "Cl", "Br", "I"):
            node_col = "#007849"  # halogens -> Miami green
        elif elem == "O":
            node_col = "#D72638"  # oxygen -> fire/blood red
        elif elem == "N":
            node_col = "#1E90FF"  # nitrogen -> cool blue
        else:
            node_col = ATOM_COLORS.get(elem, "#888")

        # Label: show element symbol for non-carbon heavy atoms.
        if elem == "C":
            label = None
        else:
            label = elem

        # Label color scheme: phosphorus & oxygen -> Miami green, others -> Miami orange
        # All labels should be white
        label_color = "white"
    # atom spheres must sit above bonds/edges
        ax.scatter(
            x, y, s=900, color=node_col,
            edgecolors="black", linewidth=2.2,
            zorder=4
        )
        # add element label for non-carbon heavy atoms
        if label:
            ax.text(x, y, label, fontsize=12, weight="bold",
                    ha="center", va="center", color=label_color, zorder=5)
    # ---------------------------------------------------------
    # 7. Place residues around ligand with automatic no-overlap
    # ---------------------------------------------------------
    residues = df_filt["Residue"].tolist()
    res_positions = auto_place_residues_guided(
        residues, residue_to_ligatom, atom_positions, center
    )

    # ---------------------------------------------------------
    # 8. Draw EDGES — MULTI-LIGATOM SUPPORT
    # ---------------------------------------------------------
    drawn_residues = set()   # collect residues that actually get edges

    for e_idx, (_, row) in enumerate(df_filt.iterrows()):

        res = row["Residue"]
        if res not in res_positions:
            continue

        xR, yR = res_positions[res]

        # Multiple ligand atoms are allowed
        lig_atoms = residue_to_ligatom.get(res, [])
        if isinstance(lig_atoms, str):
            lig_atoms = [lig_atoms]

        # ---------------------------
        # Determine interaction type
        # ---------------------------
        has_hbond = bool(row.get("HasHbond", False))
        hydro = row.get("Hydrophobic", 0) or 0
        water = row.get("Water bridge", 0) or 0
        ionic = row.get("Ionic", 0) or 0

        if has_hbond:
            col, lw = EDGE_HBOND, 2.8
        elif hydro > 0:
            col, lw = EDGE_HYDRO, 1.6
        elif water > 0:
            col, lw = EDGE_WATER, 2.0
        elif ionic > 0:
            col, lw = "#C39BD3", 2.2
        else:
            continue

        # Draw each atom → residue interaction separately
        for lig_atom in lig_atoms:

            # ---------------------------
            # Fix ligand atom coordinate — map hydrogens to heavy atoms
            # ---------------------------
            heavy = _resolve_h_to_heavy(lig_atom, parent_map, atom_positions, atom_elements, center)
            if heavy and heavy in atom_positions:
                xL, yL = atom_positions[heavy]
            else:
                # Unresolved ligand atom — skip this edge to avoid stray center edges
                print(f"⚠ Skipping unresolved ligand atom '{lig_atom}' for residue {res}")
                continue

            # ---------------------------
            # Compute arc geometry
            # ---------------------------
            sx, sy = shorten_to_circle(xL, yL, xR, yR, LIG_ATOM_RADIUS)
            ex, ey = shorten_to_circle(xR, yR, xL, yL, RESIDUE_RADIUS)

            rad = 0.18 if e_idx % 2 == 0 else -0.18

            # draw bicolor curved edge: ligand atom color -> residue node color
            lig_elem = atom_elements.get(heavy, 'C')
            lig_node_col = _atom_node_color(lig_elem)
            res_node_col = RESIDUE_COLORS.get(classify_residue(res), '#999999')
            _draw_bicolor_arc(ax, sx, sy, ex, ey, rad, lig_node_col, res_node_col, lw=lw, linestyle='--', zorder=2)

            drawn_residues.add(res)   # mark residue as having an edge

            # Label (%) — compute pct from the row BEFORE using it
            pct = row.get("Fraction", 0.0) * 100.0
            lx, ly = label_point_frac(sx, sy, ex, ey, 0.75)
            ax.text(
                lx, ly, f"{pct:.0f}%",
                fontsize=9, weight="bold",
                ha="center", va="center", color=col,
                zorder=3
            )

    # ---------------------------------------------------------
    # 8b. Re-draw ligand atom markers with higher zorder so ligand nodes are visually in front of edges.
    for nm, (x, y) in atom_positions.items():
        if nm.startswith("H"):
            continue
        elem = atom_elements.get(nm, "C")
        UM_ORANGE = "#FF6A00"
        UM_GREEN = "#007849"

        # choose node color as above (element-specific overrides)
        if elem == "C":
            node_col = "#000000"
        elif elem == "S":
            node_col = "#FFD700"
        elif elem == "P":
            node_col = "#FF69B4"
        elif elem in ("F", "Cl", "Br", "I"):
            node_col = "#007849"
        elif elem == "O":
            node_col = "#D72638"
        elif elem == "N":
            node_col = "#1E90FF"
        else:
            node_col = ATOM_COLORS.get(elem, "#888")

        if elem == "C":
            label = None
        else:
            label = elem

        # all labels white
        label_color = "white"

        ax.scatter(
            x, y, s=360, color=node_col,
            edgecolors="black", linewidth=1.8,
            zorder=5
        )
        if label:
            ax.text(x, y, label, fontsize=9, weight="bold",
                    ha="center", va="center", color=label_color, zorder=6)

    # ---------------------------------------------------------
    # 9. Draw residue circles + labels ONLY for residues with edges
    # ---------------------------------------------------------
    for res in sorted(drawn_residues, key=residue_number):
        x, y = res_positions[res]
        cat = classify_residue(res)
        col = RESIDUE_COLORS[cat]

        # residue nodes should be above everything else
        ax.scatter(
            x, y, s=3000, color=col,
            edgecolors="black", linewidth=2.0,
            zorder=6
        )
        ax.text(x, y, res,
                fontsize=11, weight="bold",
                ha="center", va="center",
                zorder=7)

def build_two_panel_figure(ligand: str, minfrac: float, output: str, no_edges: bool=False):
    """
    Build figures. 'no_edges' is True when --no-edges was specified.
    """
    # Write all results into a single directory:
    # Analysis_Results/NETWORX_Figures/
    outdir = os.path.join(ANALYSIS_DIR, "NETWORX")
    os.makedirs(outdir, exist_ok=True)

    # -----------------------
    # PANEL A + B DATA
    # -----------------------
    df_summary, total_frames = load_interaction_summary()
    residue_to_ligatom = load_residue_ligand_atom_map()

    # --------------------------------------------------
    # Ensure a single, persistent residue->ligatom mapping:
    # - load ligand layout (needed for similarity mapping / halogen fallback)
    # - apply LP mapping and similarity/halogen fixes ONCE here
    # This avoids situations where one panel uses a different map than another.
    # --------------------------------------------------
    atom_positions, atom_elements, bonds, center, parent_map, mol, atom_names = \
        get_ligand_layout(ligand)
    # first apply LP mapping (if any), then safe fallback mapping
    residue_to_ligatom_fixed = resolve_lp_mapping(ligand, residue_to_ligatom)
    residue_to_ligatom_fixed = fix_edges_to_center_safe(
        residue_to_ligatom_fixed, atom_positions, atom_elements, center
    )
    # quick sanity log
    n_fixed = sum(1 for r in residue_to_ligatom if residue_to_ligatom.get(r) != residue_to_ligatom_fixed.get(r))
    print(f"\n🔁 Applied ligand mapping fixes — {n_fixed} residues remapped (persisted for all panels)\n")

    # replace original mapping so subsequent calls use the fixed one
    residue_to_ligatom = residue_to_ligatom_fixed

    df_filt = df_summary[df_summary["Fraction"] >= minfrac].copy()
    if df_filt.empty:
        raise ValueError(f"No residues exceed minfrac={minfrac}")

    # -----------------------
    # PANEL C DATA — THE ONLY CSV WE ARE ALLOWED TO USE
    # -----------------------
    framewise_path = os.path.join(ANALYSIS_DIR, "InteractionTypes_Framewise.csv")
    df_full = pd.read_csv(framewise_path)

    # Basic validation
    required_cols = ["Time_ns", "Residue", "Type", "LigAtomName"]
    for col in required_cols:
        if col not in df_full.columns:
            raise ValueError(f"InteractionTypes_Framewise.csv missing required column: {col}")

    # Make ligand atom names consistent
    df_full["LigAtomName"] = df_full["LigAtomName"].apply(normalize_atom)

    # -----------------------------------------------------------
    #  COMPUTE FRACTION PER RESIDUE FROM FRAMEWISE DATA ONLY
    # -----------------------------------------------------------
    n_frames = df_full["Time_ns"].nunique()

    # frames where each residue appears at least once
    res_frame_counts = (
        df_full.groupby("Residue")["Time_ns"]
        .nunique()
    )

    # Fraction = (# frames with this residue contact) / (total frames)
    frac_map = (res_frame_counts / float(n_frames)).to_dict()

    df_full["Fraction"] = df_full["Residue"].map(frac_map).fillna(0.0)

    # store in global so Panel C can use it
    global MIN_FRACTION
    MIN_FRACTION = minfrac

    # -----------------------------------------------------------
    # NEW: use dedicated builders for each panel (clean, testable)
    # -----------------------------------------------------------
    # Panel A — bar
    bar_out = make_bar_figure(ligand, df_summary, total_frames, minfrac, outdir, dpi=900)

    # Panel B — collapsed network
    # pass the fixed mapping so collapsed and expanded panels use identical mapping
    net_out = make_collapsed_network_figure(ligand, df_filt, residue_to_ligatom, outdir, dpi=900)

    # Panel C — expanded network (prepare ligand layout first)
    # reuse ligand layout loaded above and call expanded builder with fixed mapping
    rings, bad_atoms = safe_kekulize_and_get_rings(mol, atom_names)
    print("Valid rings:", rings)
    print("Skipped (bad) atoms:", bad_atoms)

    exp_out = make_expanded_network_figure(
        ligand,
        df_full=df_full,
        residue_to_ligatom=residue_to_ligatom,
        atom_positions=atom_positions,
        atom_elements=atom_elements,
        bonds=bonds,
        mol=mol,
        atom_names=atom_names,
        center=center,
        parent_map=parent_map,
        outdir=outdir,
        show_edges=not no_edges,
        dpi=900
    )

    # Also save a clean version (nodes-only) for rapid figure export.
    # Use the public builder but write directly to the clean filename so
    # it does not overwrite the full expanded image produced above.
    clean_out = os.path.join(outdir, f"{ligand}_expanded_clean.png")
    clean_tmp = make_expanded_network_figure(
        ligand,
        df_full=df_full,
        residue_to_ligatom=residue_to_ligatom,
        atom_positions=atom_positions,
        atom_elements=atom_elements,
        bonds=bonds,
        mol=mol,
        atom_names=atom_names,
        center=center,
        parent_map=parent_map,
        outdir=outdir,
        show_edges="single",
        dpi=900,
        out_name=clean_out,
    )
    print(f"  ✓ Clean expanded network saved → {clean_out}")

    # =============================================================
    # Combined PDF (Panels A + Expanded C)
    # =============================================================
    combined_out = os.path.join(outdir, f"{ligand}_combined.pdf")
    figD = plt.figure(figsize=(12, 16))
    gs = gridspec.GridSpec(2, 1, height_ratios=[0.8, 1.2], hspace=0.25)

    # Panel A (bar plot)
    axD1 = figD.add_subplot(gs[0])
    axD1.imshow(plt.imread(bar_out))
    axD1.axis("off")

    # Panel C (expanded network)
    axD2 = figD.add_subplot(gs[1])
    axD2.imshow(plt.imread(exp_out))
    axD2.axis("off")

    figD.savefig(combined_out, dpi=600)
    plt.close(figD)   # make sure we close the correct figure object

    # =============================================================
    # DONE
    # =============================================================
    print("\n===================================================")
    print(f"  ✓ Bar plot saved → {bar_out}")
    print(f"  ✓ Collapsed network saved → {net_out}")
    print(f"  ✓ Expanded network saved → {exp_out}")
    print(f"  ✓ Combined PDF saved → {combined_out}")
    print("===================================================")










# --- Minimal loader + figure-builder fallbacks (keeps script self-contained) ---
def load_interaction_summary():
    """
    Simple loader for Analysis_Results/InteractionTypes_Summary.csv.
    Returns (df, total_frames)
    """
    fpath = os.path.join(ANALYSIS_DIR, "InteractionTypes_Summary.csv")
    if not os.path.exists(fpath):
        raise FileNotFoundError(f"{fpath} not found")

    df = pd.read_csv(fpath)

    if "Residue" not in df.columns:
        raise ValueError("InteractionTypes_Summary.csv missing 'Residue' column")

    # Coerce numeric columns where possible
    for c in df.columns:
        if c == "Residue":
            continue
        try:
            df[c] = pd.to_numeric(df[c], errors="coerce")
        except Exception:
            pass

    if "Total" not in df.columns:
        raise ValueError("InteractionTypes_Summary.csv must contain 'Total' column")

    total_frames = int(df["Total"].max())
    df["Fraction"] = df["Total"] / float(total_frames)

    # boolean convenience flag for H-bonds
    if "H-bond" in df.columns:
        df["HasHbond"] = df["H-bond"].fillna(0) > 0
    else:
        df["HasHbond"] = False

    return df, total_frames


def load_residue_ligand_atom_map():
    """
    Build residue -> ligand-atom map by preferring FilteredContacts over full framewise file.
    Returns mapping dict.
    """
    f_filt = os.path.join(ANALYSIS_DIR, "FilteredContacts_Framewise.csv")
    f_full = os.path.join(ANALYSIS_DIR, "InteractionTypes_Framewise.csv")

    if not os.path.exists(f_filt):
        raise FileNotFoundError(f"{f_filt} not found")
    if not os.path.exists(f_full):
        raise FileNotFoundError(f"{f_full} not found")

    df_filt = pd.read_csv(f_filt)
    df_full = pd.read_csv(f_full)

    # safe normalization
    if "LigAtomName" in df_filt.columns:
        df_filt["LigAtomName"] = df_filt["LigAtomName"].astype(str)
    else:
        df_filt["LigAtomName"] = ""

    if "LigAtomName" in df_full.columns:
        df_full["LigAtomName"] = df_full["LigAtomName"].astype(str)
    else:
        df_full["LigAtomName"] = ""

    # most frequent atom per residue
    if "Residue" in df_filt.columns and not df_filt.empty:
        filt_map = (
            df_filt.groupby("Residue")["LigAtomName"]
            .agg(lambda x: x.value_counts().idxmax())
            .apply(normalize_atom)
            .to_dict()
        )
    else:
        filt_map = {}

    if "Residue" in df_full.columns and not df_full.empty:
        full_map = (
            df_full.groupby("Residue")["LigAtomName"]
            .agg(lambda x: x.value_counts().idxmax())
            .apply(normalize_atom)
            .to_dict()
        )
    else:
        full_map = {}

    merged = full_map.copy()
    merged.update(filt_map)

    # minimal logging
    print("\n🔗 Residue → LigAtom (merged):")
    for k, v in sorted(merged.items()):
        print(f"  {k} → {v}")

    return merged


def make_bar_figure(ligand, df_summary, total_frames, minfrac, outdir, dpi=600):
    """Create and save Panel A (stacked bar). Returns path."""
    fig, ax = plt.subplots(figsize=(10, 4))
    draw_contact_bar_panel(ax, df_summary, total_frames, minfrac)
    path = os.path.join(outdir, f"{ligand}_bar.png")
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ Bar panel saved → {path}")
    return path


def make_collapsed_network_figure(ligand, df_filt, residue_to_ligatom, outdir, dpi=600):
    """Create and save Panel B (collapsed network). Returns path."""
    fig, ax = plt.subplots(figsize=(12, 10))
    draw_network_panel(ax, ligand, df_filt, residue_to_ligatom)
    path = os.path.join(outdir, f"{ligand}_network.png")
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ Collapsed network saved → {path}")
    return path


def make_expanded_network_figure(ligand, df_full, residue_to_ligatom,
    atom_positions, atom_elements, bonds,
    mol, atom_names, center, parent_map,
    outdir, show_edges=True, dpi=600, out_name=None):
    """Create and save Panel C (expanded network). Returns path.

    show_edges may be:
      - True: full expanded rendering
      - False: nodes-only (no edges)
      - "single": simplified fallback that draws one representative edge per residue
    """
    fig, ax = plt.subplots(figsize=(12, 10))

    # If a full implementation exists in the globals -> call it (preferred),
    # except when the caller specifically requested a 'single' simplified export.
    panel_fn = globals().get("draw_expanded_network_panel", None)
    if callable(panel_fn) and show_edges != "single":
        panel_fn(
            ax,
            ligand=ligand,
            df_full=df_full,
            residue_to_ligatom=residue_to_ligatom,
            atom_positions=atom_positions,
            atom_elements=atom_elements,
            bonds=bonds,
            mol=mol,
            atom_names=atom_names,
            center=center,
            parent_map=parent_map,
            show_edges=show_edges,
            show_edge_labels=show_edges,
        )
        if out_name:
            path = out_name
        else:
            path = os.path.join(outdir, f"{ligand}_expanded.png")
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  ✓ Expanded network saved → {path}")
        return path

    # Fallback branch (either no panel_fn, or show_edges == "single")
    ax.set_aspect("equal")
    ax.axis("off")

    # 1) ligand atoms
    for nm, (x, y) in atom_positions.items():
        if nm.startswith("H"):
            continue
        elem = atom_elements.get(nm, "C")
        # node color mapping (carbon -> black)
        if elem == "C":
            node_col = "#000000"
        elif elem == "S":
            node_col = "#FFD700"
        elif elem == "P":
            node_col = "#FF69B4"
        elif elem in ("F", "Cl", "Br", "I"):
            node_col = "#007849"
        elif elem == "O":
            node_col = "#D72638"
        elif elem == "N":
            node_col = "#1E90FF"
        else:
            node_col = ATOM_COLORS.get(elem, "#888")

        # label presence & color (all labels white)
        if elem == "C":
            label = None
        else:
            label = elem
        label_color = "white"

        ax.scatter(x, y, s=900, color=node_col, edgecolors="black", linewidth=2.2, zorder=3)
        if label:
            ax.text(x, y, label, fontsize=11, weight="bold", ha="center", va="center", color=label_color, zorder=4)

    # 2) ligand bonds
    for a1_idx, a2_idx, order in bonds:
        a1 = normalize_bond_endpoint(a1_idx, atom_names)
        a2 = normalize_bond_endpoint(a2_idx, atom_names)
        if a1 is None or a2 is None:
            continue
        if a1.startswith("H") or a2.startswith("H"):
            continue
        x1, y1 = atom_positions[a1]
        x2, y2 = atom_positions[a2]
        draw_bond(ax, x1, y1, x2, y2, order, atom_elements.get(a1, "C"), atom_elements.get(a2, "C"))

    # 3) aggregated edge table
    df_edges = df_full.copy()
    df_edges["LigAtomName"] = df_edges["LigAtomName"].apply(normalize_atom)
    df_edges = df_edges[(df_edges["Fraction"] >= MIN_FRACTION)]

    if df_edges.empty:
        if out_name:
            path = out_name
        else:
            path = os.path.join(outdir, f"{ligand}_expanded.png")
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  ✓ Expanded network (empty) saved → {path}")
        return path

    # collapse to unique Residue,LigAtom,Type with max fraction
    df_edges = df_edges.groupby(["Residue", "LigAtomName", "Type"], as_index=False)["Fraction"].max()

    residues_with_edges = sorted(df_edges["Residue"].unique(), key=residue_number)
    if not residues_with_edges:
        if out_name:
            path = out_name
        else:
            path = os.path.join(outdir, f"{ligand}_expanded.png")
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  ✓ Expanded network (no residues) saved → {path}")
        return path

    residue_pos, rx, ry = auto_place_residues_no_overlap(residues_with_edges, atom_positions)

    # draw residue nodes
    for res in residues_with_edges:
        x, y = residue_pos[res]
        cat = classify_residue(res)
        col = RESIDUE_COLORS[cat]
        ax.scatter(x, y, s=3000, color=col, edgecolors="black", linewidth=2.0, zorder=6)
        ax.text(x, y, res, fontsize=11, weight="bold", ha="center", va="center", zorder=7)

    # nodes-only requested -> finish here
    if show_edges is False:
        if out_name:
            path = out_name
        else:
            path = os.path.join(outdir, f"{ligand}_expanded.png")
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  ✓ Expanded network (nodes-only) saved → {path}")
        return path

    # single-edge-per-residue mode
    if show_edges == "single":
        idx = df_edges.groupby("Residue")["Fraction"].idxmax()
        best = df_edges.loc[idx].reset_index(drop=True)
        for e_idx, row in best.iterrows():
            res = row["Residue"]
            if res not in residue_pos:
                continue
            xR, yR = residue_pos[res]
            raw_atom = residue_to_ligatom.get(res, row["LigAtomName"])
            atom = normalize_atom(raw_atom)

            if atom in parent_map and parent_map.get(atom) in atom_positions:
                atom = parent_map[atom]
            if atom not in atom_positions and atom.startswith("H"):
                root = atom.rstrip("0123456789")
                for h_mol2, heavy in parent_map.items():
                    if h_mol2.startswith(root) and heavy in atom_positions:
                        atom = heavy
                        break

            if atom and atom in atom_positions:
                xL, yL = atom_positions[atom]
            else:
                print(f"⚠ Skipping unresolved ligand atom '{raw_atom}' for residue {res}")
                continue

            sx, sy = shorten_to_circle(xL, yL, xR, yR, LIG_ATOM_RADIUS)
            ex, ey = shorten_to_circle(xR, yR, xL, yL, RESIDUE_RADIUS)

            itype = str(row["Type"])
            if "H-bond" in itype:
                col, lw = EDGE_HBOND, 2.5
            elif "Hydrophobic" in itype:
                col, lw = EDGE_HYDRO, 2.0
            elif "Water" in itype:
                col, lw = EDGE_WATER, 2.0
            elif "Ionic" in itype:
                col, lw = "#C39BD3", 2.3
            else:
                col, lw = "#999999", 1.8

            rad = 0.18 if (e_idx % 2 == 0) else -0.18
            # draw curved edge colored by interaction type (match legend)
            edge_patch = FancyArrowPatch((sx, sy), (ex, ey), connectionstyle=f"arc3,rad={rad}", arrowstyle="-",
                                        linestyle="--", linewidth=lw, color=col, zorder=2)
            ax.add_patch(edge_patch)

            pct = float(row["Fraction"]) * 100.0
            lx, ly = label_point_frac(sx, sy, ex, ey, 0.75)
            ax.text(lx, ly, f"{pct:.0f}%", fontsize=8, weight="bold", ha="center", va="center", color=col, zorder=5)

        make_residue_legend(ax)
        if out_name:
            path = out_name
        else:
            path = os.path.join(outdir, f"{ligand}_expanded.png")
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  ✓ Expanded network (single-edge) saved → {path}")
        return path

    # full fallback: draw all edges
    for e_idx, row in df_edges.iterrows():
        res = row["Residue"]
        if res not in residue_pos:
            continue
        xR, yR = residue_pos[res]

        raw_atom = residue_to_ligatom.get(res, row["LigAtomName"])
        atom = normalize_atom(raw_atom)

        if atom in parent_map and parent_map.get(atom) in atom_positions:
            atom = parent_map[atom]
        if atom not in atom_positions and atom.startswith("H"):
            root = atom.rstrip("0123456789")
            for h_mol2, heavy in parent_map.items():
                if h_mol2.startswith(root) and heavy in atom_positions:
                    atom = heavy
                    break

        if atom and atom in atom_positions:
            xL, yL = atom_positions[atom]
        else:
            print(f"⚠ Skipping unresolved ligand atom '{raw_atom}' for residue {res}")
            continue

        sx, sy = shorten_to_circle(xL, yL, xR, yR, LIG_ATOM_RADIUS)
        ex, ey = shorten_to_circle(xR, yR, xL, yL, RESIDUE_RADIUS)

        itype = str(row["Type"])
        if "H-bond" in itype:
            col, lw = EDGE_HBOND, 2.5
        elif "Hydrophobic" in itype:
            col, lw = EDGE_HYDRO, 2.0
        elif "Water" in itype:
            col, lw = EDGE_WATER, 2.0
        elif "Ionic" in itype:
            col, lw = "#C39BD3", 2.3
        else:
            col, lw = "#999999", 1.8

        rad = 0.18 if (e_idx % 2 == 0) else -0.18
        edge = FancyArrowPatch((sx, sy), (ex, ey), connectionstyle=f"arc3,rad={rad}", arrowstyle="-",
                        linestyle="--", linewidth=lw, color=col, zorder=2)
        ax.add_patch(edge)

        pct = float(row["Fraction"]) * 100.0
        lx, ly = label_point_frac(sx, sy, ex, ey, 0.75)
        ax.text(lx, ly, f"{pct:.0f}%", fontsize=8, weight="bold", ha="center", va="center", color=col, zorder=5)

    make_residue_legend(ax)

    if out_name:
        path = out_name
    else:
        path = os.path.join(outdir, f"{ligand}_expanded.png")
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ Expanded network saved → {path}")
    return path
# --------------------------------------------------------------------
# Main
# --------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.ligand or args.minfrac is None or not args.output:
        args = ask_if_missing(args)

    build_two_panel_figure(args.ligand, args.minfrac, args.output, no_edges=args.no_edges)

