#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MD ANALYSIS FIGUREBOOK BUILDER
Joseph-Michael Edition — PAGINATION-SAFE VERSION

✔ All text wraps
✔ No overflow
✔ Title page guaranteed content
✔ Publication-ready PDF
"""

import os
import re
import string
import tempfile
import subprocess  # Import subprocess for Open Babel calls
from PIL import Image

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from reportlab.platypus import Image as RLImage

from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Image as RLImage, PageBreak,
    BaseDocTemplate, Frame, PageTemplate
)
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.pagesizes import letter, landscape  # Import landscape
from reportlab.lib.units import inch
from reportlab.platypus import Image as RLImage
from svglib.svglib import svg2rlg  # Import for SVG rendering
from reportlab.graphics import renderPDF
from reportlab.graphics.shapes import Drawing

TITLE_PAGE_DESCRIPTION = """
Molecular Dynamics Analysis of {LIG} in Complex with {Protein}

Trajectory Generated with GROMACS and the Customized PyMACS Analysis Framework

This report presents a comprehensive molecular dynamics (MD) analysis of a drug-bound protein complex, integrating structural, dynamic, and interaction-level metrics into a unified, publication-ready dataset. The simulation was performed using GROMACS, and all post-processing was executed using our in-house PyMACS analysis suite—a customized Python-based framework purpose-built for high-resolution drug binding trajectory evaluation.

The PyMACS pipeline couples GPU-accelerated computation with curated structural bioinformatics algorithms to quantify global stability, chain-specific motions, ligand flexibility, secondary structure evolution, and the detailed geometry of residue-level interactions. Together, these metrics allow precise characterization of how the ligand engages the protein interface, how that interface adapts over time, and which structural features contribute most to binding stability and specificity.

Across this report, each figure provides a distinct analytical lens:

Global and chain-resolved RMSD/RMSF reveal the stability, convergence, and dynamic hotspots of the complex.

Ligand RMSD/RMSF track binding mode stability and internal flexibility.

Secondary structure analyses (DSSP) highlight folding transitions, helix/sheet integrity, and ligand-induced rearrangements.

Time-resolved contact maps, frequency ranking, and interaction-type distributions capture the full chemical and temporal landscape of ligand recognition.

Distance distributions, persistence metrics, and interaction timelines decode the geometric and kinetic properties of individual contacts.

Interaction networks provide structural and quantitative visualization of binding interfaces, highlighting the key residues that anchor or modulate ligand engagement.

By integrating these orthogonal analyses, this report delivers a complete mechanistic portrait of the ligand–protein system: how stable the complex remains over time, how the ligand adapts or rigidifies within its pocket, which residues dominate affinity, and which structural rearrangements accompany sustained binding.

This document is intended as a rigorous, standalone reference for comparing drug designs, validating docking hypotheses, guiding mutagenesis or scaffold optimization, and supporting manuscript-quality mechanistic interpretation.
"""

def safe_image(img_path, max_width, max_height):
    img = RLImage(img_path)
    img._restrictSize(max_width, max_height)
    return img

# ============================================================
# IMAGE PREP
# ============================================================
def prepare_image_for_pdf(img_path, max_width_px=1200):
    if not os.path.exists(img_path):
        return None

    img = Image.open(img_path)
    w, h = img.size

    if w > max_width_px:
        scale = max_width_px / w
        img = img.resize((int(w*scale), int(h*scale)), Image.LANCZOS)

    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".jpg")
    img.convert("RGB").save(tmp.name, "JPEG", quality=90)
    return tmp.name


# ============================================================
# LIGAND DETECTION
# ============================================================
def detect_ligand_mol2():
    mol2 = [f for f in os.listdir('.') if f.endswith(".cgenff.mol2")]
    if not mol2:
        return input("Enter ligand code: ").strip().upper()

    if len(mol2) == 1:
        code = mol2[0].replace(".cgenff.mol2", "")
        if input(f"Detected {code}. Use? [Y/n]: ").lower() in ("", "y", "yes"):
            return code.upper()

    for i, f in enumerate(mol2, 1):
        print(f"{i}) {f}")
    sel = input("Pick # or ENTER: ").strip()
    if sel.isdigit():
        return mol2[int(sel)-1].replace(".cgenff.mol2", "").upper()
    return input("Enter ligand code: ").strip().upper()


# ============================================================
# ATOMINDEX PARSER
# ============================================================
def parse_atomindex():
    if not os.path.exists("atomIndex.txt"):
        return ["Protein"], [0]

    proteins = []
    with open("atomIndex.txt") as f:
        for l in f:
            if l.strip():
                proteins.append(l.split()[0])

    return proteins, list(range(len(proteins)))


# ============================================================
# CAPTION LOADER
# ============================================================
def load_main_captions(ligand, chains, proteins):
    if not os.path.exists("4_MDfigs.txt"):
        return []

    with open("4_MDfigs.txt") as f:
        blocks = f.read().split("---")

    out = []
    for b in blocks:
        lines = [x.strip() for x in b.splitlines() if x.strip()]
        if len(lines) < 2:
            continue

        fname = lines[0].replace("{LIG}", ligand)
        caption = "<br/>".join(lines[1:]).replace("{LIG}", ligand)

        if "{CHAIN}" in fname:
            for c in chains:
                out.append((fname.replace("{CHAIN}", str(c)),
                            caption.replace("{CHAIN}", str(c))))
            continue

        if "{PROTEIN_NAME}" in fname:
            for p in proteins:
                out.append((fname.replace("{PROTEIN_NAME}", p),
                            caption.replace("{PROTEIN_NAME}", p)))
            continue

        out.append((fname, caption))

    return out


# ============================================================
# FILE RESOLUTION
# ============================================================
def resolve_possible_files(target):
    base = "Analysis_Results"
    matches = []

    for root, _, files in os.walk(base):
        for f in files:
            if f == target or ("*" in target and re.fullmatch(target.replace("*", ".*"), f)):
                matches.append(os.path.join(root, f))
    return matches


# ============================================================
# HERO IMAGE
# ============================================================
def generate_rdkit_hero_image(smiles, width=900, height=900):
    """
    Generates a stylized 2D ligand image from a SMILES string.
    """
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol, addCoords=True)



    # Compute 2D coordinates for the molecule
    Chem.rdDepictor.Compute2DCoords(mol)

    # Set up the drawing options for a more stylized appearance
    drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
    options = drawer.drawOptions()
    options.bondLineWidth = 2  # Thicker bonds
    options.addAtomIndices = False  # No atom indices
    options.fixedBondLength = 25  # Standardize bond lengths
    options.dotsPerAngstrom = 100  # High resolution
    options.useBWAtomPalette()  # Black-and-white atom palette for publication-ready images

    # Draw the molecule
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    # Save the image to a temporary file
    out = tempfile.NamedTemporaryFile(delete=False, suffix=".png").name
    with open(out, "wb") as f:
        f.write(drawer.GetDrawingText())
    return out


def generate_rdkit_hero_image_from_mol2(mol2_file, width=900, height=900, use_color=True):
    """
    Generates a stylized 2D ligand image from a .mol2 file using Open Babel and RDKit.
    """
    # Convert .mol2 to SMILES using Open Babel
    try:
        result = subprocess.run(
            ["obabel", mol2_file, "-osmi", "-xk"],  # Open Babel command
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        smiles = result.stdout.strip()
        if not smiles:
            raise ValueError(f"Open Babel failed to generate SMILES for {mol2_file}")
    except FileNotFoundError:
        raise FileNotFoundError("Open Babel (obabel) is not installed or not in PATH.")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Open Babel error: {e.stderr.strip()}")

    # Generate the RDKit molecule from the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ValueError(f"RDKit failed to parse SMILES: {smiles}")

    # Compute 2D coordinates for the molecule
    Chem.rdDepictor.Compute2DCoords(mol)

    # Set up the drawing options for a more stylized appearance
    drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
    options = drawer.drawOptions()
    options.bondLineWidth = 2  # Thicker bonds
    options.addAtomIndices = False  # No atom indices
    options.fixedBondLength = 25  # Standardize bond lengths
    options.dotsPerAngstrom = 100  # High resolution

    if use_color:
        # Use the default color palette for atoms
        pass  # No need to call useBWAtomPalette for color mode
    else:
        options.useBWAtomPalette()  # Black-and-white atom palette for publication-ready images

    # Draw the molecule
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    # Save the image to a temporary file
    out = tempfile.NamedTemporaryFile(delete=False, suffix=".png").name
    with open(out, "wb") as f:
        f.write(drawer.GetDrawingText())
    return out


# ============================================================
# HEADER FUNCTION
# ============================================================
def add_header(canvas, doc):
    """
    Adds a unified header to all pages, including the first page.
    """
    canvas.saveState()
    if doc.page == 1:
        # Custom header for the first page
        header_text = "PyMacs Molecular Dynamics Simulation Report"
        canvas.setFont("Helvetica-Bold", 14)  # Larger font for the first page
        canvas.drawCentredString(doc.pagesize[0] / 2, doc.pagesize[1] - 30, header_text)  # Centered at the top
    else:
        # Header for subsequent pages with page number appended
        header_text = f"PyMacs Molecular Dynamics Simulation Report Page {doc.page}"
        canvas.setFont("Helvetica-Bold", 10)
        canvas.drawCentredString(doc.pagesize[0] / 2, doc.pagesize[1] - 30, header_text)  # Centered at the top
    canvas.restoreState()


# ============================================================
# ADD LIGAND TO PAGES FUNCTION
# ============================================================
def add_ligand_to_pages(canvas, doc, small_ligand_svg, y_offset=10):
    """
    Adds a small ligand image to every other page starting from the first page.
    The image is placed in the bottom-left corner, ignoring margin constraints.
    
    Parameters:
        canvas: The canvas object for the page.
        doc: The document object.
        small_ligand_svg: The SVG file path for the ligand image.
        y_offset: The vertical offset for the image (default is 10).
    """
    if doc.page % 2 == 1:  # Add the ligand image to odd-numbered pages
        canvas.saveState()
        drawing = svg2rlg(small_ligand_svg)
        drawing.width = 0.8 * inch  # Adjust the size of the ligand image
        drawing.height = 0.8 * inch
        renderPDF.draw(drawing, canvas, 40, y_offset)  # Position at (40, y_offset), ignoring margins
        canvas.restoreState()


# ============================================================
# BUILD PDF
# ============================================================
def build_pdf():
    ligand = detect_ligand_mol2()
    pretty = input(f"Pretty ligand name for {ligand}: ").strip()

    proteins, chains = parse_atomindex()
    protein_title = " & ".join(proteins)

    captions = load_main_captions(ligand, chains, proteins)

    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(
        name="TitleBig",
        fontSize=22,
        spaceAfter=24,
        alignment=1,
        leading=28  # Adjust line spacing for multi-line titles
    ))
    styles.add(ParagraphStyle(
        name="Body",
        fontSize=11,
        leading=14,
        spaceAfter=12
    ))
    styles.add(ParagraphStyle(
        name="Caption",
        fontSize=10,
        leading=13,
        spaceBefore=10
    ))

    # Set the document to landscape orientation with a reduced bottom margin
    doc = BaseDocTemplate(
        "MD_ANALYSIS_FIGUREBOOK.pdf",
        pagesize=landscape(letter),  # Landscape orientation
        rightMargin=50,
        leftMargin=50,
        topMargin=50,
        bottomMargin=10  # Reduced bottom margin to allow the image to go lower
    )

    # Generate the small ligand SVG for the bottom-left corner
    mol2_file = f"{ligand}.cgenff.mol2"
    small_ligand_svg = generate_rdkit_hero_image_as_svg(mol2_file, width=300, height=300, use_color=True)

    # Define the frame for content
    frame = Frame(
        doc.leftMargin,
        doc.bottomMargin,
        doc.width,
        doc.height - 40,  # Leave space for the header
        id="normal"
    )

    # Add separate page templates for the first page and subsequent pages
    doc.addPageTemplates([
        PageTemplate(
            id="all_pages",
            frames=frame,
            onPage=lambda canvas, doc: (add_header(canvas, doc), add_ligand_to_pages(canvas, doc, small_ligand_svg, y_offset=20))
        )
    ])

    flow = []

    # ---------------- Title Page ----------------
    title_text = f"MD ANALYSIS REPORT — {pretty} Bound to {protein_title}"
    flow.append(Paragraph(title_text, styles["TitleBig"]))

    # Add the TITLE_PAGE_DESCRIPTION to the title page
    description = TITLE_PAGE_DESCRIPTION.format(LIG=pretty, Protein=protein_title)
    flow.append(Paragraph(description, styles["Body"]))  # Add the description as a paragraph
    flow.append(PageBreak())

    # ---------------- Figures ----------------
    # Process main captions and figures
    for fname, caption in captions:
        direct = os.path.join("Analysis_Results", fname)
        files = [direct] if os.path.exists(direct) else resolve_possible_files(fname)

        if not files:
            print(f"⚠️ Figure not found: {fname}")
            continue

        for f in files:
            img = prepare_image_for_pdf(f)
            if not img:
                print(f"⚠️ Unable to process figure: {f}")
                continue

            # Add the figure on one page
            flow.append(RLImage(img, width=9 * inch, height=6 * inch))  # Adjusted for landscape
            flow.append(PageBreak())  # Ensure the figure is on its own page

            # Add the caption on the next page
            flow.append(Paragraph(caption, styles["Caption"]))
            flow.append(PageBreak())  # Ensure the caption is on its own page

    # ---------------- Distance Histograms ----------------
    # Handle residue-specific distance distribution plots
    DIST_DIR = os.path.join("Analysis_Results", "RESIDUE_Contact_Distance")
    DIST_CAPTION = """Figure 30{suffix}. Interaction Distance Distribution
Narrow peaks suggest stable geometry; broad peaks indicate dynamic or transient contacts.
"""

    if os.path.exists(DIST_DIR):
        dists = sorted(f for f in os.listdir(DIST_DIR) if f.startswith("dist_"))
        if not dists:
            print(f"⚠️ No distance histograms found in {DIST_DIR}")
        for i, dist_file in enumerate(dists):
            suffix = chr(97 + i)  # Generate suffixes a, b, c, etc.
            img = os.path.join(DIST_DIR, dist_file)
            prepared_img = prepare_image_for_pdf(img)

            if not prepared_img:
                print(f"⚠️ Unable to process distance histogram: {img}")
                continue

            # Add the figure on one page
            flow.append(RLImage(prepared_img, width=9 * inch, height=6 * inch))
            flow.append(PageBreak())  # Ensure the figure is on its own page

            # Add the caption on the next page
            flow.append(Paragraph(DIST_CAPTION.format(suffix=suffix), styles["Caption"]))
            flow.append(PageBreak())  # Ensure the caption is on its own page

    else:
        print(f"⚠️ Distance histogram directory not found: {DIST_DIR}")

    doc.build(flow)
    print("📘 MD_ANALYSIS_FIGUREBOOK.pdf COMPLETE")


# ============================================================
# HERO IMAGE AS SVG
# ============================================================
def generate_rdkit_hero_image_as_svg(mol2_file, width=900, height=900, use_color=True):
    result = subprocess.run(
        ["obabel", mol2_file, "-osmi", "-xk"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True
    )

    smiles = result.stdout.strip().split()[0]  # 🔑 STRIP UNK

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        raise ValueError(f"RDKit failed to parse SMILES: {smiles}")

    Chem.SanitizeMol(
        mol,
        sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES
    )

    Chem.rdDepictor.Compute2DCoords(mol)

    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    options = drawer.drawOptions()
    options.bondLineWidth = 2
    options.fixedBondLength = 25
    options.dotsPerAngstrom = 100

    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    out = tempfile.NamedTemporaryFile(delete=False, suffix=".svg").name
    with open(out, "w") as f:
        f.write(drawer.GetDrawingText())

    return out




# ============================================================
if __name__ == "__main__":
    build_pdf()
