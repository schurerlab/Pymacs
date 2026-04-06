import MDAnalysis as mda
from MDAnalysis.analysis import rms, distances
from MDAnalysis.analysis.rms import RMSD, RMSF
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis
from MDAnalysis.analysis.contacts import Contacts

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re
import glob
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

import mdtraj as md
from scipy.cluster.hierarchy import linkage, fcluster
import networkx as nx

from rdkit import Chem
from rdkit.Chem import rdMolTransforms, AllChem
from openbabel import openbabel
from sklearn.decomposition import PCA
from collections import Counter
from mdtraj import compute_dssp
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rcParams
import matplotlib.font_manager as fm




from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rcParams
import numpy as np
import pandas as pd
import MDAnalysis as mda
import MDAnalysis.analysis.distances as distances
import matplotlib.pyplot as plt
import os
import re

import MDAnalysis as mda
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import plotly.graph_objects as go



# Create output directory if it doesn't exist
OUTPUT_DIR = "Analysis_Graphs"
os.makedirs(OUTPUT_DIR, exist_ok=True)  # Ensure folder exist



# ---- CONFIGURATION ----
PDB_REFERENCE = "Complex.pdb"  # PDB file with correct chain info
TOPOLOGY_FILE = "npt.gro"  # Topology filex
TRAJECTORY_FILE = "md_0_1.xtc"  # Trajectory file
LIGAND_NAME = "PTC"  # Update if your ligand has a different name
DISTANCE_CUTOFF = 4.0  # Ligand-Protein interaction cutoff in Å

# Define distance cutoff for interactions in pLOTLY graphs (e.g., 3.5 Å)
CUTOFF = 3.5


# Assign distinct colors for each chain
CHAIN_COLORS = ["blue", "red", "green", "purple", "orange", "brown", "pink"]

# Map chains to custom labels for the legend
CHAIN_LABELS = {
    "A": "Pomalidomide",
    "B": "HIV Protease (Monomer 1)",
    "C": "HIV Protease (Monomer 2)"
}




# ---- LOAD THE UNIVERSE ----
print("🔄 Loading trajectory...")
try:
    u = mda.Universe(TOPOLOGY_FILE, TRAJECTORY_FILE)
    protein = u.select_atoms("protein")
    ligand = u.select_atoms(f"resname {LIGAND_NAME}")
    print(f"✅ Loaded system with {len(u.atoms)} atoms and {u.trajectory.n_frames} frames.")
except Exception as e:
    print(f"❌ Error loading trajectory: {e}")
    exit()




# Define chain cutoffs manually
CHAIN_A_START = 71  # Adjust based on system
CHAIN_B_START = 1   # First reset after A
CHAIN_B_END = 99    # Define where Chain B ends
CHAIN_C_START = CHAIN_B_END + 1  # Chain C starts right after B ends

def assign_chain(residue_id):
    """Assigns a chain label based on residue numbering logic."""
    if residue_id >= CHAIN_A_START:  # First chain
        return "A"
    elif residue_id >= CHAIN_B_START and residue_id <= CHAIN_B_END:  # Second chain
        return "B"
    else:  # Third chain starts after Chain B ends
        return "C"

def mda_to_rdkit(mda_ligand):
    """Convert an MDAnalysis ligand selection to an RDKit molecule"""
    pdb_block = mda_ligand.write("PDB")

    if pdb_block is None:
        raise ValueError("PDB block is None. MDAnalysis failed to write ligand PDB.")
    
    if not pdb_block.strip():
        raise ValueError("Generated PDB block is empty. Check ligand selection.")

    print("Generated PDB Block (first 500 chars):\n", pdb_block[:500])

    mol = Chem.MolFromPDBBlock(pdb_block)
    if mol is None:
        raise ValueError("RDKit failed to parse the PDB block.")
    
    if mol.GetNumConformers() == 0:
        AllChem.EmbedMolecule(mol, useRandomCoords=True)  # Ensure 3D coordinates
    
    return mol



# ---- RMSD ANALYSIS ----
print("🔄 Running RMSD analysis...")
try:
    rmsd_analysis = rms.RMSD(u, select="protein", ref_frame=0)
    rmsd_analysis.run()
    print("✅ RMSD analysis completed.")

    plt.figure(figsize=(8, 5))
    plt.plot(rmsd_analysis.results.rmsd[:, 0], rmsd_analysis.results.rmsd[:, 2], label="Protein RMSD")
    plt.xlabel("Time (Frames)")
    plt.ylabel("RMSD (Å)")
    plt.title("RMSD Over Time")
    plt.legend()
    plt.savefig(os.path.join(OUTPUT_DIR, "rmsd_plot.png"))
    plt.close()
except Exception as e:
    print(f"❌ Error in RMSD analysis: {e}")


print("🔄 Running RMSF analysis (Optimized)...")
try:
    rmsf_analysis = RMSF(protein).run()
    rmsf_values = rmsf_analysis.results.rmsf
    print("✅ RMSF analysis completed.")

    plt.figure(figsize=(8, 5))
    plt.plot(protein.resids, rmsf_values, label="Protein RMSF")
    plt.xlabel("Residue")
    plt.ylabel("RMSF (Å)")
    plt.title("Residue Flexibility (RMSF)")
    plt.legend()
    plt.savefig(os.path.join(OUTPUT_DIR,"rmsf_plot.png"))
    plt.close()
except Exception as e:
    print(f"❌ Error in RMSF analysis: {e}")


    print("✅ RMSF analysis completed.")

    plt.figure(figsize=(8, 5))
    plt.plot(protein.resids, rmsf_values, label="Protein RMSF")
    plt.xlabel("Residue")
    plt.ylabel("RMSF (Å)")
    plt.title("Residue Flexibility (RMSF)")
    plt.legend()
    plt.savefig(os.path.join(OUTPUT_DIR,"rmsf_plot.png"))
    plt.close()
except Exception as e:
    print(f"❌ Error in RMSF analysis: {e}")

# ---- RADIUS OF GYRATION ----
print("🔄 Running Radius of Gyration analysis...")
try:
    rg_values = []
    for ts in u.trajectory:
        rg_values.append(protein.radius_of_gyration())

    print("✅ Radius of Gyration analysis completed.")

    plt.figure(figsize=(8, 5))
    plt.plot(rg_values, label="Radius of Gyration")
    plt.xlabel("Time (frames)")
    plt.ylabel("Radius of Gyration (Å)")
    plt.title("Protein Compactness Over Time")
    plt.legend()
    plt.savefig(os.path.join(OUTPUT_DIR,"rg_plot.png"))
    plt.close()
except Exception as e:
    print(f"❌ Error in Radius of Gyration analysis: {e}")

# ---- HYDROGEN BOND ANALYSIS ----
print("🔄 Running Hydrogen Bond analysis (Geometric Method)...")
try:
    hbond_analysis = HydrogenBondAnalysis(u, 'protein', 'protein')
    hbond_analysis.run()

    print("✅ Hydrogen Bond analysis completed.")

    plt.figure(figsize=(8, 5))
    plt.plot(hbond_analysis.results.times, hbond_analysis.results.count_by_time(), label="H-Bonds")
    plt.xlabel("Time (ps)")
    plt.ylabel("Number of Hydrogen Bonds")
    plt.title("Hydrogen Bonds Over Time (Geometric)")
    plt.legend()
    plt.savefig(os.path.join(OUTPUT_DIR, "hbonds_plot.png"))
    plt.close()
except Exception as e:
    print(f"❌ Error in Hydrogen Bond analysis: {e}")


print("🔄 Computing Intramolecular Hydrogen Bonds (intraHB)...")

try:
    hbond_analysis = HydrogenBondAnalysis(u, 'resname PTC', 'resname PTC')  # Self H-bonds
    hbond_analysis.run()

    plt.figure(figsize=(8, 5))
    plt.plot(hbond_analysis.results.times, hbond_analysis.results.count_by_time(), label="Ligand intraHB", color="orange")
    plt.xlabel("Time (ps)")
    plt.ylabel("Number of H-Bonds")
    plt.title("Intramolecular Hydrogen Bonds Over Time")
    plt.legend()
    plt.savefig(os.path.join(OUTPUT_DIR, "ligand_intraHB_plot.png"))
    plt.close()

    print("✅ Ligand intraHB analysis completed.")
except Exception as e:
    print(f"❌ Error in Ligand intraHB analysis: {e}")


# ---- LIGAND-PROTEIN INTERACTIONS ----
print(f"🔄 Running Ligand ({LIGAND_NAME})-Protein interaction analysis...")
try:
    contact_counts = []
    for ts in u.trajectory:
        distances_matrix = distances.distance_array(ligand.positions, protein.positions)
        contacts = np.sum(distances_matrix < DISTANCE_CUTOFF)  # Contacts within 4Å
        contact_counts.append(contacts)

    print("✅ Ligand-Protein interaction analysis completed.")

    plt.figure(figsize=(8, 5))
    plt.plot(contact_counts, label=f"Ligand ({LIGAND_NAME})-Protein Contacts")
    plt.xlabel("Time (frames)")
    plt.ylabel("Number of Contacts")
    plt.title(f"Ligand ({LIGAND_NAME}) Contacts with Protein Over Time")
    plt.legend()
    plt.savefig(os.path.join(OUTPUT_DIR,"ligand_contacts_plot.png"))
    plt.close()
except Exception as e:
    print(f"❌ Error in Ligand-Protein interaction analysis: {e}")



print(f"🔄 Generating Contact Map for Ligand ({LIGAND_NAME})-Protein Interactions...")

try:
    interaction_data = []

    for ts in u.trajectory:
        frame = ts.frame
        distances_matrix = distances.distance_array(ligand.positions, protein.positions)

        for residue in protein.residues:
            if np.any(distances_matrix[:, residue.atoms.indices] < DISTANCE_CUTOFF):
                residue_label = f"{residue.resname} {residue.resid}"
                interaction_data.append([frame, residue_label])

    # Convert to DataFrame
    interaction_df = pd.DataFrame(interaction_data, columns=["Frame", "Residue"])

    # Create a pivot table (Residue vs. Time)
    pivot_table = interaction_df.pivot_table(index="Residue", columns="Frame", aggfunc="size", fill_value=0)

    # Plot heatmap with improved visibility
    plt.figure(figsize=(12, 6))
    sns.heatmap(pivot_table, cmap="coolwarm", cbar=True, annot=False)
    plt.xlabel("Time (Frame)")
    plt.ylabel("Residue")
    plt.title(f"Ligand ({LIGAND_NAME})-Protein Interaction Map")
    plt.xticks(rotation=45)
    plt.yticks(rotation=0)
    plt.savefig(os.path.join(OUTPUT_DIR, "ligand_contact_map_with_residues.png"))
    plt.close()

    print("✅ Contact Map saved as ligand_contact_map_with_residues.png")
except Exception as e:
    print(f"❌ Error in Contact Map generation: {e}")






print("🔄 Running Ligand-Protein Contact Count analysis...")

try:
    u = mda.Universe("npt.gro", "md_0_1.xtc")
    protein = u.select_atoms("protein")
    ligand = u.select_atoms("resname PTC")  # Change 'PTC' to your ligand name
    DISTANCE_CUTOFF = 4.0  # Contact threshold in Å

    contact_counts = []

    for ts in u.trajectory[::25]:  # Process every 25th frame
        distances_matrix = distances.distance_array(ligand.positions, protein.positions)
        contacts = np.sum(distances_matrix < DISTANCE_CUTOFF)
        contact_counts.append(contacts)

    plt.figure(figsize=(8, 5))
    plt.plot(contact_counts, label="Ligand-Protein Contacts")
    plt.xlabel("Time (frames)")
    plt.ylabel("Number of Contacts")
    plt.title("Ligand-Protein Contacts Over Time")
    plt.legend()
    plt.savefig(os.path.join(OUTPUT_DIR, "ligand_contacts_plot_v2.png"))
    plt.close()

    print("✅ Ligand-Protein Contact Count analysis completed.")
except Exception as e:
    print(f"❌ Error in Contact Count analysis: {e}")


print(f"🔄 Generating Ligand-Protein Interaction Network with Residues...")

try:
    G = nx.Graph()

    for ts in u.trajectory:
        frame = ts.frame
        distances_matrix = distances.distance_array(ligand.positions, protein.positions)

        for residue in protein.residues:
            if np.any(distances_matrix[:, residue.atoms.indices] < DISTANCE_CUTOFF):
                residue_label = f"{residue.resname} {residue.resid}"
                G.add_edge(LIGAND_NAME, residue_label)  # Connect ligand to interacting residue

    # Plot network with labels
    plt.figure(figsize=(10, 8))
    pos = nx.spring_layout(G, seed=42)
    nx.draw(G, pos, with_labels=True, node_color="lightblue", edge_color="gray", node_size=1200, font_size=10)

    # Add edge labels for residue interactions
    # edge_labels = {(LIGAND_NAME, res): res for res in G.nodes if res != LIGAND_NAME}
    # nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=8, font_color="black")

    plt.title(f"Ligand-Protein Interaction Network ({LIGAND_NAME})")
    plt.savefig(os.path.join(OUTPUT_DIR, "ligand_interaction_network_with_residues.png"))
    plt.close()

    print("✅ Interaction Network saved as ligand_interaction_network_with_residues.png")
except Exception as e:
    print(f"❌ Error in Interaction Network generation: {e}")


print("🔄 Running Ligand RMSD analysis...")
try:
    ligand_rmsd = rms.RMSD(u, select=f"resname {LIGAND_NAME}", ref_frame=0)
    ligand_rmsd.run()

    plt.figure(figsize=(8, 5))
    plt.plot(ligand_rmsd.results.rmsd[:, 0], ligand_rmsd.results.rmsd[:, 2], label="Ligand RMSD")
    plt.xlabel("Time (Frames)")
    plt.ylabel("RMSD (Å)")
    plt.title(f"Ligand ({LIGAND_NAME}) RMSD Over Time")
    plt.legend()
    plt.savefig(os.path.join(OUTPUT_DIR,"ligand_rmsd_plot.png"))
    plt.close()

    print("✅ Ligand RMSD analysis completed.")
except Exception as e:
    print(f"❌ Error in Ligand RMSD analysis: {e}")




print("🔄 Running Ligand RMSF analysis...")

try:
    ligand_rmsf_analysis = RMSF(ligand).run()
    ligand_rmsf_values = ligand_rmsf_analysis.results.rmsf

    # Convert NaN values to zero
    ligand_rmsf_values = np.nan_to_num(ligand_rmsf_values, nan=0.0)

    # Ensure correct indexing
    atom_indices = np.arange(len(ligand_rmsf_values))

    plt.figure(figsize=(8, 5))
    plt.plot(atom_indices, ligand_rmsf_values, label="Ligand RMSF", marker='o', linestyle='-', color='b')

    plt.xlabel("Ligand Atom Index")
    plt.ylabel("RMSF (Å)")
    plt.title(f"Ligand ({LIGAND_NAME}) RMSF Per Atom")
    plt.legend()
    plt.savefig(os.path.join(OUTPUT_DIR, "ligand_rmsf_plot_fixed.png"))
    plt.close()

    print("✅ Ligand RMSF analysis completed.")
except Exception as e:
    print(f"❌ Error in Ligand RMSF analysis: {e}")


print("🔄 Computing Ligand Radius of Gyration (rGyr)...")

try:
    rg_values = []
    times = []

    for ts in u.trajectory:
        rg_values.append(ligand.radius_of_gyration())
        times.append(u.trajectory.time)

    plt.figure(figsize=(8, 5))
    plt.plot(times, rg_values, label="Ligand Radius of Gyration", color="purple")
    plt.xlabel("Time (ps)")
    plt.ylabel("rGyr (Å)")
    plt.title("Ligand Radius of Gyration Over Time")
    plt.legend()
    plt.savefig(os.path.join(OUTPUT_DIR, "ligand_rGyr_plot.png"))
    plt.close()

    print("✅ Ligand rGyr analysis completed.")
except Exception as e:
    print(f"❌ Error in Ligand rGyr analysis: {e}")





print("🔄 Running Hydrophobic vs. Hydrophilic Interaction analysis...")
try:
    # Ensure ligand name is clean
    LIGAND_NAME = LIGAND_NAME.strip()
    print(f"Using ligand name: {LIGAND_NAME}")

    # Define hydrophobic and hydrophilic selections
    hydrophobic_selection = protein.select_atoms(
        "(resname PHE or resname LEU or resname ILE or resname VAL or resname MET or resname TRP or resname TYR or resname ALA)"
    )
    hydrophilic_selection = protein.select_atoms(
        "(resname ASP or resname GLU or resname ASN or resname GLN or resname SER or resname THR or resname HIS or resname ARG or resname LYS)"
    )

    # Properly formatted select argument for Contacts
    ligand_selection = f"resname '{LIGAND_NAME}'"
    print(f"Selection query: {ligand_selection}")

    hydrophobic = Contacts(u, select=ligand_selection, refgroup=hydrophobic_selection, radius=4.0)
    hydrophilic = Contacts(u, select=ligand_selection, refgroup=hydrophilic_selection, radius=4.0)

    hydrophobic.run()
    hydrophilic.run()

    # Plot results
    plt.figure(figsize=(8, 5))
    plt.plot(hydrophobic.results.times, hydrophobic.results.timeseries, label="Hydrophobic Contacts", color="orange")
    plt.plot(hydrophilic.results.times, hydrophilic.results.timeseries, label="Hydrophilic Contacts", color="blue")
    plt.xlabel("Time (ps)")
    plt.ylabel("Number of Contacts")
    plt.title("Hydrophobic vs. Hydrophilic Interactions Over Time")
    plt.legend()
    plt.savefig(os.path.join(OUTPUT_DIR, "hydrophobic_vs_hydrophilic.png"))
    plt.close()

    print("✅ Hydrophobic vs. Hydrophilic Interaction analysis completed.")

except Exception as e:
    print(f"❌ Error in Hydrophobic vs. Hydrophilic analysis: {e}")







#### UPDATE THIS SCRIPT TO USE OUR NEW RESIDUE ADJUSTMENT TO ENSURE THAT WE CAPTURE OUR DIMERS####


print(f"🔄 Generating Ligand-Protein Interaction Network with Chain Information...")

try:
    G = nx.Graph()

    for ts in u.trajectory:
        frame = ts.frame
        distances_matrix = distances.distance_array(ligand.positions, protein.positions)

        for residue in protein.residues:
            if np.any(distances_matrix[:, residue.atoms.indices] < DISTANCE_CUTOFF):
                chain_id = assign_chain(residue.resid)
                residue_label = f"{chain_id}:{residue.resname} {residue.resid}"  # e.g., "A:ARG 45"

                G.add_edge(LIGAND_NAME, residue_label)  # Connect ligand to interacting residue

    # Plot network with chain-separated labels
    plt.figure(figsize=(10, 8))
    pos = nx.spring_layout(G, seed=42)
    nx.draw(G, pos, with_labels=True, node_color="lightblue", edge_color="gray", node_size=1200, font_size=10)

    # Add edge labels with chains
    # edge_labels = {(LIGAND_NAME, res): res for res in G.nodes if res != LIGAND_NAME}
    # nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=8, font_color="black")

    plt.title(f"Ligand-Protein Interaction Network ({LIGAND_NAME}) with Chain IDs")
    plt.savefig(os.path.join(OUTPUT_DIR, "ligand_interaction_network_with_chains.png"))
    plt.close()

    print("✅ Interaction Network saved as ligand_interaction_network_with_chains.png")
except Exception as e:
    print(f"❌ Error in Interaction Network generation: {e}")





# Load trajectory
u = mda.Universe("npt.gro", "md_0_1.xtc")  
protein = u.select_atoms("protein and name CA")  # Use Cα atoms

# Convert trajectory to NumPy array
positions = np.array([protein.positions for ts in u.trajectory])
positions -= np.mean(positions, axis=0)  # Centering

# Apply PCA
pca = PCA(n_components=2)
pca_coords = pca.fit_transform(positions.reshape(len(positions), -1))

# Plot PCA projection
plt.figure(figsize=(7, 5))
plt.scatter(pca_coords[:, 0], pca_coords[:, 1], c=range(len(pca_coords)), cmap='viridis')
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("PCA of Protein Motions")
plt.colorbar(label="Frame")
plt.savefig(os.path.join(OUTPUT_DIR, "pca_protein_motion.png"))
plt.show()



# print("🔄 Generating Normalized Protein-Ligand Interaction Chart...")

# # ---- CLASSIFY INTERACTIONS ----
# interaction_data = []

# for ts in u.trajectory:
#     frame = ts.frame
#     distances_matrix = mda.analysis.distances.distance_array(ligand.positions, protein.positions)

#     for residue in protein.residues:
#         if np.any(distances_matrix[:, residue.atoms.indices] < DISTANCE_CUTOFF):
#             # Classify interaction type (Schrödinger-style)
#             if residue.resname in ["PHE", "LEU", "ILE", "VAL", "MET", "TRP", "TYR", "ALA"]:  # Hydrophobic
#                 interaction_type = "Hydrophobic"
#             elif residue.resname in ["ASP", "GLU", "ARG", "LYS", "HIS"]:  # Ionic
#                 interaction_type = "Ionic"
#             elif residue.resname in ["SER", "THR", "ASN", "GLN"]:  # H-bonds
#                 interaction_type = "H-bonds"
#             else:
#                 interaction_type = "Water bridges"

#             interaction_data.append([frame, residue.resid, residue.resname, interaction_type])

# # Convert to DataFrame
# interaction_df = pd.DataFrame(interaction_data, columns=["Frame", "Residue", "ResidueName", "InteractionType"])

# # Aggregate interaction counts
# interaction_counts = interaction_df.groupby(["Residue", "ResidueName", "InteractionType"]).size().unstack(fill_value=0)

# # ---- NORMALIZE THE COUNTS ----
# # Normalize interactions to range 0-1 per residue
# interaction_fractions = interaction_counts.div(interaction_counts.sum(axis=1), axis=0).fillna(0)

# # ---- PLOT STACKED BAR CHART ----
# plt.figure(figsize=(14, 6))

# # Create stacked bar chart with normalized values
# interaction_fractions.plot(
#     kind="bar",
#     stacked=True,
#     colormap="Set2",
#     width=0.8,
#     figsize=(14, 6)
# )

# plt.xlabel("Residue")
# plt.ylabel("Interactions Fraction (0-1)")
# plt.title("Normalized Protein-Ligand Interactions")
# plt.legend(title="Interaction Type", bbox_to_anchor=(1.05, 1), loc="upper left")
# plt.xticks(rotation=45, ha="right")
# plt.tight_layout()

# # Save the figure
# plt.savefig(os.path.join(OUTPUT_DIR, "protein_ligand_interactions_normalized.png"))
# plt.close()

# print("✅ Protein-Ligand Interaction Graph saved as protein_ligand_interactions_normalized.png")


############## generate geometric datasets USING ONLY THE BINDING POCKET ATOMS PREVIOUSLY EXTRACTED AS THE INPUT


# ---- CONFIGURATION ----
INPUT_PDB = "binding_pocket_only.pdb"
INPUT_XTC = "binding_pocket_only.xtc"
LIGAND_RESNAME = "PTC"

# ---- FRAME CONTROL ----
START_FRAME = 0       # Set starting frame (None = first frame)
END_FRAME = None      # Set ending frame (None = last frame)
FRAME_STEP = 50       # Process every Nth frame

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ---- FUNCTION TO INFER ELEMENT TYPE ----
def infer_element(atom_name):
    match = re.match(r"([A-Za-z]+)", atom_name)
    if match:
        element = match.group(1).upper()
        if element in ["C", "O", "N", "H", "S", "P"]:
            return element
    return None

# ---- FUNCTION TO CALCULATE ANGLE ----
def calculate_angle(a, b, c):
    ba = a - b
    bc = c - b
    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    return np.degrees(np.arccos(np.clip(cosine_angle, -1.0, 1.0)))

# ---- LOADING THE TRAJECTORY ----
print("🔄 Loading Binding Pocket trajectory...")
u = mda.Universe(INPUT_PDB, INPUT_XTC)

# Define selections
ligand = u.select_atoms(f"resname {LIGAND_RESNAME}")
protein = u.select_atoms("protein")
waters = u.select_atoms("resname HOH")

# Determine frame range
total_frames = len(u.trajectory)
start_frame = START_FRAME if START_FRAME is not None else 0
end_frame = END_FRAME if END_FRAME is not None else total_frames

print(f"✅ Loaded {len(protein.residues)} protein residues and {len(ligand.atoms)} ligand atoms.")
print(f"🔍 Processing frames {start_frame} to {end_frame} with step size {FRAME_STEP}")

# ---- CLASSIFY INTERACTIONS ----
interaction_data = []
frame_count = 0

for ts in u.trajectory[start_frame:end_frame:FRAME_STEP]:
    frame = ts.frame
    frame_count += 1
    print(f"📍 Processing Frame {frame}...")

    distances_matrix = distances.distance_array(ligand.positions, protein.positions)

    for residue in protein.residues:
        min_dist = np.min(distances_matrix[:, residue.atoms.indices])
        interactions = []  # Allow multiple interactions per residue

        # ---- HYDROGEN BONDS ----
        for donor in residue.atoms:
            donor_element = infer_element(donor.name)
            if donor_element in ["O", "N"]:
                for acceptor in ligand.atoms:
                    acceptor_element = infer_element(acceptor.name)
                    if acceptor_element in ["O", "N"]:
                        dist = np.linalg.norm(donor.position - acceptor.position)
                        if dist <= 2.5:
                            donor_angle = calculate_angle(donor.position, acceptor.position, acceptor.position + np.array([0, 0, 1]))
                            acceptor_angle = calculate_angle(acceptor.position, donor.position, donor.position + np.array([0, 0, 1]))
                            if donor_angle >= 120 and acceptor_angle >= 90:
                                subtype = "Backbone Donor" if "N" in donor.name else "Side-chain Donor"
                                interactions.append(f"H-bonds ({subtype})")

        # ---- HYDROPHOBIC INTERACTIONS ----
        if min_dist <= 3.6:
            if residue.resname in ["PHE", "TYR", "TRP"]:
                interactions.append("Hydrophobic (π-π)")
            elif residue.resname in ["ARG", "LYS", "HIS"] and any(infer_element(a.name) == "C" for a in ligand.atoms):
                if min_dist <= 4.5:
                    interactions.append("Hydrophobic (π-Cation)")
            elif residue.resname in ["LEU", "ILE", "VAL", "MET", "ALA"]:
                interactions.append("Hydrophobic (General)")

        # ---- IONIC INTERACTIONS ----
        if min_dist <= 3.7:
            if residue.resname in ["ASP", "GLU"] and any(infer_element(a.name) == "N" for a in ligand.atoms):
                interactions.append("Ionic (Negative)")
            elif residue.resname in ["ARG", "LYS", "HIS"] and any(infer_element(a.name) == "O" for a in ligand.atoms):
                interactions.append("Ionic (Positive)")

        # ---- WATER BRIDGES ----
        if len(waters) > 0:
            for water in waters.residues:
                water_atoms = water.atoms
                water_dist_ligand = distances.distance_array(ligand.positions, water_atoms.positions)
                water_dist_protein = distances.distance_array(residue.atoms.positions, water_atoms.positions)

                if np.any(water_dist_ligand < 2.8) and np.any(water_dist_protein < 2.8):
                    interactions.append("Water Bridges")
                    break

        if interactions:
            for interaction in interactions:
                interaction_data.append([frame, residue.resid, residue.resname, interaction])

print(f"🔍 Finished processing {frame_count} frames.")

# ---- CONVERT TO DATAFRAME ----
if interaction_data:
    print(f"✅ Processed {len(interaction_data)} interactions. Converting to DataFrame...")
    interaction_df = pd.DataFrame(interaction_data, columns=["Frame", "Residue", "ResidueName", "InteractionType"])

    # Compute interaction counts per residue
    interaction_counts = interaction_df.groupby(["Residue", "ResidueName", "InteractionType"]).size().unstack(fill_value=0)

    # Normalize each interaction type individually (0-1 range)
    interaction_fractions = interaction_counts.div(frame_count).fillna(0).clip(upper=1)




    # ---- PLOT STACKED BAR CHART ----
    print("📊 Generating visualization...")
    plt.figure(figsize=(14, 6))
    interaction_fractions.plot(
        kind="bar",
        stacked=True,
        colormap="Set2",
        width=0.8,
        figsize=(14, 6)
    )

    plt.xlabel("Residue")
    plt.ylabel("Interaction Persistence (0-4 scale)")
    plt.title(f"Stacked Protein-Ligand Interactions with Water Bridges (Frames {start_frame}-{end_frame}, Step {FRAME_STEP})")
    plt.legend(title="Interaction Type", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.xticks(rotation=45, ha="right")
    plt.ylim(0, 4)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "protein_ligand_interactions_geomtry.png"))
    plt.show()
else:
    print("❌ No interactions detected! Check trajectory, selection criteria, or cutoff distances.")



######PLOTLY TRAJECTORY ANALSYSIS GRPAHS######
print("🔄 Initializing trajectory analysis for perfect_ligand_interaction_map_scaled")

################################
# Store contacts across trajectory
contact_frames = []
residue_contacts = {}

for ts in u.trajectory:
    contacts = []
    for res in protein.residues:
        res_atoms = res.atoms
        distances = mda.lib.distances.distance_array(ligand.positions, res_atoms.positions)

        if np.any(distances < CUTOFF):
            res_id = f"{res.resname}{res.resid}"
            contacts.append(res_id)

            if res_id not in residue_contacts:
                residue_contacts[res_id] = 0
            residue_contacts[res_id] += 1

    contact_frames.append(contacts)
print(f"✅ Collected contacts across {len(u.trajectory)} frames.")

# Compute contact frequencies
total_frames = len(u.trajectory)
contact_frequencies = {res: (count / total_frames) * 100 for res, count in residue_contacts.items()}
print(f"✅ Computed contact frequencies for {len(contact_frequencies)} residues.")

# Remove artifacts: Filter out residues with <1% interaction frequency
filtered_contacts = {res: freq for res, freq in contact_frequencies.items() if freq >= 1}
print(f"✅ Filtered out low-frequency contacts, keeping {len(filtered_contacts)} residues.")

# Residue color coding based on type
residue_colors = {
    "HIS": "blue", "ARG": "blue", "LYS": "blue",
    "ASP": "red", "GLU": "red",
    "SER": "green", "THR": "green", "ASN": "green", "GLN": "green",
    "CYS": "yellow", "SEC": "yellow",
    "GLY": "gray", "PRO": "gray",
    "ALA": "orange", "VAL": "orange", "ILE": "orange", "LEU": "orange", "MET": "orange",
    "PHE": "purple", "TYR": "purple", "TRP": "purple"
}

def get_residue_color(resname):
    return residue_colors.get(resname[:3], "black")

# Generate Matplotlib Graph with Correct Proportional Scaling
def generate_matplotlib_graph():
    G = nx.Graph()
    
    # Add edges
    for residue, freq in filtered_contacts.items():
        G.add_edge("PTC", residue, weight=freq)

    # **Ensure PTC is at the Center**
    pos = nx.spring_layout(G, seed=42, k=50, iterations=900)
    pos["PTC"] = np.array([0, 0])  

    # **Spread Out Nodes More**
    for node in pos:
        if node != "PTC":
            pos[node] *= 14  # Push nodes further away

    plt.figure(figsize=(18, 18))  # Increase figure size
    nx.draw_networkx_edges(G, pos, alpha=0.4, width=1, edge_color="gray")

    node_sizes, node_colors, labels = [], [], {}

    # **Find Min/Max Frequencies to Normalize Sizes**
    min_freq = min(filtered_contacts.values())
    max_freq = max(filtered_contacts.values())

    # **Find the Maximum Residue Node Size**
    max_residue_size = 0

    for node in G.nodes():
        if node != "PTC":
            resname = ''.join(filter(str.isalpha, node))
            freq = filtered_contacts.get(node, 0)
            
            # Normalize size between 800 and 3000
            node_size = 800 + (freq - min_freq) / (max_freq - min_freq) * 2200  
            max_residue_size = max(max_residue_size, node_size)

    # **PTC Node should be 3x the Largest Residue Node**
    ptc_size = int(3 * max_residue_size)

    for node in G.nodes():
        if node == "PTC":
            node_sizes.append(ptc_size)  # Largest node
            node_colors.append("pink")
        else:
            resname = ''.join(filter(str.isalpha, node))
            freq = filtered_contacts.get(node, 0)
            
            # Normalize size between 800 and 3000
            node_size = 800 + (freq - min_freq) / (max_freq - min_freq) * 2200  
            
            node_colors.append(get_residue_color(resname))
            node_sizes.append(node_size)  

        labels[node] = f"{node}\n({filtered_contacts.get(node, 0):.1f}%)"

    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, edgecolors="black")

    # **Labels Inside Nodes**
    for node, (x, y) in pos.items():
        plt.text(x, y, labels[node], fontsize=14, ha='center', va='center',
                 bbox=dict(facecolor="white", alpha=0.7, edgecolor='black', boxstyle="circle"))

    plt.title("Ligand-Protein Interaction Map (PTC) - Optimized Layout", fontsize=24)
    plt.axis("off")

    # Add legend
    legend_labels = {
        "Positively Charged (His, Arg, Lys)": "blue",
        "Negatively Charged (Asp, Glu)": "red",
        "Polar (Ser, Thr, Asn, Gln)": "green",
        "Sulfur-Containing (Cys, Sec)": "yellow",
        "Special (Gly, Pro)": "gray",
        "Hydrophobic (Ala, Val, Ile, Leu, Met)": "orange",
        "Aromatic (Phe, Tyr, Trp)": "purple"
    }
    handles = [plt.Line2D([0], [0], marker='o', color='w', markersize=12, markerfacecolor=color, label=label)
               for label, color in legend_labels.items()]
    plt.legend(handles=handles, loc='lower left', fontsize=14, title="Residue Types")

    return plt

# Generate Plotly Graph (Interactive HTML with Proper PTC Sizing)
def generate_plotly_graph():
    G = nx.Graph()
    
    for residue, freq in filtered_contacts.items():
        G.add_edge("PTC", residue, weight=freq)

    pos = nx.spring_layout(G, seed=42, k=40, iterations=800)
    pos["PTC"] = np.array([0, 0])  

    for node in pos:
        if node != "PTC":
            pos[node] *= 12  # Push nodes outward

    edge_trace = go.Scatter(
        x=[], y=[], line=dict(width=1, color='#888'), hoverinfo='none', mode='lines'
    )

    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_trace.x += (x0, x1, None)
        edge_trace.y += (y0, y1, None)

    node_trace = go.Scatter(
        x=[], y=[], text=[], mode='markers+text',
        marker=dict(size=[], color=[], opacity=0.9, line=dict(width=1, color='black')),
        textposition="middle center"
    )

    node_colors, node_sizes, node_texts, x_values, y_values = [], [], [], [], []

    for node in G.nodes():
        x, y = pos[node]
        x_values.append(x)
        y_values.append(y)

        if node == "PTC":
            node_colors.append("pink")
            node_sizes.append(int(2 * max(node_sizes, default=100)))  # Ensure PTC is 2x larger
        else:
            resname = ''.join(filter(str.isalpha, node))
            freq = min(filtered_contacts.get(node, 0), 200)
            node_colors.append(get_residue_color(resname))
            node_sizes.append(80 + (freq / 1.5))  # Make residue nodes larger

        node_texts.append(f"{node}\n({filtered_contacts.get(node, 0):.1f}%)")

    node_trace.x, node_trace.y, node_trace.text = x_values, y_values, node_texts
    node_trace.marker["color"], node_trace.marker["size"] = node_colors, node_sizes

    fig = go.Figure(data=[edge_trace, node_trace])
    fig.update_layout(showlegend=False, hovermode='closest', margin=dict(b=0, l=0, r=0, t=0),
                      title="Ligand-Protein Interaction Map (PTC) - Interactive", font=dict(size=16))

    return fig# Generate Matplotlib & Plotly Figures

plt_graph = generate_matplotlib_graph()
plt_graph.savefig(os.path.join(OUTPUT_DIR, "perfect_ligand_interaction_map_scaled.png"), dpi=300, bbox_inches="tight")
plt_graph.savefig(os.path.join(OUTPUT_DIR, "perfect_ligand_interaction_map_scaled.svg"), format="svg", bbox_inches="tight")


# ✅ Generate Plotly Interactive HTML
plotly_fig = generate_plotly_graph()  # Correctly use Plotly function
plotly_fig.write_html("perfect_ligand_interaction_map_scaled.html")  # Now this works!

print("✅ Analysis complete! Check the generated plots.")





# print("✅ MDAnalysis complete! Check the generated plots.")







# #################################################################
# #################################################################
# ############## Calpha Analysis of Warheads#######################
# #################################################################
# #################################################################





def compute_rmsf_whole_protein(topology, trajectory, pdb_reference, ligand_name, distance_cutoff, save_plot=True):
    """
    Compute RMSF for the whole protein complex while handling multiple chains.
    Uses PDB reference to correct chain numbering and highlights ligand-binding residues.
    Ensures each chain has a distinct X-axis range.
    """
    print("🔄 Running RMSF analysis (Optimized for multiple chains)...")

    try:
        # Load topology and trajectory
        u = mda.Universe(topology, trajectory)

        # Load PDB reference to get correct chain-resolved numbering
        u_pdb = mda.Universe(pdb_reference)

        # Select Cα atoms only (avoiding ligand and extra atoms)
        protein = u.select_atoms("protein and name CA")

        # Compute RMSF
        from MDAnalysis.analysis.rms import RMSF
        rmsf_analysis = RMSF(protein).run()
        rmsf_values = rmsf_analysis.results.rmsf

        # Extract residue numbers and chain identifiers
        protein_residues = u_pdb.select_atoms("protein").residues  # Only protein residues
        original_resids = protein_residues.resids
        chain_ids = protein_residues.segids

        print(f"🔍 Residue count after ligand exclusion: {len(original_resids)} residues")

        # Ensure residue count matches RMSF count
        if len(original_resids) != len(rmsf_values):
            print(f"⚠️ Warning: Residue count ({len(original_resids)}) does not match RMSF values ({len(rmsf_values)})")
            min_length = min(len(original_resids), len(rmsf_values))
            original_resids = original_resids[:min_length]
            rmsf_values = rmsf_values[:min_length]

        # Assign unique X-axis positions for each chain
        corrected_resids, corrected_chain_ids, chain_start_positions = adjust_chain_residue_numbers(original_resids, chain_ids)

        # Identify ligand-binding residues
        binding_residues = find_binding_site_residues(u_pdb, ligand_name, distance_cutoff)

        # --- PLOTTING ---
        plt.figure(figsize=(10, 5))

        # ✅ **Plot full RMSF curve first (light gray for background reference)**
        plt.plot(corrected_resids, rmsf_values, color="lightgray", alpha=0.5, label="Full RMSF Reference")

        # Assign unique colors to each chain
        unique_chains = np.unique(corrected_chain_ids)
        for i, chain in enumerate(unique_chains):
            chain_mask = (corrected_chain_ids == chain)
            chain_resids = corrected_resids[chain_mask]
            chain_rmsf = rmsf_values[chain_mask]

            if len(chain_resids) > 0:
                color = CHAIN_COLORS[i % len(CHAIN_COLORS)]  # Cycle through colors
                chain_label = CHAIN_LABELS.get(chain, f"Chain {chain}")  # Use custom label if available
                plt.plot(chain_resids, chain_rmsf, label=chain_label, color=color)
                plt.scatter(chain_resids, chain_rmsf, color=color, s=10)

        # ✅ **Ensure ligand-binding residues are plotted correctly per chain**
        for chain, residues in binding_residues.items():
            if chain in chain_start_positions:  # Make sure the chain exists in our corrected numbering
                for br in residues:
                    if br in original_resids:
                        # Find correct X position for the binding residue
                        br_x_position = chain_start_positions[chain] + (br - np.min(original_resids[chain_ids == chain]))
                        rmsf_value = rmsf_values[np.where(corrected_resids == br_x_position)][0]

                        plt.vlines(x=br_x_position, ymin=0, ymax=rmsf_value, color="green", linestyle="--", alpha=0.7, linewidth=1.5)

        plt.xlabel("Residue Number (Separated by Chain)")
        plt.ylabel("RMSF (Å)")
        plt.title("Residue Flexibility (RMSF) - Whole Protein (Cα Atoms)")
        plt.legend()

        if save_plot:
            output_path = os.path.join(OUTPUT_DIR, "rmsf_whole_protein_separated.png")
            plt.savefig(output_path)
            plt.close()
            print(f"📊 RMSF plot saved as {output_path}")

        print("✅ RMSF analysis completed.")
        return corrected_resids, rmsf_values

    except Exception as e:
        print(f"❌ Error in RMSF analysis: {e}")
        return None, None


def adjust_chain_residue_numbers(resids, chain_ids):
    """
    Adjust residue numbering to create a continuous sequence across multiple chains.
    Ensures that dimer chains (e.g., Chain B and Chain C) are kept separate.
    """
    unique_chains = np.unique(np.array(chain_ids, dtype=str))  # Convert chain IDs to string type
    corrected_resids = []
    corrected_chain_ids = []
    chain_start_positions = {}  # Track where each chain starts

    offset = 0  # Offset to ensure continuous numbering

    for chain in unique_chains:
        # Select residues belonging to this specific chain
        chain_mask = np.array(chain_ids, dtype=str) == chain  # Ensure comparison is done as strings
        chain_resids = np.array(resids)[chain_mask]  # Apply mask

        if len(chain_resids) == 0:
            continue

        min_resid = np.min(chain_resids)

        # Store where this chain starts in corrected numbering
        chain_start_positions[chain] = offset  

        # Apply offset and preserve correct numbering
        corrected_resids.extend(chain_resids - min_resid + offset)
        corrected_chain_ids.extend([chain] * len(chain_resids))

        # Offset for the next unique chain
        offset += len(chain_resids)

    return np.array(corrected_resids, dtype=int), np.array(corrected_chain_ids, dtype=str), chain_start_positions


def find_binding_site_residues(pdb_universe, ligand_name, distance_cutoff):
    """
    Identify protein residues within 'distance_cutoff' Å of ligand atoms,
    while tracking chain IDs for differentiation.
    """
    binding_residues = {}  # Dictionary to store residues as {chain: [resid1, resid2, ...]}
    ligand = pdb_universe.select_atoms(f"resname {ligand_name}")

    for residue in pdb_universe.select_atoms("protein and name CA").residues:
        chain_id = residue.segid  # Get chain ID
        for atom in residue.atoms:
            distances = np.linalg.norm(atom.position - ligand.atoms.positions, axis=1)
            if np.any(distances < distance_cutoff):  # If any atom is within cutoff
                if chain_id not in binding_residues:
                    binding_residues[chain_id] = []
                binding_residues[chain_id].append(residue.resid)  # Store residue under the correct chain

    return binding_residues

compute_rmsf_whole_protein(TOPOLOGY_FILE, TRAJECTORY_FILE, PDB_REFERENCE, LIGAND_NAME, DISTANCE_CUTOFF)














#################################################################
#################################################################
#################################################################
########Calpha Warhead and ligase Comparisons and Data Gen#######
#################################################################
#################################################################
#################################################################







# ---- INPUT FILES ----
TOPOLOGY_LIGASE = "ligase_ligand.gro"
TOPOLOGY_WARHEAD = "warhead_ligand.gro"
LIGASE_LIGAND_TRAJ = "ligase_ligand.xtc"
WARHEAD_LIGAND_TRAJ = "warhead_ligand.xtc"

LIGAND_NAME = "PTC"  # Adjust to match ligand name in the structure
DISTANCE_CUTOFF = 4.0  # Ligand-Protein interaction cutoff (Å)


def adjust_residue_numbers(resids):
    """
    Adjust residue numbers when a dimer structure is detected.
    This ensures that monomer 2 continues sequentially instead of resetting.
    """
    adjusted_resids = []
    offset = 0  # Offset to add when numbering resets

    for i, resid in enumerate(resids):
        if i > 0 and resid < resids[i - 1]:  # Detect residue reset
            offset = resids[i - 1]  # Set offset to the last residue before reset
        
        adjusted_resids.append(resid + offset)

    return np.array(adjusted_resids), offset  # Return adjusted numbers and offset

def find_binding_site_residues(topology, trajectory):
    """
    Identifies protein residues that interact with the ligand within the distance cutoff.
    Adjusts for dimer numbering.
    Returns a list of corrected residue IDs involved in ligand interactions.
    """
    print(f"🔎 Identifying ligand-binding residues in {topology}...")

    try:
        u = mda.Universe(topology, trajectory)

        # Select protein Cα atoms and ligand atoms
        calpha = u.select_atoms("protein and name CA")
        ligand = u.select_atoms(f"resname {LIGAND_NAME}")

        binding_residues = set()

        # Iterate through trajectory to find binding residues
        for ts in u.trajectory:
            distances = mda.lib.distances.distance_array(calpha.positions, ligand.positions)
            interacting_residues = np.where(distances < DISTANCE_CUTOFF)[0]
            binding_residues.update(calpha.resids[interacting_residues])

        # Adjust numbering if a dimer is detected
        corrected_resids, offset = adjust_residue_numbers(list(binding_residues))
        print(f"✅ Found {len(corrected_resids)} ligand-contacting residues.")

        return corrected_resids.tolist()

    except Exception as e:
        print(f"❌ Error identifying binding site residues: {e}")
        return []

def compute_rmsf(topology, trajectory, label, color, save_individual=True):
    """
    Compute Cα RMSF for a given topology and trajectory.
    Detects dimer residue numbering resets, and marks ligand-contact residues.
    Saves individual graphs if required.
    """
    print(f"🔄 Running Cα RMSF analysis for {label}...")

    try:
        # Load the universe
        u = mda.Universe(topology, trajectory)

        # Select Cα atoms
        calpha = u.select_atoms("protein and name CA")

        # Compute RMSF
        rmsf_analysis = RMSF(calpha).run()
        rmsf_values = rmsf_analysis.results.rmsf

        # Adjust residue numbering if a reset is detected
        corrected_resids, offset = adjust_residue_numbers(calpha.resids)

        # Identify ligand-binding residues
        binding_residues = find_binding_site_residues(topology, trajectory)

        # Create a new figure for this individual protein
        if save_individual:
            plt.figure(figsize=(8, 5))
            plt.plot(corrected_resids, rmsf_values, label=f"{label} Cα RMSF", color=color)

            # Ensure vertical bars do not exceed the RMSF peak values
            for br in binding_residues:
                if br in corrected_resids:
                    rmsf_value = rmsf_values[np.where(corrected_resids == br)][0]
                    plt.vlines(x=br, ymin=0, ymax=rmsf_value, color="green", linestyle="--", alpha=0.6)

            plt.xlabel("Residue Number")
            plt.ylabel("RMSF (Å)")
            plt.title(f"{label} Cα RMSF Over Time")
            plt.legend()

            # Save individual plot
            individual_path = os.path.join(OUTPUT_DIR, f"calpha_RMSF_{label.lower()}.png")
            plt.savefig(individual_path)
            plt.close()
            print(f"📊 {label} plot saved as {individual_path}")

        print(f"✅ Completed Cα RMSF analysis for {label}.")
        return corrected_resids, rmsf_values, binding_residues

    except Exception as e:
        print(f"❌ Error in Cα RMSF analysis for {label}: {e}")
        return None, None, []

# Create combined figure
plt.figure(figsize=(8, 5))

# Compute RMSF for both ligase and warhead (saving individual plots)
resids_ligase, rmsf_ligase, binding_ligase = compute_rmsf(TOPOLOGY_LIGASE, LIGASE_LIGAND_TRAJ, "Ligase", "blue", save_individual=True)
resids_warhead, rmsf_warhead, binding_warhead = compute_rmsf(TOPOLOGY_WARHEAD, WARHEAD_LIGAND_TRAJ, "Warhead", "red", save_individual=True)

# Add both datasets to combined plot
plt.plot(resids_ligase, rmsf_ligase, label="Ligase Cα RMSF", color="blue")
plt.plot(resids_warhead, rmsf_warhead, label="Warhead Cα RMSF", color="red")

# Ensure vertical bars do not exceed RMSF peak values
all_resids = np.concatenate([resids_ligase, resids_warhead])
all_rmsf_values = np.concatenate([rmsf_ligase, rmsf_warhead])

for br in binding_ligase + binding_warhead:
    if br in all_resids:
        rmsf_value = all_rmsf_values[np.where(all_resids == br)][0]
        plt.vlines(x=br, ymin=0, ymax=rmsf_value, color="green", linestyle="--", alpha=0.6)

# Final plot settings
plt.xlabel("Residue Number")
plt.ylabel("RMSF (Å)")
plt.title("Cα RMSF Over Time (Ligase & Warhead)")
plt.legend()

# Save combined plot
combined_path = os.path.join(OUTPUT_DIR, "calpha_RMSF_combined.png")
plt.savefig(combined_path)
plt.close()

print(f"📊 Combined plot saved as {combined_path}")


#################################################################
#################################################################
#################################################################
#################################################################
####################LIGASE and WARHEAD Q VALUE ANALYSIS##########
#################################################################
#################################################################
#################################################################
#################################################################
# ---- SETUP OUTPUT DIRECTORIES ----
OUTPUT_DIR = "Analysis_Graphs"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ---- INPUT FILES ----
TOPOLOGY_LIGASE = "ligase_ligand.gro"
TOPOLOGY_WARHEAD = "warhead_ligand.gro"
LIGASE_LIGAND_TRAJ = "ligase_ligand.xtc"
WARHEAD_LIGAND_TRAJ = "warhead_ligand.xtc"
LIGAND_NAME = "PTC"
DISTANCE_CUTOFF = 4.0  # Ligand-Protein interaction cutoff (Å)
CONTACT_CUTOFF = 0.5  # Q(X) contact cutoff (nm)
NUM_CLUSTERS = 3  # Number of clusters for Q(X)

# ---- FUNCTION TO COMPUTE LIGAND-PROTEIN CONTACTS ----
def compute_ligand_contacts(topology, trajectory, output_name):
    """Computes ligand-protein contacts over time with consistent colors."""
    print(f"🔄 Computing ligand-protein contacts for {output_name}...")

    # Assign colors explicitly
    color_map = {
        "Ligase_Ligand": "blue",
        "Warhead_Ligand": "red"
    }
    color = color_map.get(output_name, "black")  # Default to black if not found

    u = mda.Universe(topology, trajectory)
    protein = u.select_atoms("protein")
    ligand = u.select_atoms(f"resname {LIGAND_NAME}")

    if len(ligand) == 0:
        raise ValueError(f"❌ No ligand atoms found with resname {LIGAND_NAME}! Check the structure file.")

    contact_counts = []
    for ts in u.trajectory:
        distances = np.linalg.norm(ligand.positions[:, None, :] - protein.positions[None, :, :], axis=-1)
        contact_count = np.sum(distances < DISTANCE_CUTOFF)
        contact_counts.append(contact_count)

    np.savetxt(os.path.join(OUTPUT_DIR, f"{output_name}_contacts.csv"), contact_counts, delimiter=",")

    plt.figure(figsize=(8, 5))
    plt.plot(contact_counts, label=f"{output_name} Contacts", color=color)
    plt.xlabel("Frame")
    plt.ylabel("Number of Contacts")
    plt.title(f"Ligand-Protein Contacts Over Time: {output_name}")
    plt.legend()
    plt.savefig(os.path.join(OUTPUT_DIR, f"{output_name}_contacts.png"))
    plt.close()

    print(f"✅ Saved {output_name} contacts plot: {output_name}_contacts.png")

    return np.array(contact_counts)


# ---- FUNCTION TO COMPUTE Q(X) VALUES ----
def compute_q_values(traj, native_contacts, native):
    """Computes Q(X) native contact fraction."""
    if len(native_contacts) == 0:
        return np.zeros(len(traj))
    
    r_traj = md.compute_distances(traj, native_contacts)
    r_native = md.compute_distances(native, native_contacts)[0]
    return np.mean(1.0 / (1 + np.exp(50 * (r_traj - 1.8 * r_native))), axis=1)

# ---- FUNCTION TO COMPARE CONTACTS WITH Q(X) ----
def overlay_contacts_qx(ligase_contacts, warhead_contacts, q_ligase, q_warhead):
    """Plots ligand contacts and Q(X) values on the same graph."""
    print("📊 Overlaying Ligand Contacts and Q(X)...")

    plt.figure(figsize=(8, 5))
    plt.plot(ligase_contacts, label="Ligase Contacts", color="blue")
    plt.plot(q_ligase * max(ligase_contacts), label="Ligase Q(X)", linestyle="dashed", color="blue")
    plt.plot(warhead_contacts, label="Warhead Contacts", color="red")
    plt.plot(q_warhead * max(warhead_contacts), label="Warhead Q(X)", linestyle="dashed", color="red")
    plt.xlabel("Frame")
    plt.ylabel("Contacts & Q(X) (Scaled)")
    plt.title("Ligand Contacts vs. Q(X) Over Time")
    plt.legend()
    plt.savefig(os.path.join(OUTPUT_DIR, "contacts_vs_qx.png"))
    plt.close()

    print("✅ Saved Ligand Contacts vs. Q(X) comparison plot.")

# ---- FUNCTION TO CLUSTER LIGAND CONFORMATIONS BASED ON Q(X) ----
def cluster_ligand_qx(q_ligase, q_warhead):
    """Clusters ligand conformations into high-Q(X) and low-Q(X) states."""
    print("🔄 Clustering ligand conformations based on Q(X)...")

    q_matrix = np.vstack([q_ligase, q_warhead]).T
    linkage_matrix = linkage(q_matrix, method="ward")
    cluster_labels = fcluster(linkage_matrix, NUM_CLUSTERS, criterion='maxclust')

    plt.figure(figsize=(8, 5))
    for cluster in range(1, NUM_CLUSTERS + 1):
        cluster_indices = np.where(cluster_labels == cluster)[0]
        plt.scatter(cluster_indices, q_ligase[cluster_indices], label=f"Cluster {cluster}")

    plt.xlabel("Time (frames)")
    plt.ylabel("Fraction of Native Contacts (Q)")
    plt.legend()
    plt.title("Q(X) Clustering")
    plt.savefig(os.path.join(OUTPUT_DIR, "qx_clustering.png"))
    plt.close()

    print("✅ Clustering analysis completed.")

# ---- RUN CONTACT ANALYSIS ----
ligase_contacts = compute_ligand_contacts(TOPOLOGY_LIGASE, LIGASE_LIGAND_TRAJ, "Ligase_Ligand")
warhead_contacts = compute_ligand_contacts(TOPOLOGY_WARHEAD, WARHEAD_LIGAND_TRAJ, "Warhead_Ligand")

# ---- COMPUTE Q(X) ----
print("🔄 Computing Q(X) values...")

traj_ligase = md.load(LIGASE_LIGAND_TRAJ, top=TOPOLOGY_LIGASE)
traj_warhead = md.load(WARHEAD_LIGAND_TRAJ, top=TOPOLOGY_WARHEAD)
native = md.load(TOPOLOGY_LIGASE)

ligase_contacts_pairs = traj_ligase.topology.select_pairs(f"resname {LIGAND_NAME}", "protein")
warhead_contacts_pairs = traj_warhead.topology.select_pairs(f"resname {LIGAND_NAME}", "protein")

q_ligase = compute_q_values(traj_ligase, ligase_contacts_pairs, native)
q_warhead = compute_q_values(traj_warhead, warhead_contacts_pairs, native)

np.savetxt(os.path.join(OUTPUT_DIR, "q_ligase.txt"), q_ligase)
np.savetxt(os.path.join(OUTPUT_DIR, "q_warhead.txt"), q_warhead)
print("✅ Q(X) values computed.")

# ---- OVERLAY CONTACTS WITH Q(X) ----
overlay_contacts_qx(ligase_contacts, warhead_contacts, q_ligase, q_warhead)

# ---- CLUSTER LIGAND STATES BASED ON Q(X) ----
cluster_ligand_qx(q_ligase, q_warhead)










# --- INPUT FILES ---
TOPOLOGY_LIGASE = "ligase_ligand.gro"
TOPOLOGY_WARHEAD = "warhead_ligand.gro"
LIGASE_LIGAND_TRAJ = "ligase_ligand.xtc"
WARHEAD_LIGAND_TRAJ = "warhead_ligand.xtc"

# --- OUTPUT DIRECTORY ---
OUTPUT_DIR = "Analysis_Graphs"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Function to adjust residue numbering for dimers
def adjust_residue_numbers(resids):
    adjusted_resids = []
    offset = 0  # Offset to add when numbering resets

    for i, resid in enumerate(resids):
        if i > 0 and resid < resids[i - 1]:  # Detect residue reset
            offset = resids[i - 1]  # Set offset to the last residue before reset
        
        adjusted_resids.append(resid + offset)

    return np.array(adjusted_resids)

# Function to analyze secondary structure
def analyze_secondary_structure(topology_file, trajectory_file):
    print(f"Processing: {trajectory_file}")
    
    # Load trajectory
    traj = md.load(trajectory_file, top=topology_file)

    # Compute secondary structure per frame
    dssp_results = compute_dssp(traj)

    # DSSP mapping
    sse_mapping = {
        'H': 'Helix', 'G': 'Helix', 'I': 'Helix',
        'E': 'Strand', 'B': 'Strand',
        'T': 'Turn', 'S': 'Turn',
        '-': 'Coil'
    }

    # Map DSSP results to SSE types
    categorized_sse = np.vectorize(sse_mapping.get)(dssp_results)

    # Count occurrences for each residue
    residue_sse_counts = [Counter(categorized_sse[:, i]) for i in range(traj.n_residues)]

    # Compute fraction of time spent in each SSE
    residue_sse_fractions = [
        {key: count / traj.n_frames for key, count in counts.items()} for counts in residue_sse_counts
    ]

    # Extract SSE fractions
    helix_fractions = [fractions.get("Helix", 0) * 100 for fractions in residue_sse_fractions]
    strand_fractions = [fractions.get("Strand", 0) * 100 for fractions in residue_sse_fractions]
    turn_fractions = [fractions.get("Turn", 0) * 100 for fractions in residue_sse_fractions]
    coil_fractions = [fractions.get("Coil", 0) * 100 for fractions in residue_sse_fractions]

    # Compute total SSE fraction over time
    total_sse_fraction = np.mean((categorized_sse == "Helix") | (categorized_sse == "Strand"), axis=1) * 100

    return helix_fractions, strand_fractions, turn_fractions, coil_fractions, total_sse_fraction, categorized_sse, traj

# Run analysis for both Ligase and Warhead trajectories
ligase_helix, ligase_strand, ligase_turn, ligase_coil, ligase_total_sse, ligase_sse, ligase_traj = analyze_secondary_structure(TOPOLOGY_LIGASE, LIGASE_LIGAND_TRAJ)
warhead_helix, warhead_strand, warhead_turn, warhead_coil, warhead_total_sse, warhead_sse, warhead_traj = analyze_secondary_structure(TOPOLOGY_WARHEAD, WARHEAD_LIGAND_TRAJ)

# Adjust residue numbers for dimers
ligase_residue_indices = np.arange(len(ligase_helix))
warhead_residue_indices = np.arange(len(warhead_helix))

adjusted_ligase_residues = adjust_residue_numbers(ligase_residue_indices)
adjusted_warhead_residues = adjust_residue_numbers(warhead_residue_indices)

# --- PLOT 1: % SSE BY RESIDUE INDEX FOR LIGASE ---
plt.figure(figsize=(12, 6))
plt.fill_between(adjusted_ligase_residues, ligase_helix, color='red', alpha=0.6, label="Helix")
plt.fill_between(adjusted_ligase_residues, ligase_strand, color='blue', alpha=0.6, label="Strand")
plt.fill_between(adjusted_ligase_residues, ligase_turn, color='green', alpha=0.6, label="Turn")
plt.fill_between(adjusted_ligase_residues, ligase_coil, color='gray', alpha=0.6, label="Coil")

plt.xlabel("Residue Index")
plt.ylabel("Res. % SSE")
plt.title("Protein Secondary Structure - Ligase")
plt.legend()
plt.grid()
plt.savefig(os.path.join(OUTPUT_DIR, "sse_ligase.png"), dpi=300)
plt.show()

# --- PLOT 2: % SSE BY RESIDUE INDEX FOR WARHEAD ---
plt.figure(figsize=(12, 6))
plt.fill_between(adjusted_warhead_residues, warhead_helix, color='red', alpha=0.6, label="Helix")
plt.fill_between(adjusted_warhead_residues, warhead_strand, color='blue', alpha=0.6, label="Strand")
plt.fill_between(adjusted_warhead_residues, warhead_turn, color='green', alpha=0.6, label="Turn")
plt.fill_between(adjusted_warhead_residues, warhead_coil, color='gray', alpha=0.6, label="Coil")

plt.xlabel("Residue Index")
plt.ylabel("Res. % SSE")
plt.title("Protein Secondary Structure - Warhead")
plt.legend()
plt.grid()
plt.savefig(os.path.join(OUTPUT_DIR, "sse_warhead.png"), dpi=300)
plt.show()

# --- PLOT 3: SSE OVER TIME (Total SSE for Ligase & Warhead) ---
plt.figure(figsize=(12, 4))
plt.plot(ligase_traj.time, ligase_total_sse, label="Ligase SSE %", color="blue")
plt.plot(warhead_traj.time, warhead_total_sse, label="Warhead SSE %", color="red")

plt.xlabel("Time (nsec)")
plt.ylabel("% SSE")
plt.title("Protein Secondary Structure Over Time")
plt.legend()
plt.grid()
plt.savefig(os.path.join(OUTPUT_DIR, "sse_over_time.png"), dpi=300)
plt.show()

# --- PLOT 4: Ligase SSE Assignment Over Time ---
plt.figure(figsize=(12, 6))
plt.imshow(ligase_sse.T == "Helix", aspect='auto', interpolation='nearest', cmap="Reds", origin='lower')
plt.imshow(ligase_sse.T == "Strand", aspect='auto', interpolation='nearest', cmap="Blues", origin='lower', alpha=0.7)
plt.xlabel("Time (nsec)")
plt.ylabel("Residue Index")
plt.title("Ligase SSE Assignment Over Time")
plt.savefig(os.path.join(OUTPUT_DIR, "ligase_sse_time.png"), dpi=300)
plt.show()

# --- PLOT 5: Warhead SSE Assignment Over Time ---
plt.figure(figsize=(12, 6))
plt.imshow(warhead_sse.T == "Helix", aspect='auto', interpolation='nearest', cmap="Reds", origin='lower')
plt.imshow(warhead_sse.T == "Strand", aspect='auto', interpolation='nearest', cmap="Blues", origin='lower', alpha=0.7)
plt.xlabel("Time (nsec)")
plt.ylabel("Residue Index")
plt.title("Warhead SSE Assignment Over Time")
plt.savefig(os.path.join(OUTPUT_DIR, "warhead_sse_time.png"), dpi=300)
plt.show()

# --- PLOT 6: SSE Over Time (Ligase) ---
plt.figure(figsize=(12, 4))
plt.plot(ligase_traj.time, ligase_total_sse, label="Ligase SSE %", color="blue")
plt.xlabel("Time (nsec)")
plt.ylabel("% SSE")
plt.title("Protein Secondary Structure Over Time - Ligase")
plt.legend()
plt.grid()
plt.savefig(os.path.join(OUTPUT_DIR, "sse_over_time_ligase.png"), dpi=300)
plt.show()

# --- PLOT 7: SSE Over Time (Warhead) ---
plt.figure(figsize=(12, 4))
plt.plot(warhead_traj.time, warhead_total_sse, label="Warhead SSE %", color="red")
plt.xlabel("Time (nsec)")
plt.ylabel("% SSE")
plt.title("Protein Secondary Structure Over Time - Warhead")
plt.legend()
plt.grid()
plt.savefig(os.path.join(OUTPUT_DIR, "sse_over_time_warhead.png"), dpi=300)
plt.show()





################################################################################################################################
##################   LIGAND GRAPH GENERATION SPECIFICALLY AROUND LIGAND TRAJECTORY IN THE PRESENCE OF WARHEAD & LIGASE##########
################################################################################################################################
################################################################################################################################
################################################################################




# === Configuration ===
BOUND_LIGAND_XTC = "ligand.xtc"  # Bound ligand trajectory
REFERENCE_GRO = "ligand.pdb"  # Reference structure
OUTPUT_DIR = "Analysis_Graphs"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# === Step 1: Extract Ligand Snapshots ===
def extract_ligand_snapshots(trajectory):
    """Extracts ligand snapshots from a trajectory and saves them as PDB files."""
    print(f"🔄 Extracting ligand snapshots from {trajectory}...")

    output_path = os.path.join(OUTPUT_DIR, "bound")
    os.makedirs(output_path, exist_ok=True)

    u = mda.Universe(REFERENCE_GRO, trajectory)
    ligand = u.select_atoms("resname PTC")  # Adjust ligand resname

    if len(ligand) == 0:
        raise ValueError(f"❌ No ligand atoms found in {trajectory}! Check resname.")

    for i, ts in enumerate(u.trajectory):
        pdb_filename = os.path.join(output_path, f"ligand_frame_{i:04d}.pdb")
        ligand.write(pdb_filename)

    print(f"✅ Extracted {i+1} ligand snapshots to {output_path}/")
    return output_path

# === Step 2: Convert PDB to MOL2 Using Open Babel ===
def convert_pdb_to_mol2(pdb_dir):
    """Converts PDB files to MOL2 format using Open Babel."""
    print(f"🔄 Converting PDB snapshots to MOL2 format in {pdb_dir}...")

    pdb_files = sorted(glob.glob(os.path.join(pdb_dir, "ligand_frame_*.pdb")))
    if not pdb_files:
        raise ValueError(f"❌ No PDB snapshots found in {pdb_dir}!")

    for pdb_file in pdb_files:
        mol2_file = pdb_file.replace(".pdb", ".mol2")
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("pdb", "mol2")

        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, pdb_file)
        obConversion.WriteFile(mol, mol2_file)

    print(f"✅ Converted all PDB files to MOL2 format in {pdb_dir}.")

# === Step 3: Compute Torsion Angles ===
def compute_torsion_angles(mol2_dir):
    """Computes torsion angles from MOL2 files using RDKit."""
    print(f"🔄 Computing torsion angles from MOL2 files in {mol2_dir}...")

    mol2_files = sorted(glob.glob(os.path.join(mol2_dir, "ligand_frame_*.mol2")))
    torsion_data = []

    if not mol2_files:
        raise ValueError(f"❌ No MOL2 files found in {mol2_dir}!")

    for mol2_file in mol2_files:
        mol = Chem.MolFromMol2File(mol2_file, removeHs=False)

        if mol is None:
            print(f"⚠️ Skipping {mol2_file}, RDKit failed to parse it.")
            continue

        conf = mol.GetConformer()
        torsions = [list(tup) for tup in mol.GetSubstructMatches(
            Chem.MolFromSmarts("[!R]~[!R]~[!R]~[!R]")
        )]

        frame_torsions = []
        for (a, b, c, d) in torsions:
            angle = rdMolTransforms.GetDihedralDeg(conf, a, b, c, d)
            frame_torsions.append(angle)

        torsion_data.append(frame_torsions)

    torsion_data = np.array(torsion_data)
    np.savetxt(os.path.join(mol2_dir, "torsion_angles.csv"), torsion_data, delimiter=",")

    print(f"✅ Computed torsion angles for {len(torsion_data)} frames.")
    return torsion_data

# === Step 4: Perform PCA on Torsion Angles ===
def perform_pca(torsion_data):
    """Performs PCA on torsion angles to identify dominant motion."""
    print("🔄 Performing PCA on torsion angles...")

    pca = PCA(n_components=2)
    reduced_data = pca.fit_transform(torsion_data)

    plt.figure(figsize=(8, 6))
    plt.scatter(reduced_data[:, 0], reduced_data[:, 1], alpha=0.5, color="blue")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("PCA of Ligand Torsion Angles")
    plt.savefig(os.path.join(OUTPUT_DIR, "torsion_pca.png"))
    plt.close()

    print(f"✅ PCA plot saved as torsion_pca.png")
    return reduced_data

# === Step 5: Perform RMSD Clustering ===
def cluster_conformations(torsion_data):
    """Clusters ligand conformations based on RMSD-like torsion similarity."""
    print("🔄 Performing RMSD-based clustering on torsion angles...")

    # Hierarchical clustering
    linkage_matrix = linkage(torsion_data, method="ward")
    num_clusters = 3  # Adjust as needed
    cluster_labels = fcluster(linkage_matrix, num_clusters, criterion="maxclust")

    plt.figure(figsize=(8, 6))
    sns.scatterplot(x=torsion_data[:, 0], y=torsion_data[:, 1], hue=cluster_labels, palette="viridis")
    plt.xlabel("Torsion 1 (°)")
    plt.ylabel("Torsion 2 (°)")
    plt.title("Ligand Conformation Clusters")
    plt.legend(title="Cluster")
    plt.savefig(os.path.join(OUTPUT_DIR, "torsion_clusters.png"))
    plt.close()

    print(f"✅ Clustering plot saved as torsion_clusters.png")

    return cluster_labels


def plot_torsion_vs_time(torsion_data):
    """Plots torsion angles as a function of time."""
    print(f"📊 Plotting torsion angles over time...")

    plt.figure(figsize=(10, 5))
    for i in range(torsion_data.shape[1]):  # Iterate over torsions
        plt.plot(range(len(torsion_data)), torsion_data[:, i], label=f"Torsion {i+1}")

    plt.xlabel("Frame")
    plt.ylabel("Torsion Angle (°)")
    plt.title("Ligand Torsion Angles Over Time")
    plt.legend()
    plt.savefig(os.path.join(OUTPUT_DIR, "torsion_vs_time.png"))
    plt.close()

    print(f"✅ Saved torsion vs. time plot: torsion_vs_time.png")


# === Run the Full Workflow ===
if __name__ == "__main__":
    bound_pdb_dir = extract_ligand_snapshots(BOUND_LIGAND_XTC)
    convert_pdb_to_mol2(bound_pdb_dir)
    bound_torsion_data = compute_torsion_angles(bound_pdb_dir)

    # Perform PCA and clustering
    pca_data = perform_pca(bound_torsion_data)
    cluster_labels = cluster_conformations(bound_torsion_data)

    print("🎉 Complete! Ligand clustering and PCA analysis done.")






print("✅ Analysis complete! Check the generated plots.")

################################################################################################
################################################################################################
################      P D F  G E N E R A T I O N   P I E C E   O F   S C R I P T $##############



import os
import glob
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rcParams, font_manager as fm
from PIL import Image, ImageDraw, ImageFont
from io import BytesIO

# Paths
OUTPUT_DIR = os.path.join(os.getcwd(), "Analysis_Graphs")
final_pdf_path = os.path.join(OUTPUT_DIR, "Final_Analysis_Report.pdf")
notes_file_path = os.path.join(os.getcwd(), "4_GraphNotes.txt")
font_path = os.path.join(os.getcwd(), "EBGaramond-VariableFont_wght.ttf")  # Optional

# Load custom Garamond font if available
if os.path.exists(font_path):
    garamond_font = fm.FontProperties(fname=font_path)
    rcParams["font.family"] = garamond_font.get_name()
    print(f"✅ Using font: {garamond_font.get_name()}")
else:
    garamond_font = None
    print("⚠️ Garamond font not found, using default font.")

# Parse GraphNotes.txt
def parse_graph_notes(file_path):
    notes_dict = {}
    current_image = None
    current_description = []

    with open(file_path, "r", encoding="utf-8") as file:
        lines = file.readlines()

    for line in lines:
        line = line.strip()
        if line.endswith(".png"):
            if current_image:
                notes_dict[current_image] = "\n".join(current_description).strip()
            current_image = line  # Ensure correct filename mapping
            current_description = []
        elif current_image:
            current_description.append(line)

    if current_image:
        notes_dict[current_image] = "\n".join(current_description).strip()

    return notes_dict

graph_notes = parse_graph_notes(notes_file_path)

# Function to render text as an image with adaptive font sizing and preserving new lines
def text_to_image(text, width=800, height=1000, max_font_size=24, min_font_size=10):
    """Creates an image from text, adjusting font size while preserving line breaks."""
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    font_size = max_font_size
    
    if garamond_font:
        font = ImageFont.truetype(font_path, font_size)
    else:
        font = ImageFont.load_default()
    
    while font_size >= min_font_size:
        if garamond_font:
            font = ImageFont.truetype(font_path, font_size)
        text_lines = text.split("\n")  # Preserve manual new lines
        wrapped_lines = []
        
        for line in text_lines:
            words = line.split()
            current_line = ""
            for word in words:
                test_line = f"{current_line} {word}".strip()
                text_width, _ = draw.textbbox((0, 0), test_line, font=font)[2:]
                if text_width < width - 40:
                    current_line = test_line
                else:
                    wrapped_lines.append(current_line)
                    current_line = word
            wrapped_lines.append(current_line)
        
        wrapped_text = "\n".join(wrapped_lines)
        text_bbox = draw.textbbox((0, 0), wrapped_text, font=font)
        text_height = text_bbox[3] - text_bbox[1]
        
        if text_height <= height * 0.9:
            break
        font_size -= 2  # Reduce font size and retry
    
    draw.multiline_text((20, 20), wrapped_text, fill="black", font=font, spacing=4)
    return img

# Read image descriptions
graph_notes = parse_graph_notes(notes_file_path)
image_files = sorted(glob.glob(os.path.join(OUTPUT_DIR, "*.png")))

if not image_files:
    print("⚠️ No images found.")
else:
    with PdfPages(final_pdf_path) as pdf:
        for image in image_files:
            # Add image
            fig, ax = plt.subplots(figsize=(8.5, 11))
            img = plt.imread(image)
            ax.imshow(img)
            ax.axis("off")
            pdf.savefig(fig)
            plt.close(fig)

            # Convert description to image with adaptive font size and wrapped text
            image_filename = os.path.basename(image).strip()
            description = graph_notes.get(image_filename, "No description available.")
            text_img = text_to_image(description)

            fig, ax = plt.subplots(figsize=(8.5, 11))
            ax.imshow(text_img)
            ax.axis("off")
            pdf.savefig(fig)
            plt.close(fig)

print(f"✅ Final PDF report generated at {final_pdf_path}")
