"""
Takes a PDB file with protein and water coordinates.
Creates a graph with vertices of water/AA residues and edges if they have an H bond.
"""

from Bio import PDB

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx

import glob
import pathlib

# CONSTANTS
AA=['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
WATER=['HOH', 'WAT', 'H2O']  
HBOND_MAX_DIST=3.5
HBOND_MIN_DIST=2.4

# Default hydrogen bond donor and acceptor definitions
# From Mingbin
DEFAULT_DONOR_ATOMS = ['N', 'ND1','ND2', 'NE', 'NE1', 'NE2', 'NH1', 'NH2', 'NZ', 'OG', 'OG1', 'OH', 'SG']
DEFAULT_ACCEPTOR_ATOMS = ['ND1', 'NE2', 'O', 'OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1', 'OH', 'SD', 'SG']

def get_filtered_atoms(pdb_structure, residue_target=None, atom_target=None):
    """
    residue_target [optional] and atom_target [optional] are lists of acceptable values
    """
    atoms = []
    if residue_target:
        residue_filter = lambda x: x.resname in residue_target
    else:
        residue_filter = lambda x: True # always return true to not filter anything out
        
    if atom_target:
        atom_filter = lambda x: x.name in atom_target
    else:
        atom_filter = lambda x: True # always return true to not filter anything out
    
    for model in pdb_structure:
        for chain in model:
            for res in filter(residue_filter, chain): # keeps only residues in the residue_target list, or all if input is None
                new_atoms = list( filter(atom_filter, res) ) # keeps only atoms in the atom_target list, or all if input is None
                atoms += new_atoms
                
    return atoms

def get_coordinates_from_atoms(atoms):
    """
    Returns the coordinates of a list of N atoms in an N x 3 ndarray
    """
    return np.array( [atom.get_coord() for atom in atoms] )

def graph_from_atom_sets(water_atoms, donor_atoms, acceptor_atoms):
    """
    Takes in coordinates for water, H-bond donors, and H-bond acceptors.
    Outputs a graph with all H bonds between waters and from water to donors or acceptors.
    TODO between donor and acceptor in the protein? would represent secondary structure 
    """
    # Create an empty graph
    Hbond_graph = nx.Graph()

    # get coordinates
    water_coordinates = get_coordinates_from_atoms(water_atoms)    
    donor_coordinates = get_coordinates_from_atoms(donor_atoms)    
    acceptor_coordinates = get_coordinates_from_atoms(acceptor_atoms)    
    
    # Add nodes to the graph with position as a node attribute
    for idx, coord in enumerate(water_coordinates):
        Hbond_graph.add_node(f"w{idx}", pos=coord[:2], res="water")
    for idx, coord in enumerate(donor_coordinates):
        Hbond_graph.add_node(f"d{idx}", pos=coord[:2], res="donor")
    for idx, coord in enumerate(acceptor_coordinates):
        Hbond_graph.add_node(f"a{idx}", pos=coord[:2], res="acceptor")

    # Add H-bond edges 
    for idx, coord in enumerate(water_coordinates): 
        water_water_distances = np.linalg.norm(coord - water_coordinates, axis=1) # water to water
        connected_edges = np.nonzero( 
            np.logical_and(water_water_distances < HBOND_MAX_DIST, water_water_distances > HBOND_MIN_DIST) 
            )[0]
        Hbond_graph.add_edges_from(
            [ (f"w{idx}", f"w{jdx}") for jdx in connected_edges ],
            color="black" )
        
        water_donor_distances = np.linalg.norm(coord - donor_coordinates, axis=1) # water to donor 
        connected_edges = np.nonzero( 
            np.logical_and(water_donor_distances < HBOND_MAX_DIST, water_donor_distances > HBOND_MIN_DIST)
            )[0]
        Hbond_graph.add_edges_from( 
            [ (f"w{idx}", f"d{jdx}") for jdx in connected_edges ],
            color="violet" )
        
        water_acceptor_distances = np.linalg.norm(coord - acceptor_coordinates, axis=1) # water to acceptor
        connected_edges = np.nonzero( 
            np.logical_and(water_acceptor_distances < HBOND_MAX_DIST, water_acceptor_distances > HBOND_MIN_DIST)
            )[0]
        Hbond_graph.add_edges_from( 
            [ (f"w{idx}", f"a{jdx}") for jdx in connected_edges ],
            color="violet" )
        
    for idx, donor in enumerate(donor_coordinates): # donor to acceptor
        donor_acceptor_distances = np.linalg.norm(donor - acceptor_coordinates, axis=1) # water to donor 
        connected_edges = np.nonzero( 
            np.logical_and(donor_acceptor_distances < HBOND_MAX_DIST, donor_acceptor_distances > HBOND_MIN_DIST)
            )[0]
        Hbond_graph.add_edges_from( 
            [ (f"d{idx}", f"a{jdx}") for jdx in connected_edges 
            if donor_atoms[idx].get_parent() != acceptor_atoms[jdx].get_parent() ],
            color="darkgreen") # do not count as H bond both acceptor and donor are in the same residue
        
    return Hbond_graph
    
    
# Set the characteristics of the graph and draw projected onto XY plane
def draw_graph(Hbond_graph, ax):
    pos = nx.get_node_attributes(Hbond_graph, 'pos')
    type_to_color = {"water": "blue", "donor": "red", "acceptor":"tomato"}
    node_colors = [type_to_color.get(nx.get_node_attributes(Hbond_graph, 'res')[n], "gray") for n in Hbond_graph.nodes()]
    edge_colors = [Hbond_graph[u][v]["color"] for u, v in Hbond_graph.edges()]

    nx.draw(Hbond_graph, pos, ax=ax, node_color=node_colors, edge_color=edge_colors, node_size=3)


""" 
Load all files ending with "_final.pdb" from data_folder.
Aggregate coordinate data into pandas dataframe.
"""

data_folder = "../PDB Data/cdk2/"
pdb_path = pathlib.Path(data_folder)
pdb_files = list(pdb_path.glob('*_final.pdb')) 
parser = PDB.PDBParser(QUIET=True)

apo_id = '1PW2'
holo_id = '1PXI'
df_atoms = pd.DataFrame(columns = ["Backbone", "Water", "H-bond donors", "H-bond acceptors"])
for file_name in pdb_files:
    prot = file_name.stem.split('_')[0]
    structure = parser.get_structure(prot, file_name)
    df_atoms.loc[prot, "Backbone"] = get_filtered_atoms(structure, residue_target = AA, atom_target = ["CA"])
    df_atoms.loc[prot, "Water"] = get_filtered_atoms(structure, residue_target = WATER, atom_target = ["O"])
    df_atoms.loc[prot, "H-bond acceptors"] = get_filtered_atoms(structure, residue_target = AA, atom_target = DEFAULT_ACCEPTOR_ATOMS)
    df_atoms.loc[prot, "H-bond donors"] = get_filtered_atoms(structure, residue_target = AA, atom_target = DEFAULT_DONOR_ATOMS)

"""
Display hydrogen bonding network as a graph.
"""

apo_Hbond_graph = graph_from_atom_sets(df_atoms.loc[apo_id, "Water"], df_atoms.loc[apo_id, "H-bond donors"], df_atoms.loc[apo_id, "H-bond acceptors"])
holo_Hbond_graph = graph_from_atom_sets(df_atoms.loc[holo_id, "Water"], df_atoms.loc[holo_id, "H-bond donors"], df_atoms.loc[holo_id, "H-bond acceptors"])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24,12))

#plot graphs
draw_graph(apo_Hbond_graph, ax1)
draw_graph(holo_Hbond_graph, ax2)

#plot protein backbones as one chain
apo_coordinates = get_coordinates_from_atoms(df_atoms.loc[apo_id, "Backbone"])
ax1.plot(apo_coordinates[:,0], apo_coordinates[:,1], color="green", alpha=0.2)
ax1.set_title("apo")

holo_coordinates = get_coordinates_from_atoms(df_atoms.loc[holo_id, "Backbone"])
ax2.plot(holo_coordinates[:,0], holo_coordinates[:,1], color="green", alpha=0.2)
ax2.set_title("holo")

plt.savefig("../Figures/H_bond_graph_test.png", dpi=144)
plt.show()
