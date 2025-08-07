"""
Takes a PDB file with protein and water coordinates.
Creates a graph with vertices of water/AA residues and edges if they have an H bond.
Committed on github under protein_ligand_informatics/water/water_protein_Hbond_network.py

"""

from Bio import PDB

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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

"""
Data loading
"""
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

"""
Graph construction
"""
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
        Hbond_graph.add_node(f"w{idx}", pos=coord, res="water")
    for idx, coord in enumerate(donor_coordinates):
        Hbond_graph.add_node(f"d{idx}", pos=coord, res="donor")
    for idx, coord in enumerate(acceptor_coordinates):
        Hbond_graph.add_node(f"a{idx}", pos=coord, res="acceptor")

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
    
"""
Graph analysis
"""
def average_degree(G):
    """ 
    Takes a graph.
    Returns a float.
    """
    return len(G.edges()) / ( 2 * len(G.nodes()) )
    
    
def degree_distribution(G):
    """
    Takes a graph.
    Returns a tuple of ndarrays. The first contains the degree values and the second the corresponding proportion of nodes with that degree.
    """
    degree_sequence = [d for n, d in G.degree()]
    degree, counts = np.unique(degree_sequence, return_counts=True)
    counts_proportion = counts/np.sum(counts)
    return degree, counts_proportion
    
def bond_proportion_by_type(Hbond_graph):
    """
    Takes an H-bond graph.
    Outputs a 3-tuple of the proportion of water-water, water-donor/acceptor, and donor-acceptor H bonds.
    """
    ww = 0
    wd = 0
    da = 0
    for u, v in Hbond_graph.edges():
        if Hbond_graph.nodes()[u]["res"] == "water" and Hbond_graph.nodes()[v]["res"] == "water":
            ww += 1
        elif Hbond_graph.nodes()[u]["res"] == "water" or Hbond_graph.nodes()[v]["res"] == "water":
            wd += 1
        else:
            da += 1
            
    total_bonds = len( Hbond_graph.edges() )
    return ww/total_bonds, wd/total_bonds, da/total_bonds

def get_graph_metrics(G_list, names_list):
    """
    Takes a list of graphs and the associated protein names.
    Outputs a dataframe with metrics.
    """
    df = pd.DataFrame(index=names_list, columns = ["Average degree", "Degree distribution", "Bond proportion by type"])
    for graph, prot in zip(G_list, names_list):
        df.at[prot, "Degree distribution"] = [degree_distribution(graph)]
        df.loc[prot, "Average degree"] = average_degree(graph)
        df.at[prot, "Bond proportion by type"] = bond_proportion_by_type(graph)

    return df


"""
Visualization
"""
# utility function to get colors from graph
def get_colors(Hbond_graph):
    type_to_color = {"water": "blue", "donor": "red", "acceptor":"tomato"}
    node_colors = [ type_to_color.get(nx.get_node_attributes(Hbond_graph, "res")[n], "gray")
        for n in Hbond_graph.nodes() ]
    edge_colors = [Hbond_graph[u][v]["color"] for u, v in Hbond_graph.edges()]
    return node_colors, edge_colors

# Set the characteristics of the graph and draw projected onto XY plane
def draw_graph(Hbond_graph, ax):
    pos = nx.get_node_attributes(Hbond_graph, 'pos')
    node_colors, edge_colors = get_colors(Hbond_graph)
    nx.draw(Hbond_graph, pos, ax=ax, node_color=node_colors, edge_color=edge_colors, node_size=3)

def draw_graph_3d(Hbond_graph, ax):
    """
    Plotting function for 3D graph. ax must be with a 3-D projection.
    """
    pos = nx.get_node_attributes(Hbond_graph, 'pos')
    node_colors, edge_colors = get_colors(Hbond_graph)
    
    # Plot nodes
    # List comprehension gives a list of (x,y,z) coordinates and zip(*([])) transposes them to the format for matplotlib
    xs, ys, zs = zip( *[pos[node] for node in Hbond_graph.nodes()] ) 
    ax.scatter(xs, ys, zs, s=5, c=node_colors)

    # Plot edges
    for idx, (u,v) in enumerate(Hbond_graph.edges()):
        x, y, z = zip(pos[u],pos[v])
        ax.plot(x, y, z, color=edge_colors[idx])



""" 
Load all files ending with "_final.pdb" from data_folder.
Aggregate coordinate data into pandas dataframe.
"""

data_folder = "../PDB Data/cdk2/"
pdb_path = pathlib.Path(data_folder)
pdb_files = list(pdb_path.glob('*_final.pdb')) 
parser = PDB.PDBParser(QUIET=True)

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

Hbond_graphs = []
for prot in df_atoms.index:
    graph = graph_from_atom_sets(df_atoms.loc[prot, "Water"], df_atoms.loc[prot, "H-bond donors"], df_atoms.loc[prot, "H-bond acceptors"])
    Hbond_graphs.append(graph)
    
df_metrics = get_graph_metrics(Hbond_graphs, df_atoms.index)

fig  = plt.figure( figsize=(14, 6) )
ax1 = fig.add_subplot(121, projection='3d')
ax2 = fig.add_subplot(122, projection='3d')

#plot graphs
draw_graph_3d(Hbond_graphs[0], ax1)
draw_graph_3d(Hbond_graphs[1], ax2)

#plot protein backbones as one chain
apo_coordinates = get_coordinates_from_atoms(df_atoms.iloc[0]["Backbone"])
ax1.plot(apo_coordinates[:,0], apo_coordinates[:,1], apo_coordinates[:,2], color="green", alpha=0.2)
ax1.set_title("apo")

holo_coordinates = get_coordinates_from_atoms(df_atoms.iloc[1]["Backbone"])
ax2.plot(holo_coordinates[:,0], holo_coordinates[:,1], holo_coordinates[:,2], color="green", alpha=0.2)
ax2.set_title("holo")

# make a legend for the plot
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='Water', markerfacecolor='blue', markersize=10),
    Line2D([0], [0], marker='o', color='w', label='Donor', markerfacecolor='red', markersize=10),
    Line2D([0], [0], marker='o', color='w', label='Acceptor', markerfacecolor='tomato', markersize=10),
    Line2D([0], [1], linewidth=1, linestyle='-', color='black', label='Water to Water H-bond'),
    Line2D([0], [1], linewidth=1, linestyle='-', color='violet', label='Water to Donor/Acceptor H-bond'),
    Line2D([0], [1], linewidth=1, linestyle='-', color='darkgreen', label='Donor to Acceptor H-bond'),
    Line2D([0], [1], linewidth=1, linestyle='-', color='green', label='Protein backbone', alpha=0.2)
]

ax2.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.2, 1), title="Graph elements")
plt.tight_layout()

# plt.savefig("../Figures/H_bond_graph_test.png", dpi=144)
plt.show()