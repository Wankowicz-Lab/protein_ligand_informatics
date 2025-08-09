"""
Takes a PDB file with protein and water coordinates.
Creates a graph with vertices of water/AA residues and edges if they have an H bond.
Committed on github under protein_ligand_informatics/water/water_protein_Hbond_network.py

"""

from Bio import PDB

import numpy as np
import pandas as pd
import networkx as nx

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from scipy.spatial.transform import Rotation

import glob
import pathlib
import sys


# UNIVERSAL CONSTANTS
AA=['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
WATER=['HOH', 'WAT', 'H2O']  

HBOND_MAX_DIST=3.5
HBOND_MIN_DIST=2.4

# Default hydrogen bond donor and acceptor definitions
# From Mingbin
DEFAULT_DONOR_ATOMS = ['N', 'ND1','ND2', 'NE', 'NE1', 'NE2', 'NH1', 'NH2', 'NZ', 'OG', 'OG1', 'OH', 'SG']
DEFAULT_ACCEPTOR_ATOMS = ['ND1', 'NE2', 'O', 'OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1', 'OH', 'SD', 'SG']

# CONSTANTS FOR THIS PROJECT (change based on protein)
DATA_FOLDER = "../PDB Data/cdk2/"
SS_POSITIONS= [ (4,12), (17,23), (29,35), (45,56), (66,72), (75,81),
        (85,94), (100,121), (129,131), (133,135), (141,143), (147,152), (170,175),
        (182,199), (207,220), (229,233), (247,252), (256,267), (276,287) ] # hardcoded based on PDB entry

for idx, pos in enumerate(SS_POSITIONS):
    start, end = pos
    if start > 40: start -= 4 # adjust for four unmodeled residues 37-40 in the PDB data
    if end > 40: end -= 4 
    SS_POSITIONS[idx] = (start, end)

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
    
# TODO align them
def mean_center(data):
    data -= data.mean(axis=0)
    return data


def get_coordinates_from_atoms(atoms, rotation=None):
    """
    Returns the coordinates of a list of N atoms in an N x 3 ndarray under a transformation.
    """
    if not rotation:
        rotation = Rotation.identity()
        
    coordinates = np.array( [atom.get_coord() for atom in atoms] )
    coordinates = mean_center(coordinates)
    coordinates = rotation.apply(coordinates)
    
    return coordinates
    
    
"""
Graph construction.
Consider using C-Graphs or Bridge2 package to ensure correctness?
"""
def graph_from_atom_sets(water_atoms, donor_atoms, acceptor_atoms, rotation=None, intramolecular_only=False):
    """
    Takes in coordinates for water, H-bond donors, and H-bond acceptors.
    Outputs a graph with all H bonds between waters and from water to donors or acceptors.
    TODO between donor and acceptor in the protein? would represent secondary structure 
    """
    # Create an empty graph
    Hbond_graph = nx.Graph()
    if not rotation:
        rotation = Rotation.identity()

    # get coordinates
    water_coordinates = get_coordinates_from_atoms(water_atoms, rotation = rotation)    
    donor_coordinates = get_coordinates_from_atoms(donor_atoms, rotation = rotation)    
    acceptor_coordinates = get_coordinates_from_atoms(acceptor_atoms, rotation = rotation)    
    
    # Add nodes to the graph with position as a node attribute
    if not intramolecular_only:
        for idx, coord in enumerate(water_coordinates):
            Hbond_graph.add_node(f"w{idx}", pos=coord, res="water")
    for idx, coord in enumerate(donor_coordinates):
        Hbond_graph.add_node(f"d{idx}", pos=coord, res="donor")
    for idx, coord in enumerate(acceptor_coordinates):
        Hbond_graph.add_node(f"a{idx}", pos=coord, res="acceptor")

    # Add H-bond edges 
    if not intramolecular_only:
        for idx, coord in enumerate(water_coordinates): 
            water_water_distances = np.linalg.norm(coord - water_coordinates, axis=1) # water to water
            connected_edges = np.nonzero( 
                np.logical_and(water_water_distances < HBOND_MAX_DIST, water_water_distances > HBOND_MIN_DIST) 
                )[0]
            Hbond_graph.add_edges_from(
                [ (f"w{idx}", f"w{jdx}") for jdx in connected_edges ],
                color="black", bond_type="ww" )
            
            water_donor_distances = np.linalg.norm(coord - donor_coordinates, axis=1) # water to donor
            connected_edges = np.nonzero( 
                np.logical_and(water_donor_distances < HBOND_MAX_DIST, water_donor_distances > HBOND_MIN_DIST)
                )[0]
            Hbond_graph.add_edges_from( 
                [ (f"w{idx}", f"d{jdx}") for jdx in connected_edges ],
                color="violet", bond_type="wp" )
                
            water_acceptor_distances = np.linalg.norm(coord - acceptor_coordinates, axis=1) # water to acceptor
            connected_edges = np.nonzero( 
                np.logical_and(water_acceptor_distances < HBOND_MAX_DIST, water_acceptor_distances > HBOND_MIN_DIST)
                )[0]
            Hbond_graph.add_edges_from( 
                [ (f"w{idx}", f"a{jdx}") for jdx in connected_edges ],
                color="violet", bond_type="wp" )
        
    for idx, donor in enumerate(donor_coordinates): # donor to acceptor, intra-molecular
        donor_acceptor_distances = np.linalg.norm(donor - acceptor_coordinates, axis=1) # water to donor 
        connected_edges = np.nonzero( 
            np.logical_and(donor_acceptor_distances < HBOND_MAX_DIST, donor_acceptor_distances > HBOND_MIN_DIST)
            )[0]
        Hbond_graph.add_edges_from( 
            [ (f"d{idx}", f"a{jdx}") for jdx in connected_edges 
            if donor_atoms[idx].get_parent() != acceptor_atoms[jdx].get_parent() ], # do not count as H bond both acceptor and donor are in the same residue
            color="green", bond_type="pp") 
        
    return Hbond_graph
    
"""
Graph analysis
"""
def get_water_nodes(Hbond_graph):
    """ 
    Takes an H-bond graph.
    Returns a list of node names for water residues.
    """
    residue_dict = nx.get_node_attributes(Hbond_graph, 'res')
    water_nodes = [n for n in Hbond_graph.nodes() if residue_dict[n] == 'water']
    return water_nodes
       
       
def bond_proportion_by_type(Hbond_graph):
    """
    Takes an H-bond graph.
    Outputs a 3-tuple of the proportion of water-water, water-donor/acceptor, and donor-acceptor H bonds.
    TODO is this a normalized measure? ww scales quadratically with number of waters, wp linearly, and pp is constant
    """
    bond_type_list = list( nx.get_edge_attributes(Hbond_graph, "bond_type").values() )
    total_bonds = len( bond_type_list )
    ww = bond_type_list.count('ww')
    wp = bond_type_list.count('wp')
    pp = bond_type_list.count('pp')
    
    return ww/total_bonds, wp/total_bonds, pp/total_bonds


def get_graph_centrality_metrics(Hbond_graph_dict):
    """
    Takes a dict of graphs with (protein name, Hbond_graph) as key, value pairs.
    Outputs a dataframe with metrics of the water.
    """
    metrics = ["Eigenvalue Centrality", "Degree Centrality", "Betweenness Centrality"]
    df = pd.DataFrame(index = Hbond_graph_dict.keys(), columns = metrics)
    for prot, graph in Hbond_graph_dict.items():
        water_nodes = get_water_nodes(graph)
        eigenvector_centrality = nx.eigenvector_centrality(graph, max_iter=1000)
        degree_centrality = nx.degree_centrality(graph)
        betweenness_centrality = nx.betweenness_centrality(graph, normalized=True)
        
        df.at[prot, "Eigenvalue Centrality"] = np.sort( [eigenvector_centrality[n] for n in water_nodes] )
        df.at[prot, "Degree Centrality"] = np.sort( [degree_centrality[n] for n in water_nodes] )
        df.at[prot, "Betweenness Centrality"] = np.sort( [betweenness_centrality[n] for n in water_nodes] )

    return df


# TODO graph intersection / difference



"""
Visualization
"""
def get_colors(Hbond_graph):
    """
    Utility function to get colors from graph
    """
    type_to_color = {"water": "blue", "donor": "red", "acceptor":"tomato"}
    residue_dict = nx.get_node_attributes(Hbond_graph, 'res')
    node_colors = [ type_to_color.get(residue_dict[n], "gray") for n in Hbond_graph.nodes() ]
    edge_colors = [Hbond_graph[u][v]["color"] for u, v in Hbond_graph.edges()]
    return node_colors, edge_colors


def draw_graph_2d(Hbond_graph, ax):
    """
    Set the characteristics of the graph and draw projected onto XY plane.
    """ 
    pos = nx.get_node_attributes(Hbond_graph, 'pos') # TODO adapt this using PCA since pos is now a 3-tuple
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
    ax.scatter(xs, ys, zs, s=5, c=node_colors, alpha=0.4)

    # Plot edges
    for idx, (u,v) in enumerate(Hbond_graph.edges()):
        x, y, z = zip(pos[u],pos[v])
        ax.plot(x, y, z, color=edge_colors[idx], alpha=0.4)
        
    legend_elements = [     # make a legend for the plot
        Line2D([0], [0], marker='o', color='w', label='Water', markerfacecolor='blue', markersize=10, alpha=0.2),
        Line2D([0], [0], marker='o', color='w', label='Donor', markerfacecolor='red', markersize=10, alpha=0.2),
        Line2D([0], [0], marker='o', color='w', label='Acceptor', markerfacecolor='tomato', markersize=10, alpha=0.2),
        Line2D([0], [1], linewidth=1, linestyle='-', color='black', label='Water to Water H-bond', alpha=0.2),
        Line2D([0], [1], linewidth=1, linestyle='-', color='violet', label='Water to Donor/Acceptor H-bond', alpha=0.2),
        Line2D([0], [1], linewidth=1, linestyle='-', color='green', label='Donor to Acceptor H-bond', alpha=0.2)   
    ]
    
    return legend_elements


def plot_backbone(prot, ax, ss_positions):
    file_name = f"{DATA_FOLDER}/{prot}_final.pdb"
    structure = parser.get_structure(prot, file_name)
    backbone = get_filtered_atoms(structure, residue_target = AA, atom_target = ["CA"])
    coordinates = get_coordinates_from_atoms(backbone)
    
    # split color based on secondary structure
    last_colored = -1
    for (start, end) in ss_positions:
        ax.plot(coordinates[last_colored+1:start,0], coordinates[last_colored+1:start,1], coordinates[last_colored+1:start,2], color="darkgreen") # color the previous segment green
        ax.plot(coordinates[start-1:end,0], coordinates[start-1:end,1], coordinates[start-1:end,2], color="darkred") # color the ss segment red
        last_colored = end-2
        
    ax.plot(coordinates[last_colored+1:-1,0], coordinates[last_colored+1:-1,1], coordinates[last_colored+1:-1,2], color="darkgreen") # color the last segment green
    
    legend_elements = [ # make a legend for the plot
        Line2D([0], [1], linewidth=1, linestyle='-', color='darkgreen', label='Protein backbone'),
        Line2D([0], [1], linewidth=1, linestyle='-', color='darkred', label='Secondary structure')
    ]
    
    return legend_elements


def plot_metric_all_proteins(data_series, title=None):
    """
    Plotting function for a metric for all proteins. Data is a df series of the metric, one for each protein.
    """
    fig, ax = plt.subplots(3, 7, figsize = (12, 14), sharex=True, sharey=True) # 21 proteins
    for idx, data in enumerate(data_series):
        row, col = idx//7, idx%7
        prot = data_series.index[idx]
        ax[row,col].scatter(data)
        ax[row,col].set_title(prot)
        
    if title:
        plt.title(title)

    plt.show()



################## MAIN ##################

""" 
Load all files ending with "_final.pdb" from DATA_FOLDER.
Extract hydrogen bonds from coordinate data and represent as a graph.
"""

pdb_path = pathlib.Path(DATA_FOLDER)
pdb_files = list(pdb_path.glob('*_final.pdb')) 
parser = PDB.PDBParser(QUIET=True)

apo_id = '1PW2'
Hbond_graphs = {}
for idx, file_name in enumerate(pdb_files):
    prot = file_name.stem.split('_')[0]
    structure = parser.get_structure(prot, file_name)
    all_atoms = get_filtered_atoms(structure, residue_target = AA, atom_target = None)
    all_coordinates = get_coordinates_from_atoms(all_atoms, rotation=None)
    rotation = None
    # TODO figure out how to handle different number of atoms
    # if idx == 0:
        # align_template = all_coordinates
        # rotation = Rotation.identity()
    # else:
        # num_atoms = min(len(align_template), len(all_coordinates)) # NOT THIS
        # rotation, _ = Rotation.align_vectors(align_template[:num_atoms], all_coordinates[:num_atoms])
       
    water = get_filtered_atoms(structure, residue_target = WATER, atom_target = ["O"])
    Hbond_donors = get_filtered_atoms(structure, residue_target = AA, atom_target = DEFAULT_DONOR_ATOMS)
    Hbond_acceptors = get_filtered_atoms(structure, residue_target = AA, atom_target = DEFAULT_ACCEPTOR_ATOMS)
    graph = graph_from_atom_sets(water, Hbond_donors, Hbond_acceptors, rotation=rotation, intramolecular_only=False)
    Hbond_graphs[prot] = graph

"""
Get graph metrics.
"""
df_metrics = get_graph_centrality_metrics(Hbond_graphs)
plot_metric_all_proteins(df_metrics["Eigenvalue Centrality"], "Eigenvalue Centrality")
plot_metric_all_proteins(df_metrics["Degree Centrality"], "Degree Centrality")
plot_metric_all_proteins(df_metrics["Betweenness Centrality", "Betweenness Centrality"])


"""
Display hydrogen bonding network for the apo and one holo protein.
"""
# prots = ['1PW2', '3QQL'] # apo, holo
# num_prots = len(prots)
# fig  = plt.figure( figsize=(num_prots * 7, 6) )
# axs = [ fig.add_subplot(1, num_prots, i+1, projection='3d') for i in range(num_prots) ]
# for prot, ax in zip(prots, axs):
    # ax.set_title(f"{prot}")  
    # legend_elements = []
    # legend_elements += draw_graph_3d(Hbond_graphs[prot], ax) # plot graph        
    # legend_elements += plot_backbone(prot, ax, SS_POSITIONS) # plot backbone

# axs[-1].legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.2, 1), title="Graph elements")
# plt.tight_layout()
# plt.show()
