"""
Takes a PDB file with protein and water coordinates.
Creates a graph with vertices of water/AA residues and edges if they have an H bond.
Committed on github under protein_ligand_informatics/water/water_protein_Hbond_network.py

"""

from Bio import PDB

import numpy as np
import pandas as pd
import networkx as nx
from scipy.spatial.transform import Rotation
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
from sklearn.decomposition import PCA

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection

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


def get_coordinates_from_atoms(atoms, rotation=None):
    """
    Returns the coordinates of a list of N atoms in an N x 3 ndarray under a transformation.
    """
    if not rotation:
        rotation = Rotation.identity()
        
    coordinates = np.array( [atom.get_coord() for atom in atoms] )
    coordinates -= coordinates.mean(axis=0)
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
    """
    # Create an empty graph
    Hbond_graph = nx.Graph()
    if not rotation:
        rotation = Rotation.identity()

    # get coordinates
    water_coordinates = get_coordinates_from_atoms(water_atoms, rotation = rotation)    
    donor_coordinates = get_coordinates_from_atoms(donor_atoms, rotation = rotation)    
    acceptor_coordinates = get_coordinates_from_atoms(acceptor_atoms, rotation = rotation)    
    
    # Add nodes to the graph
    if not intramolecular_only:
        for idx, coord in enumerate(water_coordinates):
            # freq is in how many structures this atom can be mapped by map_nodes
            Hbond_graph.add_node(water_atoms[idx].get_serial_number(), pos=coord, res="water", freq=1)
    for idx, coord in enumerate(donor_coordinates):
        Hbond_graph.add_node(donor_atoms[idx].get_serial_number(), pos=coord, res="donor")
    for idx, coord in enumerate(acceptor_coordinates):
        Hbond_graph.add_node(acceptor_atoms[idx].get_serial_number(), pos=coord, res="acceptor")

    # Add H-bond edges 
    if not intramolecular_only:
        for idx, coord in enumerate(water_coordinates): 
            water_water_distances = np.linalg.norm(coord - water_coordinates, axis=1) # water to water
            connected_edges = np.nonzero( 
                np.logical_and(water_water_distances < HBOND_MAX_DIST, water_water_distances > HBOND_MIN_DIST) 
                )[0]
            Hbond_graph.add_edges_from(
                [ (water_atoms[idx].get_serial_number(), water_atoms[idx].get_serial_number())
                for jdx in connected_edges ],
                color="black", bond_type="ww" )
            
            water_donor_distances = np.linalg.norm(coord - donor_coordinates, axis=1) # water to donor
            connected_edges = np.nonzero( 
                np.logical_and(water_donor_distances < HBOND_MAX_DIST, water_donor_distances > HBOND_MIN_DIST)
                )[0]
            Hbond_graph.add_edges_from( 
                [ (water_atoms[idx].get_serial_number(), donor_atoms[idx].get_serial_number())
                for jdx in connected_edges ],
                color="violet", bond_type="wp" )
                
            water_acceptor_distances = np.linalg.norm(coord - acceptor_coordinates, axis=1) # water to acceptor
            connected_edges = np.nonzero( 
                np.logical_and(water_acceptor_distances < HBOND_MAX_DIST, water_acceptor_distances > HBOND_MIN_DIST)
                )[0]
            Hbond_graph.add_edges_from( 
                [ (water_atoms[idx].get_serial_number(), acceptor_atoms[idx].get_serial_number())
                for jdx in connected_edges ],
                color="violet", bond_type="wp" )
        
    for idx, donor in enumerate(donor_coordinates): # donor to acceptor, intra-molecular
        donor_acceptor_distances = np.linalg.norm(donor - acceptor_coordinates, axis=1) # water to donor 
        connected_edges = np.nonzero( 
            np.logical_and(donor_acceptor_distances < HBOND_MAX_DIST, donor_acceptor_distances > HBOND_MIN_DIST)
            )[0]
        Hbond_graph.add_edges_from( 
            [ (donor_atoms[idx].get_serial_number(), acceptor_atoms[idx].get_serial_number()) 
            for jdx in connected_edges 
            if donor_atoms[idx].get_parent() != acceptor_atoms[jdx].get_parent() ], # do not count as H bond if both acceptor and donor are in the same residue
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
      
      
def get_protein_nodes(Hbond_graph):
    """ 
    Takes an H-bond graph.
    Returns a list of node names for water residues.
    """
    residue_dict = nx.get_node_attributes(Hbond_graph, 'res')
    protein_nodes = [n for n in Hbond_graph.nodes() if residue_dict[n] != 'water']
    return protein_nodes
      

def map_nodes(graph_ref, node_list_ref, graph_align, node_list_align, distance_cutoff=HBOND_MAX_DIST):
    """
    Maps nodes in node_list_align from graph_align to nodes in node_list_ref from graph_ref
    """
    ref_coordinates = np.vstack( [graph_ref.nodes()[n]['pos'] for n in node_list_ref] )
    align_coordinates = np.vstack( [graph_align.nodes()[n]['pos'] for n in node_list_align] )
    
    # calculate the square distance between two sets of nodes as an proxy for Gaussian MLE mapping between nodes
    distance_matrix = cdist(ref_coordinates, align_coordinates, metric='sqeuclidean')
    ref_nodes_idx, align_nodes_idx = linear_sum_assignment(distance_matrix)
    
    # keep assignments within a cutoff distance
    keep = np.where( distance_matrix[ref_nodes_idx, align_nodes_idx] < distance_cutoff**2 ) 
    ref_nodes_idx, align_nodes_idx = ref_nodes_idx[keep], align_nodes_idx[keep]

    node_list_ref = [ node_list_ref[i] for i in ref_nodes_idx ]
    node_list_align = [ node_list_align[i] for i in align_nodes_idx ]
    
    return node_list_ref, node_list_align

    
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
    Outputs a dataframe with centrality metrics of the water.
    """
    metrics = ["Eigenvalue Centrality", "Degree Centrality", "Betweenness Centrality"]
    df = pd.DataFrame(index = Hbond_graph_dict.keys(), columns = metrics)
    for prot, graph in Hbond_graph_dict.items():
        water_nodes = get_water_nodes(graph)
        eigenvector_centrality = nx.eigenvector_centrality(graph, max_iter=1000)
        degree_centrality = nx.degree_centrality(graph)
        betweenness_centrality = nx.betweenness_centrality(graph, normalized=True)
        
        df.at[prot, "Eigenvalue Centrality"] = np.sort( [eigenvector_centrality[n] for n in water_nodes], ascending=False )
        df.at[prot, "Degree Centrality"] = np.sort( [degree_centrality[n] for n in water_nodes], ascending=False )
        df.at[prot, "Betweenness Centrality"] = np.sort( [betweenness_centrality[n] for n in water_nodes], ascending=False )

    return df


# TODO graph intersection / difference



"""
Visualization
"""
def get_colors(Hbond_graph):
    """
    Utility function to get colors from graph.
    Input is an H-bond graph.
    Outputs are the node and edge colors for drawing as defined in the graph.
    """
    type_to_color = {"water": "blue", "donor": "red", "acceptor":"tomato"}
    residue_dict = nx.get_node_attributes(Hbond_graph, 'res')
    node_colors = [ type_to_color.get(residue_dict[n], "gray") for n in Hbond_graph.nodes() ]
    edge_colors = [Hbond_graph[u][v]["color"] for u, v in Hbond_graph.edges()]
    return node_colors, edge_colors


def draw_graph_2d(Hbond_graph, ax): 
    """
    Set the characteristics of the graph and draw projected onto XY plane.
    Input is an H-bond graph, the plot axis with 2-D projection.
    Output is a list of legend elements for the plot.
    """ 
    pos = nx.get_node_attributes(Hbond_graph, 'pos')
    xy = np.vstack( list(pos.values()) )[:, :2]
    z = np.vstack( list(pos.values()) )[:,-1:] # -1: to keep as 2d array
    xy_principal = PCA(n_components=1).fit_transform( xy )
    coordinates_2d = np.hstack(  [xy_principal, z] )
    pos_2d = dict( zip(Hbond_graph.nodes(), coordinates_2d) ) # pair new coordinate back with original node
    node_colors, edge_colors = get_colors(Hbond_graph)
    
    nx.draw(Hbond_graph, pos_2d, ax=ax, node_color=node_colors, edge_color=edge_colors, node_size=3)
    
    legend_elements = [     # make a legend for the plot
    Line2D([0], [0], marker='o', color='w', label='Water', markerfacecolor='blue', markersize=10, alpha=0.2),
    Line2D([0], [0], marker='o', color='w', label='Donor', markerfacecolor='red', markersize=10, alpha=0.2),
    Line2D([0], [0], marker='o', color='w', label='Acceptor', markerfacecolor='tomato', markersize=10, alpha=0.2),
    Line2D([0], [1], linewidth=1, linestyle='-', color='black', label='Water to Water H-bond', alpha=0.2),
    Line2D([0], [1], linewidth=1, linestyle='-', color='violet', label='Water to Donor/Acceptor H-bond', alpha=0.2),
    Line2D([0], [1], linewidth=1, linestyle='-', color='green', label='Donor to Acceptor H-bond', alpha=0.2)   
    ]
    
    return legend_elements


def plot_backbone_2d(coordinates, ax, ss_positions): 
    """
    Plotting function for 2D backbone. 
    Input are the 3D coordinates of the backbone, the plot axis with 2-D projection, and a list of tuples (start, end) of indices containing secondary structure elements.
    Output is a list of legend elements for the plot.
    """
    # split color based on secondary structure
    last_colored = -1
    xy = coordinates[:,:2]
    z = coordinates[:,-1:] # -1: to keep as 2d array
    xy_principal = PCA(n_components=1).fit_transform( xy )
    coordinates_2d = np.hstack( [xy_principal, z] )
    
    for (start, end) in ss_positions:
        ax.plot(coordinates_2d[last_colored+1:start,0], coordinates_2d[last_colored+1:start,1], color="darkgreen") # color the previous segment green
        ax.plot(coordinates_2d[start-1:end,0], coordinates_2d[start-1:end,1], color="darkred") # color the ss segment red
        last_colored = end-2
        
    ax.plot(coordinates_2d[last_colored+1:-1,0], coordinates_2d[last_colored+1:-1,1], color="darkgreen") # color the last segment green
    
    legend_elements = [ # make a legend for the plot
        Line2D([0], [1], linewidth=1, linestyle='-', color='darkgreen', label='Protein backbone'),
        Line2D([0], [1], linewidth=1, linestyle='-', color='darkred', label='Secondary structure')
    ]
    
    return legend_elements


def draw_graph_3d(Hbond_graph, ax):
    """
    Plotting function for 3D graph. 
    Input is an H-bond graph, the plot axis with 3-D projection.
    Output is a list of legend elements for the plot.
    Note: to draw a graph with a subset of nodes, use a graph.subgraph(nodes) as the input
    """
    
    pos = nx.get_node_attributes(Hbond_graph, 'pos')
    node_colors, edge_colors = get_colors(Hbond_graph)
    
    # Plot nodes
    # List comprehension gives a list of (x,y,z) coordinates and zip(*([])) transposes them to the format for matplotlib
    xs, ys, zs = zip( *[pos[n] for n in Hbond_graph.nodes()] ) 
    ax.scatter(xs, ys, zs, s=50, c=node_colors, alpha=0.4)

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


    
def plot_backbone_3d(coordinates, ax, ss_positions):
    """
    Plotting function for 3D backbone. 
    Input are the coordinates of the backbone, the plot axis with 3-D projection, and a list of tuples (start, end) of indices containing secondary structure elements.
    Output is a list of legend elements for the plot.
    """
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


def plot_metric_all_proteins(data_series, limit=None, title=None):
    """
    Plotting function for a metric for all proteins. Data is a df series of the metric, one for each protein.
    """
    
    fig, ax = plt.subplots(3, 7, figsize = (12, 14), sharex=True, sharey=True) # 21 proteins
    for idx, data in enumerate(data_series):
        row, col = idx//7, idx%7
        prot = data_series.index[idx]
        if not limit:
            limit = len(data) 
        ax[row,col].scatter(data[:limit])
        ax[row,col].set_title(prot)
        
    if title:
        plt.title(title)

    plt.show()


"""
Edit PDB b factors to view in PyMol.
input_file and output_file are string names of the PDB files.
b_factor_dict is a dict of (node_id:new_value) pairs
"""
def edit_b_factor(input_file, output_file, b_factor_dict):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("test", input_file)

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom.set_bfactor( b_factor_dict.get(atom.get_serial_number(), 0) )

    # Save the modified structure
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_file)


################## MAIN ##################

""" 
Load all files ending with "_final.pdb" from DATA_FOLDER.
Extract hydrogen bonds from coordinate data and represent as a graph.
"""

pdb_path = pathlib.Path(DATA_FOLDER)
pdb_files = list(pdb_path.glob('*_final.pdb')) 
parser = PDB.PDBParser(QUIET=True)

apo_prot = '1PW2'
test_id = '3QZH'
Hbond_graphs = {}
backbone_coordinates = {}
for idx, file_name in enumerate(pdb_files):
    prot = file_name.stem.split('_')[0]
    structure = parser.get_structure(prot, file_name)
    
    backbone = get_filtered_atoms(structure, residue_target = AA, atom_target = ["CA"])
    water = get_filtered_atoms(structure, residue_target = WATER, atom_target = ["O"])
    Hbond_donors = get_filtered_atoms(structure, residue_target = AA, atom_target = DEFAULT_DONOR_ATOMS)
    Hbond_acceptors = get_filtered_atoms(structure, residue_target = AA, atom_target = DEFAULT_ACCEPTOR_ATOMS)
    
    # find rotation to align on alpha carbon atoms
    if idx == 0:
        align_backbone = backbone.copy()
        align_template = get_coordinates_from_atoms(align_backbone, rotation=None)
        rotation = Rotation.identity()
    else:
        coordinates = get_coordinates_from_atoms(backbone, rotation=None)
        rotation, _ = Rotation.align_vectors(align_template, coordinates)
       
    graph = graph_from_atom_sets(water, Hbond_donors, Hbond_acceptors, rotation=rotation, intramolecular_only=False)
    Hbond_graphs[prot] = graph
    backbone_coordinates[prot] = get_coordinates_from_atoms(backbone, rotation=rotation)

"""
Get graph metrics.
"""
# TODO hopefully will be more useful once extraneous waters have been removed programatically
# df_metrics = get_graph_centrality_metrics(Hbond_graphs)
# plot_metric_all_proteins(df_metrics["Eigenvalue Centrality"], limit=10, title="Eigenvalue Centrality")
# plot_metric_all_proteins(df_metrics["Degree Centrality"], limit=10, title="Degree Centrality")
# plot_metric_all_proteins(df_metrics["Betweenness Centrality"], limit=10, title="Betweenness Centrality")


"""
Display hydrogen bonding network for the apo and one holo protein.
"""

# starting with the apo structure waters, identify waters conserved in all structures
apo_prot = '1PW2'
holo_prot = '3QZH'
distance_cutoff = HBOND_MAX_DIST # max distance for two waters to map together

# loop through all proteins and keep the remaining mappable waters
for ref_prot in Hbond_graphs.keys():
    ref_waters = get_water_nodes(Hbond_graphs[ref_prot])
    for test_prot in Hbond_graphs.keys():
        if test_prot == ref_prot:
            continue
        test_waters = get_water_nodes(Hbond_graphs[test_prot])
        ref_waters_mapped, _ = map_nodes(Hbond_graphs[ref_prot], ref_waters, Hbond_graphs[test_prot], test_waters, distance_cutoff)  
        for w in ref_waters_mapped:
            Hbond_graphs[ref_prot].nodes()[w]['freq'] += 1    

        edit_b_factor( input_file = f"{DATA_FOLDER}{ref_prot}_final.pdb", 
            output_file = f"{DATA_FOLDER}/all_proteins_edited/{ref_prot}_edited.pdb", 
            b_factor_dict = nx.get_node_attributes(Hbond_graphs[ref_prot], 'freq') )
    
    
# prots = [apo_prot]
# num_prots = len(prots)
# fig  = plt.figure( figsize=(num_prots * 7, 6) )
# axs = [ fig.add_subplot(1, num_prots, i+1, projection='3d') for i in range(num_prots) ]

# for prot, ax in zip(prots, axs):
    # ax.set_title(f"{prot}")  
    # legend_elements = []
    # legend_elements += draw_graph_3d(Hbond_graphs[prot].subgraph(common_waters), ax) # plot graph        
    # legend_elements += plot_backbone_3d(backbone_coordinates[prot], ax, SS_POSITIONS) # plot backbone

# axs[-1].legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.2, 1), title="Graph elements")
# plt.tight_layout()
# plt.show()

