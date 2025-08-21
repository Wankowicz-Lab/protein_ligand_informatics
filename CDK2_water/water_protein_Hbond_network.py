"""
Takes a PDB file with protein and water coordinates.
Creates a graph with vertices of water/AA residues and edges if they have an H bond.
Committed on github under protein_ligand_informatics/water/water_protein_Hbond_network.py
"""

from Bio import PDB
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

import numpy as np
import networkx as nx
from scipy.spatial.transform import Rotation
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment

from sklearn.decomposition import PCA
from sklearn.cluster import DBSCAN
from sklearn.manifold import MDS
from skbio.stats.distance import mantel

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection

import glob
import pathlib
import sys
import requests
import csv

# UNIVERSAL CONSTANTS
AA=['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
WATER=['HOH', 'WAT', 'H2O']  

HBOND_MAX_DIST=3.5
HBOND_MIN_DIST=2.4

# Default hydrogen bond donor and acceptor definitions from Mingbin
DEFAULT_DONOR_ATOMS = ['N', 'ND1', 'ND2', 'NE', 'NE1', 'NE2', 'NH1', 'NH2', 'NZ', 'OG', 'OG1', 'OH', 'SG']
DEFAULT_ACCEPTOR_ATOMS = ['ND1', 'NE2', 'O', 'OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1', 'OH', 'SD', 'SG']

# CONSTANTS FOR THIS PROJECT (change based on protein)
DATA_FOLDER = "../PDB_Data_cdk2/cdk2_original/"
NUM_RES = 298
OTHER_SOLVENT=["EDO","PO4"] # constructed by looking through all non WATER + AA or ligand labels
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
    residue_target [optional] and atom_target [optional] are lists of acceptable values. If none, return all values.
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
    

def get_coordinates_from_atoms(atoms, translation=None, rotation=None):
    """
    Returns the coordinates of a list of N atoms in an N x 3 ndarray under a transformation.
    translation [optional] is an ndarray with shape (3,)
    rotation [optional] is a scipy Rotation object
    """
    if atoms is None:
        return None
    if translation is None:
        translation = np.zeros([3])
    if rotation is None:
        rotation = Rotation.identity()
        
    coordinates = np.vstack( [atom.get_coord() for atom in atoms] )
    coordinates += translation
    coordinates = rotation.apply(coordinates)
    
    return coordinates


def get_ligand_residue(pdb_structure):
    """
    Gets the heteroatom not labeled as an amino acid or a solvent.
    """
    for model in pdb_structure:
        for chain in model:
            for res in chain:
                resname = res.get_resname()
                hetflag = res.get_id()[0]
                if hetflag != ' ' and resname not in WATER + OTHER_SOLVENT:
                    return res


def b_factor_by_residue(pdb_structure, normalize=True):
    """
    https://pubs.acs.org/doi/pdf/10.1021/acs.chemrev.8b00290?ref=article_openPDF section 2.1 for method
    """
    average_b_factor = {}
    for model in pdb_structure:
        for chain in model:
            for res in chain:
                if res.resname not in AA: continue
                avg_b = np.mean( [a.get_bfactor() for a in res] )
                average_b_factor[ res.get_id()[1] ] = avg_b
    
    if normalize:
        med = np.median( list(average_b_factor.values()) )
        average_b_factor = {res:b/med - 1 for res, b in average_b_factor.items()}
    
    return average_b_factor
    


def residue_to_mol(residue):
    """
    Written by microsoft copilot.
    Takes a residue and writes a mol object.
    """

    if residue is None:
        return ""

    ligand_id = residue.get_resname()
    url = f"https://files.rcsb.org/ligands/download/{ligand_id}_ideal.sdf"
    response = requests.get(url)
    if response.status_code == 200:
        sdf_data = response.text
        suppl = Chem.SDMolSupplier()
        suppl.SetData(sdf_data)
        mol = suppl[0]
        return mol 
    else:
        raise ValueError(f"Ligand {ligand_id} not found at RCSB")

      
    
    
def tanimoto_dissimilarity(mol1, mol2):
    """
    Takes two mol objects and outputs the tanimoto dissimilarity between them
    """

    gen = Chem.rdFingerprintGenerator.GetRDKitFPGenerator(
        maxPath=7,       # max bond path length
        fpSize=2048,     # fingerprint size
        useHs=False,     # ignore explicit hydrogens
        branchedPaths=True,  # include branched paths
        useBondOrder=True    # include bond order info
    )
    
    fp1 = gen.GetFingerprint(mol1)
    fp2 = gen.GetFingerprint(mol2)

    dissimilarity = 1 - DataStructs.TanimotoSimilarity(fp1, fp2)
    
    return dissimilarity
    
    
"""
Graph construction.
"""
def graph_from_atom_sets(water_atoms, donor_atoms, acceptor_atoms, translation=None, rotation=None):
    """
    Takes in coordinates for water, H-bond donors, and H-bond acceptors.
    Outputs a graph with all H bonds between waters and from water to donors or acceptors.
    """
    # Create an empty graph
    Hbond_graph = nx.Graph()

    # get coordinates
    water_coordinates = get_coordinates_from_atoms(water_atoms, translation=translation, rotation = rotation)     
    donor_coordinates = get_coordinates_from_atoms(donor_atoms, translation=translation, rotation = rotation)    
    acceptor_coordinates = get_coordinates_from_atoms(acceptor_atoms, translation=translation, rotation = rotation)    
    
    # Add nodes to the graph
    for idx, coord in enumerate(water_coordinates):
        w_tag = ("w", str( water_atoms[idx].get_serial_number() ) )
        # freq is in how many structures this atom can be mapped by map_nodes
        Hbond_graph.add_node(w_tag, pos=coord, res="water", freq=0)    
    for idx, coord in enumerate(donor_coordinates):
        d_tag = ("d", donor_atoms[idx].get_parent().get_id()[1], donor_atoms[idx].name) 
        Hbond_graph.add_node(d_tag, pos=coord, res="donor", freq=0)
    for idx, coord in enumerate(acceptor_coordinates):
        a_tag = ("a", acceptor_atoms[idx].get_parent().get_id()[1], acceptor_atoms[idx].name) 
        Hbond_graph.add_node(a_tag, pos=coord, res="acceptor", freq=0)

    # Add H-bond edges 
    for idx, coord in enumerate(water_coordinates): 
        w_tag = ("w", str( water_atoms[idx].get_serial_number() ) )
        
        water_water_distances = np.linalg.norm(coord - water_coordinates, axis=1) # water to water
        connected_edges = np.nonzero( 
            np.logical_and(water_water_distances < HBOND_MAX_DIST, water_water_distances > HBOND_MIN_DIST) 
            )[0]
        Hbond_graph.add_edges_from(
            [ (w_tag, ("w", str( water_atoms[jdx].get_serial_number() ) ))
            for jdx in connected_edges ],
            color="black", bond_type="ww", style="solid" )
        
        water_donor_distances = np.linalg.norm(coord - donor_coordinates, axis=1) # water to donor
        connected_edges = np.nonzero( 
            np.logical_and(water_donor_distances < HBOND_MAX_DIST, water_donor_distances > HBOND_MIN_DIST)
            )[0]
        
        Hbond_graph.add_edges_from(
            [ ( w_tag, ("d", donor_atoms[jdx].get_parent().get_id()[1], donor_atoms[jdx].name) )
            for jdx in connected_edges ],
            color="violet", bond_type="wp", style="solid" )
            
        water_acceptor_distances = np.linalg.norm(coord - acceptor_coordinates, axis=1) # water to acceptor
        connected_edges = np.nonzero( 
            np.logical_and(water_acceptor_distances < HBOND_MAX_DIST, water_acceptor_distances > HBOND_MIN_DIST)
            )[0]        
        Hbond_graph.add_edges_from( 
            [ (w_tag, ("a", acceptor_atoms[jdx].get_parent().get_id()[1], acceptor_atoms[jdx].name) )
            for jdx in connected_edges ],
            color="violet", bond_type="wp", style="solid" )
        
    # TODO maybe add geometry condition here for increased accuracy? 90 deg or 120 deg angle 
    for idx, donor in enumerate(donor_coordinates): # donor to acceptor, intra-molecular
        d_tag = ("d", donor_atoms[idx].get_parent().get_id()[1], donor_atoms[idx].name)
        
        donor_acceptor_distances = np.linalg.norm(donor - acceptor_coordinates, axis=1) # water to donor 
        connected_edges = np.nonzero( 
            np.logical_and(donor_acceptor_distances < HBOND_MAX_DIST, donor_acceptor_distances > HBOND_MIN_DIST)
            )[0]
        Hbond_graph.add_edges_from( 
            [ (d_tag, ("a", acceptor_atoms[jdx].get_parent().get_id()[1], acceptor_atoms[jdx].name) ) 
            for jdx in connected_edges 
            if donor_atoms[idx].get_parent() != acceptor_atoms[jdx].get_parent() ], # do not count as H bond if both acceptor and donor are in the same residue
            color="green", bond_type="pp", style="solid" ) 
        
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
    water_nodes = [n for n in Hbond_graph.nodes if residue_dict[n] == 'water']
    return water_nodes
      
      
def get_protein_nodes(Hbond_graph):
    """ 
    Takes an H-bond graph.
    Returns a list of node names for water residues.
    """
    residue_dict = nx.get_node_attributes(Hbond_graph, 'res')
    protein_nodes = [n for n in Hbond_graph.nodes if residue_dict[n] != 'water']
    return protein_nodes
      

def get_attr_from_nodes(G, nodes, attr):
    """
    Gets attr from nodes of graph G.
    Nodes can be 'water', 'protein', 'all', or a list of node names
    """
    if nodes == 'all':
        nodes = G.nodes
    elif nodes == 'water':
        nodes = get_water_nodes(G)
    elif nodes == 'protein':
        nodes = get_protein_nodes(G)
        
    return [G.nodes[node][attr] for node in nodes]


def map_nodes(graph_ref, node_list_ref, graph_align, node_list_align, distance_cutoff=HBOND_MAX_DIST):
    """
    Maps nodes in node_list_align from graph_align to nodes in node_list_ref from graph_ref
    """
    ref_coordinates = np.vstack( get_attr_from_nodes(graph_ref, node_list_ref, 'pos' ) )
    align_coordinates = np.vstack( get_attr_from_nodes(graph_align, node_list_align, 'pos') )
    
    # calculate the square distance between two sets of nodes as an proxy for Gaussian MLE mapping between nodes
    distance_matrix = cdist(ref_coordinates, align_coordinates, metric='sqeuclidean')
    distance_matrix[distance_matrix > distance_cutoff**2] = 1e6 # penalize mappings outside of the bound equally
    ref_nodes_idx, align_nodes_idx = linear_sum_assignment(distance_matrix)
        
    # filter out assignments outside of cutoff distance
    keep = np.where( distance_matrix[ref_nodes_idx, align_nodes_idx] < distance_cutoff**2 ) 
    ref_nodes_idx, align_nodes_idx = ref_nodes_idx[keep], align_nodes_idx[keep]
    
    node_list_ref = [ node_list_ref[i] for i in ref_nodes_idx ]
    node_list_align = [ node_list_align[i] for i in align_nodes_idx ]   
           
    return node_list_ref, node_list_align


def get_totally_conserved_waters(Hbond_graph):
    waters = get_water_nodes(Hbond_graph)
    freqs = get_attr_from_nodes(Hbond_graph, waters, 'freq')
    max_freq = max(freqs)
    return [w for w,f in zip(waters, freqs) if f==max_freq]
    
    
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
    node_colors = [ type_to_color.get(residue_dict[n], "gray") for n in Hbond_graph.nodes ]
    edge_colors = [Hbond_graph[u][v]["color"] for u, v in Hbond_graph.edges]
    edge_styles = [Hbond_graph[u][v]["style"] for u, v in Hbond_graph.edges]
    return node_colors, edge_colors, edge_styles


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
    pos_2d = dict( zip(Hbond_graph.nodes, coordinates_2d) ) # pair new coordinate back with original node
    node_colors, edge_colors, edge_styles = get_colors(Hbond_graph)
    
    nx.draw(Hbond_graph, pos_2d, ax=ax, node_color=node_colors, edge_color=edge_colors, style=edge_styles, node_size=3, alpha=0.4)
    
    legend_elements = [     # make a legend for the plot
    Line2D([0], [0], marker='o', color='w', label='Water', markerfacecolor='blue', markersize=10, alpha=0.4),
    Line2D([0], [0], marker='o', color='w', label='Donor', markerfacecolor='red', markersize=10, alpha=0.4),
    Line2D([0], [0], marker='o', color='w', label='Acceptor', markerfacecolor='tomato', markersize=10, alpha=0.4),
    Line2D([0], [1], linewidth=1, linestyle='-', color='black', label='Water to Water H-bond', alpha=0.4),
    Line2D([0], [1], linewidth=1, linestyle='-', color='violet', label='Water to Donor/Acceptor H-bond', alpha=0.4),
    Line2D([0], [1], linewidth=1, linestyle='-', color='green', label='Donor to Acceptor H-bond', alpha=0.4)   
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
    node_colors, edge_colors, edge_styles = get_colors(Hbond_graph)
    edge_styles = ['--' if s=='dashed' else '-' for s in edge_styles ]
    
    # Plot nodes
    # List comprehension gives a list of (x,y,z) coordinates and zip(*([])) transposes them to the format for matplotlib
    xs, ys, zs = zip( *[pos[n] for n in Hbond_graph.nodes] ) 
    ax.scatter(xs, ys, zs, s=5, c=node_colors, alpha=0.4)

    # Plot edges
    for idx, (u,v) in enumerate(Hbond_graph.edges):
        x, y, z = zip(pos[u],pos[v])
        ax.plot(x, y, z, edge_styles[idx], color=edge_colors[idx], alpha=0.4)
        
    legend_elements = [     # make a legend for the plot
        Line2D([0], [0], marker='o', color='w', label='Water', markerfacecolor='blue', markersize=10, alpha=0.4),
        Line2D([0], [0], marker='o', color='w', label='Donor', markerfacecolor='red', markersize=10, alpha=0.4),
        Line2D([0], [0], marker='o', color='w', label='Acceptor', markerfacecolor='tomato', markersize=10, alpha=0.4),
        Line2D([0], [1], linewidth=1, linestyle='-', color='black', label='Water to Water H-bond', alpha=0.4),
        Line2D([0], [1], linewidth=1, linestyle='-', color='violet', label='Water to Donor/Acceptor H-bond', alpha=0.4),
        Line2D([0], [1], linewidth=1, linestyle='-', color='green', label='Donor to Acceptor H-bond', alpha=0.4)   
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
        if limit is None:
            limit = len(data) 
        ax[row,col].scatter(data[:limit])
        ax[row,col].set_title(prot)
        
    if title:
        plt.title(title)

    plt.show()


def edit_b_factor(input_file, output_file, b_factor_dict):
    """
    Edit PDB b factors to view in PyMol.
    input_file and output_file are string names of the PDB files.
    coord_dict is a dict of (node_id:[x,y,z]) pairs
    b_factor_dict is a dict of (node_id:new_value) pairs
    """
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
all_prots = []
Hbond_graphs = {}
backbone_coordinates = {}
ligands = {}
    
print(f"Parsing pdb files")
for idx, file_name in enumerate(pdb_files):
    prot = file_name.stem.split('_')[0]
    all_prots.append(prot)
    structure = parser.get_structure(prot, file_name)

    ligand = get_ligand_residue(structure)
    ligands[prot] = ligand
    
    backbone = get_filtered_atoms(structure, residue_target = AA, atom_target = ["CA"])
    water = get_filtered_atoms(structure, residue_target = WATER, atom_target = ["O"])
    Hbond_donors = get_filtered_atoms(structure, residue_target = AA, atom_target = DEFAULT_DONOR_ATOMS)
    Hbond_acceptors = get_filtered_atoms(structure, residue_target = AA, atom_target = DEFAULT_ACCEPTOR_ATOMS)    
    
    # find translation and rotation to align on alpha carbon atoms
    coordinates = get_coordinates_from_atoms(backbone, translation=None, rotation=None)
    translation = -np.mean(coordinates, axis=0)
    coordinates += translation
    if idx == 0: # use first structure as template
        align_template = coordinates.copy()
        rotation = Rotation.identity()
    else: # compute translation and rotation
        rotation, _ = Rotation.align_vectors(align_template, coordinates)

    graph = graph_from_atom_sets(water, Hbond_donors, Hbond_acceptors, translation=translation, rotation=rotation)
    Hbond_graphs[prot] = graph
    backbone_coordinates[prot] = get_coordinates_from_atoms(backbone, translation=translation, rotation=rotation)
    
    
"""
Identify conserved waters and compare between ligands.
"""
distance_cutoff = HBOND_MAX_DIST/2 # max distance for two waters to map together, estimated by inflection point in number of water pairs as a function of this distance, leaves 37 waters conserved in all structures
num_prot = len(all_prots)
mapped_water_counts = np.zeros([num_prot, num_prot])

print(f"Mapping waters")
# loop through all pairs of proteins and find the mappable waters
for idx, ref_prot in enumerate(all_prots): 
    ref_waters = get_water_nodes(Hbond_graphs[ref_prot])
    
    for jdx, test_prot in enumerate(all_prots): 
        test_waters = get_water_nodes(Hbond_graphs[test_prot])
        ref_waters_mapped, test_waters_mapped = map_nodes(Hbond_graphs[ref_prot], ref_waters, 
                                                          Hbond_graphs[test_prot], test_waters, distance_cutoff)
        mapped_water_counts[idx, jdx] = len(ref_waters_mapped)
        
        # update the freq trait based on mapping
        for w in ref_waters_mapped:
            Hbond_graphs[ref_prot].nodes[w]['freq'] += 1
           
    # output a new pdb with the freq count in the b-factor column to visualize in pymol
    # edit_b_factor( input_file = f"{DATA_FOLDER}{ref_prot}_final.pdb", 
        # output_file = f"{DATA_FOLDER}/all_proteins_edited/{ref_prot}_edited.pdb", 
        # b_factor_dict = nx.get_node_attributes(Hbond_graphs[ref_prot], 'freq') )


# # calculate tanimoto dissimilarity metrics 
# lig_dissimilarities = np.zeros([num_prot-1, num_prot-1])
# water_dissimilarities = np.zeros_like(lig_dissimilarities)

# # downloading the ligands takes a while, so do all at once only once
# print("Downloading ligand molecules")
# ligands = {prot: residue_to_mol(res) for prot, res in ligands.items()}

# # +1 indexing since the ligand matrices have one less structure (apo) than the 
# for idx, ref_prot in enumerate(all_prots[1:]): 
    # mol_ref = ligands[ref_prot]
    # print(f"Comparing ligands for {ref_prot}")

    # for jdx, test_prot in enumerate(all_prots[1:]): 
        # if ref_prot == test_prot: continue # do not map to self
        # mol_test = ligands[test_prot]
        
        # chem_dissimilarity = tanimoto_dissimilarity(mol_ref, mol_test)
        # lig_dissimilarities[idx, jdx] = chem_dissimilarity
        
        # # calculate the tanimoto = |intersection| / |union| of water atoms that can be mapped to at least one other structure
        # total_ref_waters = np.sum( np.array(get_attr_from_nodes(Hbond_graphs[ref_prot], 'water', "freq")) > 1 )
        # total_test_waters = np.sum( np.array(get_attr_from_nodes(Hbond_graphs[test_prot], 'water', "freq")) > 1 )
        # mapped_waters = mapped_water_counts[idx+1, jdx+1] 
        
        # water_dissimilarity = 1 - ( mapped_waters / ( total_ref_waters + total_test_waters - mapped_waters ) )
        # water_dissimilarities[idx, jdx] = water_dissimilarity 


# Aggregate waters over all structures, get conservation level of water, compare to number of H bonds to protein
num_neighbors_dict = {i+1:[] for i in range(num_prot)}
for prot in all_prots:
    waters = get_water_nodes(Hbond_graphs[prot])
    freqs = get_attr_from_nodes(Hbond_graphs[prot], waters, 'freq')
    for w, f in zip(waters, freqs):
        num_bonds = len( [_ for _, target in Hbond_graphs[prot].edges(w) if target[0] != 'w'] )
        num_neighbors_dict[f].append(num_bonds)

means = {k:np.mean(v) for k, v in num_neighbors_dict.items()}

# Create bar plot
fig, ax = plt.subplots(figsize=(12, 6))
bars = ax.bar(means.keys(), means.values(), alpha=0.7)

# Annotate bars with values
for bar in bars:
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2, height + 0.05, f'{height:.1f}',
            ha='center', va='bottom', fontsize=10)

# Customize plot
ax.set_xlabel('Frequency of water')
ax.set_xlim([0.5,21.5])
ax.set_xticks(np.linspace(1,21,21))
ax.set_ylim([0,3])
ax.set_ylabel('Average number of H-bonds with protein')

plt.tight_layout()
plt.show()

sys.exit()

  
"""
Calculate graph difference for the set of proteins.
"""  

for other_prot in all_prots[1:]:
    w1 = get_water_nodes(Hbond_graphs[apo_prot])
    w2 = get_water_nodes(Hbond_graphs[other_prot]) 
    
    # remap water nodes to make node sets equivalent
    w1, w2 = map_nodes(Hbond_graphs[apo_prot], w1, Hbond_graphs[other_prot], w2, distance_cutoff=distance_cutoff) 
    subnodes1 = get_protein_nodes(Hbond_graphs[apo_prot]) + w1
    subnodes2 = get_protein_nodes(Hbond_graphs[other_prot]) + w2
    G1 = Hbond_graphs[apo_prot].subgraph(subnodes1)
    G2 = Hbond_graphs[other_prot].subgraph(subnodes2)
    
    # map and relabel nodes
    G2 = nx.relabel_nodes( G2, dict(zip( w2, w1 )) )        
    
    # calculate difference and restore properties
    diff = nx.symmetric_difference(G1, G2)
    
    for n in diff.nodes:
        diff.nodes[n].update(G1.nodes[n])
        
    for u, v in diff.edges:            
        if G1.has_edge(u, v):
            diff.edges[u, v].update(G1.edges[u, v])
            diff.edges[u, v]['style'] = 'dashed'
        else:
            diff.edges[u, v].update(G2.edges[u, v])


"""
Get graph metrics.
"""
# TODO hopefully will be more useful once extraneous waters have been removed programatically
# df_metrics = get_graph_centrality_metrics(Hbond_graphs)
# plot_metric_all_proteins(df_metrics["Eigenvalue Centrality"], limit=10, title="Eigenvalue Centrality")
# plot_metric_all_proteins(df_metrics["Degree Centrality"], limit=10, title="Degree Centrality")
# plot_metric_all_proteins(df_metrics["Betweenness Centrality"], limit=10, title="Betweenness Centrality")


"""
Plot in 3D
"""
# other_prot = '3QTW'
# prots = [apo_prot, other_prot]
# fig = plt.figure( figsize=(12,6) )
# ax = fig.add_subplot(1, 1, 1, projection='3d')

# legend_elements = draw_graph_3d(Hbond_graphs[apo_prot], ax) # plot graph    
# legend_elements += plot_backbone_3d(backbone_coordinates[apo_prot], ax, SS_POSITIONS) # plot backbone    
# # legend_elements += [    
        # # Line2D([0], [1], linewidth=1, linestyle='--', color='orange', label=f'Bonds broken {bonds_broken}', alpha=0.4),
        # # Line2D([0], [1], linewidth=1, linestyle='-', color='orange', label=f'Bonds formed {bonds_formed}', alpha=0.4),
# # ]

# ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.05, 1), title="Graph elements")
# plt.title(f"Symmetric difference between Hbond graphs of {apo_prot} and {other_prot}") 
# plt.tight_layout()
# plt.show()


"""
Plot in 2D
"""
# fig, ax  = plt.subplots( figsize=(12, 6) )

# legend_elements = draw_graph_2d(diff, ax) # plot graph        
# # legend_elements += plot_backbone_2d(backbone_coordinates[other_prot], ax, SS_POSITIONS) # plot backbone
# legend_elements += [    
        # Line2D([0], [1], linewidth=1, linestyle='--', color='orange', label=f'Bonds broken', alpha=0.4),
        # Line2D([0], [1], linewidth=1, linestyle='-', color='orange', label=f'Bonds formed', alpha=0.4),
# ]
# ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.05, 1), title="Graph elements")
# plt.title(f"Symmetric difference between Hbond graphs of {apo_prot} and {other_prot}") 
# plt.tight_layout()
# plt.show()