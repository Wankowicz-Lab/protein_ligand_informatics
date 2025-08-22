"""
Takes a PDB file with protein and water coordinates.
Creates a graph with vertices of water/AA residues and edges if they have an H bond.
"""

import numpy as np
import networkx as nx
from scipy.spatial.transform import Rotation

# UNIVERSAL CONSTANTS
AA=['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
WATER=['HOH', 'WAT', 'H2O']  

HBOND_MAX_DIST=3.5
HBOND_MIN_DIST=2.4

# Default hydrogen bond donor and acceptor definitions from Mingbin
DEFAULT_DONOR_ATOMS = ['N', 'ND1', 'ND2', 'NE', 'NE1', 'NE2', 'NH1', 'NH2', 'NZ', 'OG', 'OG1', 'OH', 'SG']
DEFAULT_ACCEPTOR_ATOMS = ['ND1', 'NE2', 'O', 'OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1', 'OH', 'SD', 'SG']


"""
Atom and coordinate loading
"""
def get_filtered_atoms(structure, residue_target=None, atom_target=None):
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
    
    for model in structure:
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
    

def get_transformation(align_atoms, ref_atoms):
    """
    Takes align_atoms and finds the translation, rotation pair to align to ref_atoms.
    """
    align_coords = get_coordinates_from_atoms(align_atoms)
    ref_coords = get_coordinates_from_atoms(ref_atoms)

    centered_align = align_coords - np.mean(align_coords, axis=0)
    centered_ref = ref_coords - np.mean(ref_coords, axis=0)

    translation = -np.mean(align_coords, axis=0)
    rotation, _ = Rotation.align_vectors(centered_ref, centered_align)
    
    return translation, rotation
    
    
"""
Graph construction.
"""
def graph_from_structure(new_structure, ref_structure):
    """
    Takes in the structure from Bio.PDB to analyze and a reference structure for alignment.
    Outputs a graph with all H bonds.
    """
    # Create an empty graph
    Hbond_graph = nx.Graph()

    # get atoms and coordinates
    water_atoms = get_filtered_atoms(new_structure, residue_target = WATER, atom_target = ["O"])
    donor_atoms = get_filtered_atoms(new_structure, residue_target = AA, atom_target = DEFAULT_DONOR_ATOMS)
    acceptor_atoms = get_filtered_atoms(new_structure, residue_target = AA, atom_target = DEFAULT_ACCEPTOR_ATOMS)    
    backbone_atoms = get_filtered_atoms(new_structure, residue_target = AA, atom_target = ["CA"])    
    
    ref_backbone = get_filtered_atoms(ref_structure, residue_target = AA, atom_target = ["CA"])  
    translation, rotation = get_transformation(backbone_atoms, ref_backbone)
    
    water_coordinates = get_coordinates_from_atoms(water_atoms, translation=translation, rotation=rotation)     
    donor_coordinates = get_coordinates_from_atoms(donor_atoms, translation=translation, rotation=rotation)    
    acceptor_coordinates = get_coordinates_from_atoms(acceptor_atoms, translation=translation, rotation=rotation)  
    backbone_coordinates = get_coordinates_from_atoms(backbone_atoms, translation=translation, rotation=rotation) 
    
    # Add nodes to the graph
    for idx, coord in enumerate(water_coordinates):
        w_tag = ("w", water_atoms[idx].get_parent().get_id()[1], water_atoms[idx].name)
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
        w_tag = ("w", water_atoms[idx].get_parent().get_id()[1], f"{water_atoms[idx].name}")
        
        water_water_distances = np.linalg.norm(coord - water_coordinates, axis=1) # water to water
        connected_edges = np.nonzero( 
            np.logical_and(water_water_distances < HBOND_MAX_DIST, water_water_distances > HBOND_MIN_DIST) 
            )[0]
        Hbond_graph.add_edges_from(
            [ (w_tag, ("w", water_atoms[jdx].get_parent().get_id()[1], f"{water_atoms[jdx].name}"))
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
        
    return Hbond_graph, backbone_coordinates
    
   
