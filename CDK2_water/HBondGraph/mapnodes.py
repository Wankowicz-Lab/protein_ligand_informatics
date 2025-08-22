"""
Takes two H bond graphs and maps nodes between them based on distance.
"""

import numpy as np
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
from HBondGraph import properties


# UNIVERSAL CONSTANTS
DEFAULT_MAX_MAP_DISTANCE = 1.75

def map_nodes(graph_ref, node_list_ref, graph_align, node_list_align, distance_cutoff=DEFAULT_MAX_MAP_DISTANCE):
    """
    Maps nodes in node_list_align from graph_align to nodes in node_list_ref from graph_ref
    """
    ref_coordinates = np.vstack( properties.get_attr_from_nodes(graph_ref, node_list_ref, 'pos' ) )
    align_coordinates = np.vstack( properties.get_attr_from_nodes(graph_align, node_list_align, 'pos') )
    
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

