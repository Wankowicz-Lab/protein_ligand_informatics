"""
Takes Hbond graphs and plots them in 2D or 3D projections.
"""

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from adjustText import adjust_text

import networkx as nx
import numpy as np
from sklearn.decomposition import PCA

# UNIVERSAL CONSTANTS
SS_POSITIONS= [ (4,12), (17,23), (29,35), (45,56), (66,72), (75,81),
        (85,94), (100,121), (129,131), (133,135), (141,143), (147,152), (170,175),
        (182,199), (207,220), (229,233), (247,252), (256,267), (276,287) ] # hardcoded based on PDB entry    
    

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


def plot_backbone_2d(coordinates, ax, ss_positions=SS_POSITIONS): 
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


    
def plot_backbone_3d(coordinates, ax, ss_positions=SS_POSITIONS):
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
