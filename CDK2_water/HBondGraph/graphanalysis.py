"""
Calculates specific metrics from a group of H-bond graphs
Very much in progress.
"""

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from HBondGraph import properties

    
def get_graph_centrality_metrics(Hbond_graph_dict):
    """
    Takes a dict of graphs with (protein name, Hbond_graph) as key, value pairs.
    Outputs a dataframe with centrality metrics of the water.
    """
    metrics = ["Eigenvalue Centrality", "Degree Centrality", "Betweenness Centrality"]
    df = pd.DataFrame(index = Hbond_graph_dict.keys(), columns = metrics)
    for prot, graph in Hbond_graph_dict.items():
        water_nodes = properties.get_water_nodes(graph)
        
        eigenvector_centrality = nx.eigenvector_centrality(graph, tol=1e-4, max_iter=1000)
        degree_centrality = nx.degree_centrality(graph)
        betweenness_centrality = nx.betweenness_centrality(graph, normalized=True)
        
        df.at[prot, "Eigenvalue Centrality"] = np.sort( [eigenvector_centrality[n] for n in water_nodes])[::-1]
        df.at[prot, "Degree Centrality"] = np.sort( [degree_centrality[n] for n in water_nodes])[::-1]
        df.at[prot, "Betweenness Centrality"] = np.sort( [betweenness_centrality[n] for n in water_nodes])[::-1]

    return df
    
    
def plot_metric_all_proteins(data_series, limit=None, title=None):
    """
    Plotting function for a metric for all proteins. Data is a df series of the metric, one for each protein.
    """
    
    fig, axes = plt.subplots(9, 7, figsize = (18, 12), dpi=100, sharex=True, sharey=True) # 63 proteins
    for idx, (data, ax) in enumerate( zip(data_series, axes.flat) ):
        prot = data_series.index[idx]
        if limit is None:
            limit = len(data) 
        ax.plot(data[:limit], "o")
        ax.set_title(prot, fontsize=8)
        ax.tick_params(labelsize=6)
        
    if title:
        fig.suptitle(title)

    plt.subplots_adjust(top=0.97)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.show()
