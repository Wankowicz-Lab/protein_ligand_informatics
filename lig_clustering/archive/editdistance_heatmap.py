
import torch
import glob
import os
import networkx as nx
import numpy as np
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd

ROOT_DIR = "./built_graphs"
ligand_graphs = []
ligand_names = []

def load_graph(graph_name):
    graph = torch.load(graph_name)
    dataset_name = graph_name.split("/")[-1].split("_")[3].split('\\')[-1].split(".pt")[0]
    ligand_names.append(dataset_name)
    
    nodes = graph[0]
    edges = graph[1]
    edge_features = graph[2]

    if edges.dim() == 1:
        num_edges = edges.numel() // 2
        edges = edges.view(num_edges, 2)
    elif edges.shape[0] == 2:
        edges = edges.t()
    else:
        edges = edges

    G = nx.Graph()
    num_nodes = nodes.size(0)
    if num_nodes == 0:
        G.add_node(0)
    else:
        G.add_nodes_from(range(num_nodes))
    
    edge_list = edges.tolist()
    G.add_edges_from(edge_list)
    
    node_features = nodes.numpy()
    for i in range(num_nodes):
        G.nodes[i]['feature'] = node_features[i]
    
    edge_features = edge_features.numpy()
    for idx, (u, v) in enumerate(edge_list):
        G.edges[u, v]['feature'] = edge_features[idx]
    
    ligand_graphs.append(G)

def get_files():
    return glob.glob(os.path.join(ROOT_DIR, "lig_graph_*.pt"))

if __name__ == "__main__":
    # Load all ligand graphs
    allLigFiles = get_files()
    for testFile in allLigFiles:
        print("Loading", testFile)
        load_graph(testFile)
    
    print(f"Total graphs loaded: {len(ligand_graphs)}")
    
    n_graphs = len(ligand_graphs)
    ged_matrix = np.zeros((n_graphs, n_graphs))
    
    
    for i in range(n_graphs):
        print("Nodes in graph", i, len(ligand_graphs[i].nodes))
        print("Edges in graph", i, len(ligand_graphs[i].edges))
        for j in range(n_graphs - i):
            real_j = n_graphs - j - 1
            print("Working on ", i, real_j)
            if i != real_j:
                ged = nx.graph_edit_distance(ligand_graphs[i], ligand_graphs[real_j], timeout=2.5)
                if (ged == None):
                    ged = 0
                ged_matrix[i, real_j] = ged
                ged_matrix[real_j, i] = ged
                print("Working on matrix entry", i, real_j, "GED:", ged)
    
    # export to csv
    df = pd.DataFrame(ged_matrix, columns=ligand_names[:n_graphs], index=ligand_names[:n_graphs])
    df.to_csv("ged_matrix.csv", index=True)


    # plt.figure(figsize=(10, 8))
    # sns.heatmap(ged_matrix, annot=False, fmt=".2f", 
    #             xticklabels=ligand_names, yticklabels=ligand_names, cmap="viridis")
    # plt.title("Graph Edit Distance Heatmap")
    # plt.xlabel("Ligand Graphs")
    # plt.ylabel("Ligand Graphs")
    # plt.xticks(rotation=45, ha='right')
    # plt.tight_layout()
    # plt.show()
