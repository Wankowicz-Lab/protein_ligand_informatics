# TONNAR CHECKING ISOMORPHISM
# Result: Some graphs that were not physically the same returned to
#   be isomorphic, indicating that our original graph creation scheme may be inadequate or missing something.



import torch
import glob
import os
import networkx as nx
import numpy as np
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
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
    allLigFiles = get_files()
    for testFile in allLigFiles:
        print("Loading", testFile)
        load_graph(testFile)
    
    feature_list = []
    valid_ligand_names = []
    
    print(ligand_graphs)

graph1 = ligand_graphs[0]

a = 0
j = 0

for graph in ligand_graphs:
    j = j + 1
    c = 0
    for graph2 in ligand_graphs:
        c = c + 1
        if c < 100: 
            if nx.could_be_isomorphic(graph, graph2) and graph != graph2:
                a = a + 1
                print(ligand_names[j], ligand_names[c], c, j, "is", "Graph is isomorphic")
    

print(len(ligand_graphs), "total")
print(a, "isomorphic")