# cn1264 
#  job 684076 (failed)
# job 690036
# job 690560
# job 699802 (1610)

import torch
import glob
import os
import networkx as nx
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed


ROOT_DIR = "./built_graphs"
ligand_graphs = []
ligand_names = []

def load_graph(graph_name):
    graph = torch.load(graph_name)
    dataset_name = graph_name.split("_")[-1].split(".pt")[0]
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


def calculate_ged(i, j, g1, g2):
    ged = nx.graph_edit_distance(g1, g2, timeout=30)
    print("GED between", i, j, "is", ged)
    if ged is None:
        ged = 0
    return ged


if __name__ == "__main__":
    # Load all ligand graphs
    allLigFiles = get_files()
    for testFile in allLigFiles:
        print("Loading", testFile)
        load_graph(testFile)
    
    print(f"Total graphs loaded: {len(ligand_graphs)}")
    
    n_graphs = len(ligand_graphs)
    ged_matrix = np.zeros((n_graphs, n_graphs))

    # get existing data from csv
    try:
        df = pd.read_csv("ged_matrix.csv", index_col=0)
        ged_matrix = df.values
        print ("Loaded matrix of size ", ged_matrix.shape)
    except FileNotFoundError:
        print("GED matrix not found, creating new matrix")
        pass
        
    with ProcessPoolExecutor() as executor:
        batch_size = 32
        tasks = [(i, j) for i in range(n_graphs) for j in range(i + 1, n_graphs)]

        for batch_start in range(0, len(tasks), batch_size):
        
            batch_tasks = tasks[batch_start:batch_start + batch_size]
            maybeI, maybeJ = batch_tasks[-1]
            if (ged_matrix[maybeI, maybeJ] != 0):
                print(f"Skipping batch {batch_start // batch_size + 1}")
                continue
            futures = {
                executor.submit(calculate_ged, i, j, ligand_graphs[i], ligand_graphs[j]): (i, j)
                for i, j in batch_tasks
            }
            print(f"Submitted batch {batch_start // batch_size + 1}")
            for future in as_completed(futures):
                i, j = futures[future]
                try:
                    ged = future.result()
                    ged_matrix[i, j] = ged
                    ged_matrix[j, i] = ged
                    print(f"Completed: ({i}, {j}) - GED: {ged}")
                except Exception as e:
                    print(f"Error calculating GED for ({i}, {j}): {e}")
            
            df = pd.DataFrame(ged_matrix, columns=ligand_names[:n_graphs], index=ligand_names[:n_graphs])
            df.to_csv("ged_matrix.csv", index=True)

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
