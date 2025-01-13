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
ligand_files = []
ligand_names = []

def load_graph_from_file(graph_name):
    graph = torch.load(graph_name)
    
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
    return G

def get_files():
    return glob.glob(os.path.join(ROOT_DIR, "lig_graph_*.pt"))

def calculate_ged(i, j):
    # Load the graphs from the files
    graph_file_i = ligand_files[i]
    graph_file_j = ligand_files[j]

    g1 = load_graph_from_file(graph_file_i)
    g2 = load_graph_from_file(graph_file_j)

    ged = nx.graph_edit_distance(g1, g2, timeout=30)
    print("GED between", i, j, "is", ged)
    if ged is None:
        ged = 0
    return i, j, ged

if __name__ == "__main__":
    allLigFiles = get_files()
    for testFile in allLigFiles:
        print("Processing", testFile)
        ligand_files.append(testFile)
        dataset_name = testFile.split("_")[-1].split(".pt")[0]
        ligand_names.append(dataset_name)
    
    n_graphs = len(ligand_files)
    ged_matrix = np.zeros((n_graphs, n_graphs))

    try:
        df = pd.read_csv("ged_matrix_256c.csv", index_col=0)
        ged_matrix = df.values
        print("Loaded matrix of size", ged_matrix.shape)
    except FileNotFoundError:
        print("GED matrix not found, creating new matrix")
        pass

    batch_size = 256
    with ProcessPoolExecutor(max_workers=batch_size) as executor:
        tasks = [(i, j) for i in range(n_graphs) for j in range(i + 1, n_graphs)]

        for batch_start in range(0, len(tasks), batch_size):
            batch_tasks = tasks[batch_start:batch_start + batch_size]
            maybeI, maybeJ = batch_tasks[-1]
            if ged_matrix[maybeI, maybeJ] != 0:
                print(f"Skipping batch {batch_start // batch_size + 1}")
                continue
            futures = {
                executor.submit(calculate_ged, i, j): (i, j)
                for i, j in batch_tasks
            }
            print(f"Submitted batch {batch_start // batch_size + 1}")
            for future in as_completed(futures):
                i, j = futures[future]
                try:
                    idx_i, idx_j, ged = future.result()
                    ged_matrix[idx_i, idx_j] = ged
                    ged_matrix[idx_j, idx_i] = ged
                    print(f"Completed: ({idx_i}, {idx_j}) - GED: {ged}")
                except Exception as e:
                    print(f"Error calculating GED for ({i}, {j}): {e}")
            
            df = pd.DataFrame(ged_matrix, columns=ligand_names[:n_graphs], index=ligand_names[:n_graphs])
            df.to_csv("ged_matrix_256c.csv", index=True)

    df = pd.DataFrame(ged_matrix, columns=ligand_names[:n_graphs], index=ligand_names[:n_graphs])
    df.to_csv("ged_matrix_256c.csv", index=True)
