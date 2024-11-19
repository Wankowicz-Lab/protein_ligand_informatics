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


    # Error fixing "Edge tuple {e} must be a 2-tuple or 3-tuple"
    if edges.dim() == 1:
        num_edges = edges.numel() // 2
        edges = edges.view(num_edges, 2)
    elif edges.shape[0] == 2:
        edges = edges.t()
    else:
        edges = edges

    # ignore empty graphs
    num_nodes = nodes.size(0)
    if num_nodes > 0: 
        G = nx.Graph()
        G.add_nodes_from(range(num_nodes))
        
        edge_list = edges.tolist()
        G.add_edges_from(edge_list)
        
        # importing features into nx graph
        node_features = nodes.numpy()
        for i in range(num_nodes):
            G.nodes[i]['feature'] = node_features[i]
        
        edge_features = edge_features.numpy()
        for idx, (u, v) in enumerate(edge_list):
            G.edges[u, v]['feature'] = edge_features[idx]
        
        ligand_graphs.append(G)

def get_files():
    return glob.glob(os.path.join(ROOT_DIR, "lig_graph_*.pt"))


def compute_graph_features(G):
    features = {}

    # find features/metrics about the graph to find similarities numerically.

    features['num_nodes'] = G.number_of_nodes()
    features['num_edges'] = G.number_of_edges()
    degrees = [degree for node, degree in G.degree()]
    features['avg_degree'] = np.mean(degrees) if degrees else 0
    try:
        features['degree_assortativity'] = nx.degree_assortativity_coefficient(G)
    except Exception:
        features['degree_assortativity'] = 0
    features['avg_clustering'] = nx.average_clustering(G)
    features['transitivity'] = nx.transitivity(G)
    if nx.is_connected(G):
        features['diameter'] = nx.diameter(G)
    else:
        features['diameter'] = 0
    laplacian = nx.normalized_laplacian_matrix(G).todense()
    eigenvalues = np.linalg.eigvals(laplacian)
    eigenvalues = np.sort(eigenvalues)
    eigenvalues = eigenvalues[:5] if len(eigenvalues) >= 5 else np.pad(eigenvalues, (0, 5 - len(eigenvalues)), 'constant')
    for i, eigenvalue in enumerate(eigenvalues):
        features[f'eigenvalue_{i}'] = eigenvalue.real
    
    return features


if __name__ == "__main__":
    allLigFiles = get_files()
    for testFile in allLigFiles:
        print("Loading", testFile)
        load_graph(testFile)
    
    # put all computed metrics of graph into an array
    feature_list = []
    valid_ligand_names = []
    for idx, G in enumerate(ligand_graphs):
        features = compute_graph_features(G)
        feature_list.append(features)
        valid_ligand_names.append(ligand_names[idx])
    

    # convert the array into a pandas dataframe
    df = pd.DataFrame(feature_list)
    df.fillna(0, inplace=True) # fix "The condensed distance matrix must contain only finite values" by filling NaN with 0

    # standardize the data
    scaler = StandardScaler()
    X = scaler.fit_transform(df.values)
    
    # compute the distance matrix, giving us the similarities between the ligands
    distance_matrix = pdist(X, metric='euclidean')
    linkage_matrix = linkage(distance_matrix, method='complete') # ward

    # plot the dendrogram from the linkage matrix

    plt.figure(figsize=(12, 8))
    dendrogram(
        linkage_matrix,
        labels=valid_ligand_names,
        leaf_rotation=90,
        leaf_font_size=10,
    )
    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('Ligand')
    plt.ylabel('Distance')
    plt.tight_layout()
    plt.show()