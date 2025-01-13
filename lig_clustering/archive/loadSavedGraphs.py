# For general use to perform various tests

import torch
import glob;

ROOT_DIR = "./built_graphs"

ligand_graphs = []

def load_graph(graph_name):
    graph = torch.load(graph_name)
    dataset_name = graph_name.split("/")[-1].split("_")[3].split('\\')[-1].split(".pt")[0]
    
    nodes = graph[0]
    edges = graph[1]
    edge_features = graph[2]

    ligand_graphs.append([dataset_name, nodes, edges, edge_features])




def get_files():
    return glob.glob(ROOT_DIR + "/lig_graph_*.pt")


if __name__ == "__main__":
    allLigFiles = get_files()
    for testFile in allLigFiles:
        print("Loading", testFile)
        load_graph(testFile)
    print(ligand_graphs)