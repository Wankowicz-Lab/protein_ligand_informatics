import argparse
import networkx as nx
import pandas as pd

def get_args():
    parser = argparse.ArgumentParser(description="edit distance 2 graphs")
    parser.add_argument("graph_1", type=str, help="first graph")
    parser.add_argument("graph_2", type=str, help="second graph")
    args = parser.parse_args()
    return args.graph_1, args.graph_2


def load_nx_adj_mx_from_json(graph):
    with open(graph, "r") as f:
        adj = f.read()
        adj = eval(adj)
        df = pd.DataFrame(adj)
        G = nx.from_pandas_adjacency(df)
    return G

def append_to_text_file(content):
    with open("ged_task_output_2.txt", "a") as f:
        f.write(content + "\n")
        f.close()

if __name__ == "__main__":
    args = get_args()
    graph_1 = args[0]
    graph_2 = args[1]
    G1 = load_nx_adj_mx_from_json(graph_1)
    G2 = load_nx_adj_mx_from_json(graph_2)
    ged = nx.graph_edit_distance(G1, G2, timeout=40)
    if ged is None:
        ged = 0
    graph_1 = graph_1.split("/")[-1]
    graph_2 = graph_2.split("/")[-1]
    append_to_text_file( graph_1 + " " + graph_2 + " " + str(ged))
    print(graph_1, graph_2, ged)