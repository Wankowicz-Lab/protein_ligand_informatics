# Cluster GEDs with Louvain

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import networkx as nx
import community

import matplotlib.pyplot as plt
import numpy as np


ged_ligands = pd.read_csv('./output/lig_adj_comp_no_dupes.txt', sep=' ', header=None)
ged_ligands.columns = ['ligand_1', 'ligand_2', 'similarity']

ged_ligands.replace(0, 4000, inplace=True)
ged_ligands['ligand_1'] = ged_ligands['ligand_1'].str.replace('lig_adj_', '').str.replace('.json', '')
ged_ligands['ligand_2'] = ged_ligands['ligand_2'].str.replace('lig_adj_', '').str.replace('.json', '')


################################################################
### LOUVAIN

G = nx.Graph()

for i, row in ged_ligands.iterrows():
    dist = row['similarity']
    if dist < 4000:
        weight = 1.0 / (dist + 1.0)
        G.add_edge(row['ligand_1'], row['ligand_2'], weight=weight)

partition = community.best_partition(G, weight='weight')
num_communities = len(set(partition.values()))
print(f"Detected {num_communities} communities by Louvain.")


pos = nx.spring_layout(G, k=0.1, seed=42)
communities = list(partition.values())

nx.draw_networkx_nodes(
    G, pos,
    node_color=communities,
    cmap=plt.cm.rainbow,
    node_size=100
)
nx.draw_networkx_edges(G, pos, alpha=0.5)
nx.draw_networkx_labels(G, pos, font_size=6)
plt.axis("off")
plt.title("Louvain Communities")
plt.show()
