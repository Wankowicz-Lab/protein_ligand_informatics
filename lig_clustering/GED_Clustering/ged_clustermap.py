# Visualize GEDs

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

ged_ligands = pd.read_csv('./output/lig_adj_comp_no_dupes.txt', sep=' ', header=None)
ged_ligands.columns = ['ligand_1', 'ligand_2', 'similarity']

ged_ligands.replace(0, 4000, inplace=True)
ged_ligands['ligand_1'] = ged_ligands['ligand_1'].str.replace('lig_adj_', '').str.replace('.json', '')
ged_ligands['ligand_2'] = ged_ligands['ligand_2'].str.replace('lig_adj_', '').str.replace('.json', '')

print(ged_ligands.sort_values(by='similarity', ascending=True)[ged_ligands['similarity']<300])

similarity_matrix = ged_ligands.pivot(
    index='ligand_1',
    columns='ligand_2',
    values='similarity'
).fillna(0)

sns.set_theme(context='notebook', style='white')
g = sns.clustermap(
    similarity_matrix,
    cmap='coolwarm',
    figsize=(10, 8),
    row_cluster=True,
    col_cluster=True,
    metric='euclidean',
    method='average'
)
plt.show()