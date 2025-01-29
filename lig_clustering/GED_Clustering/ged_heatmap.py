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

ligands = pd.concat([ged_ligands['ligand_1'], ged_ligands['ligand_2']]).unique()
similarity_matrix = pd.DataFrame(np.zeros((len(ligands), len(ligands))), index=ligands, columns=ligands)

for _, row in ged_ligands.iterrows():
    similarity_matrix.at[row['ligand_1'], row['ligand_2']] = row['similarity']
    similarity_matrix.at[row['ligand_2'], row['ligand_1']] = row['similarity']

plt.figure(figsize=(8, 6))
sns.heatmap(similarity_matrix)
plt.show()

plt.figure(figsize=(8, 6))
similarity_matrix = ged_ligands.pivot('ligand_2', 'ligand_1', 'similarity').fillna(0)
sns.clustermap(similarity_matrix, 
               cmap='coolwarm', 
               figsize=(10, 8), 
               col_cluster=True)
plt.show()