"""
Examples of how to use the modules to produce figures.
"""

from Bio import PDB
from rdkit import Chem

import numpy as np
import networkx as nx

from sklearn.decomposition import PCA
from sklearn.cluster import DBSCAN
from sklearn.manifold import MDS
from skbio.stats.distance import mantel

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from adjustText import adjust_text

import glob
import pathlib
import sys
import requests
import csv

from HBondGraph import create, mapnodes, ligand, plot, properties, graphanalysis

# CONSTANTS
DATA_FOLDER = "../PDB_Data_cdk2/"
HBOND_MAX_DIST = 3.5

""" 
Load all files ending with "_cleaned.pdb" from DATA_FOLDER/cdk2_expanded.
Extract hydrogen bonds from coordinate data and represent as a graph.
REQUIRED
"""

print(f"Parsing pdb files")
pdb_path = pathlib.Path(DATA_FOLDER)
pdb_files = list(pdb_path.glob('./cdk2_expanded/*_cleaned.pdb')) 
parser = PDB.PDBParser(QUIET=True)

apo_prot = '1PW2'
all_prots = []
Hbond_graphs = {}
backbone_coordinates = {}
ligs = {}
    
for idx, file_name in enumerate(pdb_files):
    prot = file_name.stem.split('_')[0]
    structure = parser.get_structure(prot, file_name)
    
    # align to first structure
    if idx == 0:
        ref_structure = structure.copy()    
    graph, backbone = create.graph_from_structure(structure, ref_structure)
    
    all_prots.append(prot)      
    Hbond_graphs[prot] = graph
    backbone_coordinates[prot] = backbone
    
    if prot != apo_prot:
        lig = ligand.get_ligand_from_structure(structure)
        ligs[prot] = lig
    

# downloading the ligands takes a while, so do all at once only once
# print("Downloading ligand molecules")
# ligs = {prot: ligand.residue_to_mol(res) for prot, res in ligs.items()}

"""
Plot in 3D
OPTIONAL
"""
# fig = plt.figure( figsize=(12,6) )
# ax = fig.add_subplot(1, 1, 1, projection='3d')

# legend_elements = plot.draw_graph_3d(Hbond_graphs[apo_prot], ax) # plot graph    
# legend_elements += plot.plot_backbone_3d(backbone_coordinates[apo_prot], ax) # plot backbone    

# ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.05, 1), title="Graph elements")
# plt.title(f"Hbonds in {apo_prot}") 
# plt.tight_layout()
# plt.show()


"""
Plot in 2D
OPTIONAL
"""
# fig, ax  = plt.subplots( figsize=(12, 6) )
# legend_elements = plot.draw_graph_2d(Hbond_graphs[apo_prot], ax) # plot graph        
# legend_elements += plot.plot_backbone_2d(backbone_coordinates[apo_prot], ax) # plot backbone
# ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.05, 1), title="Graph elements")
# plt.title(f"Hbonds in {apo_prot}") 
# plt.tight_layout()
# plt.show()


"""
Identify conserved waters.
REQUIRED
"""

print(f"Mapping waters")
distance_cutoff = HBOND_MAX_DIST/2 # max distance for two waters to map together, estimated by inflection point in number of water pairs as a function of this distance, leaves 27 waters conserved in all structures
num_prot = len(all_prots)
mapped_water_counts = np.zeros([num_prot, num_prot])

# loop through all pairs of proteins and find the mappable waters
for idx, ref_prot in enumerate(all_prots): 
    ref_waters = properties.get_water_nodes(Hbond_graphs[ref_prot])
    
    for jdx, test_prot in enumerate(all_prots): 
        test_waters = properties.get_water_nodes(Hbond_graphs[test_prot])
        ref_waters_mapped, test_waters_mapped = mapnodes.map_nodes(Hbond_graphs[ref_prot], ref_waters, 
                                                                   Hbond_graphs[test_prot], test_waters, distance_cutoff)
        mapped_water_counts[idx, jdx] = len(ref_waters_mapped)
        
        # update the freq trait based on mapping
        for w in ref_waters_mapped:
            Hbond_graphs[ref_prot].nodes[w]['freq'] += 1


# # count the number of H bonds to residues from most conserved waters
# residues_bonded = []
# for prot in all_prots:
    # most_conserved_waters = properties.get_most_conserved_waters(Hbond_graphs[prot])
    # max_freq = properties.get_attr_from_nodes(Hbond_graphs[prot], most_conserved_waters, 'freq')[0]
    # for w in most_conserved_waters:
        # residues_bonded += [v[1] for u, v in Hbond_graphs[prot].edges(w) if Hbond_graphs[prot].edges[u, v]['bond_type'] == 'wp']
        
# residues_near_conserved_waters, counts = np.unique(residues_bonded, return_counts = True)
# sorted_ind = np.argsort(counts)[::-1] # descending sort
# residues_near_conserved_waters, counts = residues_near_conserved_waters[sorted_ind], counts[sorted_ind]


# # count the number of waters in each protein rather than over all
# counts = counts/max_freq

# file = '../Figures/residues_with_H_bond_to_conserved_waters.txt'
# with open(file, 'w', newline='') as f:
    # f.write("Note: only waters that can be mapped between the all structures were counted.\n\n")
    # f.write("Residues with H-bonds to most conserved waters.\n")
    # f.write(f"{residues_near_conserved_waters}\n\n")
    # f.write("Number of H-bonds to respective residue.\n")
    # f.write(f"{np.round(counts, 2)}")    

"""
Write water conservation to PDB file.
OPTIONAL.
"""        
# # output a new pdb with the freq count in the b-factor column to visualize in pymol
# for ref_prot in all_prots: 
    # properties.edit_water_b_factor( input_file = f"{DATA_FOLDER}/cdk2_expanded/{ref_prot}_cleaned.pdb", 
        # output_file = f"{DATA_FOLDER}/cdk2_edited/{ref_prot}.pdb", 
        # b_factor_dict = nx.get_node_attributes(Hbond_graphs[ref_prot], 'freq') )

"""
Count number of H bonds per water and compare to conservation.
OPTIONAL
"""
# # Aggregate waters over all structures, get conservation level of water, compare to number of H bonds to protein
# num_neighbors_dict = {i+1:[] for i in range(num_prot)}
# for prot in all_prots:
    # waters = properties.get_water_nodes(Hbond_graphs[prot])
    # freqs = properties.get_attr_from_nodes(Hbond_graphs[prot], waters, 'freq')
    # for w, f in zip(waters, freqs):
        # num_bonds = len( [_ for _, target in Hbond_graphs[prot].edges(w) if target[0] != 'w'] )
        # num_neighbors_dict[f].append(num_bonds)

# group_size = 7
# num_steps = num_prot // group_size
# means = []
# for i in range(num_steps):
    # accum = []
    # for j in range(i*group_size, (i+1)*group_size):
        # accum += num_neighbors_dict[j+1] 
    # means.append(np.mean(accum))

# # Create bar plot
# fig, ax = plt.subplots(figsize=(12, 6))
# bars = ax.bar(range(num_steps), means, alpha=0.7)

# # Annotate bars with values
# barxs = []
# for bar in bars:
    # height = bar.get_height()
    # ax.text(bar.get_x() + bar.get_width()/2, height + 0.05, f'{height:.1f}',
            # ha='center', va='bottom', fontsize=10)
            
# ax.set_xlabel('Frequency of water')
# ax.set_xlim([-0.5, num_steps-0.5])
# ax.set_xticks(range(num_steps))
# ax.set_xticklabels( [f"{i*group_size + 1} to {(i+1)*group_size}" for i in range(num_steps)] )
# ax.set_ylim([0,3])
# ax.set_ylabel('Average number of H-bonds with protein')
# plt.tight_layout()
# plt.show()


"""
Calculate graph difference for the set of proteins using common nodes mapped with map_nodes.
Compare with ligand molecular weight. Write changes and associated residues to a CSV file.
OPTIONAL
"""  
# print("Calculating graph differences")
# bond_changes = []
# wp_bond_change_residues = []
# for other_prot in all_prots[1:]:
    # w1 = properties.get_water_nodes(Hbond_graphs[apo_prot])
    # w2 = properties.get_water_nodes(Hbond_graphs[other_prot]) 
    
    # # remap water nodes to make node sets equivalent
    # w1, w2 = mapnodes.map_nodes(Hbond_graphs[apo_prot], w1, Hbond_graphs[other_prot], w2, distance_cutoff=distance_cutoff) 
    # subnodes1 = properties.get_protein_nodes(Hbond_graphs[apo_prot]) + w1
    # subnodes2 = properties.get_protein_nodes(Hbond_graphs[other_prot]) + w2
    # G1 = Hbond_graphs[apo_prot].subgraph(subnodes1)
    # G2 = Hbond_graphs[other_prot].subgraph(subnodes2)
    
    # # relabel mapped nodes
    # G2 = nx.relabel_nodes( G2, dict(zip( w2, w1 )) )        
    
    # # calculate graph difference and restore properties. Count bond changes
    # # broken if present in apo but not in ligand bound
    # # formed if not present in apo but in ligand bound
    # diff = nx.symmetric_difference(G1, G2)
    # pp_bonds_broken = 0
    # pp_bonds_formed = 0
    # wp_bonds_broken = 0
    # wp_bonds_formed = 0
    # for n in diff.nodes:
        # diff.nodes[n].update(G1.nodes[n])
        
    # for u, v in diff.edges:            
        # if G1.has_edge(u, v):
            # diff.edges[u, v].update(G1.edges[u, v])
            # diff.edges[u, v]['style'] = 'dashed'
            # if diff.edges[u, v]['bond_type'] == 'wp': 
                # wp_bonds_broken += 1
            # elif diff.edges[u, v]['bond_type'] == 'pp': 
                # pp_bonds_broken += 1
        # else:
            # diff.edges[u, v].update(G2.edges[u, v])
            # if diff.edges[u, v]['bond_type'] == 'wp': 
                # wp_bonds_formed += 1
            # elif diff.edges[u, v]['bond_type'] == 'pp': 
                # pp_bonds_formed += 1
                
    # # get all residues that have a water-protein H bond change
    # wp_bond_change_residues += [v[1] for u, v in diff.edges if diff.edges[u, v]['bond_type'] == 'wp']

    # bonds_broken = wp_bonds_broken + pp_bonds_broken
    # bonds_formed = wp_bonds_formed + pp_bonds_formed
    # bond_changes.append( (pp_bonds_broken, pp_bonds_formed, wp_bonds_broken, wp_bonds_formed) )


# # count the residues associated with H bond changes
# residues_near_changes, counts = np.unique(wp_bond_change_residues, return_counts = True)
# sorted_ind = np.argsort(counts)[::-1] # descending sort
# residues_near_changes, counts = residues_near_changes[sorted_ind], counts[sorted_ind]

# file = '../Figures/residues_affected_by_H_bond_changes.txt'
# with open(file, 'w', newline='') as f:
    # f.write("Note: only waters that can be mapped between the apo and holo structure were counted.\n\n")
    # f.write("Residues with H-bonds formed and broken with water across all ligands.\n")
    # f.write(f"{residues_near_changes}\n\n")
    # f.write("Number of H-bonds formed and broken with water across all ligands, respectively.\n")
    # f.write(f"{counts}")    
    

# # plot number of bond changes against molecular weight of ligand
# wp_bond_changes = [ (b[2] + b[3])/mapped_water_counts[0, idx+1] for idx, b in enumerate(bond_changes) ]
# mol_weights = [ Chem.rdMolDescriptors.CalcExactMolWt(l) for l in ligs.values() ]

# plt.figure( figsize = (12,8) )
# plt.plot(mol_weights, wp_bond_changes, "o")

# texts = []
# for wp, mw, label in zip(wp_bond_changes, mol_weights, all_prots[1:]):
    # texts.append( plt.text( mw, wp, label, fontsize=8 ) )
# adjust_text(texts)

# plt.xlabel("Ligand molecular weight")
# plt.ylabel("Water-protein H bonds formed and broken per shared water")
# plt.show()


# # write bond changes to CSV
# file = '../Figures/apo_to_holo_Hbond_change_comparison.csv'
# with open(file, 'w', newline='') as f:
    # f.write("# Using the set of totally conserved waters and the protein Hbond donors/acceptors"
        # "find the difference in the Hbond graph between the apo structure and a ligand bound structure.\n")  # comment
    # writer = csv.writer(f)
    # writer.writerow(['PDB', 'protein-protein broken', 'protein-protein formed', 'water-protein broken', 'water-protein formed'])  # header
    # for idx, other_prot in enumerate(all_prots[1:]):
        # pp_bonds_broken, pp_bonds_formed, wp_bonds_broken, wp_bonds_formed = bond_changes[idx]
        # writer.writerow([other_prot, pp_bonds_broken, pp_bonds_formed, wp_bonds_broken, wp_bonds_formed])


"""
Plot graph difference
OPTIONAL. Requires calculating graph difference.
"""
# fig = plt.figure( figsize=(12, 6) )
# ax  = fig.add_subplot( projection='3d' )
# legend_elements = plot.draw_graph_3d(diff, ax) # plot graph        
# legend_elements += [    
        # Line2D([0], [1], linewidth=1, linestyle='--', color='orange', label=f'Bonds broken {bonds_broken}', alpha=0.4),
        # Line2D([0], [1], linewidth=1, linestyle='-', color='orange', label=f'Bonds formed {bonds_formed}', alpha=0.4),
# ]
# ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.05, 1), title="Graph elements")
# plt.title(f"Symmetric difference between {apo_prot} and {all_prots[-1]}") 
# plt.tight_layout()
# plt.show()


"""
Get graph metrics.
OPTIONAL.
"""
# # TODO needs to be tested and validated
# print("Calculating and plotting graph metrics")
# df_metrics = graphanalysis.get_graph_centrality_metrics(Hbond_graphs)
# graphanalysis.plot_metric_all_proteins(df_metrics["Eigenvalue Centrality"], limit=10, title="Eigenvalue Centrality")
# graphanalysis.plot_metric_all_proteins(df_metrics["Degree Centrality"], limit=10, title="Degree Centrality")
# graphanalysis.plot_metric_all_proteins(df_metrics["Betweenness Centrality"], limit=10, title="Betweenness Centrality")


"""
Compare ligands and water distances.
OPTIONAL.
"""

# # calculate tanimoto dissimilarity metrics 
# lig_dissimilarities = np.zeros([num_prot-1, num_prot-1])
# water_dissimilarities = np.zeros_like(lig_dissimilarities)

# # +1 indexing since the ligand matrices have one less structure (apo) than the proteins
# for idx, ref_prot in enumerate(all_prots[1:]): 
    # mol_ref = ligs[ref_prot]

    # for jdx, test_prot in enumerate(all_prots[1:]): 
        # if ref_prot == test_prot: continue # do not map to self
        # mol_test = ligs[test_prot]
        
        # chem_dissimilarity = ligand.tanimoto_dissimilarity(mol_ref, mol_test)
        # lig_dissimilarities[idx, jdx] = chem_dissimilarity
        
        # # calculate the tanimoto = |intersection| / |union| of water atoms that can be mapped to at least one other structure
        # total_ref_waters = np.sum( np.array(properties.get_attr_from_nodes(Hbond_graphs[ref_prot], 'water', "freq")) > 1 )
        # total_test_waters = np.sum( np.array(properties.get_attr_from_nodes(Hbond_graphs[test_prot], 'water', "freq")) > 1 )
        # mapped_waters = mapped_water_counts[idx+1, jdx+1] 
        
        # water_dissimilarity = 1 - ( mapped_waters / ( total_ref_waters + total_test_waters - mapped_waters ) )
        # water_dissimilarities[idx, jdx] = water_dissimilarity 

# # remove repeats for correlation calculation and plotting
# mask_ind = np.triu_indices_from(lig_dissimilarities, k=1)
# lig_masked = lig_dissimilarities[mask_ind]
# water_masked = water_dissimilarities[mask_ind]

# spearmanr_coeff, mantel_pvalue, n = mantel(lig_masked, water_masked, method='spearman', permutations=1000)
# plt.scatter(lig_masked, water_masked, s=10)

# plt.xlabel("Tanimoto dissimilarity between ligands")
# plt.xlim([0,1])
# plt.ylabel("Tanimoto dissimilarity between conserved waters")
# plt.ylim([0,1])
# plt.title(f"All pairs of ligands")

# custom_labels = [f"Spearman correlation: {spearmanr_coeff:.3f}", f"Mantel test p-value: {mantel_pvalue:.4f}"]
# legend_handles = [Line2D([0], [0], color='none') for _ in custom_labels]
# plt.legend(legend_handles, custom_labels, loc='upper left', frameon=False)

# plt.show()       


"""
Cluster and plot ligands according to chemical tanimoto and water tanimoto distance.
OPTIONAL. Requires comparing ligand and water distances.
"""

# print("Clustering ligands")
# db = DBSCAN(eps=0.4, min_samples=2, metric='precomputed')
# labels = db.fit_predict(lig_dissimilarities)
# colors = [f"C{l+1}" for l in labels]

# # project to 2D using MDS
# mds = MDS(n_components=2, n_init=4, dissimilarity='precomputed', random_state=42, normalized_stress='auto')
# coords_lig = mds.fit_transform(lig_dissimilarities)
# coords_water = mds.fit_transform(water_dissimilarities)

# # plot with labels
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6), sharex=True, sharey=True)
# ax1.scatter(coords_lig[:, 0], coords_lig[:, 1], color=colors, edgecolor='black', s=100)
# ax2.scatter(coords_water[:, 0], coords_water[:, 1], color=colors, edgecolor='black', s=100)

# texts1 = []
# texts2 = []
# for i, name in enumerate(all_prots[1:]):
    # texts1.append( ax1.text(coords_lig[i, 0], coords_lig[i, 1], name, fontsize=12) )
    # texts2.append( ax2.text(coords_water[i, 0], coords_water[i, 1], name, fontsize=12) )

# adjust_text(texts1, ax=ax1)
# adjust_text(texts2, ax=ax2)

# ax1.set_title("Ligand Projection via Ligand Tanimoto Distance (MDS)")
# ax2.set_title("Ligand Projection via Water Map Tanimoto Distance (MDS)")
# ax1.set_xlabel("MDS Dimension 1")
# ax2.set_xlabel("MDS Dimension 1")
# ax1.set_ylabel("MDS Dimension 2")
# plt.tight_layout()
# plt.show()
