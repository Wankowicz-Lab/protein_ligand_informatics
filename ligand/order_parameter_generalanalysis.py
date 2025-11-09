#packages
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
from Bio.PDB import *
from scipy import stats
import argparse

from sklearn.feature_selection import mutual_info_regression
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.metrics.pairwise import cosine_similarity
from analysis_functions import *

# Set font and figure params
plt.rcParams.update({
    'font.size': 24,
    'axes.titlesize': 24,
    'axes.labelsize': 24,
    'xtick.labelsize': 24,
    'ytick.labelsize': 24,
    'legend.fontsize': 24
})

boxprops = {'edgecolor': 'k', 'linewidth': 2}
lineprops = {'color': 'k', 'linewidth': 2}

boxplot_kwargs = dict({'boxprops': boxprops, 'medianprops': lineprops,
                       'whiskerprops': lineprops, 'capprops': lineprops,
                       'width': 0.75})


# Set up argument parser
parser = argparse.ArgumentParser(description="Data cleaning and general analysis for mac1 ensemble.")
parser.add_argument('--mac1_affinity', type=str, default='/Users/stephaniewanko/Downloads/vanderbilt/mac1/mac1_affinity.csv', help='Path to mac1_affinity.csv')
parser.add_argument('--apo_op', type=str, default='/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/7kqo_OP.out', help='Path to apo order parameter file')
parser.add_argument('--op_dir', type=str, default='/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/OP/', help='Directory containing *OP.out files')
parser.add_argument('--use_close_resi', action='store_true', help='Read and process close_resi files')
parser.add_argument('--close_resi', type=float, default=3.0, help='Distance cutoff (in angstroms) for defining close residues')
parser.add_argument('--far_resi', type=float, default=10.0, help='Distance cutoff (in angstroms) for defining far residues')
parser.add_argument('--resi_dist_dir', type=str, default='/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/non_normalized_OP/', help='Directory containing *_closeres.csv files')
parser.add_argument('--crystal_contacts', type=str, nargs='*', default=['48', '55', '101', '102', '159', '11', '17', '158', '166', '169', '58', '87', '156'], help='List of residue IDs for crystal contacts')
parser.add_argument('--chain', type=str, default='A', help='Specify which chain to filter for (default: A)')
parser.add_argument('--output', type=str, default='Mac1', help='Base output filename (no extension)')

args, _ = parser.parse_known_args()

# ____________IMPORT FILES__________________
mac1_affinity = pd.read_csv(args.mac1_affinity)
apo_op = pd.read_csv(args.apo_op)

# Read all *OP.out files
all_files = glob.glob(f"{args.op_dir}*_OP.out")
li = []
pdb_remove = []

for filename in all_files:
    df = pd.read_csv(filename, index_col=None, sep=',', header=0)
    # Attempt to intelligently parse PDB ID from filename
    # You may want to make the slicing user-configurable if needed
    df['PDB'] = filename[-12:-7] if len(filename) > 15 else filename[-9:-4]
    li.append(df)

order_all = pd.concat(li, axis=0, ignore_index=True)

crystal_contacts = args.crystal_contacts
order_all_A = order_all[order_all['chain'] == args.chain]
if crystal_contacts:
    order_all_A = order_all_A[~order_all_A['resi'].isin(crystal_contacts)]

#____________process, merge, differences all order parameters_____________
order_all_A = pd.merge(order_all_A, apo_op, on=['resi', 'chain'], how='left', suffixes=('', '_apo'))
order_all_A['s2calc_diff'] =  order_all_A['s2calc'] - order_all_A['s2calc_apo']
order_all_A['s2calc_diff'] = order_all_A['s2calc_diff']
order_all_A.to_csv(f'{args.output}_order_all_A.csv', index=False)

plt.figure(figsize=(20, 6))
sns.boxplot(x='resi', y='s2calc_diff', data=order_all_A, flierprops=dict(marker='o', color='gray', markersize=1))
plt.xlabel('Residue')
plt.xlabel('Order Parameter Differences')
plt.xticks(rotation=45)
plt.savefig(f'{args.output}_residue_distribution_orderparams.png')

# Create clustermap of residues by order parameters
pivot_op = order_all_A.pivot_table(index='PDB', columns='resi', values='s2calc_diff').fillna(0)
pivot_op = pivot_op.drop(columns=[150], errors='ignore')

# Generate the clustermap
plt.figure()
g = sns.clustermap(pivot_op, cmap='magma_r', center=0, figsize=(10, 8), col_cluster=True)
g.ax_heatmap.set_xticks([])
g.ax_heatmap.set_yticks([])
g.ax_heatmap.set_xlabel('Residues')
g.ax_heatmap.set_ylabel('PDBs')

row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind
pivot_op_reordered = pivot_op.iloc[row_order, col_order]
pivot_op_reordered.to_csv(args.output + '_pivot_op_reordered.csv')

plt.savefig(f'{args.output}_OP_clustermap.png', dpi=300, bbox_inches='tight')  # High resolution for publication
plt.close()  


# Calculate the average (mean) and standard deviation of s2calc_diff for each residue (resi, chain)
s2_sd_stats = order_all_A.groupby(['resi', 'chain'])['s2calc_diff'].agg(
    s2calc_mean_diff='mean',
    s2calc_std='std'
).reset_index()

pearson_corr, _ = stats.pearsonr(s2_sd_stats['s2calc_mean_diff'], s2_sd_stats['s2calc_std'])
print(f'Pearson correlation between average s2calc difference and standard deviation: {pearson_corr:.2f}')

plt.figure(figsize=(10, 8))
sns.scatterplot(x='s2calc_mean_diff', y='s2calc_std', data=s2_sd_stats, s=100, color='blue', edgecolor='w', alpha=0.8)
sns.regplot(x='s2calc_mean_diff', y='s2calc_std', data=s2_sd_stats, scatter=False, color='red', line_kws={'linewidth': 2})
plt.xlabel('Average Order Parameter Difference')
plt.ylabel('Standard Deviation of Order Parameter Difference')
plt.xticks()
plt.yticks()
plt.tight_layout()
plt.savefig(f'{args.output}_scatter_avg_s2calc_vs_stddev.png', dpi=300, bbox_inches='tight')
plt.close()

# Calculate the correlation matrix
pivot_op = order_all_A.pivot_table(index='resi', columns='PDB', values='s2calc_diff').fillna(0)
corr_matrix = pivot_op.corrwith(pivot_op, axis=1)
corr_matrix = pivot_op.T.corr()

plt.figure(figsize=(14, 12))
clustermap = sns.clustermap(
    corr_matrix,
    cmap="inferno",
    square=True,
    cbar_kws={"label": "Correlation"},
    figsize=(14, 12)
)
# Correct: reference the heatmap axis from clustermap, not undefined ax_heatmap
clustermap.ax_heatmap.set_xticks([])
clustermap.ax_heatmap.set_yticks([])
clustermap.ax_heatmap.set_xlabel('Residues')
clustermap.ax_heatmap.set_ylabel('Residues')

plt.tight_layout()
plt.savefig(f'{args.output}_OP_correlation_all_cluster.png', dpi=300)
plt.close()


# Optionally read in close_resi files
if args.use_close_resi:
    all_files = glob.glob(f"{args.resi_dist_dir}*_closeres.csv")

    # downstream processing logic should follow here (not included in the selection)
    li = []
    for filename in all_files:
        df = pd.read_csv(filename)
        li.append(df)

    res_distance = pd.concat(li, axis=0, ignore_index=True)
    res_distance = res_distance[res_distance['chain']==args.chain]

    close_resi = res_distance[res_distance['distance'] <= args.close_resi]
    far_resi = res_distance[res_distance['distance'] > args.far_resi]

    #create close and distal 
    close_OP = order_all_A.merge(
        close_resi[['resi', 'chain', 'PDB']],    # columns used to match
        on=['resi','chain','PDB'],            # merge keys
        how='inner'                            # only keep matches in both
    )

    distal_OP = order_all_A.merge(
        close_resi[['resi', 'chain', 'PDB']],    # columns used to match
        on=['resi','chain','PDB'],            # merge keys
        how='inner'                            # only keep matches in both
    )

    close_OP = close_OP.drop_duplicates()
    distal_OP = distal_OP.drop_duplicates()

    close_OP.to_csv(f'{args.output}_close_OP.csv', index=False)
    distal_OP.to_csv(f'{args.output}_distal_OP.csv', index=False)
    
    #create mean value for each subset
    avg_close_OP = close_OP.groupby('PDB')['s2calc_diff'].mean().reset_index()
    avg_distal_OP = distal_OP.groupby('PDB')['s2calc_diff'].mean().reset_index()
    avg_OP = order_all_A.groupby('PDB')['s2calc_diff'].mean().reset_index()

    # Correct label for each group
    avg_close_OP['Residue Proximity'] = 'Binding Site Residues'
    avg_distal_OP['Residue Proximity'] = 'Distal Residues'
    avg_OP['Residue Proximity'] = 'All Residues'
   
   # Concatenate the dataframes
    combined_OP = pd.concat([close_OP, distal_OP, order_all_A], ignore_index=True)
   
    plt.figure(figsize=(12, 8))
    sns.histplot(data=combined_OP, x='s2calc_diff', hue='Residue Proximity', multiple='stack', palette={'Binding Site Residues': 'blue', 'Distal Residues': 'green', 'All Residues': 'orange'}, bins=50)
    plt.xlabel('s2calc Difference', fontsize=14)
    plt.ylabel('Count', fontsize=14)
    plt.title('Histogram of s2calc Differences by Residue Proximity', fontsize=16)
    plt.legend(title='Residue Proximity', title_fontsize='13', fontsize='12')
    plt.tight_layout()
    plt.savefig(f'{args.output}_OP_histogram_close_far_all.png', dpi=300, bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(14, 10))
    sns.kdeplot(data=avg_close_OP, x='s2calc_diff', label='Close Residues', color='blue', linewidth=8)
    sns.kdeplot(data=avg_distal_OP, x='s2calc_diff', label='Distal Residues', color='green', linewidth=8)
    sns.kdeplot(data=avg_OP, x='s2calc_diff', label='All', color='orange', linewidth=8)
    plt.legend(title='Residue Location')
    plt.xlabel('Order Parameter Differences')
    plt.ylabel('Density')
    plt.tight_layout()
    plt.savefig(f'{args.output}_OP_distribution_close_far.png', dpi=300, bbox_inches='tight')
    plt.close()

