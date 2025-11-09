#packages
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import argparse


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
parser.add_argument('--apo_rmsf', type=str, default='/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/other_macro_output/ensemble/rmsf_output_7TX0.updated_refine_001_ensemble.csv', help='Path to apo order parameter file')
parser.add_argument('--rmsf_dir', type=str, default='/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/RMSF/', help='Directory containing *OP.out files')
parser.add_argument('--use_close_resi', action='store_true', help='Read and process close_resi files')
parser.add_argument('--close_resi', type=float, default=3.0, help='Distance cutoff (in angstroms) for defining close residues')
parser.add_argument('--far_resi', type=float, default=10.0, help='Distance cutoff (in angstroms) for defining far residues')
parser.add_argument('--resi_dist_dir', type=str, default='/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/non_normalized_OP/', help='Directory containing *_closeres.csv files')
parser.add_argument('--crystal_contacts', type=str, nargs='*', default=['48', '55', '101', '102', '159', '11', '17', '158', '166', '169', '58', '87', '156'], help='List of residue IDs for crystal contacts')
parser.add_argument('--chain', type=str, default='A', help='Specify which chain to filter for (default: A)')
parser.add_argument('--output', type=str, default='Mac1', help='Base output filename (no extension)')

args, _ = parser.parse_known_args()

# ____________IMPORT FILES__________________
apo_RMSF = pd.read_csv(args.apo_rmsf)

# Read all *OP.out files
all_files = glob.glob(f"{args.rmsf_dir}*ensemble.csv")
li = []
pdb_remove = []

for filename in all_files:
    df = pd.read_csv(filename, index_col=None, sep=',', header=0)
    # Attempt to intelligently parse PDB ID from filename
    # You may want to make the slicing user-configurable if needed
    df['PDB'] = filename[-42:-37] if len(filename) > 15 else filename[-9:-4]
    li.append(df)

RMSF_df = pd.concat(li, axis=0, ignore_index=True)

crystal_contacts = args.crystal_contacts
RMSF_df = RMSF_df[RMSF_df['Chain'] == args.chain]
if crystal_contacts:
    RMSF_df = RMSF_df[~RMSF_df['Resi'].isin(crystal_contacts)]


#____________process, merge, differences all order parameters_____________
delta_RMSF_A = pd.merge(RMSF_df, apo_RMSF, on=['Resi', 'Chain'], how='left', suffixes=('', '_apo'))
delta_RMSF_A['delta_RMSF'] =  delta_RMSF_A['RMSF'] - delta_RMSF_A['RMSF_apo']
delta_RMSF_A.to_csv(f'{args.output}_order_param_diff.csv', index=False)

plt.figure(figsize=(20, 6))
sns.boxplot(x='Resi', y='delta_RMSF', data=delta_RMSF_A, flierprops=dict(marker='o', color='gray', markersize=1))
plt.xlabel('Residue')
plt.xlabel('Order Parameter Differences')
plt.xticks(rotation=45)
plt.savefig(f'{args.output}_residue_distribution_RMSF.png')

# Create clustermap of residues by order parameters
pivot_RMSF = delta_RMSF_A.pivot_table(index='PDB', columns='Resi', values='delta_RMSF').fillna(0)

# Generate the clustermap
plt.figure()
g = sns.clustermap(pivot_RMSF, cmap='magma_r', center=0, figsize=(10, 8), col_cluster=True)
g.ax_heatmap.set_xticks([])
g.ax_heatmap.set_yticks([])
g.ax_heatmap.set_xlabel('Residues')
g.ax_heatmap.set_ylabel('PDBs')

row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind
pivot_RMSF_reordered = pivot_RMSF.iloc[row_order, col_order]
pivot_RMSF_reordered.to_csv(args.output + '_pivot_RMSF_reordered.csv')

plt.savefig(f'{args.output}_RMSF_clustermap.png', dpi=300, bbox_inches='tight')  # High resolution for publication
plt.close()  

# Calculate the correlation matrix
pivot_rmsf = delta_RMSF_A.pivot_table(index='Resi', columns='PDB', values='delta_RMSF').fillna(0)
corr_matrix = pivot_rmsf.corrwith(pivot_rmsf, axis=1)
corr_matrix = pivot_rmsf.T.corr()

plt.figure(figsize=(14, 12))
clustermap = sns.clustermap(
    corr_matrix,
    cmap="inferno",
    square=True,
    figsize=(14, 12)
)
# Correct: reference the heatmap axis from clustermap, not undefined ax_heatmap
clustermap.ax_heatmap.set_xticks([])
clustermap.ax_heatmap.set_yticks([])
clustermap.ax_heatmap.set_xlabel('Residues')
clustermap.ax_heatmap.set_ylabel('Residues')

plt.tight_layout()
plt.savefig(f'{args.output}_RMSF_correlation_all_cluster.png', dpi=300)
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
    close_RMSF = delta_RMSF_A.merge(
        close_resi[['resi', 'chain', 'PDB']],    # columns used to match
        on=['resi','chain','PDB'],            # merge keys
        how='inner'                            # only keep matches in both
    )

    distal_RMSF = delta_RMSF_A.merge(
        close_resi[['resi', 'chain', 'PDB']],    # columns used to match
        on=['resi','chain','PDB'],            # merge keys
        how='inner'                            # only keep matches in both
    )

    close_RMSF = close_RMSF.drop_duplicates()
    distal_RMSF = distal_RMSF.drop_duplicates()

    close_RMSF.to_csv(f'{args.output}_close_RMSF.csv', index=False)
    distal_RMSF.to_csv(f'{args.output}_distal_RMSF.csv', index=False)
    
    #create mean value for each subset
    avg_close_RMSF = close_RMSF.groupby('PDB')['s2calc_diff'].mean().reset_index()
    avg_distal_RMSF = distal_RMSF.groupby('PDB')['s2calc_diff'].mean().reset_index()
    avg_RMSF = delta_RMSF_A.groupby('PDB')['s2calc_diff'].mean().reset_index()

    # Correct label for each group
    avg_close_RMSF['Residue Proximity'] = 'Binding Site Residues'
    avg_distal_RMSF['Residue Proximity'] = 'Distal Residues'
    avg_RMSF['Residue Proximity'] = 'All Residues'
   
   # Concatenate the dataframes
    combined_RMSF = pd.concat([close_RMSF, distal_RMSF, delta_RMSF_A], ignore_index=True)
   
    plt.figure(figsize=(12, 8))
    sns.histplot(data=combined_RMSF, x='delta_RMSF', hue='Residue Proximity', multiple='stack', palette={'Binding Site Residues': 'blue', 'Distal Residues': 'green', 'All Residues': 'orange'}, bins=50)
    plt.xlabel('Δ RMSF', fontsize=14)
    plt.ylabel('Count', fontsize=14)
    plt.legend(title='Residue Proximity', title_fontsize='13', fontsize='12')
    plt.tight_layout()
    plt.savefig(f'{args.output}_RMSF_histogram_close_far_all.png', dpi=300, bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(14, 10))
    sns.kdeplot(data=avg_close_RMSF, x='s2calc_diff', label='Close Residues', color='blue', linewidth=8)
    sns.kdeplot(data=avg_distal_RMSF, x='s2calc_diff', label='Distal Residues', color='green', linewidth=8)
    sns.kdeplot(data=avg_RMSF, x='s2calc_diff', label='All', color='orange', linewidth=8)
    plt.xlabel('Δ RMSF')
    plt.ylabel('Density')
    plt.tight_layout()
    plt.savefig(f'{args.output}_RMSF_distribution_close_far.png', dpi=300, bbox_inches='tight')
    plt.close()

