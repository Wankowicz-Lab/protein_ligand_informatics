#packages
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski
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
parser = argparse.ArgumentParser(description='Generate data description figures from ligand affinity data')
parser.add_argument('affinity_csv', type=str, help='Path to the affinity data CSV file. This must incuded PDB, occ or BDC, and smiles of all ligands')
parser.add_argument('--BDC', action='store_true', help='Flag indicating whether to use BDC filtering or analysis')
args = parser.parse_args()

ligand_data = args.affinity_csv
use_bdc = args.BDC

# Function to calculate properties from SMILES
def calculate_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    properties = {
        'MW': Descriptors.MolWt(mol),
        'logP': Crippen.MolLogP(mol),
        'Hbond_donor': Lipinski.NumHDonors(mol),
        'Hbond_acceptor': Lipinski.NumHAcceptors(mol),
        'Lipinski': Lipinski.rdMolDescriptors.CalcNumLipinskiHBA(mol) + Lipinski.rdMolDescriptors.CalcNumLipinskiHBD(mol),
        'LogD': Crippen.MolLogP(mol),  # LogD is often approximated by logP
        'Polar_SA': Descriptors.TPSA(mol),
        'rotatable': Lipinski.NumRotatableBonds(mol)
    }
    return properties



if use_bdc:
    ligand_data['occ'] = ligand_data['1-BDC'] * 2

# Subset combined_data to have an occupancy > 0.2
ligand_data = ligand_data[ligand_data['occ'] > 0.2]

# Get the number of structures in the subset
num_structures = ligand_data.shape[0]
print(f'Number of structures with occupancy > 0.2: {num_structures}')
pdb_values = ligand_data['PDB'].tolist()
print("List of PDB values:", pdb_values)

# Histogram of Occupancy 
plt.figure(figsize=(10, 6))
plt.hist(ligand_data['occ'], bins=30, color='blue', edgecolor='black', alpha=0.7)
plt.xlabel('Occupancy')
plt.ylabel('Frequency')
plt.xticks()
plt.yticks()
plt.tight_layout()
plt.savefig('histogram_occ.png', dpi=300)
plt.close()

# Calculate statistics for 'occ'
occ_mean = ligand_data['occ'].mean()
occ_median = ligand_data['occ'].median()
occ_iqr = ligand_data['occ'].quantile(0.75) - ligand_data['occ'].quantile(0.25)
occ_min = ligand_data['occ'].min()
occ_max = ligand_data['occ'].max()
print(f'occ - Mean: {occ_mean}, Median: {occ_median}, IQR: {occ_iqr}, Min: {occ_min}, Max: {occ_max}')


properties_list = ligand_data['smiles'].apply(calculate_properties)
properties_df = pd.DataFrame(properties_list.tolist())
ligand_data = pd.concat([ligand_data, properties_df], axis=1)
ligand_data['Hbond_MW'] = (
    ligand_data['Hbond_donor'] + ligand_data['Hbond_acceptor']
) / ligand_data['MW']

# Create a mapping from original column names to display names
column_name_map = {
    'MW': 'Molecular Weight',
    'Hbond_donor': 'Hydrogen Bond Donors',
    'Hbond_acceptor': 'Hydrogen Bond Acceptors',
    'Polar_SA': 'Polar Surface Area',
    'rotatable': 'Rotatable Bonds',
    'Hbond_MW': 'Hydrogen Bonds/Molecular Weight'
}

# Replace column headers with new display names where appropriate
ligand_data.rename(columns=column_name_map, inplace=True)


columns_to_plot = ['Molecular Weight', 'logP', 'Hydrogen Bond Donors', 'Hydrogen Bond Acceptors',
                    'LogD', 'Polar Surface Area', 'Rotatable Bonds', 'Hydrogen Bonds/Molecular Weight']

for column in columns_to_plot:
    # Plot histogram for current property
    plt.figure(figsize=(10, 6))
    plt.hist(ligand_data[column], bins=30, color='green', edgecolor='black', alpha=0.7)
    plt.xlabel(column)
    plt.ylabel('Frequency')
    plt.xticks()
    plt.yticks()
    plt.tight_layout()
    # Save figure with a filename corresponding to the property
    plt.savefig(f'histogram_{column.replace(" ", "_").lower()}.png', dpi=300)
    plt.close()
    
    # Calculate descriptive statistics for current property
    col_mean = ligand_data[column].mean()
    col_median = ligand_data[column].median()
    col_iqr = ligand_data[column].quantile(0.75) - ligand_data[column].quantile(0.25)
    col_min = ligand_data[column].min()
    col_max = ligand_data[column].max()
    print(f'{column} - Mean: {col_mean}, Median: {col_median}, IQR: {col_iqr}, Min: {col_min}, Max: {col_max}')

out_csv = f"{args.affinity_csv}_ligand_properties.csv"
ligand_data.to_csv(out_csv, index=False)
