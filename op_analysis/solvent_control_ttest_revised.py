import pandas as pd
import numpy as np
from scipy.stats import ttest_rel
import seaborn as sns
import matplotlib.pyplot as plt
import sys

# File paths
close_residues_path = "/dors/wankowicz_lab/ellas/main_dataset/updated_single_uniprot_dataset_sasa_op_merged_binding.csv"
non_binding_path = "/dors/wankowicz_lab/ellas/main_dataset/updated_single_uniprot_dataset_sasa_op_merged_nonbinding.csv"
output_path = "/dors/wankowicz_lab/ellas/main_dataset"

print("Loading CSV files...")
sys.stdout.flush()

# Load input files
try:
    close_residues = pd.read_csv(close_residues_path)
    non_binding_residues = pd.read_csv(non_binding_path)
except FileNotFoundError as e:
    print(f"ERROR: {e}")
    sys.exit(1)

print(f"Binding residues loaded: {len(close_residues)}")
print(f"Non-binding residues loaded: {len(non_binding_residues)}")
sys.stdout.flush()

# Drop rows with missing solvent_exposure values
close_residues = close_residues.dropna(subset=['solvent_exposure'])
non_binding_residues = non_binding_residues.dropna(subset=['solvent_exposure'])

# Ensure correct types
for df in [close_residues, non_binding_residues]:
    df['pdb_id'] = df['pdb_id'].astype(str)
    df['chain'] = df['chain'].astype(str)
    df['resn'] = df['resn'].astype(str)
    df['solvent_exposure'] = df['solvent_exposure'].astype(int)

# Create merge key
close_residues['merge_key'] = (
    close_residues['pdb_id'] + "_" + close_residues['chain'] +
    "_" + close_residues['resn'] + "_" + close_residues['solvent_exposure'].astype(str)
)
non_binding_residues['merge_key'] = (
    non_binding_residues['pdb_id'] + "_" + non_binding_residues['chain'] +
    "_" + non_binding_residues['resn'] + "_" + non_binding_residues['solvent_exposure'].astype(str)
)

# Sample one non-binding per group
non_binding_sampled = non_binding_residues.groupby('merge_key').apply(
    lambda x: x.sample(n=1, random_state=42)
).reset_index(drop=True)

# Merge binding with sampled non-binding
combined_data = pd.merge(
    close_residues[['merge_key', 's2calc']],
    non_binding_sampled[['merge_key', 's2calc']],
    on='merge_key',
    suffixes=('_binding', '_non_binding')
)

print(f"Merged pairs: {len(combined_data)}")
sys.stdout.flush()

if len(combined_data) < 5:
    print("Too few matched pairs. Check merge conditions or input data.")
    sys.exit(1)

# Save merged raw values
raw_csv_path = f"{output_path}/s2calc_binding_vs_non_binding_values.csv"
combined_data.to_csv(raw_csv_path, index=False)
print(f"Saved raw merged values to: {raw_csv_path}")
sys.stdout.flush()

# Summary stats
mean_binding = combined_data['s2calc_binding'].mean()
mean_non_binding = combined_data['s2calc_non_binding'].mean()
std_binding = combined_data['s2calc_binding'].std()
std_non_binding = combined_data['s2calc_non_binding'].std()

print(f"Mean s2calc (Binding): {mean_binding:.3f}, Std: {std_binding:.3f}")
print(f"Mean s2calc (Non-Binding): {mean_non_binding:.3f}, Std: {std_non_binding:.3f}")
sys.stdout.flush()

# Paired t-test
t_stat, p_value = ttest_rel(
    combined_data['s2calc_binding'],
    combined_data['s2calc_non_binding']
)
print(f"Paired t-test: t = {t_stat:.3f}, p = {p_value:.3e}")
sys.stdout.flush()

# === Plot boxplot ===
print("Generating and saving boxplot...")
sys.stdout.flush()

# Prepare data
melted = combined_data.melt(
    value_vars=['s2calc_binding', 's2calc_non_binding'],
    var_name='Residue Type',
    value_name='s2calc'
)

# Rename residue type values for clean x-axis labels
melted['Residue Type'] = melted['Residue Type'].map({
    's2calc_binding': 'Binding Site Residues',
    's2calc_non_binding': 'Non-binding Site Residues'
})

# Set consistent style
sns.set(style="whitegrid", context="talk")
plt.figure(figsize=(6.4, 4.8))  # Consistent figure size

# Define consistent color palette
palette = {
    'Binding Site Residues': '#00BFC4',
    'Non-binding Site Residues': '#F564E3'
}

# Draw boxplot
sns.boxplot(
    data=melted,
    x='Residue Type',
    y='s2calc',
    palette=palette,
    width=0.8,
    linewidth=2.0,
    fliersize=3
)

# Customize plot
plt.title("Rigidity (s2calc) of Binding vs Non-binding Residues", fontsize=14)
plt.xlabel("")  # No x-axis label
plt.ylabel("Order Parameter", fontsize=12)
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)
plt.tight_layout()

# Save figure
plt.savefig(f"{output_path}/updated_single_uniprot_dataset_control_s2calc_binding_vs_non_binding_boxplot.png", dpi=300, bbox_inches="tight")
plt.close()

print("Boxplot saved successfully.")
print("Script finished.")

