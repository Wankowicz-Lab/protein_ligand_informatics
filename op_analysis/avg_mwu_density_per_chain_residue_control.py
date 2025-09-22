import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
import os

# === File paths ===
binding_csv = "/path/to/main_dataset_binding_op.csv"
nonbinding_csv = "/path/to/main_dataset_nonbinding_op.csv"
output_path = "/path/to/results"
os.makedirs(output_path, exist_ok=True)

# === Load data ===
df_binding = pd.read_csv(binding_csv)
df_nonbinding = pd.read_csv(nonbinding_csv)

# === Prepare list to hold per-pocket averaged values ===
pocket_data = []

# === Group by pdb_id + chain (i.e., per pocket) ===
pocket_keys = df_binding[['pdb_id', 'chain']].drop_duplicates()

for _, row in pocket_keys.iterrows():
    pdb_id = row['pdb_id']
    chain = row['chain']

    b_sub = df_binding[(df_binding['pdb_id'] == pdb_id) & (df_binding['chain'] == chain)]
    nb_sub = df_nonbinding[(df_nonbinding['pdb_id'] == pdb_id) & (df_nonbinding['chain'] == chain)]

    # If either is empty, skip
    if b_sub.empty or nb_sub.empty:
        continue

    # Match residue types
    matched_binding = []
    matched_nonbinding = []

    for resn in b_sub['resn'].unique():
        b_res = b_sub[b_sub['resn'] == resn]
        nb_res = nb_sub[nb_sub['resn'] == resn]

        n_match = min(len(b_res), len(nb_res))
        if n_match == 0:
            continue

        matched_binding.append(b_res.sample(n=n_match, random_state=42))
        matched_nonbinding.append(nb_res.sample(n=n_match, random_state=42))

    # Combine matched residues
    if matched_binding and matched_nonbinding:
        matched_b = pd.concat(matched_binding)
        matched_nb = pd.concat(matched_nonbinding)

        avg_b = matched_b['s2calc'].mean()
        avg_nb = matched_nb['s2calc'].mean()
        delta = avg_b - avg_nb

        pocket_data.append({
            'pdb_id': pdb_id,
            'chain': chain,
            'avg_s2calc_binding': avg_b,
            'avg_s2calc_nonbinding': avg_nb,
            'delta_s2calc': delta
        })

# === Convert to DataFrame ===
df_out = pd.DataFrame(pocket_data)

# === Save merged data ===
df_out.to_csv(f"{output_path}/residue_type_controlled_pocket_avg.csv", index=False)

# === Statistical Test ===
u_stat, p_value = mannwhitneyu(
    df_out['avg_s2calc_binding'],
    df_out['avg_s2calc_nonbinding'],
    alternative='two-sided'
)

# === Save Mann–Whitney U Test Results ===
mwu_results = f"=== Mann–Whitney U Test (Binding vs Non-binding, residue-type controlled) ===\n"
mwu_results += f"U statistic: {u_stat:.3f}\n"
mwu_results += f"P-value: {p_value}\n"

with open(f"{output_path}/mwu_test_results.txt", "w") as f:
    f.write(mwu_results)

# === Save Descriptive Statistics ===
desc_stats = f"=== Descriptive Statistics (Residue-type Controlled) ===\n"
desc_stats += f"Binding Site Residues:\n"
desc_stats += f"  Mean: {df_out['avg_s2calc_binding'].mean():.4f}\n"
desc_stats += f"  Median: {df_out['avg_s2calc_binding'].median():.4f}\n"
desc_stats += f"  Std: {df_out['avg_s2calc_binding'].std():.4f}\n"
desc_stats += f"  Count: {len(df_out['avg_s2calc_binding'])}\n"
desc_stats += f"\nNon-binding Site Residues:\n"
desc_stats += f"  Mean: {df_out['avg_s2calc_nonbinding'].mean():.4f}\n"
desc_stats += f"  Median: {df_out['avg_s2calc_nonbinding'].median():.4f}\n"
desc_stats += f"  Std: {df_out['avg_s2calc_nonbinding'].std():.4f}\n"
desc_stats += f"  Count: {len(df_out['avg_s2calc_nonbinding'])}\n"

with open(f"{output_path}/descriptive_statistics.txt", "w") as f:
    f.write(desc_stats)

# === Melt for plotting ===
melted = df_out.melt(id_vars=['pdb_id', 'chain'],
                     value_vars=['avg_s2calc_binding', 'avg_s2calc_nonbinding'],
                     var_name='Residue Type', value_name='s2calc')
melted['Residue Type'] = melted['Residue Type'].map({
    'avg_s2calc_binding': 'Binding Site Residues',
    'avg_s2calc_nonbinding': 'Non-binding Site Residues'
})

# === Color palette ===
palette_main = {
    'Binding Site Residues': '#00BFC4',
    'Non-binding Site Residues': '#F564E3'
}

# === Plot setup ===
sns.set(style="white", context="talk")
fig, (ax_box, ax_kde) = plt.subplots(
    nrows=2, sharex=True, figsize=(6.4, 5.5),
    gridspec_kw={"height_ratios": (1, 4)}
)

# === Boxplot ===
sns.boxplot(
    x="s2calc", y="Residue Type",
    data=melted, ax=ax_box,
    palette=palette_main, linewidth=2.0, fliersize=2, orient="h"
)
ax_box.set_xlabel("")  # No x-label here
ax_box.set_ylabel("")  # No y-label
ax_box.set_yticklabels([])  # Remove y-axis labels
ax_box.tick_params(axis='y', left=False)  # Remove y-axis ticks
ax_box.tick_params(axis='x', labelbottom=False)  # Hide x-axis labels on boxplot
ax_box.grid(True, axis='x', linestyle='--', alpha=0.4)

# === KDE Plot ===
for label, color in palette_main.items():
    sns.kdeplot(
        data=melted[melted['Residue Type'] == label],
        x="s2calc", ax=ax_kde,
        label=label, fill=True, alpha=0.4, linewidth=2, color=color
    )

# Labeling and appearance
ax_kde.set_xlabel("Average Order Parameter", fontsize=14)
ax_kde.set_ylabel("Density", fontsize=14)
ax_kde.tick_params(axis='both', labelsize=13)
ax_kde.grid(True, axis='x', linestyle='--', alpha=0.4)

# Remove y-axis ticks and labels
ax_kde.set_yticklabels([])
ax_kde.tick_params(axis='y', left=False)

# Add legend
ax_kde.legend(title="", loc="upper left", fontsize=12, frameon=False)

# === Finalize ===
plt.tight_layout(h_pad=0.2)
plt.savefig(f"{output_path}/avg_s2calc_density_boxplot_residue_type_controlled.png", dpi=300, bbox_inches="tight")
plt.close()
