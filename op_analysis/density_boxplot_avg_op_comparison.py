import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

# === File paths ===
binding_csv = "/.../main_dataset_binding_op.csv"
nonbinding_csv = "/.../main_dataset_nonbinding_op.csv"
output_path = "/.../results"

# === Load data ===
df_binding = pd.read_csv(binding_csv)
df_nonbinding = pd.read_csv(nonbinding_csv)

# === Compute average s2calc per pdb_id and chain (per pocket) ===
binding_avg = df_binding.groupby(['pdb_id', 'chain'])['s2calc'].mean().reset_index()
binding_avg.columns = ['pdb_id', 'chain', 'avg_s2calc_binding']

nonbinding_avg = df_nonbinding.groupby(['pdb_id', 'chain'])['s2calc'].mean().reset_index()
nonbinding_avg.columns = ['pdb_id', 'chain', 'avg_s2calc_nonbinding']

# === Merge for paired plot ===
merged = pd.merge(binding_avg, nonbinding_avg, on=['pdb_id', 'chain'])
merged['delta_s2calc'] = merged['avg_s2calc_binding'] - merged['avg_s2calc_nonbinding']

# === Statistical Test (Mann–Whitney U test) ===
u_stat, p_value = mannwhitneyu(
    merged['avg_s2calc_binding'],
    merged['avg_s2calc_nonbinding'],
    alternative='two-sided'
)

# === Save Mann–Whitney U Test Results ===
mwu_results = f"=== Mann–Whitney U Test (Binding vs Non-binding) ===\n"
mwu_results += f"U statistic: {u_stat:.3f}\n"
mwu_results += f"P-value: {p_value}\n"

with open(f"{output_path}/mwu_test_results.txt", "w") as f:
    f.write(mwu_results)

print(f"=== Mann–Whitney U Test (Binding vs Non-binding) ===")
print(f"U statistic: {u_stat:.3f}")
print(f"P-value: {p_value}")

# === Save Descriptive Statistics ===
desc_stats = f"=== Descriptive Statistics ===\n"
desc_stats += f"Binding Site Residues:\n"
desc_stats += f"  Mean: {merged['avg_s2calc_binding'].mean():.4f}\n"
desc_stats += f"  Median: {merged['avg_s2calc_binding'].median():.4f}\n"
desc_stats += f"  Std: {merged['avg_s2calc_binding'].std():.4f}\n"
desc_stats += f"  Count: {len(merged['avg_s2calc_binding'])}\n"
desc_stats += f"\nNon-binding Site Residues:\n"
desc_stats += f"  Mean: {merged['avg_s2calc_nonbinding'].mean():.4f}\n"
desc_stats += f"  Median: {merged['avg_s2calc_nonbinding'].median():.4f}\n"
desc_stats += f"  Std: {merged['avg_s2calc_nonbinding'].std():.4f}\n"
desc_stats += f"  Count: {len(merged['avg_s2calc_nonbinding'])}\n"

with open(f"{output_path}/descriptive_statistics.txt", "w") as f:
    f.write(desc_stats)

print(f"\n=== Descriptive Statistics ===")
print(f"Binding Site Residues:")
print(f"  Mean: {merged['avg_s2calc_binding'].mean():.4f}")
print(f"  Median: {merged['avg_s2calc_binding'].median():.4f}")
print(f"  Std: {merged['avg_s2calc_binding'].std():.4f}")
print(f"  Count: {len(merged['avg_s2calc_binding'])}")

print(f"\nNon-binding Site Residues:")
print(f"  Mean: {merged['avg_s2calc_nonbinding'].mean():.4f}")
print(f"  Median: {merged['avg_s2calc_nonbinding'].median():.4f}")
print(f"  Std: {merged['avg_s2calc_nonbinding'].std():.4f}")
print(f"  Count: {len(merged['avg_s2calc_nonbinding'])}")

# Melt for plotting
melted = merged.melt(id_vars=['pdb_id', 'chain'],
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
plt.savefig(f"{output_path}/avg_s2calc_density_boxplot_single_uniprot.png", dpi=300, bbox_inches="tight")
plt.close()
