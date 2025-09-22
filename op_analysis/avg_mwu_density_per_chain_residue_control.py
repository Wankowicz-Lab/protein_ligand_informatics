import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

# === File paths ===
binding_csv = "/.../main_dataset_binding_op.csv"
nonbinding_csv = "/.../main_dataset_nonbinding_op.csv"
output_path = "/.../"
output_figure = f"{output_path}/....png"
output_csv = f"{output_path}/....csv"
output_stats_txt = f"{output_path}/....txt"

# === Load and clean ===
df_binding = pd.read_csv(binding_csv)[['pdb_id', 's2calc', 'resn']].dropna().copy()
df_nonbinding = pd.read_csv(nonbinding_csv)[['pdb_id', 's2calc', 'resn']].dropna().copy()

df_binding['resn'] = df_binding['resn'].astype(str).str.strip().str.upper()
df_nonbinding['resn'] = df_nonbinding['resn'].astype(str).str.strip().str.upper()

# === Residue-type-controlled sampling ===
sampled_binding = []
sampled_nonbinding = []

print("Residue-type sampling summary:")
shared_residues = sorted(set(df_binding['resn']) & set(df_nonbinding['resn']))

for resn in shared_residues:
    group_bind = df_binding[df_binding['resn'] == resn]
    group_nonbind = df_nonbinding[df_nonbinding['resn'] == resn]
    n = min(len(group_bind), len(group_nonbind))

    if n >= 10:
        sampled_b = group_bind.sample(n=n, random_state=42)
        sampled_nb = group_nonbind.sample(n=n, random_state=42)
        sampled_binding.append(sampled_b)
        sampled_nonbinding.append(sampled_nb)
        print(f"{resn}: sampled {n} from each group")

df_binding_bal = pd.concat(sampled_binding)
df_nonbinding_bal = pd.concat(sampled_nonbinding)

# === Compute average s2calc per pdb ===
binding_avg = df_binding_bal.groupby('pdb_id')['s2calc'].mean().reset_index()
binding_avg.columns = ['pdb_id', 'avg_s2calc_binding']

nonbinding_avg = df_nonbinding_bal.groupby('pdb_id')['s2calc'].mean().reset_index()
nonbinding_avg.columns = ['pdb_id', 'avg_s2calc_nonbinding']

# === Merge & analyze ===
merged = pd.merge(binding_avg, nonbinding_avg, on='pdb_id')
merged['delta_s2calc'] = merged['avg_s2calc_binding'] - merged['avg_s2calc_nonbinding']
merged.to_csv(output_csv, index=False)
print(f"\nSaved merged per-PDB averages to: {output_csv}")

# === Descriptive stats & save to txt ===
with open(output_stats_txt, 'w') as f:
    f.write("Binding Site Residues (per structure average s2calc):\n")
    f.write(f"{binding_avg['avg_s2calc_binding'].describe()}\n\n")
    f.write("Non-binding Site Residues (per structure average s2calc):\n")
    f.write(f"{nonbinding_avg['avg_s2calc_nonbinding'].describe()}\n\n")

    stat, p = mannwhitneyu(binding_avg['avg_s2calc_binding'],
                           nonbinding_avg['avg_s2calc_nonbinding'],
                           alternative='two-sided')

    f.write("Mann–Whitney U Test (Residue-Type Controlled):\n")
    f.write(f"U = {stat:.2f}, p = {p:.4f}\n")

print(f"Descriptive statistics and test results saved to: {output_stats_txt}")

# === Plot (horizontal boxplot + KDE) ===
melted = merged.melt(id_vars='pdb_id',
                     value_vars=['avg_s2calc_binding', 'avg_s2calc_nonbinding'],
                     var_name='Residue Type', value_name='s2calc')

melted['Residue Type'] = melted['Residue Type'].map({
    'avg_s2calc_binding': 'Binding Site Residues',
    'avg_s2calc_nonbinding': 'Non-binding Site Residues'
})

palette_main = {
    'Binding Site Residues': '#00BFC4',
    'Non-binding Site Residues': '#F564E3'
}

sns.set(style="white", context="talk")
fig, (ax_box, ax_kde) = plt.subplots(
    nrows=2, sharex=True, figsize=(6.4, 5.5),
    gridspec_kw={"height_ratios": (1, 4)}
)

# Boxplot
sns.boxplot(
    x="s2calc", y="Residue Type",
    data=melted, ax=ax_box,
    palette=palette_main, linewidth=2.0, fliersize=2, orient="h"
)
ax_box.set_xlabel("")
ax_box.set_ylabel("")
ax_box.set_yticklabels([])
ax_box.tick_params(axis='y', left=False)
ax_box.tick_params(axis='x', labelbottom=False)
ax_box.grid(True, axis='x', linestyle='--', alpha=0.4)

# KDE Plot
for label, color in palette_main.items():
    sns.kdeplot(
        data=melted[melted['Residue Type'] == label],
        x="s2calc", ax=ax_kde,
        label=label, fill=True, alpha=0.4, linewidth=2, color=color
    )

ax_kde.set_xlabel("Average Order Parameter", fontsize=14)
ax_kde.set_ylabel("Density", fontsize=14)
ax_kde.tick_params(axis='both', labelsize=13)
ax_kde.set_yticklabels([])
ax_kde.tick_params(axis='y', left=False)
ax_kde.grid(True, axis='x', linestyle='--', alpha=0.4)
ax_kde.legend(title="", loc="upper left", fontsize=12, frameon=False)

plt.tight_layout(h_pad=0.2)
plt.savefig(output_figure, dpi=300, bbox_inches="tight")
plt.close()
print(f"Figure saved to {output_figure}")
