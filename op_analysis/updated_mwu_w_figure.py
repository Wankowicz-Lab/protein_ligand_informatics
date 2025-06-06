import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

# === File paths ===
binding_csv = "/dors/wankowicz_lab/ellas/main_dataset/updated_single_uniprot_dataset_binding_op.csv"
nonbinding_csv = "/dors/wankowicz_lab/ellas/main_dataset/updated_single_uniprot_dataset_nonbinding_op.csv"
output_figure = "/dors/wankowicz_lab/ellas/main_dataset/updated_single_uniprot_dataset_s2calc_binding_vs_nonbinding.png"

# === Load data ===
df_binding = pd.read_csv(binding_csv)
df_nonbinding = pd.read_csv(nonbinding_csv)

# === Drop NA values in s2calc ===
df_binding = df_binding[['s2calc']].dropna().copy()
df_binding['group'] = 'Binding Site Residues'

df_nonbinding = df_nonbinding[['s2calc']].dropna().copy()
df_nonbinding['group'] = 'Non-binding Site Residues'

# === Combine ===
combined = pd.concat([df_binding, df_nonbinding], ignore_index=True)

# === Mann-Whitney U Test ===
stat, p = mannwhitneyu(df_binding['s2calc'], df_nonbinding['s2calc'], alternative='two-sided')
print(f"Mann-Whitney U: U={stat:.2f}, p={p:.3e}")

# === Plot ===
sns.set(style="whitegrid", context="talk")  
plt.figure(figsize=(6.4, 4.8))              # Reasonable size 

# Define consistent color palette
palette = {
    "Binding Site Residues": "#00BFC4",
    "Non-binding Site Residues": "#F564E3"
}

# Ensure category order
combined['group'] = pd.Categorical(
    combined['group'],
    categories=["Binding Site Residues", "Non-binding Site Residues"],
    ordered=True
)

# Draw boxplot
sns.boxplot(
    x="group",
    y="s2calc",
    data=combined,
    palette=palette,
    width=0.8,
    linewidth=2.0,
    fliersize=3
)

# === Customize labels and title with smaller font sizes ===
plt.title("Rigidity (s2calc) of Binding vs Non-binding Residues", fontsize=14)
plt.xlabel("")  # No x-label
plt.ylabel("Order Parameter", fontsize=12)
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)

plt.tight_layout()

# === Save ===
plt.savefig(output_figure, dpi=300, bbox_inches="tight")
plt.close()

