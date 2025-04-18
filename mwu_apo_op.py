import pandas as pd
import os
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt

# === Paths ===
all_path = "/dors/wankowicz_lab/ellas/apo_op/all_op_combined.csv"
binding_path = "/dors/wankowicz_lab/ellas/apo_op/close_resi_filtered/binding_combined.csv"
nonbinding_output = "/dors/wankowicz_lab/ellas/apo_op/close_resi_filtered/nonbinding_combined.csv"

# === Load data ===
df_all = pd.read_csv(all_path)
df_binding = pd.read_csv(binding_path)

# === Ensure consistent types for matching ===
for df in [df_all, df_binding]:
    df['resi'] = df['resi'].astype(str)
    df['chain'] = df['chain'].astype(str)
    df['resn'] = df['resn'].astype(str)

# === Create unique match key (you can include more columns if needed) ===
df_all['match_id'] = df_all['resi'] + "_" + df_all['chain'] + "_" + df_all['resn']
df_binding['match_id'] = df_binding['resi'] + "_" + df_binding['chain'] + "_" + df_binding['resn']

# === Identify non-binding residues ===
df_nonbinding = df_all[~df_all['match_id'].isin(df_binding['match_id'])].copy()

# === Save non-binding combined CSV ===
df_nonbinding.to_csv(nonbinding_output, index=False)
print(f"Saved non-binding combined CSV with {len(df_nonbinding)} residues.")

# === Mann-Whitney U test ===
u_stat, p_value = mannwhitneyu(df_nonbinding['s2calc'].dropna(), df_binding['s2calc'].dropna(), alternative='two-sided')

print(f"Non-binding residues: {len(df_nonbinding)}")
print(f"Binding residues: {len(df_binding)}")
print(f"Mann-Whitney U statistic: {u_stat}")
print(f"P-value: {p_value}")

# === Optional plot ===
plt.figure(figsize=(10, 6))
plt.boxplot([df_nonbinding['s2calc'].dropna(), df_binding['s2calc'].dropna()],
            labels=['Non-binding', 'Binding'],
            patch_artist=True,
            boxprops=dict(facecolor='lightblue'),
            medianprops=dict(color='red'))
plt.title("s2calc: Non-binding vs Binding")
plt.ylabel("s2calc")
plt.grid(True, linestyle='--', alpha=0.3)
plt.tight_layout()
plt.savefig("/dors/wankowicz_lab/ellas/apo_op/close_resi_filtered/nonbinding_vs_binding_boxplot.png", dpi=300)
plt.show()
