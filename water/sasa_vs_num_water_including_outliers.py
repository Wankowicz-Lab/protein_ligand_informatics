import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

# === Load input files ===
sasa_path = "/.../binding_site_sasa_per_residue.csv"
water_path = "/.../binding_site_water_distribution.csv"
output_dir = "/.../"  # Make sure this exists

sasa_df = pd.read_csv(sasa_path)
water_df = pd.read_csv(water_path)

# === Clean column names ===
sasa_df.columns = sasa_df.columns.str.strip()
water_df.columns = water_df.columns.str.strip()

# === Rename columns for consistency ===
water_df = water_df.rename(columns={
    'chainID': 'chain',
    'residueNumber': 'resi',
    'residueName': 'resname'
})
sasa_df['resi'] = sasa_df['resi'].astype(str)
water_df['resi'] = water_df['resi'].astype(str)

# === Count waters per residue ===
water_counts = (
    water_df.groupby(['pdb_id', 'chain', 'resi', 'resname'])
    .size()
    .reset_index(name='water_count')
)

# === Merge with SASA data ===
merged_df = pd.merge(
    sasa_df,
    water_counts,
    on=['pdb_id', 'chain', 'resi', 'resname'],
    how='inner'
)

print(f"Total merged data points: {len(merged_df)}\n")

# === Plot and save ===
plot_configs = [
    ("total_sasa", "total_sasa_vs_water.png"),
    ("polar_sasa", "polar_sasa_vs_water.png"),
    ("apolar_sasa", "apolar_sasa_vs_water.png"),
]

for sasa_type, filename in plot_configs:
    x = merged_df["water_count"]
    y = merged_df[sasa_type]

    # === Print correlations ===
    pearson_corr, _ = pearsonr(x, y)
    spearman_corr, _ = spearmanr(x, y)
    print(f"{sasa_type.upper()}")
    print(f"  Pearson correlation:  {pearson_corr:.4f}")
    print(f"  Spearman correlation: {spearman_corr:.4f}\n")

    # === Plot ===
    plt.figure(figsize=(8, 6))
    sns.regplot(
        x=x,
        y=y,
        scatter_kws={'alpha': 0.4},
        line_kws={'color': 'red'},
        ci=None
    )
    plt.xlabel("Number of Nearby Waters", fontsize=18)
    plt.ylabel(sasa_type.replace("_", " ").title(), fontsize=18)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{filename}", dpi=300)
    plt.close()
    print(f"Saved plot: {output_dir}/{filename}")
