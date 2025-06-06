import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr
import numpy as np

# === Paths ===
pdb_list_file = "/dors/wankowicz_lab/ellas/main_dataset/updated_main_dataset.txt"
fpocket_base = "/dors/wankowicz_lab/all_pdb/fpocket"
avg_op_file = "/dors/wankowicz_lab/ellas/main_dataset/updated_main_dataset_avg_op_per_pdb.csv"
output_csv = "/dors/wankowicz_lab/ellas/main_dataset/pdb_volume_op_merged.csv"
plot_file = "/dors/wankowicz_lab/ellas/main_dataset/volume_vs_op_spearman.png"

# === Load PDB list ===
with open(pdb_list_file, "r") as f:
    pdb_ids = [line.strip().upper() for line in f if line.strip()]

# === Extract volume from fpocket CSVs ===
volume_data = []
for pdb in pdb_ids:
    fpocket_file = os.path.join(fpocket_base, f"{pdb}_fpocket_pocket1.csv")
    if not os.path.exists(fpocket_file):
        continue
    df = pd.read_csv(fpocket_file)
    match = df[df["Attribute"] == "Volume"]
    if not match.empty:
        volume = match["Value"].values[0]
        volume_data.append((pdb, float(volume)))

volume_df = pd.DataFrame(volume_data, columns=["pdb_id", "volume"])

# === Load OP dataset ===
op_df = pd.read_csv(avg_op_file)
op_df["pdb_id"] = op_df["pdb_id"].str.upper()

# === Merge datasets ===
merged = pd.merge(op_df, volume_df, on="pdb_id")
merged.to_csv(output_csv, index=False)

# === Spearman correlation ===
rho, pval = spearmanr(merged["volume"], merged["avg_s2calc"])
print(f"[Spearman] ρ = {rho:.3f}, p = {pval:.2e}")

# === Plot ===
plt.figure(figsize=(6, 5))
sns.scatterplot(x="volume", y="avg_s2calc", data=merged, alpha=0.7)
sns.regplot(x="volume", y="avg_s2calc", data=merged, scatter=False, lowess=True, color='red', label="LOESS fit")

plt.title("Binding Site Rigidity vs Pocket Volume (Spearman)")
plt.xlabel("Pocket Volume")
plt.ylabel("Average s2calc (Binding Site Rigidity)")

# === Annotate Spearman stats ===
plt.text(0.05, 0.95, f"ρ = {rho:.2f}\np = {pval:.2e}", transform=plt.gca().transAxes,
         fontsize=10, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

plt.legend()
plt.tight_layout()
plt.savefig(plot_file)
plt.close()
