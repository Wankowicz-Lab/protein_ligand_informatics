import pandas as pd
import matplotlib.pyplot as plt

# === Input & Output paths ===
input_csv = "/.../main_dataset_binding_site_residues_with_resn.csv"   # <-- update with your file
output_plot = "/.../main_dataset_binding_site_residue_distribution.png"
output_stats = ".../main_dataset_binding_site_residue_stats.txt"

# === Load data ===
df = pd.read_csv(input_csv)

# Group by pdb_id + chain, count residues
binding_counts = (
    df.groupby(["pdb_id", "chain"])
    .size()
    .reset_index(name="n_residues")
)

# === Descriptive stats ===
desc_stats = binding_counts["n_residues"].describe()

# Save stats to .txt
with open(output_stats, "w") as f:
    f.write("Descriptive statistics for # of residues in binding sites\n")
    f.write(desc_stats.to_string())
    f.write("\n")
    f.write(f"Median: {binding_counts['n_residues'].median():.2f}\n")
    f.write("\n")
    f.write(f"Total unique pdb_id + chain pairs: {binding_counts.shape[0]}\n")
    
# === Plot distribution ===
plt.figure(figsize=(8,6))
plt.hist(binding_counts["n_residues"], bins=30, edgecolor="black")

plt.xlabel("Number of Residues in Binding Pocket", fontsize=18)
plt.ylabel("Count of Binding Pockets", fontsize=18)

# Increase Y-axis tick font size
plt.tick_params(axis='y', labelsize=16)
plt.tick_params(axis='x', labelsize=16)
plt.tight_layout()
plt.savefig(output_plot, dpi=300)
plt.close()

