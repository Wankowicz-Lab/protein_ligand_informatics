import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# === INPUT FILES ===
pdb_avg_op_csv = "/dors/wankowicz_lab/ellas/main_dataset/updated_main_dataset_avg_op_per_pdb.csv"
residue_csv = "/dors/wankowicz_lab/ellas/main_dataset/all_binding_residues.csv"

# === LOAD DATA ===
avg_op_df = pd.read_csv(pdb_avg_op_csv)
residue_df = pd.read_csv(residue_csv)

# === RIGIDITY CATEGORY FUNCTION ===
def rigidity_category(op):
    if op > 0.873404:
        return "very rigid"
    elif op > 0.822762:
        return "rigid"
    elif op > 0.756902:
        return "flexible"
    else:
        return "very flexible"

# Apply rigidity category
avg_op_df["rigidity_category"] = avg_op_df["avg_s2calc"].apply(rigidity_category)

# === CLEAN COLUMNS ===
residue_df["residue_name"] = residue_df["residue_name"].astype(str).str.strip().str.upper()
residue_df["pdb_id"] = residue_df["pdb_id"].astype(str).str.upper()

# === AMINO ACID ORDER BY POLAR SASA ===
ordered_aas = [
    'GLN', 'LYS', 'ASN', 'ARG', 'GLU', 'TYR', 'SER', 'ASP', 'THR', 'HIS',
    'MET', 'GLY', 'TRP', 'PRO', 'CYS', 'ALA', 'LEU', 'PHE', 'VAL', 'ILE'
]

# === COMPUTE COMPOSITION PER RIGIDITY CATEGORY ===
composition_data = []

for category in ["very flexible", "flexible", "rigid", "very rigid"]:
    pdbs = avg_op_df[avg_op_df["rigidity_category"] == category]["pdb_id"].unique()
    subset = residue_df[residue_df["pdb_id"].isin(pdbs)]
    total = subset["residue_name"].isin(ordered_aas).sum()
    aa_freq = subset["residue_name"].value_counts(normalize=True) * 100

    print(f"\n--- {category.upper()} ---")
    print(f"Total PDBs: {len(pdbs)}")
    print(f"Total binding residues considered: {total}")
    print("Top 5 most frequent AAs:")
    print(aa_freq.loc[ordered_aas].sort_values(ascending=False).head(5).round(2))

    for aa in ordered_aas:
        composition_data.append({
            "resn": aa,
            "rigidity_category": category,
            "percentage": aa_freq.get(aa, 0.0)
        })

# === CREATE FINAL DF FOR PLOTTING ===
composition_df = pd.DataFrame(composition_data)
composition_df["resn"] = pd.Categorical(composition_df["resn"], categories=ordered_aas, ordered=True)
composition_df["rigidity_category"] = pd.Categorical(
    composition_df["rigidity_category"],
    categories=["very flexible", "flexible", "rigid", "very rigid"],
    ordered=True
)

# === SAVE STATS ===
composition_df.to_csv("/dors/wankowicz_lab/ellas/aa_composition_by_rigidity.csv", index=False)

# === PLOT ===
sns.set_style("whitegrid")
sns.set_context("talk", font_scale=1.5)

plt.figure(figsize=(18, 6))
ax = sns.barplot(
    data=composition_df,
    x="resn",
    y="percentage",
    hue="rigidity_category",
    hue_order=["very flexible", "flexible", "rigid", "very rigid"],
    palette=sns.color_palette("viridis", n_colors=4)
)

plt.title("Amino Acid Composition (%) per Rigidity Category", fontsize=22)
plt.ylabel("Percentage (%)", fontsize=20)
plt.xlabel("Amino Acid (Polar → Apolar)", fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.legend(title="Rigidity Category", fontsize=16, title_fontsize=18, bbox_to_anchor=(1.02, 1), loc='upper left')
ax.yaxis.grid(True, linestyle='--', linewidth=0.5)
ax.xaxis.grid(True, linestyle=':', linewidth=0.4)
ax.set_axisbelow(True)

# === SAVE PLOT ===
output_path = "/dors/wankowicz_lab/ellas/aa_composition_by_pdb_rigidity.png"
plt.tight_layout()
plt.savefig(output_path, dpi=300)

print(f"\nPlot saved to: {output_path}")
print("Composition table saved to: /dors/wankowicz_lab/ellas/aa_composition_by_rigidity.csv")
