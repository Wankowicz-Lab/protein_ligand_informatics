import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# === Load CSVs ===
bonds_df = pd.read_csv("/.../water_bond_clash_summary.csv")
rigidity_df = pd.read_csv("/.../updated_main_dataset_avg_op_per_pdb.csv")

# === Prepare and merge ===
bonds_df['pdb'] = bonds_df['pdb'].astype(str).str.upper()
rigidity_df['pdb_id'] = rigidity_df['pdb_id'].astype(str).str.upper()
merged_df = pd.merge(bonds_df, rigidity_df, left_on='pdb', right_on='pdb_id')

# === Remove outliers within each rigidity category using IQR ===
def remove_outliers(group):
    q1 = group['total_bonds'].quantile(0.25)
    q3 = group['total_bonds'].quantile(0.75)
    iqr = q3 - q1
    return group[(group['total_bonds'] >= q1 - 1.5 * iqr) & (group['total_bonds'] <= q3 + 1.5 * iqr)]

filtered_df = merged_df.groupby('rigidity_quartile', group_keys=False).apply(remove_outliers)

# === Count per rigidity category (after outlier removal) ===
counts = filtered_df['rigidity_quartile'].value_counts().to_dict()

# === Descriptive stats ===
stats_df = filtered_df.groupby('rigidity_quartile')['total_bonds'].describe()
print("=== Total Bonds per Rigidity Category (Outliers Removed) ===")
print(stats_df)

# Save stats
stats_df.to_csv("total_bonds_stats_by_rigidity_no_outliers.csv_qFit")

# === Label with PDB counts ===
filtered_df['rigidity_labeled'] = filtered_df['rigidity_quartile'].map(
    lambda x: f"{x} (n={counts.get(x, 0)})"
)

# === Define order and labeled order ===
order = ['very flexible', 'flexible', 'rigid', 'very rigid']
labeled_order = [f"{cat} (n={counts.get(cat, 0)})" for cat in order]

# === Set category color palette ===
color_map = {
    "very flexible": "#483D8B",  # darkslateblue
    "flexible": "#6495ED",       # cornflowerblue
    "rigid": "#20B2AA",          # lightseagreen
    "very rigid": "#90EE90"      # lightgreen
}
palette = {f"{cat} (n={counts.get(cat, 0)})": color_map[cat] for cat in order}

# === Boxplot ===
sns.set(style="whitegrid")
plt.figure(figsize=(8, 6))
sns.boxplot(
    x='rigidity_labeled',
    y='total_bonds',
    data=filtered_df,
    order=labeled_order,
    palette=palette
)

plt.xlabel("Rigidity Category")
plt.ylabel("Total Bonds per PDB")
plt.title("Total Bonds by Rigidity Category (Outliers Removed)")
plt.xticks(rotation=20)
plt.tight_layout()
plt.savefig("total_bonds_by_rigidity_boxplot_no_outliers_qFit.png")
