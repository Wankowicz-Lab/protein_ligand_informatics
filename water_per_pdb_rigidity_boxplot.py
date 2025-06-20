import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pandas.api.types import CategoricalDtype

# === Load your CSVs here ===
# Replace with your actual file paths
water_csv = "/dors/wankowicz_lab/ellas/binding_site_water_distribution.csv"
rigidity_csv = "/dors/wankowicz_lab/ellas/main_dataset/updated_main_dataset_avg_op_per_pdb.csv"

# === Load data ===
water_df = pd.read_csv(water_csv)
rigid_df = pd.read_csv(rigidity_csv)

# === Count water rows per PDB ===
water_counts = water_df.groupby("pdb_id").size().reset_index(name="num_waters")

# === Merge with rigidity data ===
merged_df = pd.merge(water_counts, rigid_df, on="pdb_id")

# === Remove outliers using IQR per group ===
def remove_outliers(group):
    q1 = group["num_waters"].quantile(0.25)
    q3 = group["num_waters"].quantile(0.75)
    iqr = q3 - q1
    lower = q1 - 1.5 * iqr
    upper = q3 + 1.5 * iqr
    return group[(group["num_waters"] >= lower) & (group["num_waters"] <= upper)]

filtered_df = merged_df.groupby("rigidity_quartile", group_keys=False).apply(remove_outliers)

# === Set category color palette and order ===
color_map = {
    "very flexible": "#483D8B",  # darkslateblue
    "flexible": "#6495ED",      # cornflowerblue
    "rigid": "#20B2AA",         # lightseagreen
    "very rigid": "#90EE90"     # lightgreen
}

category_order = ["very flexible", "flexible", "rigid", "very rigid"]
rigidity_dtype = CategoricalDtype(categories=category_order, ordered=True)
filtered_df["rigidity_quartile"] = filtered_df["rigidity_quartile"].astype(rigidity_dtype)

# === Plot boxplot without outliers ===
plt.figure(figsize=(8, 6))
sns.boxplot(data=filtered_df, x="rigidity_quartile", y="num_waters", palette=color_map, order=category_order)

plt.xlabel("Rigidity Category")
plt.ylabel("Number of Waters per PDB")
plt.title("Water Counts per PDB by Rigidity Category (Outliers Removed)")
plt.xticks(rotation=15)
plt.tight_layout()

# === Save plot ===
plt.savefig("water_counts_by_rigidity_category_no_outliers.png", dpi=300)
