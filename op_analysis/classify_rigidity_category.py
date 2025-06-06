import pandas as pd
import os

# === CONFIG ===
input_csv = "/dors/wankowicz_lab/ellas/main_dataset/updated_main_dataset_binding_op.csv"
output_dir = "/dors/wankowicz_lab/ellas/main_dataset"
os.makedirs(output_dir, exist_ok=True)

# === LOAD DATA ===
df = pd.read_csv(input_csv)

# === CATEGORIZE RIGIDITY ===
def rigidity_category(op):
    if op > 0.873404:
        return "very rigid"
    elif op > 0.822762:
        return "rigid"
    elif op > 0.756902:
        return "flexible"
    else:
        return "very flexible"

df["rigidity_category"] = df["s2calc"].apply(rigidity_category)

# === AA COUNTS AND PERCENTAGES PER CATEGORY ===
aa_counts = df.groupby(["rigidity_category", "resn"]).size().reset_index(name="count")

# Get total residues per category
total_counts = df.groupby("rigidity_category").size().reset_index(name="total")

# Merge to calculate percentage
aa_counts = aa_counts.merge(total_counts, on="rigidity_category")
aa_counts["percentage"] = (aa_counts["count"] / aa_counts["total"]) * 100

# Save result
aa_counts[["rigidity_category", "resn", "count", "percentage"]].to_csv(
    f"{output_dir}/aa_counts_percent_by_category.csv", index=False
)

print(f"Saved counts and percentages to {output_dir}/aa_counts_percent_by_category.csv")
