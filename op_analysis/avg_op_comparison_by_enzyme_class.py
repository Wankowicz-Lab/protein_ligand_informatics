import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

# === Paths ===
cluster_files = {
    "30": ".../cluster_baseline_summary_with_distributions_30.csv",
    "60": ".../cluster_baseline_summary_with_distributions_60.csv",
    "90": ".../cluster_baseline_summary_with_distributions_90.csv",
}
uniprot_ec_path = ".../apo_merged_uniprot_data.csv"
output_dir = "..."
os.makedirs(output_dir, exist_ok=True)

# === Load UniProt to EC mapping ===
uniprot_ec = pd.read_csv(uniprot_ec_path)
uniprot_ec["uniprot_id"] = uniprot_ec["uniprot_id"].str.strip()
uniprot_ec["enzyme_number"] = uniprot_ec["enzyme_number"].astype(str).str.strip()
uniprot_to_ec = dict(zip(uniprot_ec["uniprot_id"], uniprot_ec["enzyme_number"]))

# === Main loop for EC class distributions ===
for id_level, path in cluster_files.items():
    df = pd.read_csv(path)
    df["avg_s2calc_distribution"] = df["avg_s2calc_distribution"].astype(str)
    df["uniprot_ids"] = df["uniprot_ids"].astype(str)

    records = []
    for _, row in df.iterrows():
        cluster_id = row["cluster_id"]
        s2_vals = [float(x.strip()) for x in row["avg_s2calc_distribution"].split(";") if x.strip()]
        uniprots = [u.strip() for u in row["uniprot_ids"].split(",") if u.strip()]

        ec_classes = []
        for u in uniprots:
            ec = uniprot_to_ec.get(u)
            if ec and ec != "nan":
                try:
                    ec_class = int(ec.split(".")[0])
                    ec_classes.append(ec_class)
                except:
                    continue

        for ec_class in ec_classes:
            for s2 in s2_vals:
                records.append({
                    "cluster_id": cluster_id,
                    "s2calc": s2,
                    "ec_class": f"EC {ec_class}"
                })

    long_df = pd.DataFrame(records)

    # === Save CSV of records ===
    csv_out_path = os.path.join(output_dir, f"cluster_s2calc_by_ec_class_{id_level}.csv")
    long_df.to_csv(csv_out_path, index=False)

    # === Boxplot ===
    plt.figure(figsize=(10, 6))
    ec_order = sorted(long_df["ec_class"].unique(), key=lambda x: int(x.split(" ")[1]))
    sns.boxplot(data=long_df, x="ec_class", y="s2calc", order=ec_order)
    plt.ylabel("Order Parameter", fontsize=16)
    plt.xlabel("EC Class", fontsize=16)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"s2calc_by_ec_class_{id_level}.png"), dpi=300)
    plt.close()
