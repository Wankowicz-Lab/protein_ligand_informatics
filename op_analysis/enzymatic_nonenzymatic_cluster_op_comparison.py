import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import scipy.stats as stats

# === Paths ===
cluster_files = {
    "30": "/.../cluster_baseline_summary_with_distributions_30.csv",
    "60": "/.../cluster_baseline_summary_with_distributions_60.csv",
    "90": "/.../cluster_baseline_summary_with_distributions_90.csv",
}
uniprot_ec_path = "/.../apo_merged_uniprot_data.csv"
output_dir = "/.../enzymatic_cluster_analysis"
os.makedirs(output_dir, exist_ok=True)

# === Load UniProt to EC mapping ===
uniprot_ec = pd.read_csv(uniprot_ec_path)
uniprot_ec["uniprot_id"] = uniprot_ec["uniprot_id"].str.strip()
uniprot_ec["enzyme_number"] = uniprot_ec["enzyme_number"].astype(str).str.strip()

uniprot_to_ec = dict(zip(uniprot_ec["uniprot_id"], uniprot_ec["enzyme_number"]))

# === Color settings ===
palette = {"Enzymatic": "#4682B4", "Non-enzymatic": "#B0B0B0"}

# === Helper ===
def is_enzymatic(uniprot_list):
    for u in uniprot_list:
        if u in uniprot_to_ec:
            if uniprot_to_ec[u] and uniprot_to_ec[u] != "nan":
                return True
        else:
            print(f"[WARNING] UniProt ID not found in EC mapping: {u}")
    return False

# === Main loop ===
for id_level, path in cluster_files.items():
    df = pd.read_csv(path)
    
    # Clean and parse lists
    df["avg_s2calc_distribution"] = df["avg_s2calc_distribution"].astype(str)
    df["uniprot_ids"] = df["uniprot_ids"].astype(str)

    records = []
    for _, row in df.iterrows():
        cluster_id = row["cluster_id"]
        s2_vals = [float(x.strip()) for x in row["avg_s2calc_distribution"].split(";") if x.strip()]
        uniprots = [u.strip() for u in row["uniprot_ids"].split(",") if u.strip()]
        label = "Enzymatic" if is_enzymatic(uniprots) else "Non-enzymatic"
        for s2 in s2_vals:
            records.append({"cluster_id": cluster_id, "s2calc": s2, "enzymatic_status": label})

    long_df = pd.DataFrame(records)

    # === Statistical test ===
    enz = long_df[long_df["enzymatic_status"] == "Enzymatic"]["s2calc"]
    non_enz = long_df[long_df["enzymatic_status"] == "Non-enzymatic"]["s2calc"]
    stat, p = stats.mannwhitneyu(enz, non_enz, alternative="two-sided")
    print(f"\n===== Identity {id_level}% =====")
    print(f"Mann-Whitney U test: U = {stat:.3f}, p = {p:.4e}")
    print(f"Median s²calc - Enzymatic: {enz.median():.4f}, Non-enzymatic: {non_enz.median():.4f}")
    print(f"Count - Enzymatic: {len(enz)}, Non-enzymatic: {len(non_enz)}")

    # === Plot ===
    plt.figure(figsize=(6, 5))
    sns.boxplot(
        data=long_df,
        x="enzymatic_status",
        y="s2calc",
        palette=palette
    )
    plt.xlabel("")
    plt.ylabel("Order Parameter", fontsize=16)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"s2calc_enzymatic_comparison_{id_level}.png"), dpi=300)
    plt.close()
