import pandas as pd
import glob
import os
import matplotlib.pyplot as plt
import seaborn as sns

# === Load top 3 scores for each group ===

def load_scores_from_list(pdb_list_file):
    base_path = '/dors/wankowicz_lab/ellas/p2rank_2.5/p2rank_results'
    all_scores = []

    with open(pdb_list_file, 'r') as f:
        pdb_list = [line.strip() for line in f if line.strip()]

    for pdb in pdb_list:
        pattern = os.path.join(base_path, '**', f'{pdb}_clean.pdb_predictions.csv')
        matches = glob.glob(pattern, recursive=True)

        if not matches:
            continue

        for file in matches:
            try:
                df = pd.read_csv(file)
                df.columns = df.columns.str.strip()
                if 'score' in df.columns:
                    scores = df['score'].dropna().tolist()
                    if len(scores) > 3:
                        scores = sorted(scores, reverse=True)[:3]
                    all_scores.extend(scores)
            except Exception:
                continue
    return pd.Series(all_scores)

def load_all_scores():
    all_scores = []
    all_files = glob.glob("/dors/wankowicz_lab/ellas/p2rank_2.5/p2rank_results/**/*_clean.pdb_predictions.csv", recursive=True)

    for file in all_files:
        try:
            df = pd.read_csv(file)
            df.columns = df.columns.str.strip()
            if 'score' in df.columns:
                scores = df['score'].dropna().tolist()
                if len(scores) > 3:
                    scores = sorted(scores, reverse=True)[:3]
                all_scores.extend(scores)
        except Exception:
            continue
    return pd.Series(all_scores)

# === Load data for all groups ===

apo_scores = load_scores_from_list('/dors/wankowicz_lab/ellas/apo_structures_pdb_list.txt')
holo_scores = load_scores_from_list('/dors/wankowicz_lab/ellas/num_atoms/pdbs_15_to_100_atoms.txt')
all_scores = load_all_scores()

# === Combine into one DataFrame ===

df = pd.DataFrame({
    'score': pd.concat([apo_scores, holo_scores, all_scores], ignore_index=True),
    'group': (['apo'] * len(apo_scores)) +
             (['holo'] * len(holo_scores)) +
             (['all'] * len(all_scores))
})

# === Set consistent color palette ===
palette = {'apo': '#1f77b4', 'holo': '#2ca02c', 'all': '#ff7f0e'}  # blue, green, orange

# === 1. Full Histogram ===
plt.figure(figsize=(10, 6))
sns.histplot(data=df, x='score', hue='group', bins=50, kde=True, palette=palette, multiple='stack')
plt.xlabel("Score")
plt.ylabel("Frequency")
plt.title("Comparison of Top 3 P2Rank Scores (Full Data)")
plt.grid(True)
plt.savefig('/dors/wankowicz_lab/ellas/p2rank_2.5/p2rank_results/combined_histogram_full.png')
plt.close()

# === 2. 3σ Cutoff Histogram ===
cutoff_3sigma = df.groupby('group')['score'].transform(lambda x: x.mean() + 3 * x.std())
df_3sigma = df[df['score'] <= cutoff_3sigma]

plt.figure(figsize=(10, 6))
sns.histplot(data=df_3sigma, x='score', hue='group', bins=50, kde=True, palette=palette, multiple='stack')
plt.xlabel("Score")
plt.ylabel("Frequency")
plt.title("Comparison of Top 3 P2Rank Scores (≤ μ + 3σ)")
plt.grid(True)
plt.savefig('/dors/wankowicz_lab/ellas/p2rank_2.5/p2rank_results/combined_histogram_3sigma.png')
plt.close()

# === 3. 95th Percentile Cutoff Histogram ===
cutoff_95 = df.groupby('group')['score'].transform(lambda x: x.quantile(0.95))
df_95 = df[df['score'] <= cutoff_95]

plt.figure(figsize=(10, 6))
sns.histplot(data=df_95, x='score', hue='group', bins=50, kde=True, palette=palette, multiple='stack')
plt.xlabel("Score")
plt.ylabel("Frequency")
plt.title("Comparison of Top 3 P2Rank Scores (≤ 95th Percentile)")
plt.grid(True)
plt.savefig('/dors/wankowicz_lab/ellas/p2rank_2.5/p2rank_results/combined_histogram_95percentile.png')
plt.close()

print("Histograms generated.")
