import pandas as pd
import glob
import os
from joblib import Parallel, delayed

# === Load top 3 scores for each group ===

def process_file(file, group_name):
    data = []
    try:
        pdb_id = os.path.basename(file).split('_clean.pdb')[0]
        df = pd.read_csv(file)
        df.columns = df.columns.str.strip()
        if 'score' in df.columns:
            scores = df['score'].dropna().tolist()
            if len(scores) > 3:
                scores = sorted(scores, reverse=True)[:3]
            for score in scores:
                data.append((pdb_id, score, group_name))
    except Exception:
        pass
    return data

def load_scores_from_list(pdb_list_file, group_name):
    base_path = '/dors/wankowicz_lab/ellas/p2rank_2.5/p2rank_results'
    with open(pdb_list_file, 'r') as f:
        pdb_list = [line.strip() for line in f if line.strip()]

    patterns = [os.path.join(base_path, '**', f'{pdb}_clean.pdb_predictions.csv') for pdb in pdb_list]
    matches = [file for pattern in patterns for file in glob.glob(pattern, recursive=True)]

    # Parallel processing using Joblib
    results = Parallel(n_jobs=-1)(delayed(process_file)(file, group_name) for file in matches)
    data = [item for sublist in results for item in sublist]  # Flatten the list of lists

    return pd.DataFrame(data, columns=['pdb_id', 'score', 'group'])

def load_all_scores():
    all_files = glob.glob("/dors/wankowicz_lab/ellas/p2rank_2.5/p2rank_results/**/*_clean.pdb_predictions.csv", recursive=True)

    # Parallel processing for "all" group
    results = Parallel(n_jobs=-1)(delayed(process_file)(file, 'all') for file in all_files)
    data = [item for sublist in results for item in sublist]  # Flatten the list of lists

    return pd.DataFrame(data, columns=['pdb_id', 'score', 'group'])

# === Load data for all groups ===

apo_df = load_scores_from_list('/dors/wankowicz_lab/ellas/apo_structures_pdb_list.txt', 'apo')
holo_df = load_scores_from_list('/dors/wankowicz_lab/ellas/num_atoms/pdbs_15_to_100_atoms.txt', 'holo')
all_df = load_all_scores()

# === Save CSV files ===
apo_df.to_csv('/dors/wankowicz_lab/ellas/p2rank_2.5/p2rank_results/apo_scores.csv', index=False)
holo_df.to_csv('/dors/wankowicz_lab/ellas/p2rank_2.5/p2rank_results/holo_scores.csv', index=False)
all_df.to_csv('/dors/wankowicz_lab/ellas/p2rank_2.5/p2rank_results/all_scores.csv', index=False)

print("CSV files generated.")

# === Generate Statistics ===
def save_stats(df, group_name):
    stats = df['score'].describe()
    mean = stats['mean']
    std = stats['std']
    percentile_95 = df['score'].quantile(0.95)

    # Calculate 3σ cutoff range
    lower_3sigma = mean - 3 * std
    upper_3sigma = mean + 3 * std

    output_path = f'/dors/wankowicz_lab/ellas/p2rank_2.5/p2rank_results/{group_name}_stats.txt'
    with open(output_path, 'w') as f:
        f.write(f"Statistics for {group_name} group\n")
        f.write(f"Total PDBs: {df['pdb_id'].nunique()}\n")
        f.write(f"Total Scores: {len(df)}\n\n")
        f.write(f"Mean: {mean}\n")
        f.write(f"Median: {stats['50%']}\n")
        f.write(f"Standard Deviation: {std}\n")
        f.write(f"Min: {stats['min']}\n")
        f.write(f"25th Percentile: {stats['25%']}\n")
        f.write(f"50th Percentile (Median): {stats['50%']}\n")
        f.write(f"75th Percentile: {stats['75%']}\n")
        f.write(f"Max: {stats['max']}\n")
        f.write(f"Range: {stats['max'] - stats['min']}\n")
        f.write(f"95th Percentile: {percentile_95}\n")
        f.write(f"3σ Lower Bound: {lower_3sigma}\n")
        f.write(f"3σ Upper Bound: {upper_3sigma}\n")
        f.write("\n")

    print(f"Statistics saved to: {output_path}")

# === Save Stats for Each Group ===
save_stats(apo_df, 'apo')
save_stats(holo_df, 'holo')
save_stats(all_df, 'all')

print("Statistics text files generated.")
