import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# === File paths ===
noncontrolled_csv = "/dors/wankowicz_lab/ellas/main_dataset/updated_single_uniprot_dataset_binding_op.csv"
noncontrolled_nonbind_csv = "/dors/wankowicz_lab/ellas/main_dataset/updated_single_uniprot_dataset_nonbinding_op.csv"
controlled_csv = "/dors/wankowicz_lab/ellas/main_dataset/updated_single_uniprot_dataset_sasa_op_merged_binding.csv"
controlled_nonbind_csv = "/dors/wankowicz_lab/ellas/main_dataset/updated_single_uniprot_dataset_sasa_op_merged_nonbinding.csv"
output_path = "/dors/wankowicz_lab/ellas"

# === Load data ===
df_binding_nc = pd.read_csv(noncontrolled_csv)
df_nonbinding_nc = pd.read_csv(noncontrolled_nonbind_csv)
df_binding_c = pd.read_csv(controlled_csv)
df_nonbinding_c = pd.read_csv(controlled_nonbind_csv)

# === Non-controlled analysis ===
bind_avg_nc = df_binding_nc.groupby('pdb_id')['s2calc'].mean().reset_index()
bind_avg_nc.columns = ['pdb_id', 'avg_s2calc_binding']
nonbind_avg_nc = df_nonbinding_nc.groupby('pdb_id')['s2calc'].mean().reset_index()
nonbind_avg_nc.columns = ['pdb_id', 'avg_s2calc_nonbinding']
merged_nc = pd.merge(bind_avg_nc, nonbind_avg_nc, on='pdb_id')
merged_nc['delta_s2calc'] = merged_nc['avg_s2calc_binding'] - merged_nc['avg_s2calc_nonbinding']
merged_nc['group'] = 'Non-controlled'

# === Controlled analysis ===
for df in [df_binding_c, df_nonbinding_c]:
    df['pdb_id'] = df['pdb_id'].astype(str)
    df['chain'] = df['chain'].astype(str)
    df['resn'] = df['resn'].astype(str)
    df['solvent_exposure'] = pd.to_numeric(df['solvent_exposure'], errors='coerce')
    df.dropna(subset=['solvent_exposure'], inplace=True)
    df['solvent_exposure'] = df['solvent_exposure'].astype(int)

df_binding_c = df_binding_c[df_binding_c['solvent_exposure'].isin([0, 1])]
df_nonbinding_c = df_nonbinding_c[df_nonbinding_c['solvent_exposure'].isin([0, 1])]

df_binding_c['merge_key'] = df_binding_c['pdb_id'] + "_" + df_binding_c['chain'] + "_" + df_binding_c['resn'] + "_" + df_binding_c['solvent_exposure'].astype(str)
df_nonbinding_c['merge_key'] = df_nonbinding_c['pdb_id'] + "_" + df_nonbinding_c['chain'] + "_" + df_nonbinding_c['resn'] + "_" + df_nonbinding_c['solvent_exposure'].astype(str)

df_nonbinding_sampled = df_nonbinding_c.groupby('merge_key').apply(lambda x: x.sample(n=1, random_state=42)).reset_index(drop=True)
merged_c = pd.merge(
    df_binding_c[['merge_key', 's2calc']],
    df_nonbinding_sampled[['merge_key', 's2calc']],
    on='merge_key',
    suffixes=('_binding', '_non_binding')
)
merged_c['delta_s2calc'] = merged_c['s2calc_binding'] - merged_c['s2calc_non_binding']
merged_c['group'] = 'Controlled'

# === Combine datasets ===
combined_delta = pd.concat([
    merged_nc[['delta_s2calc', 'group']],
    merged_c[['delta_s2calc', 'group']]
], ignore_index=True)

# === Plot ===
sns.set(style="whitegrid", context="talk")
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(6.4, 8), sharex=True)

sns.boxplot(
    ax=axes[0],
    x=combined_delta[combined_delta['group'] == 'Non-controlled']['delta_s2calc'],
    color="#999933", linewidth=2.0, fliersize=3, orient='h'
)
axes[0].set_title("Box Plot of Delta s2calc (Binding - Non-binding)", fontsize=14)
axes[0].set_ylabel("Non-controlled", fontsize=12)
axes[0].set_xlabel("")

sns.boxplot(
    ax=axes[1],
    x=combined_delta[combined_delta['group'] == 'Controlled']['delta_s2calc'],
    color="#882255", linewidth=2.0, fliersize=3, orient='h'
)
axes[1].set_ylabel("Controlled", fontsize=12)
axes[1].set_xlabel("Delta s2calc", fontsize=12)

for ax in axes:
    ax.tick_params(labelsize=11)

plt.tight_layout()
plt.savefig(f"{output_path}/updated_single_uniprot_dataset_delta_s2calc_boxplot_controlled_vs_noncontrolled.png", dpi=300, bbox_inches="tight")
plt.close()
print("Boxplot figure saved.")

# === Print Descriptive Statistics ===
def describe_group(name, data):
    stats = data['delta_s2calc'].describe()
    print(f"\n=== {name} Statistics ===")
    print(stats)
    print(f"Median: {data['delta_s2calc'].median():.4f}")
    print(f"Standard Deviation: {data['delta_s2calc'].std():.4f}")
    print(f"Min: {data['delta_s2calc'].min():.4f}")
    print(f"Max: {data['delta_s2calc'].max():.4f}")

describe_group("Non-controlled", combined_delta[combined_delta['group'] == 'Non-controlled'])
describe_group("Controlled", combined_delta[combined_delta['group'] == 'Controlled'])
