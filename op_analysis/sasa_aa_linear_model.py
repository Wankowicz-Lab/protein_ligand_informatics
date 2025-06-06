import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt
import seaborn as sns

# === INPUT FILES ===
s2_avg_file = "/dors/wankowicz_lab/ellas/main_dataset/pdb_avg_s2calc_z.csv"
sasa_file = "/dors/wankowicz_lab/ellas/main_dataset/relative_exposure_with_zscore.csv"
aa_raw_file = "/dors/wankowicz_lab/ellas/main_dataset/all_binding_residues.csv"

# === LOAD AA COMPOSITION FILE ===
aa_df = pd.read_csv(aa_raw_file, dtype=str)
aa_df['residue_name'] = aa_df['residue_name'].str.upper().str.strip()
aa_df['pdb_id'] = aa_df['pdb_id'].str.upper() 

# Keep only standard amino acids
canonical_aas = [
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
]
aa_df = aa_df[aa_df['residue_name'].isin(canonical_aas)]

# === 20 AA COUNT PER PDB ===
aa_counts = aa_df.groupby(['pdb_id', 'residue_name']).size().unstack(fill_value=0)

# === LOAD S2CALC (averaged) ===
s2_avg = pd.read_csv(s2_avg_file).rename(columns={'avg_s2calc_z': 's2calc'})
s2_avg['pdb_id'] = s2_avg['pdb_id'].str.upper()

# === LOAD SASA FILE ===
sasa_df = pd.read_csv(sasa_file)[['pdb_id', 'relative_exposure_z']]
sasa_df['pdb_id'] = sasa_df['pdb_id'].str.upper()

# === MERGE ALL DATA ===
merged_df = s2_avg.merge(sasa_df, on='pdb_id').merge(aa_counts, on='pdb_id')

# === DROP MISSING DATA AND PRINT WHAT'S LOST ===
before = merged_df.shape[0]
dropped_pdbs = merged_df[merged_df.isnull().any(axis=1)]['pdb_id'].tolist()
merged_df = merged_df.dropna()
after = merged_df.shape[0]

print(f"\nDropped {before - after} PDBs due to missing data.")
if dropped_pdbs:
    print("Dropped PDB IDs:", dropped_pdbs)

# === SET PREDICTORS ===
aa_cols = aa_counts.columns.tolist()
X = merged_df[aa_cols + ['relative_exposure_z']]
X = sm.add_constant(X)

# === RESPONSE VARIABLE: S2CALC AVERAGE ===
y = merged_df['s2calc'].values.ravel()

# === DEBUG SHAPES ===
print("Final AA columns in X:", aa_counts.columns.tolist())
print("Merged dataframe shape:", merged_df.shape)
print(f"\nX shape: {X.shape}")
print(f"y shape: {y.shape}")
print(f"X columns: {list(X.columns)}")
assert y.shape[0] == X.shape[0], "Mismatch between X and y row counts"

# === FIT MODEL ===
model = sm.OLS(y, X).fit()

# === OUTPUT MODEL SUMMARY ===
with open("/dors/wankowicz_lab/ellas/main_dataset/lm_model_summary.txt", "w") as f:
    f.write(model.summary().as_text())

# === COEFFICIENT BARPLOT ===
coef_df = model.params.drop('const').sort_values().reset_index()
coef_df.columns = ['feature', 'coefficient']

plt.figure(figsize=(10, 6))
sns.barplot(data=coef_df, x='coefficient', y='feature', orient='h')
plt.title('Linear Model Coefficients: s2calc ~ SASA + AA Composition')
plt.tight_layout()
plt.savefig("/dors/wankowicz_lab/ellas/main_dataset/lm_coefficients_plot.png")

# === SAVE MERGED DATA AND COEFFICIENTS ===
merged_df.to_csv("/dors/wankowicz_lab/ellas/main_dataset/linear_model_input_data.csv", index=False)
coef_df.to_csv("/dors/wankowicz_lab/ellas/main_dataset/lm_coefficients.csv", index=False)

print("\nLinear model complete. Output files:")
print("- Model summary: lm_model_summary.txt")
print("- Coefficient plot: lm_coefficients_plot.png")
print("- Coefficients CSV: lm_coefficients.csv")
print("- Input data CSV: linear_model_input_data.csv")
