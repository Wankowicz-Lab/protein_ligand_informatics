import pandas as pd
import matplotlib.pyplot as plt

# Load your input CSV
df = pd.read_csv('/.../apo_merged_uniprot_data.csv')  # replace with actual path

# Group by uniprot_id and aggregate
grouped = df.groupby('uniprot_id')['pdb_id'].agg(['count', lambda x: ','.join(sorted(set(x)))])
grouped.columns = ['frequency', 'pdb_ids']
grouped = grouped.reset_index()

# Save full CSV with frequency and pdb list
grouped.to_csv('/.../uniprot_summary.csv', index=False)

# Plot top 10 most frequent uniprot_ids
top10 = grouped.sort_values('frequency', ascending=False).head(20)
plt.figure(figsize=(12, 6))
plt.bar(top10['uniprot_id'], top10['frequency'])
plt.xlabel('UniProt ID')
plt.ylabel('Frequency (# of PDBs)')
plt.title('Top 10 Most Frequent UniProt IDs by PDB Count')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig('/.../top20_uniprot_barplot.png')

print("CSV and bar plot generated successfully.")
