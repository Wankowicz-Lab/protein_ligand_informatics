import os
from Bio.PDB import PDBParser
import pandas as pd
from collections import defaultdict

# === SETTINGS ===
binding_csv = '/dors/wankowicz_lab/ellas/main_dataset/updated_main_dataset_binding_op.csv'
pdb_dir = '/dors/wankowicz_lab/ellas/dsasa_results'
output_csv = '/dors/wankowicz_lab/ellas/binding_site_sasa_per_residue.csv'

# === ELEMENT CLASSIFICATION ===
polar_elements = {'O', 'N', 'S'}
apolar_elements = {'C'}

# === Load binding site residue list ===
binding_df = pd.read_csv(binding_csv)
binding_df['resi'] = binding_df['resi'].astype(str)

# === Initialize data collection ===
residue_data = []

# === Process each unique PDB ===
for pdb_id in binding_df['pdb_id'].unique():
    print(f"Processing {pdb_id}...")

    pdb_file = os.path.join(pdb_dir, f"{pdb_id}_asa_nowater.pdb")
    if not os.path.isfile(pdb_file):
        print(f"Missing PDB file for {pdb_id}, skipping.")
        continue

    # Parse PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, pdb_file)

    bs_residues = binding_df[binding_df['pdb_id'] == pdb_id]

    for model in structure:
        for chain in model:
            for residue in chain:
                res_id = residue.get_id()[1]
                chain_id = chain.id
                if not ((bs_residues['chain'] == chain_id) & (bs_residues['resi'] == str(res_id))).any():
                    continue

                polar_sasa = 0.0
                apolar_sasa = 0.0
                total_sasa = 0.0
                polar_atoms = 0
                apolar_atoms = 0

                for atom in residue:
                    sasa = atom.bfactor
                    element = atom.element.strip().upper()
                    total_sasa += sasa
                    if element in polar_elements:
                        polar_sasa += sasa
                        polar_atoms += 1
                    elif element in apolar_elements:
                        apolar_sasa += sasa
                        apolar_atoms += 1

                residue_data.append({
                    'pdb_id': pdb_id,
                    'chain': chain_id,
                    'resi': res_id,
                    'resname': residue.get_resname(),
                    'total_sasa': total_sasa,
                    'polar_sasa': polar_sasa,
                    'apolar_sasa': apolar_sasa,
                    'polar_atoms': polar_atoms,
                    'apolar_atoms': apolar_atoms
                })

# === Save full residue-level results ===
df = pd.DataFrame(residue_data)
df.to_csv(output_csv, index=False)
print(f"Full residue data saved to: {output_csv}")

# === Rank amino acids by average polar SASA ===
aa_group = df.groupby('resname').agg({
    'polar_sasa': 'mean',
    'apolar_sasa': 'mean',
    'total_sasa': 'mean'
}).sort_values('polar_sasa', ascending=False)

print("\n=== Amino Acid Ranking by Polar SASA (Binding Sites Only) ===")
print(aa_group[['polar_sasa']])

