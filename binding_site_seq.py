import os
import pandas as pd
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from pathlib import Path

# === Inputs ===
binding_site_csv = "/dors/wankowicz_lab/ellas/main_dataset/all_binding_residues.csv"
pdb_list_path = "/dors/wankowicz_lab/ellas/main_dataset/updated_main_dataset.txt"
base_path = "/dors/wankowicz_lab/all_pdb"
pdb_root_dirs = [
    "1_10000", "10001_20000", "20001_30000", "30001_40000",
    "40001_50000", "50001_60000", "60001_70000", "70001_80000", "80001_end"
]
output_fasta = "/dors/wankowicz_lab/ellas/main_dataset/binding_site_seqs.fasta"

# === Load binding site info ===
df = pd.read_csv(binding_site_csv, dtype=str, low_memory=False)
df['pdb_id'] = df['pdb_id'].str.lower()
df['residue_label'] = df['residue_label'].str.strip()

# Extract numeric part and insertion code
df['resnum'] = df['residue_label'].str.extract(r"(\d+)")[0].astype(int)
df['insertion'] = df['residue_label'].str.extract(r"([A-Za-z]?)")[0].fillna('').apply(lambda x: x if x else ' ')

# === Set up parser ===
parser = PDBParser(QUIET=True)

# === Group by PDB ID and chain ===
grouped = df.groupby(['pdb_id', 'chain'])

# === Load all pdb ids from your list ===
with open(pdb_list_path, "r") as f:
    valid_pdb_ids = set(line.strip().lower() for line in f if line.strip())

with open(output_fasta, "w") as fasta_out:
    for (pdb_id, chain_id), group in grouped:
        if pdb_id not in valid_pdb_ids:
            continue

        pdb_id_upper = pdb_id.upper()
        final_path = None

        for subfolder in pdb_root_dirs:
            pdb_folder = os.path.join(base_path, subfolder, pdb_id_upper)
            pdb_path = os.path.join(pdb_folder, f"{pdb_id_upper}.pdb")
            pdb_updated_path = os.path.join(pdb_folder, f"{pdb_id_upper}_updated.pdb")

            if os.path.isfile(pdb_path):
                final_path = pdb_path
                break
            elif os.path.isfile(pdb_updated_path):
                final_path = pdb_updated_path
                break

        if final_path is None:
            print(f"PDB file for {pdb_id} not found.")
            continue

        # Parse PDB
        try:
            structure = parser.get_structure(pdb_id, final_path)
            chain = structure[0][chain_id]
        except Exception as e:
            print(f"Error parsing {pdb_id} chain {chain_id}: {e}")
            continue

        # Extract binding residues
        binding_residues = []
        for _, row in group.iterrows():
            res_num = row['resnum']
            ins_code = row['insertion']
            try:
                residue = chain[(' ', res_num, ins_code)]
                resname = residue.get_resname().strip().upper()
                aa = seq1(resname)
                binding_residues.append(aa)
            except Exception as e:
                print(f"Skipping {pdb_id} {chain_id} {res_num}{ins_code}: {e}")

        # Output FASTA
        if binding_residues:
            sequence = ''.join(binding_residues)
            fasta_out.write(f">{pdb_id}_{chain_id}\n{sequence}\n")
