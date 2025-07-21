from Bio.PDB import PDBParser, PPBuilder
from Bio import SeqIO
import os

pdb_list_path = "/dors/wankowicz_lab/ellas/main_dataset/updated_main_dataset.txt"
output_fasta = "/dors/wankowicz_lab/ellas/main_dataset/all_pdbs.fasta"
pdb_root_dirs = [
    "1_10000", "10001_20000", "20001_30000", "30001_40000",
    "40001_50000", "50001_60000", "60001_70000", "70001_80000", "80001_end"
]
base_path = "/dors/wankowicz_lab/all_pdb"

parser = PDBParser(QUIET=True)
ppb = PPBuilder()

with open(pdb_list_path) as f:
    pdb_ids = [line.strip().lower() for line in f if line.strip()]

with open(output_fasta, "w") as fasta_out:
    for pdb_id in pdb_ids:
        pdb_id_upper = pdb_id.upper()
        found = False

        for subfolder in pdb_root_dirs:
            pdb_folder = os.path.join(base_path, subfolder, pdb_id_upper)
            pdb_path = os.path.join(pdb_folder, f"{pdb_id_upper}.pdb")
            pdb_updated_path = os.path.join(pdb_folder, f"{pdb_id_upper}_updated.pdb")

            if os.path.isfile(pdb_path):
                final_path = pdb_path
                found = True
                break
            elif os.path.isfile(pdb_updated_path):
                final_path = pdb_updated_path
                found = True
                break

        if not found:
            print(f"Missing: {pdb_id} in any allowed subfolder")
            continue

        try:
            structure = parser.get_structure(pdb_id, final_path)
            for model in structure:
                for chain in model:
                    seqs = ppb.build_peptides(chain)
                    if seqs:
                        seq = ''.join(str(peptide.get_sequence()) for peptide in seqs)
                        fasta_out.write(f">{pdb_id}_{chain.id}\n{seq}\n")
        except Exception as e:
            print(f"Error parsing {pdb_id}: {e}")
