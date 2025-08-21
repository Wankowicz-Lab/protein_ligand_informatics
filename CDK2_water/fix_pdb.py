"""
Written by microsoft copilot.
"""
from openmm.app import PDBFile
from pdbfixer import PDBFixer

from glob import glob
from pathlib import Path

import matplotlib.pyplot as plt
from Bio import PDB

# pip install pdb-tools
import subprocess
import sys


# utility functions to trim residues before methionine that get added from pdb fixer
def get_first_met_residue_number(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() == "MET":
                    return residue.get_id()[1]  # Residue number
    return None


def trim_pdb_from_met(pdb_file, met_resnum, output_file):
    cmd_selchain = f"pdb_selchain -A {pdb_file}"
    cmd_selres = f"pdb_selres -{met_resnum}: "
    cmd_keep_het = f"pdb_selhetatm {pdb_file}"
    
    p1 = subprocess.Popen(cmd_selchain, stdout=subprocess.PIPE)
    prot_lines = subprocess.run(cmd_selres, stdin=p1.stdout, capture_output=True, text=True)
    lines = prot_lines.stdout.splitlines()[:-1]
    p1.stdout.close()

    het_lines = subprocess.run(cmd_keep_het, capture_output=True, text=True)
    lines += het_lines.stdout.splitlines()
    lines += ["END"]
    
    with open(output_file, "w") as out:
        file_str = "\n".join(lines)
        out.write(file_str)
        

def renumber_pdb(pdb_file, output_file):
    cmd_reres = f"pdb_reres -1 {pdb_file}"
    cmd_reatom = "pdb_reatom"

    # Pipe reres → reatom
    with open(output_file, "w") as out:
        p1 = subprocess.Popen(cmd_reres, stdout=subprocess.PIPE)
        subprocess.run(cmd_reatom, stdin=p1.stdout, stdout=out)
        p1.stdout.close()


def standardize_pdb(pdb_path, output_path):
    met_resnum = get_first_met_residue_number(pdb_path)    
    if met_resnum is None:
        raise ValueError("No methionine found in the structure.")
    
    trimmed_path = 'temp.pdb'
    trim_pdb_from_met(pdb_path, met_resnum, trimmed_path)
    renumber_pdb(trimmed_path, output_path)
    Path.unlink(trimmed_path)


DATA_FOLDER = "../PDB_Data_cdk2/cdk2_expanded/"
pdb_path = Path(DATA_FOLDER)
pdb_files = list( pdb_path.glob('*_final.pdb') )
parser = PDB.PDBParser(QUIET=True)
all_prots = []
  
  
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16,8), sharey = True)
for idx, file_name in enumerate(pdb_files):
    prot = file_name.stem.split('_')[0]
    all_prots.append( prot )
    
    fixer = PDBFixer( filename = f"{file_name}" )
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    
    new_file = f"{DATA_FOLDER}/{prot}_cleaned.pdb"
    PDBFile.writeFile( fixer.topology, fixer.positions, open(new_file, 'w') )
    standardize_pdb(new_file, new_file)
    print(f"Done fixing {prot}")
    
    # display the new numbering to ensure evenness
    resnums = []
    solnums = []
    structure = parser.get_structure(prot, new_file)
    for model in structure:
        for chain in model:
            resnums += [res.get_id()[1] for res in chain if res.get_id()[0] == ' ']
            solnums += [res.get_id()[1] for res in chain if res.get_id()[0] != ' ']
                    
    ax1.plot(resnums, [idx]*len(resnums), "o", ms=2)
    ax2.plot(solnums, [idx]*len(solnums), "o", ms=2)
    
    
ax1.set_yticks( range(len(all_prots)), all_prots, fontsize = 8 )
plt.show()
