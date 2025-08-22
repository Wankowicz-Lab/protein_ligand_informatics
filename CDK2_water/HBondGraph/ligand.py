"""
Gets ligands from PDB and calculates chemical properties.
"""

from Bio import PDB
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

import requests

# UNIVERSAL CONSTANTS
WATER=['HOH', 'WAT', 'H2O']  
OTHER_SOLVENT=["EDO", "PO4"] # constructed by looking through all non WATER + AA or ligand labels


def get_ligand_from_structure(structure):
    """
    Gets the heteroatom not labeled as an amino acid or a solvent.
    """
    for model in structure:
        for chain in model:
            for res in chain:
                resname = res.get_resname()
                hetflag = res.get_id()[0]
                if hetflag != ' ' and resname not in WATER + OTHER_SOLVENT:
                    return res 


def residue_to_mol(residue):
    """
    Written by microsoft copilot.
    Takes a residue and writes a mol object.
    """

    if residue is None:
        return ""

    ligand_id = residue.get_resname()
    url = f"https://files.rcsb.org/ligands/download/{ligand_id}_ideal.sdf"
    response = requests.get(url)
    if response.status_code == 200:
        sdf_data = response.text
        suppl = Chem.SDMolSupplier()
        suppl.SetData(sdf_data)
        mol = suppl[0]
        return mol 
    else:
        raise ValueError(f"Ligand {ligand_id} not found at RCSB")
        
    
def tanimoto_dissimilarity(mol1, mol2):
    """
    Takes two mol objects and outputs the tanimoto dissimilarity between them
    """

    gen = Chem.rdFingerprintGenerator.GetRDKitFPGenerator(
        maxPath=7,       # max bond path length
        fpSize=2048,     # fingerprint size
        useHs=False,     # ignore explicit hydrogens
        branchedPaths=True,  # include branched paths
        useBondOrder=True    # include bond order info
    )
    
    fp1 = gen.GetFingerprint(mol1)
    fp2 = gen.GetFingerprint(mol2)

    dissimilarity = 1 - DataStructs.TanimotoSimilarity(fp1, fp2)
    
    return dissimilarity
    