
ELEMENT_MAP = {
    'H': 1,
    'C': 6,
    'N': 7,
    'O': 8,
    'F': 9,
    'P': 15,
    'S': 16,
    'Cl': 17,
    'Br': 35,
    'I': 53
}

def getElementForTensor(element):
    if element in ELEMENT_MAP:
        return ELEMENT_MAP[element]
    print("[tensormap] Element not found: ", element)
    return 0

RESIDUE_NAME_MAP = {
    'LIG': 1,

    'ALA': 2,
    'ARG': 3,
    'ASN': 4,
    'ASP': 5,
    'CYS': 6,
    'GLN': 7,
    'GLU': 8,
    'GLY': 9,
    'HIS': 10,
    'ILE': 11,
    'LEU': 12,
    'LYS': 13,
    'MET': 14,
    'PHE': 15,
    'PRO': 16,
    'SER': 17,
    'THR': 18,
    'TRP': 19,
    'TYR': 20,
    'VAL': 21
}


def getResidueForTensor(residue): 
    if residue in RESIDUE_NAME_MAP:
        return RESIDUE_NAME_MAP[residue]
    print("[tensormap] Residue not found: ", residue)
    return 0


ATOM_NAME_MAP = {
    'N': 1,
    'CA': 2,
    'C': 3,
    'O': 4,
    'CB': 5,
    'CG': 6,
    'CG1': 7,
    'CG2': 8,
    'CD': 9,
    'CD1': 10,
    'CD2': 11,
    'CE': 12,
    'CE1': 13,
    'CE2': 14,
    'CE3': 15,
    'CZ': 16,
    'CZ2': 17,
    'CZ3': 18,
    'CH2': 19,
    'ND1': 20,
    'ND2': 21,
    'NE': 22,
    'NE1': 23,
    'NE2': 24,
    'NZ': 25,
    'NH1': 26,
    'NH2': 27,
    'SD': 28,
    'SG': 29,
    'OG': 30,
    'OG1': 31,
    'OD1': 32,
    'OD2': 33,
    'OE1': 34,
    'OE2': 35,

    'OH': 36,
    'OXT': 37,
    'C01': 38,
    'C02': 39,
    'C03': 40,
    'C04': 41,
    'C05': 42,
    'C06': 43,
    'C07': 44,
    'C08': 45,
    'C09': 46,
    'C10': 47,
    'C11': 48,
    'C12': 49,
    'C13': 50,
    'C14': 51,
    'C15': 52,
    'C16': 53,
    'C17': 54,
    'C18': 55,
    'C19': 56,
    'C20': 57,
    'C21': 58,
    'C22': 59,
    'C23': 60,
    
}


def getAtomForTensor(atom):
    if atom in ATOM_NAME_MAP:
        return ATOM_NAME_MAP[atom]
    print("[tensormap] Atom not found: ", atom)
    return 0


CHAIN_MAP = {
    'A': 1,
    'B': 2,
    'C': 3,
    'D': 4,
    'E': 5,
    'F': 6,
    'G': 7,
    'H': 8,
    'I': 9,
    'J': 10,
    'K': 11,
    'L': 12,
    'M': 13,
    'N': 14,
    'O': 15,
    'P': 16,
    'Q': 17,
    'R': 18,
    'S': 19,
    'T': 20,
    'U': 21,
    'V': 22,
    'W': 23,
    'X': 24,
    'Y': 25,
    'Z': 26
}

def getChainForTensor(chain):
    if chain in CHAIN_MAP:
        return CHAIN_MAP[chain]
    print("[tensormap] Chain not found: ", chain)
    return 0