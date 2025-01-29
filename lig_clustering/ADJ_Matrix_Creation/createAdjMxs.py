# Create ./built_adj_mxs/*
# Uses alt locs

# Create adjacency matrices for all ligands

import torch
import glob
import ADJ_Matrix_Creation.keywordmap as tm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import networkx as nx


ROOT_DIR = "./ligands"

def get_files():
    return glob.glob(ROOT_DIR + "/*_lig.pdb")

def parse_pdb_line(pdbline):
    atom = pdbline[0:6].strip()
    atom_num = pdbline[6:11].strip()
    atom_name = pdbline[12:16].strip()
    residue_name = pdbline[17:20].strip()
    chain = pdbline[21].strip()
    residue_num = pdbline[22:26].strip()
    x = pdbline[30:38].strip()
    y = pdbline[38:46].strip()
    z = pdbline[46:54].strip()
    occupancy = pdbline[54:60].strip()
    b_factor = pdbline[60:66].strip()
    atom_type = pdbline[77:78].strip()
    return {
        "atom": atom,
        "atom_num": atom_num,
        "atom_name": atom_name,
        "atom_element": atom_type,
        "residue_name": residue_name,
        "chain": chain,
        "residue_num": residue_num,
        "x": x,
        "y": y,
        "z": z,
        "occupancy": occupancy,
        "b_factor": b_factor
    }

def read_pdb(file):
    dataset_name = file.split("/")[-1].split("_")[0].split('\\')[-1]
    with open(file, "r") as f:
        lines = f.readlines()
        protein_atoms = []
        ligand_atoms = []
        for line in lines:
            if line.startswith("ATOM"):
                decodedLine = parse_pdb_line(line)
                protein_atoms.append(decodedLine)
            elif line.startswith("HETATM"):
                decodedLine = parse_pdb_line(line)
                if 1 == 1 or decodedLine["residue_name"] == "LIG": # all of these are lig and not LIG labeled
                    ligand_atoms.append(decodedLine)
        return dataset_name, protein_atoms, ligand_atoms



# =========== INTERACTIONS ============

def is_hydrogen_bond(atomElement1, atomElement2, distance):
    bondable = ["H", "O", "N", "F"]
    if atomElement1 in bondable and atomElement2 in bondable:
        if distance <= 3.2 and distance >= 2.5:
            return True
    return False

def is_covalent_bond(atomElement1, atomElement2, distance):
    if distance <= 1.8: # TODO use table for covalent bonds
        return True

def is_pi_pi_interaction(atom1, atom2, distance):
    aromatic_atoms = {'C', 'N'}
    return (atom1["atom_element"] in aromatic_atoms and atom2["atom_element"] in aromatic_atoms and
            distance < 5.0 and distance > 3.0)



def create_graph(protein_atoms, ligand_atoms):
    atoms = protein_atoms + ligand_atoms
    num_atoms = len(atoms)
    nodes = torch.tensor([[float(atom["x"]), float(atom["y"]), float(atom["z"]), 
                           tm.getElementForTensor(atom["atom_element"]), tm.getAtomForTensor(atom["atom_name"]), 
                           tm.getResidueForTensor(atom["residue_name"]), int(atom["residue_num"]), 
                           tm.getChainForTensor(atom["chain"]), 
                           float(atom["occupancy"]), float(atom["b_factor"])] for atom in atoms])
    edges = []
    edge_features = []
    print("Building edges for", len(atoms), "atoms")
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            distance = ((nodes[i][:3] - nodes[j][:3]) ** 2).sum().sqrt().item()
            if distance < 5.0:
                atom1 = atoms[i]
                atom2 = atoms[j]


                bond_order = 1 # figure out bond order later
                covalent_bond = is_covalent_bond(atom1["atom_element"], atom2["atom_element"], distance)
                hydrogen_bond = is_hydrogen_bond(atom1["atom_element"], atom2["atom_element"], distance)
                pi_pi_interaction = is_pi_pi_interaction(atom1, atom2, distance)

                weight = (float(atom1['occupancy']) + float(atom2['occupancy'])) / 2

                edge_feature = [distance, bond_order, 
                                covalent_bond and 1 or 0, 
                                hydrogen_bond and 1 or 0,
                                pi_pi_interaction and 1 or 0,
                                weight]

                edges.append([i, j])
                edges.append([j, i])
                edge_features.append(edge_feature)
                edge_features.append(edge_feature)
    edges = torch.tensor(edges).t().contiguous()
    edge_features = torch.tensor(edge_features)
    return edges, nodes, edge_features


def make_nx_graph(nodes, edges, edge_features):
 
    if edges.dim() == 1:
        num_edges = edges.numel() // 2
        edges = edges.view(num_edges, 2)
    elif edges.shape[0] == 2:
        edges = edges.t()
    else:
        edges = edges

    G = nx.Graph()
    num_nodes = nodes.size(0)
    if num_nodes == 0:
        G.add_node(0)
    else:
        G.add_nodes_from(range(num_nodes))
    
    edge_list = edges.tolist()
    G.add_edges_from(edge_list)
    
    node_features = nodes.numpy()
    for i in range(num_nodes):
        G.nodes[i]['feature'] = node_features[i]
    
    edge_features = edge_features.numpy()
    for idx, (u, v) in enumerate(edge_list):
        G.edges[u, v]['feature'] = edge_features[idx]
    return G


def save_adj_mx_to_json(adj, dsName):
    adj = adj.todense()
    adj = adj.tolist()
    with open("./built_adj_mxs/lig_adj_" + dsName + ".json", "w") as f:
        f.write(str(adj))

if __name__ == "__main__":
    allLigFiles = get_files()
    for testFile in allLigFiles:
        dataset_name, protein_atoms, ligand_atoms = read_pdb(testFile)
        if (len(protein_atoms) > 100):
            print("More than 100 found in", testFile, ", skipping for now")
            continue
        print("protein atom count:", len(protein_atoms))
        print("ligand atom count:", len(ligand_atoms))
        edges, nodes, edge_features = create_graph(protein_atoms, ligand_atoms)
        G = make_nx_graph(nodes, edges, edge_features)
        adjacency_matrix = nx.adjacency_matrix(G)

        save_adj_mx_to_json(adjacency_matrix, dataset_name)
        print("Saved adj mx for", dataset_name)
