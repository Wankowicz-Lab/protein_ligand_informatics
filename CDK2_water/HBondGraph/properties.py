"""
Get and set graph properties.
"""

from Bio import PDB
import networkx as nx

# Default hydrogen bond donor and acceptor definitions from Mingbin
DEFAULT_DONOR_ATOMS = ['N', 'ND1', 'ND2', 'NE', 'NE1', 'NE2', 'NH1', 'NH2', 'NZ', 'OG', 'OG1', 'OH', 'SG']
DEFAULT_ACCEPTOR_ATOMS = ['ND1', 'NE2', 'O', 'OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1', 'OH', 'SD', 'SG']

def get_water_nodes(Hbond_graph):
    """ 
    Takes an H-bond graph.
    Returns a list of node names for water residues.
    """
    residue_dict = nx.get_node_attributes(Hbond_graph, 'res')
    water_nodes = [n for n in Hbond_graph.nodes if residue_dict[n] == 'water']
    return water_nodes
      
      
def get_protein_nodes(Hbond_graph):
    """ 
    Takes an H-bond graph.
    Returns a list of node names for water residues.
    """
    residue_dict = nx.get_node_attributes(Hbond_graph, 'res')
    protein_nodes = [n for n in Hbond_graph.nodes if residue_dict[n] != 'water']
    return protein_nodes
      

def get_attr_from_nodes(G, nodes, attr):
    """
    Gets attr from nodes of graph G.
    Nodes can be 'water', 'protein', 'all', or a list of node names
    """
    if nodes == 'all':
        nodes = G.nodes
    elif nodes == 'water':
        nodes = get_water_nodes(G)
    elif nodes == 'protein':
        nodes = get_protein_nodes(G)
        
    return [G.nodes[node][attr] for node in nodes]


def get_node_from_atom(atom):
    name = atom.name
    parent = atom.get_parent()
    hetflag, resseq, _ = parent.get_id()
    if hetflag != " ":
        ch = 'w'
    elif name in DEFAULT_DONOR_ATOMS:
        ch = 'd'
    elif name in DEFAULT_ACCEPTOR_ATOMS:
        ch = 'a'
    else:
        return None
    
    return (ch, resseq, name)


def get_most_conserved_waters(Hbond_graph):
    waters = get_water_nodes(Hbond_graph)
    freqs = get_attr_from_nodes(Hbond_graph, waters, 'freq')
    max_freq = max(freqs)
    if max_freq == 0:
        raise AttributeError("Nodes have not been mapped yet. Please run map_nodes.")
    
    return [w for w,f in zip(waters, freqs) if f==max_freq]
    
    
def edit_water_b_factor(input_file, output_file, b_factor_dict):
    """
    Edit PDB b factors to view in PyMol. Replace with node "freq" to visualize water conservation.
    input_file and output_file are string names of the PDB files.
    coord_dict is a dict of (node_id:[x,y,z]) pairs
    b_factor_dict is a dict of (node_id:new_value) pairs
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("test", input_file)
    
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    n = get_node_from_atom(atom)
                    atom.set_bfactor( b_factor_dict.get(n, 0) ) 

    # Save the modified structure
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_file)