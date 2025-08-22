This folder contains scripts for downloading, processing, and analyzing structures of CDK2 bound to different ligands.

./cdk2\_expanded contains the full dataset of 63 structures from PDB-redo with high resolution (<2A), specified space group ("P 21 21 21"), and unit cell (52 < a < 54, 70 < b < 72, 71 < c < 73, alpha = 90, beta = 90, gamma = 90) characteristics.

./HBondGraph contains a set of python files with functions to construct graphs from PDB file structures, perform graph analysis, and collect and analyze data about the water structure and ligand properties in those structures.

./analysis_examples.py has routines that demonstrate how to use the functions in ./HBondGraph to produce the figures

fix\_pdb.py runs openmm pdbfixer to fill in missing residues with energy minimization and uses the command line tool pdb-tools to regularize all PDBs to have the same residues and atoms.

./pdbredo_download.sh is a script to grab pdb files following specific constraints from the PDB-redo data bank based on the entries in ./pdb_ids
