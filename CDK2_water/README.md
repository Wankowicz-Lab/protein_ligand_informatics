This folder contains scripts for downloading, processing, and analyzing structures of CDK2 bound to different ligands.

./cdk2\_expanded contains the full dataset of 63 structures from PDB-redo with high resolution (<2A), specified space group ("P 21 21 21"), and unit cell (52 < a < 54, 70 < b < 72, 71 < c < 73, alpha = 90, beta = 90, gamma = 90) characteristics.

fix\_pdb.py runs openmm pdbfixer to fill in missing residues with energy minimization and uses the command line tool pdb-tools to regularize all PDBs to have the same residues and atoms.


