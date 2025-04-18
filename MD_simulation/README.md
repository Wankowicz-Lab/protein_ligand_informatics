MD Workflow

source /programs/sbgrid.shrc
phenix.pdbtools 7kqo.pdb keep="chain 'A'" output.file_name=clean_chainA.pdb
phenix.pdbtools clean_chainA.pdb remove_alt_confs=True
pdb4amber -i clean_chainA_modified.pdb -o clean_chainA_final.pdb
tleap -f tleap_protein.in

Add MD_simulation folder with setup commands
