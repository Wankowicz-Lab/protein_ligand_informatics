source /programs/sbgrid.shrc
# APO
Prepare pdb file  
phenix.pdbtools 7kqo.pdb keep="chain 'A'" output.file_name=clean_chainA.pdb  
phenix.pdbtools clean_chainA.pdb remove_alt_confs=True  
pdb4amber -i clean_chainA_modified.pdb -o clean_chainA_final.pdb  

Prepare topology and coordinate file  
tleap -f tleap_protein.in  

Relaxation (change all file names)  
sander -O -i 1min.in -o 1min_2.out -p 7kqo_clean.prmtop -c 7kqo_clean.inpcrd -r 1min_2.rst -ref 7kqo_clean.inpcrd -x 1min_2.nc  
sander -O -i 2mdheat.in -o 2mdheat_3.out -p 7kqo_clean_1.prmtop -c 1min_3.rst -r 2mdheat_3.rst -ref 7kqo_clean_1.inpcrd -x 2mdheat_3.nc  
sander -O -i 3md.in -o 3md_3.out -p 7kqo_clean_1.prmtop -c 2mdheat_3.rst -r 3md_3.rst -x 3md_3.nc  
sander -O -i 4md.in -o 4md_3.out -p 7kqo_clean_1.prmtop -c 3md_3.rst -r 4md_3.rst -ref 7kqo_clean_1.inpcrd -x 4md_3.nc  

Run MD  
pmemd.cuda -O -i prod.in -o prod.out -p 7kqo.prmtop -c 4md.rst -r prod.rst -x prod.nc  

# Ligand-protein
Ligand parameterization  
obabel -ipdb file.pdb -O file.com  

Gaussian input (edit file.com)  
%mem=3000MB  
%nproc=16  
%chk=FAD_geom.chk  
#P B3LYP/6-31+G* opt scf=(tight, xqc) FAD.pdb-2 1  
P     -3.92600    -7.72400    -13.10300  
P     -4.33200    -7.21400    -15.75400  
O     -2.54100    -8.22900    -13.23200  

Run Gaussian  
g16 file.com  

Generate esp.com from log  
obabel -ig16 FAD.log -O esp.com  

ESP Gaussian input (edit esp.com)  
%mem=3000MB  
%nproc=16  
%chk=resp.chk  
#P HF/6-31+G* SCF=Tight Pop=MK IOp(6/33=2,6/41=10,6/42=6) FAD-2 1  
P ............  

g16 esp.com  

AMBER tools  
antechamber -i lig.log -fi gout -o lig.ac -fo ac -at gaff2 -c resp -pf y  
espgen -i lig.log -o lig.esp  
respgen -i lig.ac -o lig.respin1 -f resp1  
respgen -i lig.ac -o lig.respin2 -f resp2  
resp -O -i lig.respin1 -o lig.respout1 -e lig.esp -t qout_stage1  
resp -O -i lig.respin2 -o lig.respout2 -e lig.esp -q qout_stage1 -t qout_stage2  
antechamber -i lig.pdb -fi pdb -o lig.mol2 -fo mol2 -at gaff2 -c rc -cf qout_stage2 -pf y  
parmchk2 -i lig.mol2 -f mol2 -o lig.frcmod  
