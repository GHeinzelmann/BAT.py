#!/bin/bash

export MGL_ROOT=/home/germano/CELPP-challenge/mgltools_x86_64Linux2_1.5.6

# Remove all non-receptor information from the PDB file
grep ATOM 5uf0.pdb > stripped_protein.pdb

# Generate charges and other small tweaks for receptor preparation
chimera --nogui --script "chimeraPrep.py stripped_protein.pdb prepared_protein.mol2"

# And use an AutoDockTool to convert the mol2 to pdbqt format
$MGL_ROOT/bin/pythonsh $MGL_ROOT/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r prepared_protein.mol2

# Use Chimera’s DockPrep to calculate charges
babel -i pdb ligand.pdb -o mol2 ligand.mol2 -h
chimera --nogui --script "chimeraPrep.py ligand.mol2 prepared_ligand.mol2" 

# And use an AutoDockTool to convert the mol2 to pdbqt format
$MGL_ROOT/bin/pythonsh $MGL_ROOT/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l prepared_ligand.mol2

## Perform vina docking
vina --receptor prepared_protein.pdbqt --ligand prepared_ligand.pdbqt --center_x \-14.0 --center_y 12.0 --center_z 1.5 --size_x 10 --size_y 10 --size_z 10 --seed -26579 

## Separate poses from the 9 output structures
awk '/MODEL /{close("file"f);f++}{if (f > 9) next}{print $0 > "score"f".pdbqt"}' prepared_ligand_out.pdbqt
sed -e '/ENDMDL/,$d' prepared_ligand_out.pdbqt > top_pose.pdbqt
echo ENDMDL >> top_pose.pdbqt

## Convert ligand poses to pdb
for i in {1..9}
do
j=$(($i-1))
$MGL_ROOT/bin/pythonsh $MGL_ROOT/MGLToolsPckgs/AutoDockTools/Utilities24/pdbqt_to_pdb.py -f score$i\.pdbqt -o pose$j\.pdb
done
$MGL_ROOT/bin/pythonsh $MGL_ROOT/MGLToolsPckgs/AutoDockTools/Utilities24/pdbqt_to_pdb.py -f top_pose.pdbqt -o top_pose.pdb

## Convert the protein to pdb (technically we didn’t allow for protein flexibility here so we could 
## just hand back stripped_protein.pdb, but I’ll write this for the more general case)
$MGL_ROOT/bin/pythonsh $MGL_ROOT/MGLToolsPckgs/AutoDockTools/Utilities24/pdbqt_to_pdb.py -f prepared_protein.pdbqt -o 5uf0_docked.pdb

