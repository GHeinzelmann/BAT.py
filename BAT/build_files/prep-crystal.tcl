mol load pdb CCCC.pdb
set a [atomselect top "chain A and protein and noh"]
set b [atomselect top "chain A and resname MMM and noh"]
$a writepdb rec_file.pdb
$b writepdb mmm.pdb
exit

