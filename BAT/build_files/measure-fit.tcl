mol load pdb aligned-nc.pdb
mol load pdb complex.pdb
set sel1 [atomselect 0 "protein and backbone"]
set sel2 [atomselect 1 "protein and backbone"]
set all [atomselect 1 all]
$all move [measure fit $sel2 $sel1]
$all writepdb aligned.pdb
exit
