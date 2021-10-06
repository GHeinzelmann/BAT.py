mol load pdb reference.pdb
mol load pdb complex.pdb
set sel1 [atomselect 0 "backbone"]
set sel2 [atomselect 1 "backbone"]
set all [atomselect 1 all]
$all move [measure fit $sel2 $sel1]
set all [atomselect 1 all]
$all writepdb aligned.pdb
exit