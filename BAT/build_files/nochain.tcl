mol load pdb reference_amber.pdb
mol load pdb complex.pdb
set a [atomselect 0 "protein"]
set b [atomselect 1 "protein"]
$a set chain X
$b set chain X
$a writepdb reference_amber-nc.pdb
$b writepdb complex-nc.pdb
exit
