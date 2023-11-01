mol load pdb rec_file.pdb
set dum [atomselect 0 "resname DUM"]
set prot [atomselect 0 "protein"]
set othrs [atomselect 0 "resname OTHRS and not water and same residue as within SHLL of (protein or resname MMM)"]
set wat [atomselect 0 "water and same residue as within SHLL of (protein or resname MMM)"]
set lig [atomselect 0 "resname MMM"]
$wat set resname WAT 
$dum writepdb dummy.pdb
$prot writepdb protein.pdb
$othrs writepdb others.pdb
$wat writepdb crystalwat.pdb
$lig writepdb mmm.pdb
exit
