mol load pdb aligned_amber.pdb
set filini STAGE-mmm-ini.pdb
set filpdb STAGE-mmm.pdb
set filpsf STAGE-mmm.psf

set all [atomselect top all]
$all set chain A
set a [atomselect top "resname MMM and noh"]
set tot [$a get name]
set n 0
set z 100
set zmax 100
set xd XDIS 
set yd YDIS 
set ran RANG
set dmax DMAX
set dmin DMIN
set dini ZDIS
set dzmx ZMAX
set mat {}

foreach i $tot {
set t [atomselect top "resname MMM and name $i"]
set p [atomselect top "protein and resid P1A and name NN"]
set d1 [measure center $t weight mass]
set d2 [measure center $p weight mass]
set diff [vecsub $d1 $d2]
foreach {x1 y1 z1} $diff {break}
set xl $x1
set yl $y1
set zl $z1
if {[expr $zl > $dini] && [expr $zl < $dzmx] && [expr abs([expr $x1 - $xd]) < $ran] && [expr abs([expr $y1 - $yd]) < $ran ]} {
lappend mat $i
}
}

foreach i $mat {
set t [atomselect top "resname MMM and name $i"]
set m [measure center $t weight mass]
foreach {x1 y1 z1} $m {break}
if [expr $z1 < $z] {
set z $z1
set aa1 $i }
}


set exist [info exists aa1]
if {[expr $exist == 0]} {
set data ""
set filename "anchors.txt"
set fileId [open $filename "w"]
puts -nonewline $fileId $data
close $fileId
puts "Ligand first anchor not found"
exit
}

set m [atomselect top "resname MMM and name $aa1"]
set p [atomselect top "protein resid P1A and name NN"]
set d1 [measure center $m weight mass]
set d2 [measure center $p weight mass]
set diff [vecsub $d1 $d2]
foreach {x1 y1 z1} $diff {break}
set zl $z1
puts "anchor 1 is" 
puts $aa1
puts $zl



set pr [atomselect top "(protein and resid FIRST to LAST and not water and not resname MMM and noh) or (resname MMM)"]
$pr moveby [vecinvert [measure center $pr]] 
$pr writepdb $filini
mol delete all
mol load pdb $filini
set all [atomselect top all]
$all writepdb $filpdb
mol delete all
mol load pdb $filpdb
set a [atomselect top "resname MMM and name $aa1"]
$a moveby {0 0 -5}
$a writepdb dum1.pdb
mol delete all
mol load pdb dum1.pdb
mol load pdb $filpdb
set a [atomselect 3 all]
set b [atomselect 4 "protein and resid P1A and name NN"]
set c [atomselect 4 "resname MMM and name $aa1"]
set d1 [vecsub [measure center $b weight mass] [measure center $c weight mass]]
set d2 [lindex $d1 2]
set dlis [list 0 0 [expr -$d2-5]]
$b moveby $dlis
$b writepdb dum2.pdb
mol delete all
mol load pdb dum1.pdb
mol load pdb dum2.pdb
mol load pdb $filpdb
set a [atomselect 5 all]
set b [atomselect 6 all]
set c [atomselect 7 all]
set d1 [measure center $a weight mass]
set d2 [measure center $b weight mass]
set diff [vecsub $d1 $d2]
set d [veclength $diff]
foreach {x1 y1 z1} $diff {break}
set x $x1
set y $y1
set z $z1
$c move [transvecinv "$x $y $z"]
$c move [transaxis z -90]
$c moveby [vecinvert [measure center $c weight mass]]
set lig [atomselect 7 "resname MMM"]
set ligh [atomselect 7 "resname MMM and noh"]
$lig set chain S
$lig set resid 1
$lig writepdb mmm.pdb
$ligh writepdb mmm-noh.pdb
$c writepdb $filpdb
mol delete all
mol load pdb $filpdb
set a [atomselect top "resname MMM and name $aa1"]
$a moveby {0 0 -5}
$a set name Pb
$a set resname DUM
$a set resid 1
$a set chain D
$a set beta 0.00
$a writepdb dum1.pdb
mol delete all
mol load pdb dum1.pdb
mol load pdb $filpdb
set a [atomselect 9 all]
set b [atomselect 10 "protein and resid P1A and name NN"]
set c [atomselect 10 "resname MMM and name $aa1"]
set d1 [vecsub [measure center $b weight mass] [measure center $c weight mass]]
set d2 [lindex $d1 2]
set dlis [list 0 0 [expr -$d2-5]]
$b moveby $dlis
$b set name Pb
$b set resname DUM
$b set resid 2
$b set chain D
$b set beta 0.00
$b writepdb dum2.pdb
mol delete all
mol load pdb dum1.pdb
mol load pdb dum2.pdb
mol load pdb $filpdb
set a [atomselect 11 all]
set b [atomselect 12 all]
set d3 [veclength [vecsub [measure center $a weight mass] [measure center $b weight mass]]]
set dlis2 [list 0 0 $d3]
$b moveby $dlis2
$b set name Pb
$b set resname DUM
$b set resid 3
$b set chain D
$b set beta 0.00
$b writepdb dum3.pdb
mol delete all

mol load pdb dum1.pdb
mol load pdb dum2.pdb
mol load pdb dum3.pdb
mol load pdb $filpdb

set amat {}
foreach i $tot {
set alis {}
set angle1 {}
set angle2 {}
set angle3 {}
set angle {}
set t [atomselect 17 "resname MMM and name $i"]
set p [atomselect 17 "resname MMM and name $aa1"]
set d [atomselect 14 all]
if {$i ne $aa1} { set a1 [$d get index]
set d1 [measure center $t weight mass]
set d2 [measure center $p weight mass]
set leng [veclength [vecsub $d1 $d2]]
lappend angle1 $a1
lappend angle1 "14"
lappend angle $angle1
lappend alis [$d get name]
set a2 [$p get index]
lappend angle2 $a2
lappend angle2 "17"
lappend angle $angle2
lappend alis [$p get name]
set a3 [$t get index]
lappend angle3 $a3
lappend angle3 "17"
lappend angle $angle3
lappend alis [$t get name]
set ang [measure angle $angle]
if {[expr $leng > $dmin] && [expr $leng < $dmax]} {
lappend amat $i}
}
}


set amx 90
foreach i $amat {
set angle1 {}
set angle2 {}
set angle3 {}
set angle {}
set t [atomselect 17 "resname MMM and name $i"]
set p [atomselect 17 "resname MMM and name $aa1"]
set d [atomselect 14 all]
set d1 [measure center $t weight mass]
set d2 [measure center $p weight mass]
set a1 [$d get index]
lappend angle1 $a1
lappend angle1 "14"
lappend angle $angle1
lappend alis [$d get name]
set a2 [$p get index]
lappend angle2 $a2
lappend angle2 "17"
lappend angle $angle2
lappend alis [$p get name]
set a3 [$t get index]
lappend angle3 $a3
lappend angle3 "17"
lappend angle $angle3
lappend alis [$t get name]
set ang [measure angle $angle]
if {[expr abs([expr $ang - 90.0])] < $amx} {
set amx [expr abs([expr $ang - 90.0])]
set angl $ang
set aa2 $i
set leng [veclength [vecsub $d1 $d2]]
}
}


set exist [info exists aa2]
if {[expr $exist == 0]} {
set data "$aa1\n"
set filename "anchors.txt"
set fileId [open $filename "w"]
puts -nonewline $fileId $data
close $fileId
puts "Ligand second anchor not found"
exit
}

puts "anchor 2 is" 
puts $aa2
puts $angl
puts $leng

mol delete all             
mol load pdb dum1.pdb
mol load pdb dum2.pdb
mol load pdb dum3.pdb
mol load pdb $filpdb

set amat {}
foreach i $tot {
set alis {}
set angle1 {}
set angle2 {}
set angle3 {}
set angle {}
set t [atomselect 21 "resname MMM and name $i"]
set p [atomselect 21 "resname MMM and name $aa2"]
set d [atomselect 21 "resname MMM and name $aa1"]
if {$i ne $aa1 && $i ne $aa2} { set a1 [$d get index]
set d1 [measure center $t weight mass]
set d2 [measure center $p weight mass]
set leng [veclength [vecsub $d1 $d2]]
lappend angle1 $a1
lappend angle1 "21"
lappend angle $angle1
lappend alis [$d get name]
set a2 [$p get index]
lappend angle2 $a2
lappend angle2 "21"
lappend angle $angle2
lappend alis [$p get name]
set a3 [$t get index]
lappend angle3 $a3
lappend angle3 "21"
lappend angle $angle3
lappend alis [$t get name]
set ang [measure angle $angle]
if {[expr $leng > $dmin] && [expr $leng < $dmax]} {
lappend amat $i}
}
}


set adf 90
foreach i $amat {
set angle1 {}
set angle2 {}
set angle3 {}
set angle {}
set t [atomselect 21 "resname MMM and name $i"]
set p [atomselect 21 "resname MMM and name $aa2"]
set d [atomselect 21 "resname MMM and name $aa1"]
set d1 [measure center $t weight mass]
set d2 [measure center $p weight mass]
set a1 [$d get index]
lappend angle1 $a1
lappend angle1 "21"
lappend angle $angle1
lappend alis [$d get name]
set a2 [$p get index]
lappend angle2 $a2
lappend angle2 "21"
lappend angle $angle2
lappend alis [$p get name]
set a3 [$t get index]
lappend angle3 $a3
lappend angle3 "21"
lappend angle $angle3
lappend alis [$t get name]
set ang [measure angle $angle]
if {[expr abs([expr $ang - 90.0])] < $adf} {
set adf [expr abs([expr $ang - 90.0])]
set angf $ang
set aa3 $i
set leng [veclength [vecsub $d1 $d2]]
}
}

set exist [info exists aa3]
if {[expr $exist == 0]} {
set data "$aa1 $aa2\n"
set filename "anchors.txt"
set fileId [open $filename "w"]
puts -nonewline $fileId $data
close $fileId
puts "Ligand third anchor not found"
exit
}

puts "anchor 3 is"
puts $aa3
puts $angf
puts $leng

puts "The three anchors are"

set data "$aa1 $aa2 $aa3\n"
set filename "anchors.txt"
set fileId [open $filename "w"]
puts -nonewline $fileId $data
close $fileId

exit
