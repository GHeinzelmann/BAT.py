mol load pdb aligned_amber.pdb
set filini STAGE-mmm-ini.pdb
set filpdb STAGE-mmm.pdb
set n 0
set rad 100
set xd XDIS 
set yd YDIS 
set zd ZDIS
set ran RANG
set dmax DMAX
set dmin DMIN
set mat {}
set sdr_dist SDRD

set pr [atomselect 0 "(not resname MMM) and (resid FIRST to LAST and name CA C N O)"]
set all [atomselect 0 "(resid FIRST to LAST and not water and not resname MMM and noh) or (resname MMM) or (resname OTHRS WAT)"]
$all moveby [vecinvert [measure center $pr weight mass]]
$all writepdb $filini
mol delete all
mol load pdb $filini
set all [atomselect 1 all]
$all writepdb $filpdb

set lig [atomselect 1 "resname MMM"]
set ligh [atomselect 1 "resname MMM and noh"]
$lig set chain S
$lig set resid 1
$lig writepdb mmm.pdb
$ligh writepdb mmm-noh.pdb

set a [atomselect 1 "resname MMM and noh"]
set tot [$a get name]

foreach i $tot {
set t [atomselect 1 "resname MMM and name $i"]
set p [atomselect 1 "(not resname MMM) and (resid P1A and name NN)"]
set d1 [measure center $t weight mass]
set d2 [measure center $p weight mass]
foreach {x2 y2 z2} $d2 {break}
set xl [expr $x2+$xd]
set yl [expr $y2+$yd]
set zl [expr $z2+$zd]
foreach {x1 y1 z1} $d1 {break}
set xc [expr abs($x1-$xl)]
set yc [expr abs($y1-$yl)]
set zc [expr abs($z1-$zl)]
set dist [expr sqrt([expr pow($xc,2) + pow($yc,2) + pow($zc,2)])]
puts $i
puts $dist
if {[expr $dist < $ran ]} {
lappend mat $i
}
}

puts $mat

foreach i $mat {
set t [atomselect 1 "resname MMM and name $i"]
set p [atomselect 1 "(not resname MMM) and (resid P1A and name NN)"]
set d1 [measure center $t weight mass]
set d2 [measure center $p weight mass]
foreach {x2 y2 z2} $d2 {break}
set xl [expr $x2+$xd]
set yl [expr $y2+$yd]
set zl [expr $z2+$zd]
foreach {x1 y1 z1} $d1 {break}
set xc [expr abs($x1-$xl)]
set yc [expr abs($y1-$yl)]
set zc [expr abs($z1-$zl)]
set diff [expr sqrt([expr pow($xc,2) + pow($yc,2) + pow($zc,2)])]
if [expr $diff < $rad] {
set rad $diff
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

puts "anchor 1 is" 
puts $aa1
puts $rad
puts ""

set amat {}
foreach i $tot {
set alis {}
set angle1 {}
set angle2 {}
set angle3 {}
set angle {}
set t [atomselect 1 "resname MMM and name $i"]
set p [atomselect 1 "resname MMM and name $aa1"]
set d [atomselect 1 "(not resname MMM) and (resid P1A and name NN)"]
if {$i ne $aa1} { set a1 [$d get index]
set d1 [measure center $t weight mass]
set d2 [measure center $p weight mass]
set leng [veclength [vecsub $d1 $d2]]
lappend angle1 $a1
lappend angle1 "1"
lappend angle $angle1
lappend alis [$d get name]
set a2 [$p get index]
lappend angle2 $a2
lappend angle2 "1"
lappend angle $angle2
lappend alis [$p get name]
set a3 [$t get index]
lappend angle3 $a3
lappend angle3 "1"
lappend angle $angle3
lappend alis [$t get name]
set ang [measure angle $angle]
puts $i
puts $ang
puts $leng
if {[expr $leng > $dmin] && [expr $leng < $dmax]} {
lappend amat $i}
}
}

puts $amat

set amx 90
foreach i $amat {
set angle1 {}
set angle2 {}
set angle3 {}
set angle {}
set t [atomselect 1 "resname MMM and name $i"]
set p [atomselect 1 "resname MMM and name $aa1"]
set d [atomselect 1 "(not resname MMM) and (resid P1A and name NN)"]
set d1 [measure center $t weight mass]
set d2 [measure center $p weight mass]
set a1 [$d get index]
lappend angle1 $a1
lappend angle1 "1"
lappend angle $angle1
lappend alis [$d get name]
set a2 [$p get index]
lappend angle2 $a2
lappend angle2 "1"
lappend angle $angle2
lappend alis [$p get name]
set a3 [$t get index]
lappend angle3 $a3
lappend angle3 "1"
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
puts ""

set amat {}
foreach i $tot {
set alis {}
set angle1 {}
set angle2 {}
set angle3 {}
set angle {}
set t [atomselect 1 "resname MMM and name $i"]
set p [atomselect 1 "resname MMM and name $aa2"]
set d [atomselect 1 "resname MMM and name $aa1"]
if {$i ne $aa1 && $i ne $aa2} { set a1 [$d get index]
set d1 [measure center $t weight mass]
set d2 [measure center $p weight mass]
set leng [veclength [vecsub $d1 $d2]]
lappend angle1 $a1
lappend angle1 "1"
lappend angle $angle1
lappend alis [$d get name]
set a2 [$p get index]
lappend angle2 $a2
lappend angle2 "1"
lappend angle $angle2
lappend alis [$p get name]
set a3 [$t get index]
lappend angle3 $a3
lappend angle3 "1"
lappend angle $angle3
lappend alis [$t get name]
set ang [measure angle $angle]
puts $i
puts $ang
puts $leng
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
set t [atomselect 1 "resname MMM and name $i"]
set p [atomselect 1 "resname MMM and name $aa2"]
set d [atomselect 1 "resname MMM and name $aa1"]
set d1 [measure center $t weight mass]
set d2 [measure center $p weight mass]
set a1 [$d get index]
lappend angle1 $a1
lappend angle1 "1"
lappend angle $angle1
lappend alis [$d get name]
set a2 [$p get index]
lappend angle2 $a2
lappend angle2 "1"
lappend angle $angle2
lappend alis [$p get name]
set a3 [$t get index]
lappend angle3 $a3
lappend angle3 "1"
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

mol load pdb dum.pdb

set a [atomselect 1 "(not resname MMM) and (resid FIRST to LAST and name CA C N O)"]
set b [atomselect 1 "resname MMM and noh"]
set c [atomselect 2 all]
$c moveby [vecsub [measure center $a weight mass] [measure center $c weight mass]]
$c writepdb dum1.pdb
if {[expr $sdr_dist != 0]} {
set dlis [list 0 0 [expr $sdr_dist]]
$b moveby $dlis
$c moveby [vecsub [measure center $b weight mass] [measure center $c weight mass]]
$c set resid 2
$c writepdb dum2.pdb
set dlis2 [list 0 0 [expr -1*$sdr_dist]]
$b moveby $dlis2
}

exit
