

x=0
while [  $x -lt 7 ]; do
cd e0$x
qsub PBS-run
cd ../
cd v0$x
qsub PBS-run
cd ../
let x=x+1
done
