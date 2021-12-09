


x=0
while [  $x -lt 10 ]; do
cd m0$x
qsub PBS-run
cd ../
cd n0$x
qsub PBS-run
cd ../
let x=x+1
done
