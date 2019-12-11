

x=0
while [  $x -lt 10 ]; do
cd f0$x
qsub PBS-run
cd ../
cd w0$x
qsub PBS-run
cd ../
let x=x+1
done

if [ $x -ge 10 ]; then
while [  $x -lt 12 ]; do
cd f$x
qsub PBS-run
cd ../
cd w$x
qsub PBS-run
cd ../
let x=x+1
done
fi
