


x=0
while [  $x -lt 10 ]; do
cd t0$x
qsub PBS-run
cd ../
cd l0$x
qsub PBS-run
cd ../
cd a0$x
qsub PBS-run
cd ../
cd r0$x
qsub PBS-run
cd ../
let x=x+1
done

if [ $x -ge 10 ]; then
while [  $x -lt 16 ]; do
cd t$x
qsub PBS-run
cd ../
cd l$x
qsub PBS-run
cd ../
cd a$x
qsub PBS-run
cd ../
cd r$x
qsub PBS-run
cd ../
let x=x+1
done
fi


x=0
while [  $x -lt 10 ]; do
cd c0$x
qsub PBS-run
cd ../
let x=x+1
done

if [ $x -ge 10 ]; then
while [  $x -lt 16 ]; do
cd c$x
qsub PBS-run
cd ../
let x=x+1
done
fi

