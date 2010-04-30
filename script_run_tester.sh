<<<<<<< .mine
#!/bin/bash

EMAIL=wko@seas.upenn.edu

#qsub -m a -M $EMAIL -cwd -b y -j y -q core2-quad.q -pe threaded 4 ./script_test.sh 50 65536 4 1 serial
#qsub -m a -M $EMAIL -cwd -b y -j y -q core2-quad.q -pe threaded 4 ./script_test.sh 50 65536 4 1 openmp

for (( ii = 1; ii <= 1048576; ii = (ii * 2) ))
do
    echo qsub -m a -M $EMAIL -cwd -b y -j y -q core2-quad.q -pe threaded 4 ./script_test.sh 50 $ii 4 1 openmp
    qsub -m a -M $EMAIL -cwd -b y -j y -q core2-quad.q -pe threaded 4 ./script_test.sh 50 $ii 4 1 openmp
done
