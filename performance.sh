#!/bin/bash

for run in {1..10}; do
    #./ballAlg-omp 50 1000000 0 > /dev/null  
    #./ballAlg-omp 4 20000000 0 > /dev/null
    #./ballAlg 20 1000000 0 > /dev/null
    #./ballAlg 3 5000000 0 > /dev/null
    #./ballAlg 4 10000000 0 > /dev/null
    #./ballAlg 3 20000000 0 > /dev/null
    srun -x lab6p[1-9] -n 8 --cpus-per-task 4 ./ballAlg-mpi 2 20000000 0 > /dev/null
    #srun -x lab6p[1-9] -n 4 --cpus-per-task 4 ./ballAlg-mpi 2 8 0 > /dev/null
done