#!/bin/bash

for run in {1..10}; do
    #./ballAlg-omp 50 1000000 0 > /dev/null  
    #./ballAlg-omp 4 20000000 0 > /dev/null
    #srun -n 1 ./ballAlg 20 1000000 0 > /dev/null
    #srun -n 1./ballAlg 3 5000000 0 > /dev/null
    #srun -n 1./ballAlg 4 10000000 0 > /dev/null
    #srun -n 1./ballAlg 3 20000000 0 > /dev/null
    #srun -n 1 --cpus-per-task 1 ./ballAlg-mpi 20 1000000 0 > /dev/null
    #srun -n 1 --cpus-per-task 1 ./ballAlg-mpi 3 5000000 0 > /dev/null
    #srun -n 1 --cpus-per-task 1 ./ballAlg-mpi 4 10000000 0 > /dev/null
    #srun -n 1 --cpus-per-task 1 ./ballAlg-mpi 3 20000000 0 > /dev/null
    #srun -n 4 --cpus-per-task 4 ./ballAlg-mpi 2 8 0 > /dev/null
done
