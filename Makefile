OBJS      = ballAlg.o ballQuery.o ballAlg-omp.o ballAlg-mpi.o gen_points.o quickSelect.o ballAlg-large-mpi.o
SOURCE    = ballAlg.c ballAlg-omp.c ballAlg-mpi.c gen_points.c quickSelect.c ballAlg-large-mpi.c
HEADER    = gen_points.h quick_sort.h ballAlg-large-mpi.h
CC        = mpicc
LFLAGS    = -lm -fopenmp -O3 -Wextra -Wall

all: $(OBJS)
	rm -f $(OBJS)

ballAlg.o: ballAlg.c gen_points.o quickSelect.o
	$(CC) ballAlg.c gen_points.o quickSelect.o -o ballAlg $(LFLAGS)

ballAlg-omp.o: ballAlg-omp.c
	$(CC) ballAlg-omp.c gen_points.o quickSelect.o -o ballAlg-omp $(LFLAGS)

ballAlg-mpi.o: ballAlg-mpi.c
	$(CC) ballAlg-mpi.c gen_points.o quickSelect.o ballAlg-large-mpi.o -o ballAlg-mpi $(LFLAGS)

ballQuery.o: ballQuery.c
	$(CC) ballQuery.c -o ballQuery $(LFLAGS)

gen_points.o: gen_points.c
	$(CC) -c $(FLAGS) gen_points.c

quickSelect.o: quickSelect.c
	$(CC) -c $(FLAGS) quickSelect.c

ballAlg-large-mpi.o: ballAlg-large-mpi.c
	$(CC) -c $(FLAGS) ballAlg-large-mpi.c

# clean house
clean:
	rm -f $(OBJS) ballAlg ballAlg-omp ballAlg-mpi ballQuery
