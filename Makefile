OBJS      = ballAlg.o ballQuery.o ballAlg-omp.o gen_points.o quickSelect.o
SOURCE    = ballAlg.c ballAlg-omp.c gen_points.c quickSelect.c
HEADER    = gen_points.h quick_sort.h
CC        = gcc
LFLAGS    = -lm -fopenmp -O3 -Wextra -Wall

all: $(OBJS)
	rm -f $(OBJS)

ballAlg.o: ballAlg.c gen_points.o quickSelect.o
	$(CC) ballAlg.c gen_points.o quickSelect.o -o ballAlg $(LFLAGS)

ballAlg-omp.o: ballAlg-omp.c
	$(CC) ballAlg-omp.c gen_points.o quickSelect.o -o ballAlg-omp $(LFLAGS)

ballQuery.o: ballQuery.c
	$(CC) ballQuery.c -o ballQuery $(LFLAGS)

gen_points.o: gen_points.c
	$(CC) -c $(FLAGS) gen_points.c

quickSelect.o: quickSelect.c
	$(CC) -c $(FLAGS) quickSelect.c

# clean house
clean:
	rm -f $(OBJS) ballAlg ballAlg-omp ballQuery
