#
# In order to execute this "Makefile" just type "make"
#    A. Delis (ad@di.uoa.gr)
#

OBJS      = ballAlg-omp.o gen_points.o quickSelect.o
SOURCE    = main.c ballAlg-omp.c gen_points.c quickSelect.c
HEADER    = gen_points.h quick_sort.h
OUT       = ballAlg
CC        = gcc
FLAGS     = -c -Wextra -Wall -g3 -O3 -fopenmp
LFLAGS    = -lm -fopenmp -g3 -O3 -Wextra -Wall
# -g option enables debugging mode 
# -c flag generates object code for separate files


all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) $(LFLAGS)
	rm -f $(OBJS)


main.o: main.c gen_points.o quickSelect.o
	$(CC) -g main.c gen_points.o quickSelect.o -o $(OUT) $(LFLAGS)
	rm -f $(OBJS)

# create/compile the individual files >>separately<<
ballAlg-omp.o: ballAlg-omp.c
	$(CC) $(FLAGS) ballAlg-omp.c

ballQuery.o: ballQuery.c
	$(CC) $(FLAGS) ballQuery.c

gen_points.o: gen_points.c
	$(CC) $(FLAGS) gen_points.c


# clean house
clean:
	rm -f $(OBJS) $(OUT)

# run the program
run: $(OUT)
	./$(OUT) $(ARGS)

# compile program with debugging information
debug: $(OUT)
	gdb ./$(OUT) $(ARGS)

# run program with valgrind for errors
valgrind: $(OUT)
	valgrind ./$(OUT) $(ARGS)

# run program with valgrind for leak checks
valgrind_leakcheck: $(OUT)
	valgrind --leak-check=full ./$(OUT) $(ARGS)

# run program with valgrind for leak checks (extreme)
valgrind_extreme: $(OUT)
	valgrind --leak-check=full --show-leak-kinds=all --leak-resolution=high --track-origins=yes --vgdb=yes ./$(OUT) $(ARGS)
