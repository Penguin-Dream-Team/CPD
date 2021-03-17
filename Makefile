#
# In order to execute this "Makefile" just type "make"
#    A. Delis (ad@di.uoa.gr)
#

OBJS      = main.o gen_points.o quick_sort.o
SOURCE    = main.c gen_points.c quick_sort.c
HEADER    = gen_points.h quick_sort.h
OUT       = ballAlg
CC        = gcc
FLAGS     = -g3 -c -Wextra -Wall -O3
LFLAGS    = -lm -fopenmp -O3 -Wextra -Wall
# -g option enables debugging mode 
# -c flag generates object code for separate files


all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) $(LFLAGS)
	rm -f $(OBJS)


# create/compile the individual files >>separately<<
main.o: main.c
	$(CC) $(FLAGS) main.c

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
