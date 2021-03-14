#
# In order to execute this "Makefile" just type "make"
#    A. Delis (ad@di.uoa.gr)
#

OBJS      = main.o gen_points.o
SOURCE    = main.c gen_points.c
HEADER    = gen_points.h
OUT       = ballAlg
CC        = gcc
FLAGS     = -g3 -c -Wall
LFLAGS    = -lm
# -g option enables debugging mode 
# -c flag generates object code for separate files


all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) $(LFLAGS)


# create/compile the individual files >>separately<<
main.o: main.c
	$(CC) $(FLAGS) main.c -std=gnu99

ballQuery.o: ballQuery.c
	$(CC) $(FLAGS) ballQuery.c -std=gnu99

gen_points.o: gen_points.c
	$(CC) $(FLAGS) gen_points.c -std=gnu99


# clean house
clean:
	rm -f $(OBJS) $(OUT)

# run the program
run: $(OUT)
	./$(OUT)

# compile program with debugging information
debug: $(OUT)
	valgrind $(OUT)

# run program with valgrind for errors
valgrind: $(OUT)
	valgrind $(OUT)

# run program with valgrind for leak checks
valgrind_leakcheck: $(OUT)
	valgrind --leak-check=full ./$(OUT) $(ARGS)

# run program with valgrind for leak checks (extreme)
valgrind_extreme: $(OUT)
	valgrind --leak-check=full --show-leak-kinds=all --leak-resolution=high --track-origins=yes --vgdb=yes $(OUT)