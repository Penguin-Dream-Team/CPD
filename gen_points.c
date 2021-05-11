#include <stdio.h>
#include <stdlib.h>

#define RANGE 10

extern void print_point(double *, int);

double **create_array_pts(int n_dims, long np)
{
    double *_p_arr;
    double **p_arr;

    _p_arr = (double *) malloc(n_dims * np * sizeof(double));
    p_arr = (double **) malloc(np * sizeof(double *));
    if((_p_arr == NULL) || (p_arr == NULL)){
        printf("Error allocating array of points, exiting.\n");
        exit(4);
    }

    for(long i = 0; i < np; i++)
        p_arr[i] = &_p_arr[i * n_dims];

    return p_arr;
}


double **get_points(int argc, char *argv[], int n_dims, long np)
{
    double **pt_arr;
    unsigned seed;
    long i;
    int j;

    seed = atoi(argv[3]);
    srandom(seed);

    pt_arr = (double **) create_array_pts(n_dims, np);

    for(i = 0; i < np; i++)
        for(j = 0; j < n_dims; j++)
            pt_arr[i][j] = RANGE * ((double) random()) / RAND_MAX;

    #ifdef DEBUG
    for(i = 0; i < *np; i++)
        print_point(pt_arr[i], *n_dims);
    #endif

    return pt_arr;
}

double **get_points_large(int argc, char *argv[], int n_dims, long np, long start, long end)
{
    double **pt_arr;
    unsigned seed;
    long i, rand;
    int j;

    seed = atoi(argv[3]);
    srandom(seed);

    pt_arr = (double **) create_array_pts(n_dims, np);

    for(i = 0; i < np; i++)
        for(j = 0; j < n_dims; j++) {
            rand = random();
            if (start <= i && i < end) {
                pt_arr[i - start][j] = RANGE * ((double) rand) / RAND_MAX;
            }
        }

    #ifdef DEBUG
    for(i = 0; i < *np; i++)
        print_point(pt_arr[i], *n_dims);
    #endif

    return pt_arr;
}
