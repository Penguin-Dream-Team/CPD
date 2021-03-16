#include <stdio.h>
#include <stdlib.h>
#include "gen_points.h"

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


double **get_points(int argc, char *argv[], int *n_dims, long *np)
{
    double **pt_arr;
    unsigned seed;
    long i;
    int j;

    if(argc != 4){
        printf("Usage: %s <n_dims> <n_points> <seed>\n", argv[0]);
        exit(1);
    }

    *n_dims = atoi(argv[1]);
    if(*n_dims < 2){
        printf("Illegal number of dimensions (%d), must be above 1.\n", *n_dims);
        exit(2);
    }

    *np = atol(argv[2]);
    if(*np < 1){
        printf("Illegal number of points (%ld), must be above 0.\n", *np);
        exit(3);
    }

    seed = atoi(argv[3]);
    srandom(seed);

    pt_arr = (double **) create_array_pts(*n_dims, *np);

    for(i = 0; i < *np; i++)
        for(j = 0; j < *n_dims; j++)
            pt_arr[i][j] = RANGE * ((double) random()) / RAND_MAX;

#ifdef DEBUG
    for(i = 0; i < *np; i++)
        print_point(pt_arr[i], *n_dims);
#endif

    return pt_arr;
}

void print_point(double* point, int n_dims) {
    print_point_indent(point, n_dims, 1);
}

void print_point_indent(double* point, int n_dims, int indent) {
    for(int i = 0; i < n_dims; i++) {
        for(int t = 0; t < indent; t++) {
            printf("\t");
        }
        printf("Dimension: %d = %f\n", i, point[i]);
    }
}

void print_points(double** points, int n_dims, long n_samples) {
    for(int i = 0; i < n_samples; i++) {
        printf("Point %d\n", i);
        print_point(points[i], n_dims);
    }
}
