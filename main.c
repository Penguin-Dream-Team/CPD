#include <stdio.h>
#include <math.h>
#include "gen_points.h"

int main(int argc, char *argv[]) { 
    int n_dims;
    long n_samples;
    
    double** points = get_points(argc, argv, &n_dims, &n_samples);

    print_points(points, n_dims, n_samples);

    int a = get_furthest_point(points, 0, &n_dims, &n_samples); 
    int b = get_furthest_point(points, a, &n_dims, &n_samples);

    printf("A %d\nB %d\n", a, b);
}

int get_furthest_point(double** points, int first_point, int* n_dims, long* n_samples) {

    double distances[*n_samples];
    for (int i = 0; i < *n_dims; i++) {
        for (int j = 0; j < *n_samples; j++) {
            distances[j] += (points[i][first_point] - points[i][j]) * (points[i][first_point] - points[i][j]);
        }
    }

    int max = 0;
    double max_distance = 0;
    for (int i = 0; i < *n_samples; i++) {
        double new_distance = sqrt(distances[i]);
        printf("Distance %d %d = %f\n", first_point, i, new_distance);
        if (max_distance < new_distance) {
            max = i;
            max_distance = new_distance;
        }
    }

    return max;
}
