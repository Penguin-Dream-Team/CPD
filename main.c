#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gen_points.h"

int get_furthest_point(double** points, int first_point, int* n_dims, long* n_samples) {
    double distances[*n_samples];

    for (int i = 0; i < *n_samples; i++) {
        distances[i] = 0;
    }

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

void calc_ortho_projection(double **points, int n_dims, int a, int b, int p,
        int out_points_index, double **out_points) {

    double top_inner_product = 0;
    double bot_inner_product = 0;
    for (int i = 0; i < n_dims; i++) {
        double *value_a = &points[i][a];
        double *value_b = &points[i][b];
        double *value_p = &points[i][p];
        top_inner_product += (*value_p - *value_a) * (*value_b - *value_a);
        bot_inner_product += (*value_b - *value_a) * (*value_b - *value_a);
    }
    double inner_product = top_inner_product / bot_inner_product;

    for (int i = 0; i < n_dims; i++) {
        double *value_a = &points[i][a];
        double *value_b = &points[i][b];
        out_points[i][out_points_index] = inner_product * (*value_b - *value_a);
        out_points[i][out_points_index] += *value_a;
    }
}

void swap(double** points, int n_dims, int a, int b) {
    for (int i = 0; i < n_dims; i++) {
        double t = points[i][a];
        points[i][a] = points[i][b];
        points[i][b] = t;
    }
}
 
int partition (double** arr, int low, int high, int dim, int n_dims)
{ 
    double pivot = arr[dim][high]; // pivot 
    int i = (low - 1); // Index of smaller element and indicates the right position of pivot found so far
 
    for (int j = low; j <= high - 1; j++) 
    { 
        // If current element is smaller than the pivot 
        if (arr[dim][j] < pivot) 
        { 
            i++; // increment index of smaller element 
            swap(arr, n_dims, i, j);
        } 
    } 

    swap(arr, n_dims, i + 1, high);
    return (i + 1); 
} 
 
/* The main function that implements QuickSort 
arr[] --> Array to be sorted, 
low --> Starting index, 
high --> Ending index */
void quickSort(double** arr, int low, int high, int dim, int n_dims)
{ 
    if (low < high) 
    { 
        /* pi is partitioning index, arr[p] is now 
        at right place */
        int pi = partition(arr, low, high, dim, n_dims);
 
        // Separately sort elements before 
        // partition and after partition 
        quickSort(arr, low, pi - 1, dim, n_dims);
        quickSort(arr, pi + 1, high, dim, n_dims);
    } 
}

int main(int argc, char *argv[]) { 
    int n_dims;
    long n_samples;

    double** points = get_points(argc, argv, &n_dims, &n_samples);

    print_points(points, n_dims, n_samples);

    /*
     * Get furthest points
     */
    int a = get_furthest_point(points, 0, &n_dims, &n_samples); 
    int b = get_furthest_point(points, a, &n_dims, &n_samples);

    /*
     * Get ortho projection of points in line ab
     */
    int n_set = n_samples;
    double** ortho_points = malloc(sizeof(double*) * n_dims);
    for (int i = 0; i < n_dims; i++) {
        ortho_points[i] = malloc(sizeof(double) * n_set);
    }

    for (int j = 0, c = 0; j < n_set; j++, c++) {
        calc_ortho_projection(points, n_dims, a, b, j, c, ortho_points);
    }
    print_points(ortho_points, n_dims, n_set);

    /*
     * Sort ortho projection points
     */
    for (int i = 0; i < n_dims; i++) {
        quickSort(ortho_points, 0, n_set - 1, i, n_dims);
    }
    print_points(ortho_points, n_dims, n_set);

    /*
     * Get median point which will be the center of the ball
     */
    double* median_point = malloc(sizeof(double) * n_dims);
    if (n_set % 2 == 0) {
        int mid_point_i = n_set / 2 - 1;
        for (int i = 0; i < n_dims; i++) {
            median_point[i] = (ortho_points[i][mid_point_i] + ortho_points[i][mid_point_i + 1]) / 2;
        }
    } else {
        int mid_point_i = n_set / 2;
        printf("%d\n", mid_point_i);
        for (int i = 0; i < n_dims; i++) {
            median_point[i] = ortho_points[i][mid_point_i];
        }
    }

    printf("Median point:\n");
    for(int i = 0; i < n_dims; i++) {
        printf("\tDimension: %d = %f\n", i, median_point[i]);
    }

    /*
     * Get the radius of the ball (largest distance)
     */
    // TODO: THIS IS TEMPORARY
    double distances[2] = {0.0, 0.0};

    for (int i = 0; i < n_dims; i++) {
        distances[0] += (points[i][a] - median_point[i]) * (points[i][a] - median_point[i]);
        distances[1] += (points[i][b] - median_point[i]) * (points[i][b] - median_point[i]);
    }

    distances[0] = sqrt(distances[0]);
    distances[1] = sqrt(distances[1]);

    printf("Distance 1: %f\n", distances[0]);
    printf("Distance 2: %f\n", distances[1]);

    double radius = (distances[0] > distances[1]) ? distances[0] : distances[1];
    printf("RADIUS: %f\n", radius);

    free(median_point);
    for (int i = 0; i < n_dims; i++) {
        free(ortho_points[i]);
    }
    free(ortho_points);
}
