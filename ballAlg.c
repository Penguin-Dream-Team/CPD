// Serial version

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "gen_points.h"
#include "quickSelect.h"

#define ELEM_SWAP(a,b) { register node_t t = (a); (a) = (b); (b) = t; }

double **points;
node_t *ortho_points;
int n_dims;

node_t *create_node(double *point, long id, double radius) {
    node_t *node = malloc(sizeof(node_t));

    node->point_id = id;
    node->center = point;
    node->radius = radius;
    node->L = NULL;
    node->R = NULL;

    return node;
}

void free_node(node_t *node) {
    if (node->radius != 0) {
        free(node->center);
    }
    if (node->L) {
        free_node(node->L);
    }
    if (node->R) {
        free_node(node->R);
    }
    free(node);
}

double distance_sqrd(double *pt1, double *pt2) {
    double dist = 0.0;

    for (int d = 0; d < n_dims; d++) {
        double tmp = pt1[d] - pt2[d];
        dist += tmp * tmp;
    }

    return dist;
}

double distance(double *pt1, double *median) {
    return sqrt(distance_sqrd(pt1, median));
}

int get_furthest_point(long point, long start, long end) {
    long max = point;
    double *point_point = points[ortho_points[point].point_id];
    double distance, max_distance = 0.0;

    for (long i = start; i < end; i++) {
        distance = distance_sqrd(points[ortho_points[i].point_id], point_point);
        if (max_distance < distance) {
            max = i;
            max_distance = distance;
        }
    }
    return max;
}

double get_furthest_distance(double *point, long start, long end) {
    double distance, max_distance = 0.0;

    for (long i = start; i < end; i++) {
        distance = distance_sqrd(points[ortho_points[i].point_id], point);
        if ((max_distance - distance) < 0) {
            max_distance = distance;
        }
    }
    return max_distance;
}

void calc_ortho_projection(double *point_a, double *point_b, double *p1, double *p2, node_t *ortho_points, int p1_index, int p2_index) {
    double top_inner_product1 = 0;
    double top_inner_product2 = 0;
    double bot_inner_product = 0;

    fprintf(stderr, "point a inside calc =");
    for (int i = 0; i < n_dims; i++) {
        fprintf(stderr, " %f", point_a[i]);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "point b inside calc =");
    for (int i = 0; i < n_dims; i++) {
        fprintf(stderr, " %f", point_b[i]);
    }
    fprintf(stderr, "\n");

    for (int i = 0; i < n_dims; i++) {
        double b_minus_a = point_b[i] - point_a[i];
        top_inner_product1 += (p1[i] - point_a[i]) * b_minus_a;
        top_inner_product2 += (p2[i] - point_a[i]) * b_minus_a;
        bot_inner_product += b_minus_a * b_minus_a;
    }

    double inner_product1 = top_inner_product1 / bot_inner_product;
    double inner_product2 = top_inner_product2 / bot_inner_product;

    for (int i = 0; i < n_dims; i++) {
        ortho_points[p1_index].center[i] = inner_product1 * (point_b[i] - point_a[i]) + point_a[i];
        ortho_points[p2_index].center[i] = inner_product2 * (point_b[i] - point_a[i]) + point_a[i];
    }

    fprintf(stderr, "median 1 inside calc =");
    for (int i = 0; i < n_dims; i++) {
        fprintf(stderr, " %f", ortho_points[p1_index].center[i]);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "median 2 inside calc =");
    for (int i = 0; i < n_dims; i++) {
        fprintf(stderr, " %f", ortho_points[p2_index].center[i]);
    }
    fprintf(stderr, "\n");
}

// not inclusive
node_t *build_tree(long start, long end) {  
    fprintf(stderr, "BUILD TREE BEGIN\n");
    for (int i = start; i < end; i++) {
        fprintf(stderr, "POINTS %d is", i);
        for (int j = 0; j < n_dims; j++) {
            fprintf(stderr, " %f", points[ortho_points[i].point_id][j]);
        }
        fprintf(stderr, "\n");
    }
    if (start == end - 1) {  // 1 point
        return create_node(points[ortho_points[start].point_id], ortho_points[start].point_id, 0);

    } else if (start == end - 2) {  // 2 points
        double *median_point = malloc(sizeof(double) * n_dims);
        double *point1 = points[ortho_points[start].point_id];
        double *point2 = points[ortho_points[end - 1].point_id];

        for (int d = 0; d < n_dims; d++) {
            median_point[d] = (point1[d] + point2[d]) / 2;
        }
        double dist = distance(point1, median_point);

        node_t *tree = create_node(median_point, -1, dist);

        if ((point1[0] - point2[0]) < 0) {
            tree->L = create_node(point1, ortho_points[start].point_id, 0);
            tree->R = create_node(point2, ortho_points[end - 1].point_id, 0);
        } else {
            tree->L = create_node(point2, ortho_points[end - 1].point_id, 0);
            tree->R = create_node(point1, ortho_points[start].point_id, 0);
        }

        return tree;
    }

    /*
     * Get furthest points
     */
    long a = get_furthest_point(start, start, end);
    long b = get_furthest_point(a, start, end);

    double *point_a = points[ortho_points[a].point_id];
    double *point_b = points[ortho_points[b].point_id];

    fprintf(stderr, "Point a is");
    for (int j = 0; j < n_dims; j++) {
        fprintf(stderr, " %f", point_a[j]);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "Point b is");
    for (int j = 0; j < n_dims; j++) {
        fprintf(stderr, " %f", point_b[j]);
    }
    fprintf(stderr, "\n");

    /*
     * Get projections to allow median calc
     */
    for (int i = start; i < end; i++) {
        double projection = 0.0;
        double *point = points[ortho_points[i].point_id];
        for (int d = 0; d < n_dims; d++) {
            projection += (point[d] - point_a[d]) * (point_b[d] - point_a[d]);
        }
        ortho_points[i].center[0] = projection * (point_b[0] - point_a[0]);
    }

    for (int i = start; i < end; i++) {
        fprintf(stderr, "Projection at ortho_points %d is", i);
        for (int j = 0; j < n_dims; j++) {
            fprintf(stderr, " %f", ortho_points[i].center[j]);
        }
        fprintf(stderr, "\n");
    }
    
    /*
     * Get median point which will be the center of the ball
     */
    medianValues median_ids = quickSelect(ortho_points, start, end);

    node_t point_median_1 = ortho_points[median_ids.first];
    node_t point_median_2 = ortho_points[median_ids.second];

    /*
     * Separate L and R sets
     */
    double *median_point = malloc(sizeof(double) * n_dims);
    node_t *tree = create_node(median_point, -1, -1);

    /* Calc ortho projection of median points */
    double *p1 = points[point_median_1.point_id];
    double *p2 = points[point_median_2.point_id];
    fprintf(stderr, "median 1 ");
    for (int j = 0; j < n_dims; j++) {
        fprintf(stderr, " %f", p1[j]);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "median 2 ");
    for (int j = 0; j < n_dims; j++) {
        fprintf(stderr, " %f", p2[j]);
    }
    fprintf(stderr, "\n");

    calc_ortho_projection(point_a, point_b, p1, p2, ortho_points, median_ids.first, median_ids.second);

    /*
     * Get the radius of the ball (largest distance)
     */
    if ((end - start) % 2 != 0) {
        for (int d = 0; d < n_dims; d++) {
            median_point[d] = point_median_1.center[d];
        }
    } else {
        for (int d = 0; d < n_dims; d++) {
            median_point[d] = (point_median_1.center[d] + point_median_2.center[d]) / 2;
        }
    }
    tree->radius = sqrt(get_furthest_distance(median_point, start, end));

    tree->L = build_tree(start, median_ids.second);

    tree->R = build_tree(median_ids.second, end);

    return tree;
}

void count_nodes(node_t *tree, long *n_count) {
    n_count[0]++;

    if (tree->L) {
        count_nodes(tree->L, n_count);
    }
    if (tree->R) {
        count_nodes(tree->R, n_count);
    }
}

long aux_print_tree(node_t *tree, int n_dims, double **points,
        long n_count, long count) {
    long my_id, left_id = -1, right_id = -1;

    my_id = count;
    if (tree->L) {
        left_id = count + 1;
        count = aux_print_tree(tree->L, n_dims, points, n_count, left_id);
    }
    if (tree->R) {
        right_id = count + 1;
        count = aux_print_tree(tree->R, n_dims, points, n_count, right_id);
    }

    printf("%ld %ld %ld %f", my_id, left_id, right_id, tree->radius);
    for (int i = 0; i < n_dims; i++) {
        printf(" %f", tree->center[i]);
    }
    printf("\n");

    return count;
}

void print_tree(node_t *tree, int n_dims, double **points) {
    long n_count = 0;
    count_nodes(tree, &n_count);
    printf("%d %ld\n", n_dims, n_count);
    aux_print_tree(tree, n_dims, points, n_count, 0);
}

int main(int argc, char *argv[]) {
    double exec_time;
    long n_samples;

    if(argc != 4){
        printf("Usage: %s <n_dims> <n_points> <seed>\n", argv[0]);
        exit(1);
    }

    n_dims = atoi(argv[1]);
    if (n_dims < 2){
        printf("Illegal number of dimensions (%d), must be above 1.\n", n_dims);
        exit(2);
    }

    n_samples = atol(argv[2]);
    if (n_samples < 1){
        printf("Illegal number of points (%ld), must be above 0.\n", n_samples);
        exit(3);
    }

    exec_time = -omp_get_wtime();
    points = get_points(argc, argv, n_dims, n_samples);

    for (int i = 0; i < n_samples; i++) {
        fprintf(stderr, "Points at %d is", i);
        for (int j = 0; j < n_dims; j++) {
            fprintf(stderr, " %f", points[i][j]);
        }
        fprintf(stderr, "\n");
    }

    /*
     * Get ortho projection of points in line ab
     */
    ortho_points = malloc(sizeof(node_t) * n_samples);
    for (long i = 0; i < n_samples; i++) {
        ortho_points[i].center = malloc(sizeof(double) * n_dims);
        ortho_points[i].point_id = i;
    }

    node_t * tree;
    tree = build_tree(0, n_samples);

    exec_time += omp_get_wtime();
    fprintf(stderr, "%lf\n", exec_time);

    print_tree(tree, n_dims, points);

    free(ortho_points);
    free_node(tree);
    free(points[0]);
    free(points);
}
