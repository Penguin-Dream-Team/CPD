// OpenMP version

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <mpi.h>

#include "gen_points.h"
#include "quickSelect.h"

#define ELEM_SWAP(a,b) { register node_t t = (a); (a) = (b); (b) = t; }
#define WORLD MPI_COMM_WORLD
#define POINT_TAG 100
#define ARGS_TAG 101
#define COUNT_TAG 102
#define PRINT_TAG 103
#define FOR_TAG 104

typedef struct {
    double max_distance;
    long max;
} furthest_point;

double **points;
node_t *ortho_points;
int n_dims;
long max_size;
int max_threads;
int sent = 0;

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

int get_furthest_point_parallel(long point, long start, long end, int threads) {
    double *point_point = points[ortho_points[point].point_id];
    furthest_point *furthest_points = malloc((sizeof(furthest_point) + 2048) * threads);

    #pragma omp parallel num_threads(threads)
    {
        furthest_point fp = furthest_points[omp_get_thread_num()];
        #pragma omp for schedule(static) 
        for (long i = start; i < end; i++) {
            double distance = distance_sqrd(points[ortho_points[i].point_id], point_point);

            if (fp.max_distance < distance) {
                fp.max = i;
                fp.max_distance = distance;
            }
        }
        furthest_points[omp_get_thread_num()] = fp;
    }

    long max = point;
    double max_distance = 0.0;
    for (int i = 0; i < threads; i++) {
        if (max_distance < furthest_points[i].max_distance) {
            max = furthest_points[i].max;
            max_distance = furthest_points[i].max_distance;
        }
    }
    return max;
}

double get_furthest_distance_parallel(double *point, long start, long end, int threads) {
    double *furthest_distances = malloc((sizeof(double) + 2048) * threads);

    #pragma omp parallel num_threads(threads)
    {
        double fd = furthest_distances[omp_get_thread_num()];
        #pragma omp for schedule(static) 
        for (long i = start; i < end; i++) {
            double distance = distance_sqrd(points[ortho_points[i].point_id], point);
            if ((fd - distance) < 0) {
                fd = distance;
            }
        }
        furthest_distances[omp_get_thread_num()] = fd;
    }

    double max_distance = 0.0;
    for (int i = 0; i < threads; i++) {
        if (max_distance < furthest_distances[i]) {
            max_distance = furthest_distances[i];
        }
    }
    return max_distance;
}

void calc_ortho_projection(double *point_a, double *point_b, double *p1, double *p2, node_t *ortho_points, int p1_index, int p2_index) {
    double top_inner_product1 = 0;
    double top_inner_product2 = 0;
    double bot_inner_product = 0;

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
}

// not inclusive
node_t *build_tree(long start, long end) {  
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

// not inclusive
node_t *build_tree_parallel_omp(long start, long end, int threads) {  
    printf("FUNCTION OMP Threads: %d, Start: %ld, End: %ld\n", threads, start, end);

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
    long a = get_furthest_point_parallel(start, start, end, threads);
    long b = get_furthest_point_parallel(a, start, end, threads);

    double *point_a = points[ortho_points[a].point_id];
    double *point_b = points[ortho_points[b].point_id];

    /*
     * Get projections to allow median calc
     */
    #pragma omp parallel for num_threads(threads)
    for (int i = start; i < end; i++) {
        double projection = 0.0;
        double *point = points[ortho_points[i].point_id];
        for (int d = 0; d < n_dims; d++) {
            projection += (point[d] - point_a[d]) * (point_b[d] - point_a[d]);
        }
        ortho_points[i].center[0] = projection * (point_b[0] - point_a[0]);
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
    tree->radius = sqrt(get_furthest_distance_parallel(median_point, start, end, threads));

    #pragma omp parallel
    {
        #pragma omp single 
        {
            if (threads > 2) {
                #pragma omp task
                tree->L = build_tree_parallel_omp(start, median_ids.second, (threads + 1) / 2);
            } else {
                #pragma omp task
                tree->L = build_tree(start, median_ids.second);
            }

            if (threads > 2) {
                if (threads != 3) {
                    #pragma omp task
                    tree->R = build_tree_parallel_omp(median_ids.second, end, threads / 2);
                } else {
                    #pragma omp task
                    tree->R = build_tree(median_ids.second, end);
                }
            } else {
                #pragma omp task
                tree->R = build_tree(median_ids.second, end);
            }    
        }
        #pragma omp taskwait
    }
    return tree;
}

// not inclusive
node_t *build_tree_parallel_mpi(long start, long end, int process, int max_processes, int threads) {  
    printf("FUNCTION MPI Process: %d, Max_Process: %d, Threads: %d, Start: %ld, End: %ld\n", process, max_processes, threads, start, end);

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
    long a = get_furthest_point_parallel(start, start, end, threads);
    long b = get_furthest_point_parallel(a, start, end, threads);

    double *point_a = points[ortho_points[a].point_id];
    double *point_b = points[ortho_points[b].point_id];

    /*
     * Get projections to allow median calc
     */
    #pragma omp parallel for num_threads(threads)
    for (int i = start; i < end; i++) {
        double projection = 0.0;
        double *point = points[ortho_points[i].point_id];
        for (int d = 0; d < n_dims; d++) {
            projection += (point[d] - point_a[d]) * (point_b[d] - point_a[d]);
        }
        ortho_points[i].center[0] = projection * (point_b[0] - point_a[0]);
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
    tree->radius = sqrt(get_furthest_distance_parallel(median_point, start, end, threads));

    long diff = (max_processes - process) / 2;
    if (process + 1 < max_processes) {
        //SEND THE RIGHT
        long args[] = {end - median_ids.second, max_processes};
        MPI_Send(args, 2, MPI_LONG, process + diff, ARGS_TAG, WORLD);
        printf("Node %d, has sent points\n", process);
        long points_index[end - median_ids.second];
        for (int i = median_ids.second, j = 0; i < end; i++, j++){
            points_index[j] = ortho_points[i].point_id;
        }

        MPI_Send(points_index, end - median_ids.second, MPI_LONG, process + diff, POINT_TAG, WORLD);

        tree->L = build_tree_parallel_mpi(start, median_ids.second, 0, diff, threads);
    }
    else {
        #pragma omp parallel
        {
            #pragma omp single 
            {
                if (threads > 2) {
                    #pragma omp task
                    tree->L = build_tree_parallel_omp(start, median_ids.second, (threads + 1) / 2);
                } else {
                    #pragma omp task
                    tree->L = build_tree(start, median_ids.second);
                }

                if (threads > 2) {
                    if (threads != 3) {
                        #pragma omp task
                        tree->R = build_tree_parallel_omp(median_ids.second, end, threads / 2);
                    } else {
                        #pragma omp task
                        tree->R = build_tree(median_ids.second, end);
                    }
                } else {
                    #pragma omp task
                    tree->R = build_tree(median_ids.second, end);
                }    
            }
            #pragma omp taskwait
        }
    }

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
    count++;
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

    return count + 1;       
}

void print_tree(node_t *tree, int n_dims, double **points, int prev_count) {
    long n_count = 0;
    count_nodes(tree, &n_count);
    n_count = n_count + prev_count;
    printf("%d %ld\n", n_dims, n_count);
    aux_print_tree(tree, n_dims, points, n_count, 0);
}

void wait_mpi(int me) {
    node_t *new_ortho_points = malloc(sizeof(node_t) * max_size);
    MPI_Status statuses[3];

    // Receive Args
    long args[2];
    MPI_Recv(args, 2, MPI_LONG, MPI_ANY_SOURCE, ARGS_TAG, WORLD, &statuses[0]);

    // Receive point indexes
    long rec_index[args[0]];
    MPI_Recv(rec_index, args[0], MPI_LONG, MPI_ANY_SOURCE, POINT_TAG, WORLD, &statuses[1]);

    for(int i = 0; i < args[0]; i++){
        new_ortho_points[i].center = malloc(sizeof(double) * n_dims);
        *new_ortho_points[i].center = *points[rec_index[i]];
        new_ortho_points[i].point_id = rec_index[i];
        new_ortho_points[i].center[0] = 69.69;
    }
    
    printf("Node %d received points\n", me);
    ortho_points = new_ortho_points;
    node_t * tree = build_tree_parallel_mpi(0, args[0], me, args[1], max_threads);

    // Send node count
    long n_count = 0;
    count_nodes(tree, &n_count);
    long node_count[1] = {n_count};
    MPI_Send(node_count, 1, MPI_LONG, 0, COUNT_TAG, WORLD);

    // Recieve print command
    int print[1];
    MPI_Recv(print, 1, MPI_INT, 0, PRINT_TAG, WORLD, &statuses[2]);

    //aux_print_tree(tree, n_dims, points, n_count, 0);

    free(new_ortho_points);

    MPI_Finalize();
    exit(0);
}

int main(int argc, char *argv[]) {
    double exec_time;
    long n_samples;

    omp_set_nested(1);
    omp_set_dynamic(1);
    max_threads = omp_get_max_threads();


    int me, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);

    exec_time = -omp_get_wtime();
    points = get_points(argc, argv, &n_dims, &n_samples);
    max_size = n_samples;

    /*
     * Get ortho projection of points in line ab
     */
    ortho_points = malloc(sizeof(node_t) * n_samples);
    node_t * tree;

    #pragma omp parallel for 
    for (long i = 0; i < n_samples; i++) {
        ortho_points[i].center = malloc(sizeof(double) * n_dims);
        ortho_points[i].point_id = i;
    }

    if (me != 0)
        wait_mpi(me);

    tree = build_tree_parallel_mpi(0, n_samples, 0, nprocs, max_threads);



    // Receive node count
    int node_count = 0;
    int count[1];
    MPI_Status statuses[nprocs];
    for (int i = 1; i < nprocs; i++){
        count[0] = 0;
        MPI_Recv(count, 1, MPI_LONG, MPI_ANY_SOURCE, COUNT_TAG, WORLD, &statuses[i]);
        node_count += count[0];
    }

    exec_time += omp_get_wtime();
    fprintf(stderr, "%lf\n", exec_time);

    // Send print command
    int print[1] = {1};
    for (int i = 1; i < nprocs; i++){
        MPI_Send(print, 1, MPI_INT, i, PRINT_TAG, WORLD);
    }

    //print_tree(tree, n_dims, points, node_count);

    free(ortho_points);
    free_node(tree);
    free(points[0]);
    free(points);

    MPI_Finalize();

}