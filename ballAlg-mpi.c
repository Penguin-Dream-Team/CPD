// MPI version

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <omp.h>
#include <mpi.h>
#include <math.h>

#include "ballAlg-large-mpi.h"
#include "gen_points.h"
#include "quickSelect.h"

#define ELEM_SWAP(a,b) { register node_t t = (a); (a) = (b); (b) = t; }
#define WORLD MPI_COMM_WORLD
#define POINT_TAG        100
#define ARGS_TAG         101
#define CONFIRMATION_TAG 102
#define COUNT_TAG        103
#define NODE_TAG         104
#define PRINT_TAG        105
#define PRINT_2_TAG      106
#define FOR_TAG          107
#define FOR_RESPONSE_TAG 108
#define EARLY_END_TAG    109

typedef struct {
    double max_distance;
    long max;
} furthest_point;

double **points;
node_t *ortho_points;
int n_dims, max_threads, sent = 0;
long max_size;
int for_it = 0;
int current_print_proc = 0;
int nprocs = 0;

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
    furthest_point *furthest_points = malloc(sizeof(furthest_point) * threads);

    #pragma omp parallel num_threads(threads)
    {
        furthest_point fp;
	    fp.max = point;
        fp.max_distance = 0.0;
        #pragma omp for schedule(static) 
        for (long i = start; i < end; i++) {
            double distance = distance_sqrd(points[ortho_points[i].point_id], point_point);

            if ((fp.max_distance - distance) < 0) {
                fp.max = i;
                fp.max_distance = distance;
            }
        }
	
        furthest_points[omp_get_thread_num()] = fp;
    }

    long max = point;
    double max_distance = 0.0;
    for (int i = 0; i < threads; i++) {
        if ((max_distance - furthest_points[i].max_distance) < 0) {
            max = furthest_points[i].max;
            max_distance = furthest_points[i].max_distance;
        }
    }
    return max;
}

double get_furthest_distance_parallel(double *point, long start, long end, int threads) {
    double *furthest_distances = malloc(sizeof(double) * threads);

    #pragma omp parallel num_threads(threads)
    {
        double fd = 0.0;
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
        if ((max_distance - furthest_distances[i]) < 0) {
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
void calc_projections(long start, long end, int threads, double *point_a, double *point_b){

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
        //printf("Acessing memory position %d with value %f\n", i, projection * (point_b[0] - point_a[0]));
        ortho_points[i].center[0] = projection * (point_b[0] - point_a[0]);
    }
}

// not inclusive
void calc_projections_mpi(long start, long end, int threads, long interval, long processor, long a, long b, MPI_Comm comm){
    
    /*
     * Get furthest points
     */

    double *point_a = points[a];
    double *point_b = points[b];
    
    double *response = malloc(sizeof(double) * interval);
    
    /*
     * Get projections to allow median calc
     */
    #pragma omp parallel for num_threads(threads)
    for (long i = start; i < end; i++) {
        double projection = 0.0;
        double *point = points[ortho_points[i].point_id];
        for (int d = 0; d < n_dims; d++) {
            projection += (point[d] - point_a[d]) * (point_b[d] - point_a[d]);
        }
        response[i - start] = projection * (point_b[0] - point_a[0]);
    }

    
    double recieved[1];
    MPI_Gather(response, interval, MPI_DOUBLE, recieved, interval, MPI_DOUBLE, processor, comm);
}


// not inclusive
node_t *build_tree_parallel_omp(long start, long end, int threads, int me) {  
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
                tree->L = build_tree_parallel_omp(start, median_ids.second, (threads + 1) / 2, me);
            } else {
                #pragma omp task
                tree->L = build_tree(start, median_ids.second);
            }

            if (threads > 2) {
                if (threads != 3) {
                    #pragma omp task
                    tree->R = build_tree_parallel_omp(median_ids.second, end, threads / 2, me);
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
    

    int nprocs = max_processes - process;
    long interval = (end - start) / nprocs;

    /*
     * Get furthest points
     */

    long a = get_furthest_point_parallel(start, start, end, threads);
    long b = get_furthest_point_parallel(a, start, end, threads);

    double *point_a = points[ortho_points[a].point_id];
    double *point_b = points[ortho_points[b].point_id];
   
    long a1 = ortho_points[a].point_id;
    long b1 = ortho_points[b].point_id;

    int group_size = max_processes - process;
    int ranks[group_size];
    ranks[0] = process;

    double *recieved = malloc(sizeof(double) * (end - start));
    double *response = malloc(sizeof(double) * interval);
    
    // No need to send points on the first round
    if (for_it == 0){
        // Sending some arguments
        long first_for[] = {interval, process, a1, b1};
        for (int i = process + 1, d = 1; i < max_processes; i++, d++) {
            MPI_Send(first_for, 4, MPI_LONG, i, FOR_TAG, WORLD);
        }
        for_it = 1;
        calc_projections(start, start + interval, threads, point_a, point_b);

        MPI_Gather(response, interval, MPI_DOUBLE, recieved, interval, MPI_DOUBLE, process, WORLD);

    }
    // Second for onwards
    else {
        for_it++;
        long *points_indexes = malloc(sizeof(long) * (end - start) );
        long *fakes = malloc(sizeof(long) * interval );

        long for_args[] = {interval, process, a1, b1};
        // Sending some arguments
        for (int i = process + 1, d = 1; i < max_processes; i++, d++) {
            MPI_Send(for_args, 4, MPI_LONG, i, FOR_TAG, WORLD);
            ranks[d] = i;
        }
            
        for (int j = start, k = 0; j < end; j++, k++){
            points_indexes[k] = ortho_points[j].point_id;
        }

        MPI_Group world_group;
        MPI_Comm_group(WORLD, &world_group);

        MPI_Group prime_group;
        MPI_Group_incl(world_group, group_size, ranks, &prime_group);

        // Create a new communicator based on the group
        MPI_Comm prime_comm;
        MPI_Comm_create_group(MPI_COMM_WORLD, prime_group, 0, &prime_comm);

        int root_rank;
        MPI_Comm_rank(prime_comm, &root_rank);

        // Scattering indexes of points
        MPI_Scatter(points_indexes, interval, MPI_LONG, fakes, interval, MPI_LONG, 0, prime_comm);

        // Calculating its own work
        calc_projections(start, start + interval, threads, point_a, point_b);

        // Calculate Rest
        calc_projections(start + interval * nprocs, end, threads, point_a, point_b);

        // Gathering responses
        MPI_Gather(response, interval, MPI_DOUBLE, recieved, interval, MPI_DOUBLE, 0, prime_comm);

        MPI_Group_free(&prime_group);
        MPI_Comm_free(&prime_comm);
    }
    
    // Recording the responses in desidered structure
    for (int i = start + interval, d = interval; i < end; i++, d++) {
        ortho_points[i].center[0] = recieved[d];
    }
    
    free(recieved);
    free(response);

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

    int diff = (max_processes - process) / 2;
    if (process + 1 < max_processes) {
        long size = end - median_ids.second;
        long new_max_processes = process + diff;

        // Send the right subtree to another processor

        // Send some arguments
        long args[] = {size, max_processes};
        MPI_Send(args, 2, MPI_LONG, new_max_processes, ARGS_TAG, WORLD);

        long *points_index = malloc(sizeof(long) * size);
        for (int i = median_ids.second, j = 0; i < end; i++, j++){
            points_index[j] = ortho_points[i].point_id;
        }

        // Sending the indexes of poitns
        MPI_Send(points_index, size, MPI_LONG, new_max_processes, POINT_TAG, WORLD);

        tree->L = build_tree_parallel_mpi(start, median_ids.second, process, new_max_processes, threads);
        

    } else {
        #pragma omp parallel
        {
            #pragma omp single 
            {
                if (threads > 2) {
                    #pragma omp task
                    tree->L = build_tree_parallel_omp(start, median_ids.second, (threads + 1) / 2, process);
                } else {
                    #pragma omp task
                    tree->L = build_tree(start, median_ids.second);
                }

                if (threads > 2) {
                    if (threads != 3) {
                        #pragma omp task
                        tree->R = build_tree_parallel_omp(median_ids.second, end, threads / 2, process);
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

long aux_print_tree_main_proc(node_t *tree, int n_dims, double **points,
        long n_count, long count, int me) {
    long my_id, left_id = -1, right_id = -1;
    bool left = false;

    my_id = count;
    if (tree->L) {
        left = true;
        left_id = count + 1;
        count = aux_print_tree_main_proc(tree->L, n_dims, points, n_count, left_id, me);
    }
    if (tree->R) {
        right_id = count + 1;
        count = aux_print_tree_main_proc(tree->R, n_dims, points, n_count, right_id, me);
    
    } else if (left) {
        right_id = count + 1;

        // Send node id
        long print_type = 0;
        int process_id = nprocs;
        while (process_id > 2) {
            process_id /= 2;
            if (process_id + me == current_print_proc) {
                print_type = 1;
                break;
            }
        }

        long print[] = {right_id, print_type, me};
        MPI_Send(print, 3, MPI_LONG, current_print_proc, NODE_TAG, WORLD);

        long count_received[1];
        MPI_Recv(count_received, 1, MPI_LONG, current_print_proc, PRINT_TAG, WORLD, MPI_STATUS_IGNORE);
        count = count_received[0];

        if (print_type == 1) {
            int next_print_proc[1];
            MPI_Recv(next_print_proc, 1, MPI_INT, current_print_proc, PRINT_2_TAG, WORLD, MPI_STATUS_IGNORE);
            current_print_proc = next_print_proc[0] - 1;
        }        

        current_print_proc++;
    }

    printf("%ld %ld %ld %f", my_id, left_id, right_id, tree->radius);
    for (int i = 0; i < n_dims; i++) {
        printf(" %f", tree->center[i]);
    }
    printf("\n");
    
    return count;       
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

void print_tree_mpi(node_t *tree, int n_dims, double **points, int prev_count, int me) {
    long n_count = 0;

    count_nodes(tree, &n_count);
    n_count = n_count + prev_count;
    printf("%d %ld\n", n_dims, n_count);
    fflush(stdout);
    MPI_Barrier(WORLD);

    current_print_proc = me + 1;
    aux_print_tree_main_proc(tree, n_dims, points, n_count, 0, me);
}

void print_tree(node_t *tree, int n_dims, double **points, int prev_count) {
    long n_count = 0;

    count_nodes(tree, &n_count);
    n_count = n_count + prev_count;
    printf("%d %ld\n", n_dims, n_count);

    aux_print_tree(tree, n_dims, points, n_count, 0);
}

int custom_power(int base, int exponent){
    int result = 1;
    for (int i = 0 ; i < exponent; i++) {
        result = result * base;
    }
    return result;
}

void wait_mpi(int me, int start, int end, int threads) {

    int diff = (end - start) / 2 + start;
    
    long for_args[4];
    MPI_Recv(for_args, 4, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, WORLD, MPI_STATUS_IGNORE);

    long interval = for_args[0];
    int res_process = for_args[1];
    int a = for_args[2];
    int b = for_args[3];


    if (for_it == 0){
        calc_projections_mpi(interval * me, interval * me + interval, threads, interval, res_process, a, b, WORLD);
        for_it = 1;
    }
    else {

        int group_size = nprocs / custom_power(2, for_it);
        int delta = res_process + group_size;
        int ranks[group_size];
        for (int i = delta - group_size, d = 0; d < group_size ; i++, d++) {
            ranks[d] = i;
        }

        MPI_Group world_group;
        MPI_Comm_group(MPI_COMM_WORLD, &world_group);

        MPI_Group prime_group;
        MPI_Group_incl(world_group, group_size, ranks, &prime_group);

        // Create a new communicator based on the group
        MPI_Comm prime_comm;
        MPI_Comm_create_group(MPI_COMM_WORLD, prime_group, 0, &prime_comm);

        //Receiving points
        long *for_indexes = malloc(sizeof(long) * interval);
        MPI_Scatter(for_indexes, interval, MPI_LONG, for_indexes, interval, MPI_LONG, 0, prime_comm);

        for(int i = 0; i < interval; i++){
            memset(ortho_points[i].center, 0, n_dims);
            ortho_points[i].point_id = for_indexes[i];
        }

        calc_projections_mpi(0, interval, threads, interval, 0, a, b, prime_comm);
        for_it++;
	
        MPI_Group_free(&prime_group);
        MPI_Comm_free(&prime_comm);

    }
    
    if (me < diff){
        wait_mpi(me, start, diff, threads);
    }

    else if (me > diff){
        wait_mpi(me, diff, end, threads);
    }
    

    // Receive Args
    long args[2];
    MPI_Recv(args, 2, MPI_LONG, MPI_ANY_SOURCE, ARGS_TAG, WORLD, MPI_STATUS_IGNORE);

    // Receive point indexes
    long *rec_index = malloc(sizeof(long) * args[0]);
    MPI_Recv(rec_index, args[0], MPI_LONG, MPI_ANY_SOURCE, POINT_TAG, WORLD, MPI_STATUS_IGNORE);

    for(int i = 0; i < args[0]; i++){
        memset(ortho_points[i].center, 0, n_dims);
        ortho_points[i].point_id = rec_index[i];
    }
    
    node_t * tree = build_tree_parallel_mpi(0, args[0], me, args[1], max_threads);

    // Send end confirmation
    int confirmation[1] = {1};
    MPI_Send(confirmation, 1, MPI_INT, 0, CONFIRMATION_TAG, WORLD);

    MPI_Barrier(WORLD);

    // Send node count
    long n_count = 0;
    count_nodes(tree, &n_count);
    long node_count[1] = {n_count};
    MPI_Send(node_count, 1, MPI_LONG, 0, COUNT_TAG, WORLD);

    MPI_Barrier(WORLD);

    // Receive node id and print type 
    long node_id[3];
    MPI_Recv(node_id, 3, MPI_LONG, MPI_ANY_SOURCE, NODE_TAG, WORLD, MPI_STATUS_IGNORE);

    long print_result[1];
    if (node_id[1]) {
        current_print_proc = me + 1;
        print_result[0] = aux_print_tree_main_proc(tree, n_dims, points, n_count, node_id[0], me);
    } else {  
        print_result[0] = aux_print_tree(tree, n_dims, points, n_count, node_id[0]);
    }

    MPI_Send(print_result, 1, MPI_LONG, node_id[2], PRINT_TAG, WORLD);

    if (node_id[1]) {
        MPI_Send(&current_print_proc, 1, MPI_INT, node_id[2], PRINT_2_TAG, WORLD);
    }

    MPI_Finalize();
    exit(0);
}

int ballAlg_mpi(int argc, char *argv[], long n_samples) {
    double exec_time;

    omp_set_nested(1);
    omp_set_dynamic(1);
    max_threads = omp_get_max_threads();

    exec_time = -omp_get_wtime();

    int me;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    
    points = get_points(argc, argv, n_dims, n_samples);
    max_size = n_samples / 2;

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

    // If the problem is too simple to distribute, simple omp version is used
    if (nprocs / 2 * 3 > n_samples){
        if (me == 0) {
            tree = build_tree_parallel_omp(0, n_samples, max_threads, me);

            exec_time += omp_get_wtime();
            fprintf(stderr, "%lf\n", exec_time);
            print_tree(tree, n_dims, points, 0);

            free(ortho_points);
            free_node(tree);
            free(points[0]);
            free(points);
        }
        MPI_Finalize();
        exit(0);
    }

    if (me != 0)
        wait_mpi(me, 0, nprocs, max_threads);

    tree = build_tree_parallel_mpi(0, n_samples, 0, nprocs, max_threads);
    
    // Receive confirmation
    int confirmation[1];
    for (int i = 1; i < nprocs; i++){
        MPI_Recv(confirmation, 1, MPI_LONG, MPI_ANY_SOURCE, CONFIRMATION_TAG, WORLD, MPI_STATUS_IGNORE);
    }

    exec_time += omp_get_wtime();
    fprintf(stderr, "%lf\n", exec_time);
    fflush(stderr);
    MPI_Barrier(WORLD);

    // Receive node count
    int node_count = 0;
    int count[1];
    for (int i = 1; i < nprocs; i++){
        count[0] = 0;
        MPI_Recv(count, 1, MPI_LONG, MPI_ANY_SOURCE, COUNT_TAG, WORLD, MPI_STATUS_IGNORE);
        node_count += count[0];
    }

    print_tree_mpi(tree, n_dims, points, node_count, me);

    free(ortho_points);
    free_node(tree);

    free(points[0]);
    free(points);

    MPI_Finalize();
    exit(0);
}

int main(int argc, char *argv[]) {
    long np;

    if(argc != 4){
        printf("Usage: %s <n_dims> <n_points> <seed>\n", argv[0]);
        exit(1);
    }

    n_dims = atoi(argv[1]);
    if (n_dims < 2){
        printf("Illegal number of dimensions (%d), must be above 1.\n", n_dims);
        exit(2);
    }

    np = atol(argv[2]);
    if (np < 1){
        printf("Illegal number of points (%ld), must be above 0.\n", np);
        exit(3);
    }

    if (np * n_dims < 100000000) {
        ballAlg_mpi(argc, argv, np);
    } else {
        ballAlg_large_mpi(argc, argv, np);
    }
}
