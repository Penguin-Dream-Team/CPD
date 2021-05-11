// MPI version

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <omp.h>
#include <mpi.h>

#include "gen_points.h"
#include "quickSelect.h"

#define ELEM_SWAP(a,b) { register node_t t = (a); (a) = (b); (b) = t; }
#define WORLD MPI_COMM_WORLD
#define POINT_TAG 100
#define ARGS_TAG 101
#define CONFIRMATION_TAG 102
#define COUNT_TAG 103
#define NODE_TAG 104
#define PRINT_TAG 105
#define FOR_TAG 106
#define FOR_RESPONSE_TAG 107
#define EARLY_END_TAG 108

typedef struct {
    double max_distance;
    long max;
} furthest_point;

double **points;
double *recieved, *response;
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
    //printf("calc projections mpi start %ld, end %ld, interval %ld\n", start, end, interval);
    
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
    //printf("Gathering on slave to %d \n", processor);
    //MPI_Send(response, interval, MPI_DOUBLE, processor, FOR_RESPONSE_TAG, comm);
    MPI_Gather(response, interval, MPI_DOUBLE, recieved, interval, MPI_DOUBLE, processor, comm);
    //printf("Sent for response to processor %d \n", processor);
}


// not inclusive
node_t *build_tree_parallel_omp(long start, long end, int threads, int me) {  
    //fprintf(stderr, "FUNCTION OMP on Processor: %d, Threads: %d, Start: %ld, End: %ld\n", me, threads, start, end);

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
                //printf("Process %d going serial \n", me);
                tree->L = build_tree(start, median_ids.second);
            }

            if (threads > 2) {
                if (threads != 3) {
                    #pragma omp task
                    tree->R = build_tree_parallel_omp(median_ids.second, end, threads / 2, me);
                } else {
                    #pragma omp task
                    //printf("Process %d going serial \n", me);
                    tree->R = build_tree(median_ids.second, end);
                }
            } else {
                #pragma omp task
                //printf("Process %d going serial \n", me);
                tree->R = build_tree(median_ids.second, end);
            }    
        }
        #pragma omp taskwait
    }
    return tree;
}

// not inclusive
node_t *build_tree_parallel_mpi(long start, long end, int process, int max_processes, int threads) {  
    //fprintf(stderr, "FUNCTION MPI Process: %d, Max_Process: %d, Threads: %d, Start: %ld, End: %ld\n", process, max_processes, threads, start, end);

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

    memset(response, 0, sizeof(double) * interval);
    memset(recieved, 0, sizeof(double) * (end - start));
    
    // No need to send points on the first round
    if (for_it == 0){
        //printf("++++Processor %ld Starting first for\n", process);
        long first_for[] = {interval, process, a1, b1};
        for (int i = process + 1, d = 1; i < max_processes; i++, d++) {
            MPI_Send(first_for, 4, MPI_LONG, i, FOR_TAG, WORLD);
        }
        for_it = 1;
        calc_projections(start, start + interval, threads, point_a, point_b);
        //printf("*Processor %d CALCULATED PROJECTIONS\n", process);

        //printf("Gathering on root\n");
        MPI_Gather(response, interval, MPI_DOUBLE, recieved, interval, MPI_DOUBLE, process, WORLD);

        //printf("\n Gathered:\n");
        //for (int i = start, d = 0; i < end; i++, d++) {
            //printf("%f\n", recieved[d]);
        //}
        //printf("\n");

    }
    // Second for onwards
    else {
        for_it++;
        //printf("----Processor %ld Starting second for\n", process);
        long *points_indexes = malloc(sizeof(long) * (end - start) );
        long *fakes = malloc(sizeof(long) * interval );

        long for_args[] = {interval, process, a1, b1};
        for (int i = process + 1, d = 1; i < max_processes; i++, d++) {
            MPI_Send(for_args, 4, MPI_LONG, i, FOR_TAG, WORLD);
            //printf("Processor %d Sent args to %d\n",process, i);
            ranks[d] = i;
            //printf("Adding rank %d\n", i);
        }
            
        for (int j = start, k = 0; j < end; j++, k++){
            points_indexes[k] = ortho_points[j].point_id;
        }

        //printf("Create world group\n");
        MPI_Group world_group;
        MPI_Comm_group(WORLD, &world_group);

        //printf("Create prime group\n");
        MPI_Group prime_group;
        MPI_Group_incl(world_group, group_size, ranks, &prime_group);

        // Create a new communicator based on the group
        //printf("Create prime comm\n");
        MPI_Comm prime_comm;
        MPI_Comm_create_group(MPI_COMM_WORLD, prime_group, 0, &prime_comm);

        int root_rank;
        MPI_Comm_rank(prime_comm, &root_rank);
        //printf("Root rank is: %d\n", root_rank);

        //printf("-*/-*/-*/-*/Processor %d, scattering\n", process);
        MPI_Scatter(points_indexes, interval, MPI_LONG, fakes, interval, MPI_LONG, 0, prime_comm);
        //MPI_Send(points_indexes[i], interval, MPI_LONG, i, POINT_TAG, WORLD);
        //printf("Processor %d Sent points\n",process);

        //printf("\nSCATETING Mastre\n");
        //for(int i = 0; i < interval; i++){
            //printf("*** Master process %d Received point %ld and working on %ld\n", process, fakes[i], ortho_points[start+i].point_id);
        //}
        //printf("\n");
        
        //printf("*Processor %d CALCULATING PROJECTIONS with start %ld and end %ld\n", process, start, start+interval);
        calc_projections(start, start + interval, threads, point_a, point_b);
        //printf("*Processor %d CALCULATED PROJECTIONS\n", process);

        // Calculate Rest
        //printf("Calculating rest from: %ld to %ld on process %d\n", start + interval * nprocs, end, process);
        calc_projections(start + interval * nprocs, end, threads, point_a, point_b);

        //printf("Gathering on root %d\n", process);
        MPI_Gather(response, interval, MPI_DOUBLE, recieved, interval, MPI_DOUBLE, 0, prime_comm);
    }
    

    for (int i = start + interval, d = interval; i < end; i++, d++) {
        //printf("Process %d accesing memory %d with value %f on indice %d\n", process, i, recieved[d], d);
        ortho_points[i].center[0] = recieved[d];
    }
    

    /*
    double *projections_mpi = malloc(sizeof(double) * interval);
    for (int i = process + 1, d = 1; i < max_processes; i++, d++) {
        //printf("Processor %d Recieving projections from processor %d\n", process, i);
        //MPI_Recv(projections_mpi, interval, MPI_DOUBLE, i, FOR_RESPONSE_TAG, WORLD, &status);

        //printf("Processor %d Recieved projections from processor %d\n", process, i);
        for (int j = start + interval * d, p = 0; j < start + interval * d + interval; j++, p++) {
            //printf("----------Processor %d accessing memory %d from received number %d\n", process, j, p);
            ortho_points[j].center[0] = projections_mpi[p];
        }
    }
    */

    //fprintf(stderr, "PROCESS %d - Getting median points\n", process);
    /*
     * Get median point which will be the center of the ball
     */
    medianValues median_ids = quickSelect(ortho_points, start, end);

    node_t point_median_1 = ortho_points[median_ids.first];
    node_t point_median_2 = ortho_points[median_ids.second];

    //fprintf(stderr, "PROCESS %d - Separate L and R\n", process);
    /*
     * Separate L and R sets
     */
    double *median_point = malloc(sizeof(double) * n_dims);
    node_t *tree = create_node(median_point, -1, -1);

    //fprintf(stderr, "PROCESS %d - Calc ortho projection of median points\n", process);
    /* Calc ortho projection of median points */
    double *p1 = points[point_median_1.point_id];
    double *p2 = points[point_median_2.point_id];
    calc_ortho_projection(point_a, point_b, p1, p2, ortho_points, median_ids.first, median_ids.second);

    //fprintf(stderr, "PROCESS %d - Getting radius of the ball\n", process);
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

    //fprintf(stderr, "PROCESS %d - Calculating furthest distance from start %ld to end %ld and median %f\n", process, start, end, median_point[0]);
    // ISTO CRASHA EM PARALELO
    //tree->radius = sqrt(get_furthest_distance_parallel(median_point, start, end));
    tree->radius = sqrt(get_furthest_distance(median_point, start, end));
    //fprintf(stderr, "PROCESS %d - Calculated furthest distance\n", process);

    int diff = (max_processes - process) / 2;
    //fprintf(stderr, "PROCESS %d - Calculating diff %d\n", process, diff);
    if (process + 1 < max_processes) {
        long size = end - median_ids.second;
        long new_max_processes = process + diff;

        //SEND THE RIGHT
        //fprintf(stderr, "Node %d, is going to send args to process %ld when diff is %ld\n", process, new_max_processes, diff);
        long args[] = {size, max_processes};
        MPI_Send(args, 2, MPI_LONG, new_max_processes, ARGS_TAG, WORLD);
        //fprintf(stderr, "Node %d, has sent args to process %ld with size %ld\n", process, new_max_processes, size);

        long *points_index = malloc(sizeof(long) * size);
        //printf("\nSending TO THE RIGHT\n");
        for (int i = median_ids.second, j = 0; i < end; i++, j++){
            points_index[j] = ortho_points[i].point_id;
            //printf("**Sending point %ld to %d\n", ortho_points[i].point_id, new_max_processes);
        }
        //printf("\n");

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
        if ((nprocs - me) / 2 == current_print_proc) print_type = 1;

        long print[] = {right_id, print_type};
        //fprintf(stderr, "I am process %d sending print to process %d with value %ld\n", me, current_print_proc, print_type);
        MPI_Send(print, 2, MPI_LONG, current_print_proc, NODE_TAG, WORLD);

        long count_received[1];
        MPI_Status statuses[1];
        MPI_Recv(count_received, 1, MPI_LONG, current_print_proc, PRINT_TAG, WORLD, &statuses[0]);
        //fprintf(stderr, "I am process %d | received print from process %d\n", me, current_print_proc);

        /*if (current_print_proc == 2) {
            MPI_Finalize();
            exit(0);
        }*/
        count = count_received[0];
        
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

    current_print_proc = 1;
    aux_print_tree_main_proc(tree, n_dims, points, n_count, 0, me);
}

void print_tree(node_t *tree, int n_dims, double **points, int prev_count) {
    long n_count = 0;

    count_nodes(tree, &n_count);
    n_count = n_count + prev_count;
    printf("%d %ld\n", n_dims, n_count);

    current_print_proc = 1;
    aux_print_tree(tree, n_dims, points, n_count, 0);
}

void finish_early_mpi(){
    // Send end confirmation
    int confirmation[1] = {1};
    MPI_Send(confirmation, 1, MPI_INT, 0, CONFIRMATION_TAG, WORLD);


    // Send node count
    long node_count[1] = {0};
    MPI_Send(node_count, 1, MPI_LONG, 0, COUNT_TAG, WORLD);


    MPI_Finalize();
    exit(0);
}

void wait_mpi(int me, int start, int end, int threads) {
    MPI_Status statuses[5];

    int diff = (end - start) / 2 + start;
    
    long for_args[4];
    MPI_Recv(for_args, 4, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, WORLD, &statuses[4]);

    long interval = for_args[0];
    int res_process = for_args[1];
    int a = for_args[2];
    int b = for_args[3];


    if (for_it == 0){
        //printf("Process %d received first for with interval %ld\n", me, interval);
        calc_projections_mpi(interval * me, interval * me + interval, threads, interval, res_process, a, b, WORLD);
        for_it = 1;
    }
    else {
        //printf("Process: %d, receiving second for with interval %ld\n", me, interval);

        int group_size = nprocs / (for_it+1);
        int delta = res_process + group_size;
        //printf("Process %d with delta %d and group size %d\n", me, delta, group_size );
        int ranks[group_size];
        for (int i = me, d = 0; i <= delta; i++, d++){
            ranks[d] = i-1;
            //printf("*****Process %d adding rank %d\n", me, i-1);
        }

        //printf("Process %d creating group\n ", me );
        MPI_Group world_group;
        MPI_Comm_group(MPI_COMM_WORLD, &world_group);

        MPI_Group prime_group;
        MPI_Group_incl(world_group, group_size, ranks, &prime_group);

        // Create a new communicator based on the group
        MPI_Comm prime_comm;
        MPI_Comm_create_group(MPI_COMM_WORLD, prime_group, 0, &prime_comm);

        //Receiving points
        node_t *for_ortho_points = malloc(sizeof(node_t) * interval);
        long *for_indexes = malloc(sizeof(long) * interval);
        //MPI_Recv(for_indexes, interval, MPI_LONG, MPI_ANY_SOURCE, POINT_TAG, WORLD, &statuses[1]);

        //printf("------------------Process: %d, receiving points\n", me);
        MPI_Scatter(for_indexes, interval, MPI_LONG, for_indexes, interval, MPI_LONG, 0, prime_comm);
        //printf("------------------Process: %d, received points\n", me);

        //printf("\nSCATETING SLAVE\n");
        for(int i = 0; i < interval; i++){
            for_ortho_points[i].center = malloc(sizeof(double) * n_dims);
            for_ortho_points[i].point_id = for_indexes[i];
            //printf("*** Slave process %d Received point %ld\n", me, for_indexes[i]);
        }
        //printf("\n");
        
        ortho_points = for_ortho_points;

        //printf("Process: %d, starting calculations\n", me);
        calc_projections_mpi(0, interval, threads, interval, 0, a, b, prime_comm);
        for_it++;

    }
    
    if (me < diff){
        wait_mpi(me, start, diff, threads);
    }

    else if (me > diff){
        wait_mpi(me, diff, end, threads);
    }
    
    //printf("Process: %d, receiving points\n", me);


    // Receive Args
    long args[2];
    //fprintf(stderr, "Node %d, is waiting for args\n", me);
    MPI_Recv(args, 2, MPI_LONG, MPI_ANY_SOURCE, ARGS_TAG, WORLD, &statuses[0]);

    node_t *new_ortho_points = malloc(sizeof(node_t) * args[0]);
    //fprintf(stderr, "Node %d, after malloc\n", me);

    // Receive point indexes
    long *rec_index = malloc(sizeof(long) * args[0]);
    //fprintf(stderr, "Node %d, is waiting for points\n", me);
    MPI_Recv(rec_index, args[0], MPI_LONG, MPI_ANY_SOURCE, POINT_TAG, WORLD, &statuses[1]);

    for(int i = 0; i < args[0]; i++){
        new_ortho_points[i].center = malloc(sizeof(double) * n_dims);
        new_ortho_points[i].point_id = rec_index[i];
        //printf("------/-/-/-/-/-/-Process: %d, received pomt %ld\n", me, rec_index[i]);
    }
    
    //fprintf(stderr, "Node %d received points\n", me);
    ortho_points = new_ortho_points;
    node_t * tree = build_tree_parallel_mpi(0, args[0], me, args[1], max_threads);

    // Send end confirmation
    int confirmation[1] = {1};
    MPI_Send(confirmation, 1, MPI_INT, 0, CONFIRMATION_TAG, WORLD);
    //printf("//////////////%d sent confirmation\n", me);

    if (max_size < 1000) {
        // Send node count
        long n_count = 0;
        count_nodes(tree, &n_count);
        long node_count[1] = {n_count};
        MPI_Send(node_count, 1, MPI_LONG, 0, COUNT_TAG, WORLD);

        // Receive node id and print type 
        int node_sending = 0;
        int proc_number = 0, aux = nprocs / 2;
        while (aux != 0) {
            proc_number += aux;
            if (proc_number == me) {
                node_sending = proc_number - aux;
                aux = 0;
            }
            aux /= 2;
        }
        //fprintf(stderr, "I am process %d waiting for %d to contact me\n", me, node_sending);
        long node_id[2];
        MPI_Recv(node_id, 2, MPI_LONG, node_sending, NODE_TAG, WORLD, &statuses[2]);
        //fprintf(stderr, "I am process %d | received a message from process %d\n", me, node_sending);

        long print_result[1];
        if (node_id[1]) {
            //fprintf(stderr, "Process %d starting print as main proc\n", me);
            current_print_proc = me + 1;
            print_result[0] = aux_print_tree_main_proc(tree, n_dims, points, n_count, node_id[0], me);
        } else {  
            print_result[0] = aux_print_tree(tree, n_dims, points, n_count, node_id[0]);
        }

        //fprintf(stderr, "Process %d return print to proc %d\n", me, node_sending);
        MPI_Send(print_result, 1, MPI_LONG, node_sending, PRINT_TAG, WORLD);
    }

    free(new_ortho_points);

    //fprintf(stderr, "Process %d is ENDING\n", me);
    MPI_Finalize();
    exit(0);
}

int main(int argc, char *argv[]) {
    double exec_time;
    long n_samples;

    omp_set_nested(1);
    omp_set_dynamic(1);
    max_threads = omp_get_max_threads();

    int me;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);

    exec_time = -omp_get_wtime();
    points = get_points(argc, argv, &n_dims, &n_samples);
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
    /*
     * Init some variables for MPI
     */

    recieved = malloc(sizeof(double) * n_samples );
    response = malloc(sizeof(double) * n_samples / nprocs);

    if (nprocs / 2 * 3 > n_samples){
        if (me == 0) {
            //printf("NO MPI \n");

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
    MPI_Status statuses[nprocs];
    for (int i = 1; i < nprocs; i++){
        MPI_Recv(confirmation, 1, MPI_LONG, MPI_ANY_SOURCE, CONFIRMATION_TAG, WORLD, &statuses[i]);
    }

    exec_time += omp_get_wtime();
    fprintf(stderr, "%lf\n", exec_time);

    if (n_samples < 1000) {
        // Receive node count
        int node_count = 0;
        int count[1];
        for (int i = 1; i < nprocs; i++){
            count[0] = 0;
            MPI_Recv(count, 1, MPI_LONG, MPI_ANY_SOURCE, COUNT_TAG, WORLD, &statuses[i]);
            node_count += count[0];
        }

        print_tree_mpi(tree, n_dims, points, node_count, me);
    }

    free(ortho_points);
    free_node(tree);
    
    free(recieved);
    free(response);

    free(points[0]);
    free(points);

    //fprintf(stderr, "Process %d is ENDING\n", me);
    MPI_Finalize();
    exit(0);
}
