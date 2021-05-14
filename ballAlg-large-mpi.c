// Large MPI version

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <omp.h>
#include <mpi.h>

#include "gen_points.h"
#include "quickSelect.h"
#include "ballAlg-large-mpi.h"

#define ELEM_SWAP(a,b) { register node_t t = (a); (a) = (b); (b) = t; }
#define WORLD MPI_COMM_WORLD
#define A_TAG               100
#define B_TAG               101
#define PROCESS_TAG         102
#define POINT_A_TAG         103
#define POINT_B_TAG         104
#define PROCESS_TYPE_TAG    105
#define REGULAR_SAMPLES_TAG 106
#define PIVOTS_TAG          107
#define CONFIRMATION_TAG    108
#define COUNT_TAG           109
#define NODE_TAG            110
#define PRINT_TAG           111
#define FOR_TAG             112
#define FOR_RESPONSE_TAG    113

typedef struct {
    int process;
    long max;
    double max_distance;
} furthest_point;

static double **points;
static node_t *ortho_points;
static int n_dims, max_threads, sent = 0;
static long max_size;
static int first = 0;
static int current_print_proc = 0;
static int nprocs = 0;
MPI_Datatype furthest_point_mpi;
MPI_Datatype ortho_point_mpi;

static node_t *create_node(double *point, long id, double radius) {
    node_t *node = malloc(sizeof(node_t));

    node->point_id = id;
    node->center = point;
    node->radius = radius;
    node->L = NULL;
    node->R = NULL;

    return node;
}

static void create_furthest_point_mpi(MPI_Datatype *data_type) {
    const int nitems = 1; 
    int blocklengths = 1;
    MPI_Datatype types[3] = { MPI_INT, MPI_LONG, MPI_DOUBLE }; 
    MPI_Aint offsets[3]; 

    offsets[0] = offsetof(furthest_point, process);
    offsets[1] = offsetof(furthest_point, max);
    offsets[2] = offsetof(furthest_point, max_distance);

    MPI_Type_create_struct(nitems, &blocklengths, offsets, types, data_type);
    MPI_Type_commit(data_type);
}

static void free_node(node_t *node) {
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

static double distance_sqrd(double *pt1, double *pt2) {
    double dist = 0.0;

    for (int d = 0; d < n_dims; d++) {
        double tmp = pt1[d] - pt2[d];
        dist += tmp * tmp;
    }

    return dist;
}

static double distance(double *pt1, double *median) {
    return sqrt(distance_sqrd(pt1, median));
}

static int get_furthest_point(long point, long start, long end) {
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

static double get_furthest_distance(double *point, long start, long end) {
    double distance, max_distance = 0.0;

    for (long i = start; i < end; i++) {
        distance = distance_sqrd(points[ortho_points[i].point_id], point);
        if ((max_distance - distance) < 0) {
            max_distance = distance;
        }
    }
    return max_distance;
}

static furthest_point get_furthest_point_parallel(double *point_point, long start, long end, int threads) {
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

    furthest_point max;
    double max_distance = 0.0;
    for (int i = 0; i < threads; i++) {
        if ((max_distance - furthest_points[i].max_distance) < 0) {
            max = furthest_points[i];
            max_distance = furthest_points[i].max_distance;
        }
    }
    return max;
}

static double get_furthest_distance_parallel(double *point, long start, long end, int threads) {
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

static void calc_ortho_projection(double *point_a, double *point_b, double *p1, double *p2, node_t *ortho_points, int p1_index, int p2_index) {
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
static node_t *build_tree(long start, long end) {  
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
static void calc_projections(long start, long end, int threads, double *point_a, double *point_b){

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
}

// not inclusive
static void calc_projections_mpi(long start, long end, int threads, long interval, long processor, long a, long b){
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

    //printf("Sending for response to processor %d \n", processor);
    MPI_Send(response, interval, MPI_DOUBLE, processor, FOR_RESPONSE_TAG, WORLD);
    //printf("Sent for response to processor %d \n", processor);
}

// not inclusive
static node_t *build_tree_parallel_omp(long start, long end, int threads, int me) {  
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
    long a = get_furthest_point_parallel(ortho_points[start].center, start, end, threads).max;
    long b = get_furthest_point_parallel(ortho_points[a].center, start, end, threads).max;

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

static int cmpfunc_ortho (const void *pa, const void *pb) {
  const node_t *a = (const node_t *) pa;
  const node_t *b = (const node_t *) pb;

  if ((a->center[0] - b->center[0]) > 0) return 1;
  else if ((a->center[0] - b->center[0]) < 0) return -1;
  else return 0;
}

static int cmpfunc_double (const void *pa, const void *pb) {
  const double *a = (const double *) pa;
  const double *b = (const double *) pb;

  if ((*a - *b) > 0) return 1;
  else if ((*a - *b) < 0) return -1;
  else return 0;
}

// not inclusive
static node_t *build_tree_parallel_mpi(double *first_point, long start, long end, int me, int max_processes, int threads) {  
    //fprintf(stderr, "FUNCTION MPI Process: %d, Max_Process: %d, Threads: %d, Start: %ld, End: %ld\n", process, max_processes, threads, start, end);    
    double *point_a, *point_b;
    int process, sender;

    /*
     * Get furthest point a
     */
    furthest_point a = get_furthest_point_parallel(first_point, start, end, threads);

    // Receive value a
    furthest_point value_a[max_processes - 1];
    for (int i = 1; i < max_processes; i++) {
        MPI_Recv(&value_a[i - 1], sizeof(furthest_point_mpi), furthest_point_mpi, MPI_ANY_SOURCE, A_TAG, WORLD, MPI_STATUS_IGNORE);
    }

    furthest_point max_furthest = a;
    for (int i = 0; i < max_processes - 1; i++) {
        if ((max_furthest.max_distance - value_a[i].max_distance) < 0) {
            max_furthest = value_a[i];
        }
    }

    for (int i = 1; i < max_processes; i++) {
        if (max_furthest.process == 0) sender = 0;
        else if (max_furthest.process == i) sender = 1;
        else sender = 0;

        MPI_Send(&sender, 1, MPI_INT, i, PROCESS_TYPE_TAG, WORLD);
    }

    if (max_furthest.process != 0) {
        // Receive a from other process
        point_a = malloc(sizeof(double) * n_dims);
        MPI_Recv(point_a, n_dims, MPI_DOUBLE, MPI_ANY_SOURCE, POINT_A_TAG, WORLD, MPI_STATUS_IGNORE); 
   
    } else {
        // Send a to everyone
        point_a = points[ortho_points[a.max].point_id];
        for (int i = 1; i < max_processes; i++) {
            MPI_Send(point_a, n_dims, MPI_DOUBLE, i, POINT_A_TAG, WORLD);
        }
    }


    /*
     * Get furthest point b
     */
    furthest_point b = get_furthest_point_parallel(point_a, start, end, threads);

    // Receive value b
    furthest_point value_b[max_processes - 1];
    for (int i = 1; i < max_processes; i++) {
        MPI_Recv(&value_b[i - 1], sizeof(furthest_point_mpi), furthest_point_mpi, MPI_ANY_SOURCE, B_TAG, WORLD, MPI_STATUS_IGNORE);
    }

    max_furthest = b;
    for (int i = 0; i < max_processes - 1; i++) {
        if ((max_furthest.max_distance - value_b[i].max_distance) < 0) {
            max_furthest = value_b[i];
        }
    }

    for (int i = 1; i < max_processes; i++) {
        if (max_furthest.process == 0) sender = 0;
        else if (max_furthest.process == i) sender = 1;
        else sender = 0;

        MPI_Send(&sender, 1, MPI_INT, i, PROCESS_TYPE_TAG, WORLD);
    }
    
    if (max_furthest.process != 0) {
        // Receive b from other process
        point_b = malloc(sizeof(double) * n_dims);
        MPI_Recv(point_b, n_dims, MPI_DOUBLE, MPI_ANY_SOURCE, POINT_B_TAG, WORLD, MPI_STATUS_IGNORE); 

    } else {
        // Send b to everyone
        point_b = points[ortho_points[b.max].point_id];
        for (int i = 1; i < max_processes; i++) {
            MPI_Send(point_b, n_dims, MPI_DOUBLE, i, POINT_B_TAG, WORLD);
        }
    }
    
    fprintf(stderr, "[PROCESS %d]: point a is", 0);
    for (int j = 0; j < n_dims; j++) {
        fprintf(stderr, " %f", point_a[j]);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "[PROCESS %d]: point b is", 0);
    for (int j = 0; j < n_dims; j++) {
        fprintf(stderr, " %f", point_b[j]);
    }
    fprintf(stderr, "\n");

    calc_projections(start, end, threads, point_a, point_b);

    for (int i = start; i < end; i++) {
        fprintf(stderr, "[PROCESS %d]: projection at ortho_points %d is", 0, i);
        for (int j = 0; j < n_dims; j++) {
            fprintf(stderr, " %f", ortho_points[i].center[j]);
        }
        fprintf(stderr, "\n");
    }
    
    /*
     * Get median point which will be the center of the ball
     */

    // Sort ortho points
    qsort(ortho_points, end, sizeof(node_t), cmpfunc_ortho);
    
    fprintf(stderr, "[PROCESS %d]: ortho_points sorted ", me);
    for (int i = 0; i < end; i++) {
        fprintf(stderr, " %f", ortho_points[i].center[0]);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "[PROCESS %d]: ortho_points point id ", me);
    for (int i = 0; i < end; i++) {
        fprintf(stderr, " %ld", ortho_points[i].point_id);
    }
    fprintf(stderr, "\n");

    double regular_samples[max_processes * max_processes];
    for (int i = 0; i < max_processes; i++) {
        regular_samples[i] = ortho_points[(end / nprocs) * i].center[0];
    }  

    for (int i = 1; i < max_processes; i++) {
        MPI_Recv(&regular_samples[i * max_processes], max_processes, MPI_DOUBLE, MPI_ANY_SOURCE, REGULAR_SAMPLES_TAG, WORLD, MPI_STATUS_IGNORE);
    }

    // Sort regular samples
    qsort(regular_samples, max_processes * max_processes, sizeof(double), cmpfunc_double);

    fprintf(stderr, "[PROCESS %d]: regular_samples sorted ", me);
    for (int i = 0; i < max_processes * max_processes; i++) {
        fprintf(stderr, " %f", regular_samples[i]);
    }
    fprintf(stderr, "\n");

    double pivots[max_processes - 1];
    double first_pivot = max_processes + (max_processes / 2) - 1;
    double last_pivot = (max_processes - 1) * max_processes + (max_processes / 2) - 1;
    for (int i = 0, pivot = first_pivot; i <= last_pivot; i += max_processes) {
        pivots[i++] = regular_samples[pivot];
    }

    //Send pivots
    for (int i = 1; i < max_processes; i++) {
        MPI_Send(pivots, max_processes - 1, MPI_DOUBLE, i, PIVOTS_TAG, WORLD);
    }

    fprintf(stderr, "[PROCESS %d]: pivots ", me);
    for (int i = 0; i < max_processes - 1; i++) {
        fprintf(stderr, " %f", pivots[i]);
    }
    fprintf(stderr, "\n");

    int send_counts[max_processes];
    send_counts[0] = 0; 
    for (int i = 0, pi = 0; i < end; i++) {
        fprintf(stderr, "PROCESS %d => Ortho Points[%d] %f and Pivots[%d] %f\n", me, i, ortho_points[i].center[0], pi, pivots[pi]);
        if ((ortho_points[i].center[0] - pivots[pi]) > 0) {
            pi++;
            send_counts[pi] = 0;
        }
        if (pi == max_processes - 1) {
            send_counts[max_processes - 1] = (end - i) * n_dims;
            break;
        }
        send_counts[pi] += n_dims; 
    }

    int send_displs[max_processes];
    int recv_counts[max_processes];
    int recv_displs[max_processes];

    MPI_Barrier(WORLD);
    MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, WORLD);
    
    fprintf(stderr, "PROCESS %d SEND COUNTS", me);
    for (int i = 0; i < max_processes; i++) {
        fprintf(stderr, " %d", send_counts[i]);
    }
    fprintf(stderr, "\n");

    send_displs[0] = 0;
    recv_displs[0] = 0;
    for (int i = 1; i < max_processes; i++) {
        send_displs[i] = send_counts[i - 1] + send_displs[i - 1];
        recv_displs[i] = recv_counts[i - 1] + recv_displs[i - 1];
    }

    double **points_to_send = (double **) create_array_pts(n_dims, end);
    for (int i = 0; i < end; i++) {
        points_to_send[i] = points[ortho_points[i].point_id];
    }
    
    for (int i = 0; i < end; i++) {
        fprintf(stderr, "[PROCESS %d]: points to be sent: ", me);
        for (int j = 0; j < n_dims; j++) {
            fprintf(stderr, " %f", points_to_send[i][j]);
        }
        fprintf(stderr, "\n");
    }

    MPI_Barrier(WORLD);
    MPI_Alltoallv(&(points_to_send[0][0]), send_counts, send_displs, MPI_DOUBLE, &(points[0][0]), recv_counts, recv_displs, MPI_DOUBLE, WORLD);

    //free(points_to_send[0]);
    //free(points_to_send);

    for (int i = 0; i < end; i++) {
        fprintf(stderr, "[PROCESS %d]: points after Alltoall: ", me);
        for (int j = 0; j < n_dims; j++) {
            fprintf(stderr, " %f", points[i][j]);
        }
        fprintf(stderr, "\n");
    }

    MPI_Finalize();
    exit(0);

    medianValues median_ids = quickSelect(ortho_points, start, end);

    // Send median to process 0
    // Wait for new set of points
    node_t point_median_1 = ortho_points[median_ids.first];
    node_t point_median_2 = ortho_points[median_ids.second];

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

    long size = end - median_ids.second;
    if ((size > 1000)) {
        tree->L = build_tree_parallel_mpi(first_point, start, median_ids.second, process, max_processes, threads);
        
        tree->L = build_tree_parallel_mpi(first_point, median_ids.second, end, process, max_processes, threads);

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


static void count_nodes(node_t *tree, long *n_count) {
    n_count[0]++;

    if (tree->L) {
        count_nodes(tree->L, n_count);
    }
    if (tree->R) {
        count_nodes(tree->R, n_count);
    }
}

static long aux_print_tree_main_proc(node_t *tree, int n_dims, double **points,
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
        MPI_Recv(count_received, 1, MPI_LONG, current_print_proc, PRINT_TAG, WORLD, MPI_STATUS_IGNORE);
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

static long aux_print_tree(node_t *tree, int n_dims, double **points,
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

static void print_tree_mpi(node_t *tree, int n_dims, double **points, int prev_count, int me) {
    long n_count = 0;

    count_nodes(tree, &n_count);
    n_count = n_count + prev_count;
    printf("%d %ld\n", n_dims, n_count);

    current_print_proc = 1;
    aux_print_tree_main_proc(tree, n_dims, points, n_count, 0, me);
}

static void print_tree(node_t *tree, int n_dims, double **points, int prev_count) {
    long n_count = 0;

    count_nodes(tree, &n_count);
    n_count = n_count + prev_count;
    printf("%d %ld\n", n_dims, n_count);

    current_print_proc = 1;
    aux_print_tree(tree, n_dims, points, n_count, 0);
}

static void wait_mpi(double *first_point, long start, long end, int me, int max_processes, int threads) {
    //fprintf(stderr, "FUNCTION MPI Process: %d, Max_Process: %d, Threads: %d, Start: %ld, End: %ld\n", process, max_processes, threads, start, end);    
    double *point_a, *point_b;
    int sender;

    /*
     * Get furthest point a
     */
    furthest_point a = get_furthest_point_parallel(first_point, start, end, threads);
    a.process = me;

    // Send value a
    MPI_Send(&a, sizeof(furthest_point_mpi), furthest_point_mpi, 0, A_TAG, WORLD);
    MPI_Recv(&sender, 1, MPI_INT, 0, PROCESS_TYPE_TAG, WORLD, MPI_STATUS_IGNORE);

    if (sender) {
        point_a = points[ortho_points[a.max].point_id];
        for (int i = 0; i < max_processes; i++) {
            if (i == me) continue;
            MPI_Send(point_a, n_dims, MPI_DOUBLE, i, POINT_A_TAG, WORLD);
        }
    } else {
        // Receive final a
        point_a = malloc(sizeof(double) * n_dims);
        MPI_Recv(point_a, n_dims, MPI_DOUBLE, MPI_ANY_SOURCE, POINT_A_TAG, WORLD, MPI_STATUS_IGNORE);
    }

    /*
     * Get furthest point b
     */
    furthest_point b = get_furthest_point_parallel(point_a, start, end, threads);
    b.process = me;

    // Send value b
    MPI_Send(&b, sizeof(furthest_point_mpi), furthest_point_mpi, 0, B_TAG, WORLD);
    MPI_Recv(&sender, 1, MPI_INT, 0, PROCESS_TYPE_TAG, WORLD, MPI_STATUS_IGNORE);

    if (sender) {
        point_b = points[ortho_points[b.max].point_id];
        for (int i = 0; i < max_processes; i++) {
            if (i == me) continue;
            MPI_Send(point_b, n_dims, MPI_DOUBLE, i, POINT_B_TAG, WORLD);
        }
    } else {
        // Receive final b
        point_b = malloc(sizeof(double) * n_dims);
        MPI_Recv(point_b, n_dims, MPI_DOUBLE, MPI_ANY_SOURCE, POINT_B_TAG, WORLD, MPI_STATUS_IGNORE);
    }

    fprintf(stderr, "[PROCESS %d]: point a is", me);
    for (int j = 0; j < n_dims; j++) {
        fprintf(stderr, " %f", point_a[j]);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "[PROCESS %d]: point b is", me);
    for (int j = 0; j < n_dims; j++) {
        fprintf(stderr, " %f", point_b[j]);
    }
    fprintf(stderr, "\n");

    /*
     * Calc projections
     */
    calc_projections(start, end, threads, point_a, point_b);

    for (int i = start; i < end; i++) {
        fprintf(stderr, "[PROCESS %d]: projection at ortho_points %d is", me, i);
        for (int j = 0; j < n_dims; j++) {
            fprintf(stderr, " %f", ortho_points[i].center[j]);
        }
        fprintf(stderr, "\n");
    }

    /*
     * Get median point which will be the center of the ball
     */

    // Sort ortho points
    qsort(ortho_points, end, sizeof(node_t), cmpfunc_ortho);
    
    fprintf(stderr, "[PROCESS %d]: ortho_points sorted ", me);
    for (int i = 0; i < end; i++) {
        fprintf(stderr, " %f", ortho_points[i].center[0]);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "[PROCESS %d]: ortho_points point id ", me);
    for (int i = 0; i < end; i++) {
        fprintf(stderr, " %ld", ortho_points[i].point_id);
    }
    fprintf(stderr, "\n");

    double regular_samples[max_processes];
    for (int i = 0; i < max_processes; i++) {
        regular_samples[i] = ortho_points[(end / nprocs) * i].center[0];
    }  

    MPI_Send(regular_samples, max_processes, MPI_DOUBLE, 0, REGULAR_SAMPLES_TAG, WORLD);

    //Receive pivots
    double pivots[max_processes - 1];
    MPI_Recv(pivots, max_processes - 1, MPI_DOUBLE, 0, PIVOTS_TAG, WORLD, MPI_STATUS_IGNORE);

    fprintf(stderr, "[PROCESS %d]: pivots ", me);
    for (int i = 0; i < max_processes - 1; i++) {
        fprintf(stderr, " %f", pivots[i]);
    }
    fprintf(stderr, "\n");

    int send_counts[max_processes];
    send_counts[0] = 0; 
    for (int i = 0, pi = 0; i < end; i++) {
        fprintf(stderr, "PROCESS %d => Ortho Points[%d] %f and Pivots[%d] %f\n", me, i, ortho_points[i].center[0], pi, pivots[pi]);
        if ((ortho_points[i].center[0] - pivots[pi]) > 0) {
            pi++;
            send_counts[pi] = 0;
        }
        if (pi == max_processes - 1) {
            send_counts[max_processes - 1] = (end - i) * n_dims;
            break;
        }
        send_counts[pi] += n_dims; 
    }

    int send_displs[max_processes];
    int recv_counts[max_processes];
    int recv_displs[max_processes];
    
    MPI_Barrier(WORLD);
    MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, WORLD);

    fprintf(stderr, "PROCESS %d SEND COUNTS", me);
    for (int i = 0; i < max_processes; i++) {
        fprintf(stderr, " %d", send_counts[i]);
    }
    fprintf(stderr, "\n");
    
    send_displs[0] = 0;
    recv_displs[0] = 0;
    for (int i = 1; i < max_processes; i++) {
        send_displs[i] = send_counts[i - 1] + send_displs[i - 1];
        recv_displs[i] = recv_counts[i - 1] + recv_displs[i - 1];
    }

    double **points_to_send = (double **) create_array_pts(n_dims, end);
    for (int i = 0; i < end; i++) {
        points_to_send[i] = points[ortho_points[i].point_id];
    }

    MPI_Barrier(WORLD);
    MPI_Alltoallv(&(points_to_send[0][0]), send_counts, send_displs, MPI_DOUBLE, &(points[0][0]), recv_counts, recv_displs, MPI_DOUBLE, WORLD);

    //free(points_to_send[0]);
    //free(points_to_send);

    for (int i = 0; i < end; i++) {
        fprintf(stderr, "[PROCESS %d]: points after Alltoall: ", me);
        for (int j = 0; j < n_dims; j++) {
            fprintf(stderr, " %f", points[i][j]);
        }
        fprintf(stderr, "\n");
    }

    MPI_Finalize();
    exit(0);

    // Send end confirmation
    int confirmation = 1;
    MPI_Send(&confirmation, 1, MPI_INT, 0, CONFIRMATION_TAG, WORLD);
    
    // Send node count
    long n_count = 0;
    //count_nodes(tree, &n_count);
    long node_count = n_count;
    MPI_Send(&node_count, 1, MPI_LONG, 0, COUNT_TAG, WORLD);

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
    MPI_Recv(node_id, 2, MPI_LONG, node_sending, NODE_TAG, WORLD, MPI_STATUS_IGNORE);
    //fprintf(stderr, "I am process %d | received a message from process %d\n", me, node_sending);

    long print_result;
    if (node_id[1]) {
        //fprintf(stderr, "Process %d starting print as main proc\n", me);
        current_print_proc = me + 1;
        //print_result = aux_print_tree_main_proc(tree, n_dims, points, n_count, node_id[0], me);
    } else {  
        //print_result = aux_print_tree(tree, n_dims, points, n_count, node_id[0]);
    }

    //fprintf(stderr, "Process %d return print to proc %d\n", me, node_sending);
    MPI_Send(&print_result, 1, MPI_LONG, node_sending, PRINT_TAG, WORLD);

    //free(new_ortho_points);

    //fprintf(stderr, "Process %d is ENDING\n", me);
    MPI_Finalize();
    exit(0);
}

int ballAlg_large_mpi(int argc, char *argv[], long n_samples) {
    double exec_time;
    double *first_point;
    long start, end;

    n_dims = atoi(argv[1]);

    omp_set_nested(1);
    omp_set_dynamic(1);
    max_threads = omp_get_max_threads();

    exec_time = -omp_get_wtime();

    int me;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);

    start = (n_samples / nprocs) * me;
    end = start + (n_samples / nprocs);

    first_point = malloc(sizeof(double) * n_dims);
    points = get_points_large(argc, argv, n_dims, n_samples, start, end, first_point);
    n_samples = end - start;

    fprintf(stderr, "[PROCESS %d] First points", me);
    for (int j = 0; j < n_dims; j++) {
        fprintf(stderr, " %f", first_point[j]);
    }
    fprintf(stderr, "\n");

    for (int i = 0; i < n_samples; i++) {
        fprintf(stderr, "[PROCESS %d] Points at %d is", me, i);
        for (int j = 0; j < n_dims; j++) {
            fprintf(stderr, " %f", points[i][j]);
        }
        fprintf(stderr, "\n");
    }

    create_furthest_point_mpi(&furthest_point_mpi);

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
        wait_mpi(first_point, 0, n_samples, me, nprocs, max_threads);

    tree = build_tree_parallel_mpi(first_point, 0, n_samples, 0, nprocs, max_threads);
    
    // Receive confirmation
    int confirmation;
    for (int i = 1; i < nprocs; i++){
        MPI_Recv(&confirmation, 1, MPI_LONG, MPI_ANY_SOURCE, CONFIRMATION_TAG, WORLD, MPI_STATUS_IGNORE);
    }

    exec_time += omp_get_wtime();

    // Receive node count
    int node_count = 0;
    int count;
    for (int i = 1; i < nprocs; i++){
        count = 0;
        MPI_Recv(&count, 1, MPI_LONG, MPI_ANY_SOURCE, COUNT_TAG, WORLD, MPI_STATUS_IGNORE);
        node_count += count;
    }

    fprintf(stderr, "%lf\n", exec_time);
    print_tree_mpi(tree, n_dims, points, node_count, me);

    free(ortho_points);
    free_node(tree);
    
    free(points[0]);
    free(points);

    //fprintf(stderr, "Process %d is ENDING\n", me);
    MPI_Finalize();
    exit(0);
}