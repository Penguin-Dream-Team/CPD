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
#define A_TAG                 100
#define B_TAG                 101
#define PROCESS_TAG           102
#define FIRST_POINT_TAG       103
#define POINT_A_TAG           104
#define POINT_B_TAG           105
#define PROCESS_TYPE_TAG      106
#define REGULAR_SAMPLES_TAG   107
#define REQUEST_POINTS_TAG    108
#define PIVOTS_TAG            109
#define MEDIAN_POINT_TAG      110
#define FURTHEST_DISTANCE_TAG 111
#define CONFIRMATION_TAG      112
#define COUNT_TAG             113
#define NODE_TAG              114
#define PRINT_TAG             115
#define PRINT_2_TAG           116
#define FOR_TAG               117
#define FOR_RESPONSE_TAG      118

typedef struct {
    int process;
    long max;
    double max_distance;
} furthest_point;

static double **points;
static node_t *ortho_points;
static int n_dims, max_threads, sent = 0;
static long initial_size;
static int first = 0;
static int current_print_proc = 0;
static int initial_procs, nprocs = 0;
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

static furthest_point get_furthest_point_parallel(int me, long index, double *point_point, long start, long end, int threads) {
    furthest_point *furthest_points = malloc((sizeof(furthest_point) + 2048) * threads);
    for (int i = 0; i < threads; i++) {
        furthest_points[i].max = index;
        furthest_points[i].max_distance = 0.0;
        furthest_points[i].process = me;
    }

    #pragma omp parallel num_threads(threads)
    {
        furthest_point fp;
        fp.max = index;
        fp.max_distance = 0.0;
        fp.process = me;
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

    furthest_point max = furthest_points[0];
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
    double *furthest_distances = malloc(sizeof(double) * threads);
    for (int i = 0; i < threads; i++) furthest_distances[i] = 0.0;

    /*fprintf(stderr, "Median Point ");
    for (int i = start; i < end; i++) fprintf(stderr, " %f", point[i]);
    fprintf(stderr, "\n");*/
    #pragma omp parallel num_threads(threads)
    {
        double fd = 0.0;
        #pragma omp for schedule(static) 
        for (long i = start; i < end; i++) {
            double distance = distance_sqrd(points[ortho_points[i].point_id], point);
            //fprintf(stderr, "Distance %f\n", distance);
            if ((fd - distance) < 0) {
                fd = distance;
            }
        }
        furthest_distances[omp_get_thread_num()] = fd;
        //fprintf(stderr, "Furthest_point %f\n", furthest_distances[omp_get_thread_num()]);
    }

    double max_distance = 0.0;
    for (int i = 0; i < threads; i++) {
        //fprintf(stderr, "Distance vector %f\n", furthest_distances[i]);
        if ((max_distance - furthest_distances[i]) < 0) {
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

static void calc_ortho_projection_mpi(double *point_a, double *point_b, double *p1, double *p2) {
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
        p1[i] = inner_product1 * (point_b[i] - point_a[i]) + point_a[i];
        p2[i] = inner_product2 * (point_b[i] - point_a[i]) + point_a[i];
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
        //if (i == start + 3) fprintf(stderr, "ITS MEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE => %f %f\n", points[ortho_points[i].point_id][0], points[ortho_points[i].point_id][1]);
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
    /*fprintf(stderr, "BUILD TREE BEGIN PROCESS %d\n", me);
    for (int i = start; i < end; i++) {
        fprintf(stderr, "PROCESS %d POINTS %d is", me, i);
        for (int j = 0; j < n_dims; j++) {
            fprintf(stderr, " %f", points[ortho_points[i].point_id][j]);
        }
        fprintf(stderr, "\n");
    }*/
    if (start == end - 1) {  // 1 point
        return create_node(points[ortho_points[start].point_id], ortho_points[start].point_id, 0);

    } else if (start == end - 2) {  // 2 points
        //fprintf(stderr, "I am process %d and im on OMP with start %ld and end %ld\n", me, start, end);
        double *median_point = malloc(sizeof(double) * n_dims);
        double *point1 = points[ortho_points[start].point_id];
        double *point2 = points[ortho_points[end - 1].point_id];

        ////fprintf(stderr, "I am process %d and im on OMP still alive\n", me);
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

        ////fprintf(stderr, "I am process %d and im returning my tree\n", me);
        return tree;
    }

    /*
     * Get furthest points
     */
    long a = get_furthest_point_parallel(0, start, points[ortho_points[start].point_id], start, end, threads).max;
    long b = get_furthest_point_parallel(0, a, points[ortho_points[a].point_id], start, end, threads).max;

    double *point_a = points[ortho_points[a].point_id];
    double *point_b = points[ortho_points[b].point_id];

    /*fprintf(stderr, "[PROCESS %d]: point a is", 0);
    for (int j = 0; j < n_dims; j++) {
        fprintf(stderr, " %f", point_a[j]);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "[PROCESS %d]: point b is", 0);
    for (int j = 0; j < n_dims; j++) {
        fprintf(stderr, " %f", point_b[j]);
    }
    fprintf(stderr, "\n");*/

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
static node_t *build_tree_parallel_mpi(long start, long end, int me, int max_processes, MPI_Comm communicator, int threads) {  
    //fprintf(stderr, "FUNCTION MPI Process: %d, Max_Process: %d, Threads: %d, Start: %ld, End: %ld\n", process, max_processes, threads, start, end);    
    double *point_a = malloc(sizeof(double) * n_dims);
    double *point_b = malloc(sizeof(double) * n_dims);
    int process, sender;

    //Send first point to everyone
    double *first_point = malloc(sizeof(double) * n_dims);
    for (int i = 0; i < n_dims; i++)        first_point[i] = points[0][i];
    for (int i = 1; i < max_processes; i++) MPI_Send(first_point, n_dims, MPI_DOUBLE, i, FIRST_POINT_TAG, communicator);

    /*
     * Get furthest point a
     */
    furthest_point a = get_furthest_point_parallel(me, 0, first_point, start, end, threads);

    // Receive value a
    furthest_point value_a[max_processes - 1];
    for (int i = 1; i < max_processes; i++) {
        MPI_Recv(&value_a[i - 1], sizeof(furthest_point_mpi), furthest_point_mpi, MPI_ANY_SOURCE, A_TAG, communicator, MPI_STATUS_IGNORE);
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

        MPI_Send(&sender, 1, MPI_INT, i, PROCESS_TYPE_TAG, communicator);
    }

    if (max_furthest.process != 0) {
        // Receive a from other process
        MPI_Recv(point_a, n_dims, MPI_DOUBLE, MPI_ANY_SOURCE, POINT_A_TAG, communicator, MPI_STATUS_IGNORE); 
   
    } else {
        // Send a to everyone
        for (int i = 0; i < n_dims; i++) point_a[i] = points[ortho_points[a.max].point_id][i];
        for (int i = 1; i < max_processes; i++) {
            MPI_Send(point_a, n_dims, MPI_DOUBLE, i, POINT_A_TAG, communicator);
        }
    }

    /*
     * Get furthest point b
     */
    furthest_point b = get_furthest_point_parallel(me, 0, point_a, start, end, threads);

    // Receive value b
    furthest_point value_b[max_processes - 1];
    for (int i = 1; i < max_processes; i++) {
        MPI_Recv(&value_b[i - 1], sizeof(furthest_point_mpi), furthest_point_mpi, MPI_ANY_SOURCE, B_TAG, communicator, MPI_STATUS_IGNORE);
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

        MPI_Send(&sender, 1, MPI_INT, i, PROCESS_TYPE_TAG, communicator);
    }
    
    if (max_furthest.process != 0) {
        // Receive b from other process
        MPI_Recv(point_b, n_dims, MPI_DOUBLE, MPI_ANY_SOURCE, POINT_B_TAG, communicator, MPI_STATUS_IGNORE); 

    } else {
        // Send b to everyone
        for (int i = 0; i < n_dims; i++) point_b[i] = points[ortho_points[b.max].point_id][i];
        for (int i = 1; i < max_processes; i++) {
            MPI_Send(point_b, n_dims, MPI_DOUBLE, i, POINT_B_TAG, communicator);
        }
    }
    
    /*fprintf(stderr, "[PROCESS %d]: point a is", rank_on_world);
    for (int j = 0; j < n_dims; j++) {
        fprintf(stderr, " %f", point_a[j]);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "[PROCESS %d]: point b is", rank_on_world);
    for (int j = 0; j < n_dims; j++) {
        fprintf(stderr, " %f", point_b[j]);
    }
    fprintf(stderr, "\n");*/

    calc_projections(start, end, threads, point_a, point_b);

    /*for (int i = start; i < end; i++) {
        fprintf(stderr, "[PROCESS %d]: projection at ortho_points %d is", 0, i);
        for (int j = 0; j < n_dims; j++) {
            fprintf(stderr, " %f", ortho_points[i].center[j]);
        }
        fprintf(stderr, "\n");
    }*/
    
    /*
     * Get median point which will be the center of the ball
     */

    // Sort ortho points
    qsort(ortho_points, end, sizeof(node_t), cmpfunc_ortho);
    
    /*fprintf(stderr, "[PROCESS %d]: ortho_points sorted ", me);
    for (int i = 0; i < end; i++) {
        fprintf(stderr, " %f", ortho_points[i].center[0]);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "[PROCESS %d]: ortho_points point id ", me);
    for (int i = 0; i < end; i++) {
        fprintf(stderr, " %ld", ortho_points[i].point_id);
    }
    fprintf(stderr, "\n");*/

    double *regular_samples = malloc(sizeof(double) * (max_processes * max_processes));
    for (int i = 0; i < max_processes; i++) {
        regular_samples[i] = ortho_points[(end / max_processes) * i].center[0];
    }  

    for (int i = 1; i < max_processes; i++) {
        MPI_Recv(&regular_samples[i * max_processes], max_processes, MPI_DOUBLE, MPI_ANY_SOURCE, REGULAR_SAMPLES_TAG, communicator, MPI_STATUS_IGNORE);
    }

    // Sort regular samples
    qsort(regular_samples, max_processes * max_processes, sizeof(double), cmpfunc_double);

    /*fprintf(stderr, "[PROCESS %d]: regular_samples sorted ", me);
    for (int i = 0; i < max_processes * max_processes; i++) {
        fprintf(stderr, " %f", regular_samples[i]);
    }
    fprintf(stderr, "\n");*/

    double *pivots = malloc(sizeof(double) * max_processes - 1);
    double first_pivot = max_processes + (max_processes / 2) - 1;
    double last_pivot = (max_processes - 1) * max_processes + (max_processes / 2) - 1;
    for (int i = 0, pivot = first_pivot; pivot <= last_pivot; pivot += max_processes) {
        pivots[i++] = regular_samples[pivot];
    }

    //Send pivots
    for (int i = 1; i < max_processes; i++) {
        MPI_Send(pivots, max_processes - 1, MPI_DOUBLE, i, PIVOTS_TAG, communicator);
    }

    /*fprintf(stderr, "[PROCESS %d]: pivots ", me);
    for (int i = 0; i < max_processes - 1; i++) {
        fprintf(stderr, " %f", pivots[i]);
    }
    fprintf(stderr, "\n");*/

    int *send_counts = malloc(sizeof(int) * max_processes);
    send_counts[0] = 0; 
    for (int i = 0, pi = 0; i < end;) {
        //fprintf(stderr, "PROCESS %d => Ortho Points[%d] %f and Pivots[%d] %f\n", me, i, ortho_points[i].center[0], pi, pivots[pi]);
        if ((ortho_points[i].center[0] - pivots[pi]) > 0) {
            pi++;
            send_counts[pi] = 0;
        }
        if (pi == max_processes - 1) {
            send_counts[max_processes - 1] = end - i;
            break;
        }
        if ((ortho_points[i].center[0] - pivots[pi]) <= 0) {
            i++;
            send_counts[pi]++; 
        }
    }

    int *send_displs = malloc(sizeof(int) * max_processes);
    int *recv_counts = malloc(sizeof(int) * max_processes);
    int *recv_displs = malloc(sizeof(int) * max_processes);
    
    MPI_Barrier(communicator);
    MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, communicator);
    
    int recv_count = 0;
    for (int i = 0; i < max_processes; i++) {
        recv_count += recv_counts[i];
    }

    /*fprintf(stderr, "PROCESS %d SEND COUNTS", me);
    for (int i = 0; i < max_processes; i++) {
        fprintf(stderr, " %d", send_counts[i]);
    }
    fprintf(stderr, "\n");*/
    
    send_displs[0] = 0;
    recv_displs[0] = 0;
    for (int i = 1; i < max_processes; i++) {
        send_displs[i] = send_counts[i - 1] + send_displs[i - 1];
        recv_displs[i] = recv_counts[i - 1] + recv_displs[i - 1];
    }

    /*for (int i = 0; i < max_processes; i++) {
        fprintf(stderr, "PROCESS %d PARTITION %d =", me, i);
        for (int j = 0; j < send_counts[i]; j++) {
            if (i == 0) fprintf(stderr, " %f", ortho_points[j].center[0]);
            else        fprintf(stderr, " %f", ortho_points[send_displs[i] + j].center[0]);
        }
        if (send_counts[i] == 0) fprintf(stderr, " Nothing");
        fprintf(stderr, "\n");
    }

    fprintf(stderr, "PROCESS %d RECV COUNTS", me);
    for (int i = 0; i < max_processes; i++) {
        fprintf(stderr, " %d", recv_counts[i]);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "PROCESS %d SEND DISPLS", me);
    for (int i = 0; i < max_processes; i++) {
        fprintf(stderr, " %d", send_displs[i]);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "PROCESS %d RECV DISPLS", me);
    for (int i = 0; i < max_processes; i++) {
        fprintf(stderr, " %d", recv_displs[i]);
    }
    fprintf(stderr, "\n");*/

    double *points_to_send = malloc(sizeof(double)* end);
    double *points_to_recv = malloc(sizeof(double)* recv_count);
    double **points_received = malloc(sizeof(double)* recv_count);
    for (int i = 0; i < recv_count; i++) {
        points_received[i] = malloc(sizeof(double) * n_dims);
    }
    
    for (int i = 0; i < n_dims; i++) {
        for (int j = 0; j < end; j++) {
            memcpy(&(points_to_send[j]), &(points[ortho_points[j].point_id][i]), sizeof(double));
        }

        MPI_Barrier(communicator);
        MPI_Alltoallv(&(points_to_send[0]), send_counts, send_displs, MPI_DOUBLE, &(points_to_recv[0]), recv_counts, recv_displs, MPI_DOUBLE, communicator);
        MPI_Barrier(communicator);

        for (int j = 0; j < recv_count; j++) {
            memcpy(&(points_received[j][i]), &(points_to_recv[j]), sizeof(double));
        }
    }

    for (int i = 0; i < end; i++) {
        memcpy(&(points_to_send[i]), &(ortho_points[i].center[0]), sizeof(double));
    }

    MPI_Barrier(communicator);
    MPI_Alltoallv(&(points_to_send[0]), send_counts, send_displs, MPI_DOUBLE, &(points_to_recv[0]), recv_counts, recv_displs, MPI_DOUBLE, communicator);
    MPI_Barrier(communicator);

    node_t *ortho_points_to_recv = malloc(sizeof(node_t)* recv_count);
    for (int i = 0; i < recv_count; i++) {
        ortho_points_to_recv[i].center = malloc(sizeof(double) * n_dims);
        memcpy(&(ortho_points_to_recv[i].center[0]), &(points_to_recv[i]), sizeof(double));
        ortho_points_to_recv[i].point_id = i;
    }

    // Sorting points after exchange
    qsort(ortho_points_to_recv, recv_count, sizeof(node_t), cmpfunc_ortho);

    //free(points_to_send[0]);
    //free(points_to_send);

    /*for (int i = 0; i < recv_count; i++) {
        fprintf(stderr, "[PROCESS %d]: points received after Alltoall: ", rank_on_world);
        for (int j = 0; j < n_dims; j++) {
            fprintf(stderr, " %f", points_received[ortho_points_to_recv[i].point_id][j]);
        }
        fprintf(stderr, "\n");
    }

    fprintf(stderr, "[PROCESS %d]: ortho points after Alltoall: ", rank_on_world);
    for (int i = 0; i < recv_count; i++) {
        fprintf(stderr, " %f", ortho_points_to_recv[i].center[0]);
    }
    fprintf(stderr, "\n");*/

    // Equally divide the ammount of points to all processes
    int requestPoints[2] = { 0, 0 };
    if (recv_count < end) {
        int request_size = end - recv_count;
        //fprintf(stderr, "I AM PROCESS %d and I need %d points from %d\n", me, request_size, me + 1);
        requestPoints[0] = 1;
        requestPoints[1] = request_size;

        MPI_Send(requestPoints, 2, MPI_INT, me + 1, REQUEST_POINTS_TAG, communicator);
        double *points_from_next_process = malloc(sizeof(double) * n_dims);
        for (int i = 0; i < request_size; i++) {
            MPI_Recv(&points_from_next_process[0], n_dims, MPI_DOUBLE, me + 1, REQUEST_POINTS_TAG, communicator, MPI_STATUS_IGNORE);
            for (int j = 0; j < n_dims; j++) {
                points[i + recv_count][j] = points_from_next_process[j];
            }
        }

        // Merge all to points
        //fprintf(stderr, "I AM PROCESS %d and I am merging my points\n", me);
        for (int i = 0; i < recv_count; i++) {
            for (int j = 0; j < n_dims; j++) {
                points[i][j] = points_received[ortho_points_to_recv[i].point_id][j];
            }
        }

    } else if (recv_count > end) {
        int request_size = recv_count - end;
        //fprintf(stderr, "I AM PROCESS %d and I am sending %d points to %d\n", me, request_size, me + 1);
        requestPoints[0] = 2;
        requestPoints[1] = request_size;
        MPI_Send(requestPoints, 2, MPI_INT, me + 1, REQUEST_POINTS_TAG, communicator);

        double **points_for_next_process = malloc(sizeof(double) * requestPoints[1]);
        for (int i = 0; i < requestPoints[1]; i++) {
            points_for_next_process[i] = malloc(sizeof(double) * n_dims);
            for (int j = 0; j < n_dims; j++) { 
                points_for_next_process[i][j] = points_received[ortho_points_to_recv[recv_count - requestPoints[1] + i].point_id][j];
            }
            MPI_Send(points_for_next_process[i], n_dims, MPI_DOUBLE, me + 1, REQUEST_POINTS_TAG, communicator);  
        }

        // Merge all to points
        //fprintf(stderr, "I AM PROCESS %d and I am merging my points\n", me);
        for (int i = 0; i < end; i++) {
            for (int j = 0; j < n_dims; j++) {
                points[i][j] = points_received[ortho_points_to_recv[i].point_id][j];
            }
        }

    } else {
        requestPoints[0] = 3;
        requestPoints[1] = 0;
        MPI_Send(requestPoints, 2, MPI_INT, me + 1, REQUEST_POINTS_TAG, communicator);

        // Merge all to points
        //fprintf(stderr, "I AM PROCESS %d and I am merging my points\n", me);
        for (int i = 0; i < end; i++) {
            for (int j = 0; j < n_dims; j++) {
                points[i][j] = points_received[ortho_points_to_recv[i].point_id][j];
            }
        }
    }

    // Sorting and exchange finished
    //fprintf(stderr, "I AM PROCESS %d and I finished my merge\n", me);
    for (int i = 0; i < end; i++) {
        ortho_points[i].point_id = i;
    }

    /*for (int i = 0; i < end; i++) {
        fprintf(stderr, "[PROCESS %d]: points after Alltoall: ", me);
        for (int j = 0; j < n_dims; j++) {
            fprintf(stderr, " %f", points[i][j]);
        }
        fprintf(stderr, "\n");
    }*/

    free(points_to_send);
    free(points_to_recv);
    free(ortho_points_to_recv);
    free(pivots);
    free(send_counts);
    free(recv_counts);
    free(send_displs);
    free(recv_displs);
    
    double *p1 = malloc(sizeof(double) * n_dims);
    double *p2 = malloc(sizeof(double) * n_dims);
    //fprintf(stderr, "[PROCESS %d]: Receiving median points\n", me);
    if (initial_size % 2 == 0) {
        if (max_processes > 2) {
            MPI_Recv(p1, n_dims, MPI_DOUBLE, (max_processes / 2) - 1, MEDIAN_POINT_TAG, communicator, MPI_STATUS_IGNORE);
            MPI_Recv(p2, n_dims, MPI_DOUBLE, (max_processes / 2)    , MEDIAN_POINT_TAG, communicator, MPI_STATUS_IGNORE);
        } else {
            for (int i = 0; i < n_dims; i++) {
                p1[i] = points[end - 1][i];
            }
            MPI_Recv(p2, n_dims, MPI_DOUBLE, me + 1, MEDIAN_POINT_TAG, communicator, MPI_STATUS_IGNORE);
        }
    } 
    
    /*
     * Separate L and R sets
     */
    double *median_point = malloc(sizeof(double) * n_dims);
    node_t *tree = create_node(median_point, -1, -1);
    
    /* Calc ortho projection of median points */
    /*fprintf(stderr, "[PROCESS %d]: Starting ortho projections\n", me);
    fprintf(stderr, "[PROCESS %d]: median 1: ", me);
    for (int j = 0; j < n_dims; j++) {
        fprintf(stderr, " %f", p1[j]);
    }
    fprintf(stderr, "\n");
    fprintf(stderr, "[PROCESS %d]: median 2: ", me);
    for (int j = 0; j < n_dims; j++) {
        fprintf(stderr, " %f", p2[j]);
    }
    fprintf(stderr, "\n");*/

    calc_ortho_projection_mpi(point_a, point_b, p1, p2);

    /*
     * Get the radius of the ball (largest distance)
     */
    if ((end - start) % 2 != 0) {
        for (int d = 0; d < n_dims; d++) {
            median_point[d] = p1[d];
        }
    } else {
        for (int d = 0; d < n_dims; d++) {
            median_point[d] = (p1[d] + p2[d]) / 2;
        }
    }

    /*fprintf(stderr, "Sending median ");
    for (int j = 0; j < n_dims; j++) {
        fprintf(stderr, " %f", median_point[j]);
    }
    fprintf(stderr, "\n");*/
    // Send median to other processes
    for (int i = 1; i < max_processes; i++) {
        MPI_Send(median_point, n_dims, MPI_DOUBLE, i, MEDIAN_POINT_TAG, communicator);
    }

    //fprintf(stderr, "Receiving distances\n");
    // Receive other process distances
    double furthest_distance = get_furthest_distance_parallel(median_point, start, end, threads);
    double *others_furthest_distance = malloc(sizeof(double) * (max_processes - 1));
    for (int i = 1; i < max_processes; i++) {
        MPI_Recv(&others_furthest_distance[i - 1], 1, MPI_DOUBLE, i, FURTHEST_DISTANCE_TAG, communicator, MPI_STATUS_IGNORE);
    }

    for (int i = 1; i < max_processes; i++) {
        //fprintf(stderr, "I am process %d and my furthest distance is %f vs the process %d furthest distance %f\n", me, furthest_distance, i, others_furthest_distance[i - 1]);
        if ((furthest_distance - others_furthest_distance[i - 1]) < 0) {
            furthest_distance = others_furthest_distance[i - 1];
        }
    }

    tree->radius = sqrt(furthest_distance);
    //fprintf(stderr, "Tree radius is %f\n", tree->radius);

    /*
     * Split processes in smaller groups or start serial version
     */
    //points = points_to_recv;

    /*for (int i = 0; i < end; i++) {
        fprintf(stderr, "Process %d point %d =", rank_on_world, i);
        for (int j = 0; j < n_dims; j++) {
            fprintf(stderr, " %f", points[i][j]);
        }
        fprintf(stderr, "\n");
    }*/

    if (max_processes > 2) {
        int half_max_processes = max_processes / 2;
        int ranks[half_max_processes];
        if (me < half_max_processes) for (int i = 0; i < half_max_processes; i++) ranks[i] = i;
        else                         for (int i = half_max_processes; i < max_processes; i++) ranks[i - half_max_processes] = i;
        
        MPI_Group world_group;
        MPI_Comm_group(WORLD, &world_group);

        MPI_Group new_group;
        MPI_Group_incl(world_group, half_max_processes, ranks, &new_group);

        // Create a new communicator based on the new group
        MPI_Comm new_communicator;
        MPI_Comm_create_group(MPI_COMM_WORLD, new_group, 0, &new_communicator);

        MPI_Comm_rank(new_communicator, &me);
    
        tree->L = build_tree_parallel_mpi(start, end, me, half_max_processes, new_communicator, threads);

    } else {
        #pragma omp parallel
        {
            #pragma omp single 
            {  
                if (threads > 2) {
                    #pragma omp task
                    tree->L = build_tree_parallel_omp(start, end, threads, process);
                } else {
                    #pragma omp task
                    tree->L = build_tree(start, end);
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
    //fprintf(stderr, "current process 0 count %ld\n", n_count[0]);
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
        int process_id = nprocs;
        while (process_id > 2) {
            process_id /= 2;
            if (process_id + me == current_print_proc) {
                print_type = 1;
                break;
            }
        }

        long print[] = {right_id, print_type, me};
        //fprintf(stderr, "[PROCESS %d]: sendind print order to %d\n", me, current_print_proc);
        MPI_Send(print, 3, MPI_LONG, current_print_proc, NODE_TAG, WORLD);

        long count_received;
        MPI_Recv(&count_received, 1, MPI_LONG, current_print_proc, PRINT_TAG, WORLD, MPI_STATUS_IGNORE);
        count = count_received;

        if (print_type == 1) {
            int next_print_proc;
            MPI_Recv(&next_print_proc, 1, MPI_INT, current_print_proc, PRINT_2_TAG, WORLD, MPI_STATUS_IGNORE);
            current_print_proc = next_print_proc - 1;
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
    fflush(stdout);
    MPI_Barrier(WORLD);

    current_print_proc = 1;
    //fprintf(stderr, "starting aux print - process %d with count = %ld\n", me, n_count);
    aux_print_tree_main_proc(tree, n_dims, points, n_count, 0, me);
}

static void print_tree(node_t *tree, int n_dims, double **points, int prev_count) {
    long n_count = 0;

    count_nodes(tree, &n_count);
    n_count = n_count + prev_count;
    printf("%d %ld\n", n_dims, n_count);

    aux_print_tree(tree, n_dims, points, n_count, 0);
}

static void wait_mpi(long start, long end, int me, int max_processes, MPI_Comm communicator, int threads) {
    //fprintf(stderr, "FUNCTION MPI Process: %d, Max_Process: %d, Threads: %d, Start: %ld, End: %ld\n", process, max_processes, threads, start, end);    
    double *point_a = malloc(sizeof(double) * n_dims);
    double *point_b = malloc(sizeof(double) * n_dims);
    int sender;

    //Receive first point
    double *first_point = malloc(sizeof(double) * n_dims);
    MPI_Recv(&first_point[0], n_dims, MPI_DOUBLE, MPI_ANY_SOURCE, FIRST_POINT_TAG, communicator, MPI_STATUS_IGNORE);

    /*
     * Get furthest point a
     */
    furthest_point a = get_furthest_point_parallel(me, 0, first_point, start, end, threads);

    // Send value a
    MPI_Send(&a, sizeof(furthest_point_mpi), furthest_point_mpi, 0, A_TAG, communicator);
    MPI_Recv(&sender, 1, MPI_INT, 0, PROCESS_TYPE_TAG, communicator, MPI_STATUS_IGNORE);

    if (sender) {
        for (int i = 0; i < n_dims; i++) point_a[i] = points[ortho_points[a.max].point_id][i];
        for (int i = 0; i < max_processes; i++) {
            if (i == me) continue;
            MPI_Send(point_a, n_dims, MPI_DOUBLE, i, POINT_A_TAG, communicator);
        }
    } else {
        // Receive final a
        MPI_Recv(point_a, n_dims, MPI_DOUBLE, MPI_ANY_SOURCE, POINT_A_TAG, communicator, MPI_STATUS_IGNORE);
    }

    /*
     * Get furthest point b
     */
    furthest_point b = get_furthest_point_parallel(me, 0, point_a, start, end, threads);

    // Send value b
    MPI_Send(&b, sizeof(furthest_point_mpi), furthest_point_mpi, 0, B_TAG, communicator);
    MPI_Recv(&sender, 1, MPI_INT, 0, PROCESS_TYPE_TAG, communicator, MPI_STATUS_IGNORE);

    if (sender) {
        for (int i = 0; i < n_dims; i++) point_b[i] = points[ortho_points[b.max].point_id][i];
        for (int i = 0; i < max_processes; i++) {
            if (i == me) continue;
            MPI_Send(point_b, n_dims, MPI_DOUBLE, i, POINT_B_TAG, communicator);
        }
    } else {
        // Receive final b
        MPI_Recv(point_b, n_dims, MPI_DOUBLE, MPI_ANY_SOURCE, POINT_B_TAG, communicator, MPI_STATUS_IGNORE);
    }

    /*fprintf(stderr, "[PROCESS %d]: point a is", me);
    for (int j = 0; j < n_dims; j++) {
        fprintf(stderr, " %f", point_a[j]);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "[PROCESS %d]: point b is", me);
    for (int j = 0; j < n_dims; j++) {
        fprintf(stderr, " %f", point_b[j]);
    }
    fprintf(stderr, "\n");*/

    /*
     * Calc projections
     */
    calc_projections(start, end, threads, point_a, point_b);

    /*for (int i = start; i < end; i++) {
        fprintf(stderr, "[PROCESS %d]: projection at ortho_points %d is", me, i);
        for (int j = 0; j < n_dims; j++) {
            fprintf(stderr, " %f", ortho_points[i].center[j]);
        }
        fprintf(stderr, "\n");
    }*/

    /*
     * Get median point which will be the center of the ball
     */

    // Sort ortho points
    qsort(ortho_points, end, sizeof(node_t), cmpfunc_ortho);
    
    /*fprintf(stderr, "[PROCESS %d]: ortho_points sorted ", me);
    for (int i = 0; i < end; i++) {
        fprintf(stderr, " %f", ortho_points[i].center[0]);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "[PROCESS %d]: ortho_points point id ", me);
    for (int i = 0; i < end; i++) {
        fprintf(stderr, " %ld", ortho_points[i].point_id);
    }
    fprintf(stderr, "\n");*/

    double *regular_samples = malloc(sizeof(double) * max_processes);
    for (int i = 0; i < max_processes; i++) {
        regular_samples[i] = ortho_points[(end / max_processes) * i].center[0];
    }  

    MPI_Send(regular_samples, max_processes, MPI_DOUBLE, 0, REGULAR_SAMPLES_TAG, communicator);

    //Receive pivots
    double *pivots = malloc(sizeof(double) * max_processes - 1);
    MPI_Recv(pivots, max_processes - 1, MPI_DOUBLE, 0, PIVOTS_TAG, communicator, MPI_STATUS_IGNORE);

    /*fprintf(stderr, "[PROCESS %d]: pivots ", me);
    for (int i = 0; i < max_processes - 1; i++) {
        fprintf(stderr, " %f", pivots[i]);
    }
    fprintf(stderr, "\n");*/

    int *send_counts = malloc(sizeof(int) * max_processes);
    send_counts[0] = 0; 
    for (int i = 0, pi = 0; i < end;) {
        //fprintf(stderr, "PROCESS %d => Ortho Points[%d] %f and Pivots[%d] %f\n", me, i, ortho_points[i].center[0], pi, pivots[pi]);
        if ((ortho_points[i].center[0] - pivots[pi]) > 0) {
            pi++;
            send_counts[pi] = 0;
        }
        if (pi == max_processes - 1) {
            send_counts[max_processes - 1] = end - i;
            break;
        }
        if ((ortho_points[i].center[0] - pivots[pi]) <= 0) {
            i++;
            send_counts[pi]++; 
        }
    }

    int *send_displs = malloc(sizeof(int) * max_processes);
    int *recv_counts = malloc(sizeof(int) * max_processes);
    int *recv_displs = malloc(sizeof(int) * max_processes);
    
    MPI_Barrier(communicator);
    MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, communicator);

    int recv_count = 0;
    for (int i = 0; i < max_processes; i++) {
        recv_count += recv_counts[i];
    }

    /*fprintf(stderr, "PROCESS %d SEND COUNTS", me);
    for (int i = 0; i < max_processes; i++) {
        fprintf(stderr, " %d", send_counts[i]);
    }
    fprintf(stderr, "\n");*/
    
    send_displs[0] = 0;
    recv_displs[0] = 0;
    for (int i = 1; i < max_processes; i++) {
        send_displs[i] = send_counts[i - 1] + send_displs[i - 1];
        recv_displs[i] = recv_counts[i - 1] + recv_displs[i - 1];
    }

    /*for (int i = 0; i < max_processes; i++) {
        fprintf(stderr, "PROCESS %d PARTITION %d =", me, i);
        for (int j = 0; j < send_counts[i]; j++) {
            if (i == 0) fprintf(stderr, " %f", ortho_points[j].center[0]);
            else        fprintf(stderr, " %f", ortho_points[send_displs[i] + j].center[0]);
        }
        if (send_counts[i] == 0) fprintf(stderr, " Nothing");
        fprintf(stderr, "\n");
    }

    fprintf(stderr, "PROCESS %d RECV COUNTS", me);
    for (int i = 0; i < max_processes; i++) {
        fprintf(stderr, " %d", recv_counts[i]);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "PROCESS %d SEND DISPLS", me);
    for (int i = 0; i < max_processes; i++) {
        fprintf(stderr, " %d", send_displs[i]);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "PROCESS %d RECV DISPLS", me);
    for (int i = 0; i < max_processes; i++) {
        fprintf(stderr, " %d", recv_displs[i]);
    }
    fprintf(stderr, "\n");*/

    double *points_to_send = malloc(sizeof(double) * end);
    double *points_to_recv = malloc(sizeof(double) * recv_count);
    double **points_received = malloc(sizeof(double) * recv_count);
    for (int i = 0; i < recv_count; i++) {
        points_received[i] = malloc(sizeof(double) * n_dims);
    }

    for (int i = 0; i < n_dims; i++) {
        for (int j = 0; j < end; j++) {
            memcpy(&(points_to_send[j]), &(points[ortho_points[j].point_id][i]), sizeof(double));
        }

        MPI_Barrier(communicator);
        MPI_Alltoallv(&(points_to_send[0]), send_counts, send_displs, MPI_DOUBLE, &(points_to_recv[0]), recv_counts, recv_displs, MPI_DOUBLE, communicator);
        MPI_Barrier(communicator);

        for (int j = 0; j < recv_count; j++) {
            memcpy(&(points_received[j][i]), &(points_to_recv[j]), sizeof(double));
        }
    }

    for (int i = 0; i < end; i++) {
        memcpy(&(points_to_send[i]), &(ortho_points[i].center[0]), sizeof(double));
    }

    MPI_Barrier(communicator);
    MPI_Alltoallv(&(points_to_send[0]), send_counts, send_displs, MPI_DOUBLE, &(points_to_recv[0]), recv_counts, recv_displs, MPI_DOUBLE, communicator);
    MPI_Barrier(communicator);

    node_t *ortho_points_to_recv = malloc(sizeof(node_t)* recv_count);
    for (int i = 0; i < recv_count; i++) {
        ortho_points_to_recv[i].center = malloc(sizeof(double) * n_dims);
        memcpy(&(ortho_points_to_recv[i].center[0]), &(points_to_recv[i]), sizeof(double));
        ortho_points_to_recv[i].point_id = i;
    }

    // Sorting points after exchange
    qsort(ortho_points_to_recv, recv_count, sizeof(node_t), cmpfunc_ortho);

    //free(points_to_send[0]);
    //free(points_to_send);

    /*for (int i = 0; i < recv_count; i++) {
        fprintf(stderr, "[PROCESS %d]: points received after Alltoall: ", rank_on_world);
        for (int j = 0; j < n_dims; j++) {
            fprintf(stderr, " %f", points_received[ortho_points_to_recv[i].point_id][j]);
        }
        fprintf(stderr, "\n");
    }

    fprintf(stderr, "[PROCESS %d]: ortho points after Alltoall: ", rank_on_world);
    for (int i = 0; i < recv_count; i++) {
        fprintf(stderr, " %f", ortho_points_to_recv[i].center[0]);
    }
    fprintf(stderr, "\n");*/
    
    // Equally divide the ammount of points to all processes
    int requestPoints[2] = { 0, 0 };
    double **points_from_previous_process;
    MPI_Recv(&requestPoints, 2, MPI_INT, me - 1, REQUEST_POINTS_TAG, communicator, MPI_STATUS_IGNORE);
    if (requestPoints[0] == 1) {
        //fprintf(stderr, "I AM PROCESS %d and process %d is asking me %d points\n", me, me - 1, requestPoints[1]);
        double *points_for_previous_process = malloc(sizeof(double) * n_dims);
        for (int i = 0; i < requestPoints[1]; i++) {
            for (int j = 0; j < n_dims; j++) {
                points_for_previous_process[j] = points_received[ortho_points_to_recv[i].point_id][j];
            }
            MPI_Send(points_for_previous_process, n_dims, MPI_DOUBLE, me - 1, REQUEST_POINTS_TAG, communicator);
        }           

    } else if (requestPoints[0] == 2) {
        //fprintf(stderr, "I AM PROCESS %d and process %d is sending me %d points\n", me, me - 1, requestPoints[1]);
        points_from_previous_process = malloc(sizeof(double) * requestPoints[1]);
        for (int i = 0; i < requestPoints[1]; i++) {
            points_from_previous_process[i] = malloc(sizeof(double) * n_dims);
            MPI_Recv(&(points_from_previous_process[i][0]), n_dims, MPI_DOUBLE, me - 1, REQUEST_POINTS_TAG, communicator, MPI_STATUS_IGNORE);
        }
    } 

    if (me != max_processes - 1) {
        int requestPointsNextProcesses[2] = { 0, 0 };
        int request_size_with_previous_process = 0;
        if (requestPoints[0] == 1)      request_size_with_previous_process = recv_count - requestPoints[1];
        else if (requestPoints[0] == 2) request_size_with_previous_process = recv_count + requestPoints[1];

        if (request_size_with_previous_process < end) {
            int request_size = end - request_size_with_previous_process;
            requestPointsNextProcesses[0] = 1;
            requestPointsNextProcesses[1] = request_size;

            MPI_Send(requestPointsNextProcesses, 2, MPI_INT, me + 1, REQUEST_POINTS_TAG, communicator);
            double **points_from_next_process = malloc(sizeof(double) * request_size);
            for (int i = 0; i < request_size; i++) {
                points_from_next_process[i] = malloc(sizeof(double) * n_dims);
                MPI_Recv(&(points_from_next_process[i][0]), n_dims, MPI_DOUBLE, me + 1, REQUEST_POINTS_TAG, communicator, MPI_STATUS_IGNORE);
            }

            // Merge all to points
            //fprintf(stderr, "I AM PROCESS %d and I am merging my points\n", me);
            if (requestPoints[0] == 1) {
                //fprintf(stderr, "I AM PROCESS %d and I am merging my points - I SENT AND RECEIVED\n", me);
                for (int i = 0; i < request_size_with_previous_process; i++) {
                    for (int j = 0; j < n_dims; j++) {
                        points[i][j] = points_received[ortho_points_to_recv[i + requestPoints[1]].point_id][j];
                    }
                }

                for (int i = request_size_with_previous_process; i < end; i++) {
                    for (int j = 0; j < n_dims; j++) {
                        points[i][j] = points_from_next_process[i - request_size_with_previous_process][j];
                    }
                }

            } else if (requestPoints[0] == 2) {
                //fprintf(stderr, "I AM PROCESS %d and I am merging my points - I RECEIVED AND RECEIVED\n", me);
                for (int i = 0; i < requestPoints[1]; i++) {
                    for (int j = 0; j < n_dims; j++) {
                        points[i][j] = points_from_previous_process[i][j];
                    }
                }

                for (int i = requestPoints[1]; i < request_size_with_previous_process; i++) {
                    for (int j = 0; j < n_dims; j++) {
                        points[i][j] = points_received[ortho_points_to_recv[i - requestPoints[1]].point_id][j];
                    }
                }

                for (int i = request_size_with_previous_process; i < end; i++) {
                    for (int j = 0; j < n_dims; j++) {
                        points[i][j] = points_from_next_process[i - request_size_with_previous_process][j];
                    }
                }

            } else {
                //fprintf(stderr, "I AM PROCESS %d and I am merging my points - I X AND RECEIVED\n", me);
                for (int i = 0; i < request_size_with_previous_process; i++) {
                    for (int j = 0; j < n_dims; j++) {
                        points[i][j] = points_received[ortho_points_to_recv[i].point_id][j];
                    }
                }

                for (int i = request_size_with_previous_process; i < end; i++) {
                    for (int j = 0; j < n_dims; j++) {
                        points[i][j] = points_from_next_process[i - request_size_with_previous_process][j];
                    }
                }
            }

        } else if (request_size_with_previous_process > end) {
            int request_size = request_size_with_previous_process - end;
            requestPointsNextProcesses[0] = 2;
            requestPointsNextProcesses[1] = request_size;
            MPI_Send(requestPointsNextProcesses, 2, MPI_INT, me + 1, REQUEST_POINTS_TAG, communicator);

            double **points_for_next_process = malloc(sizeof(double) * requestPointsNextProcesses[1]);
            for (int i = 0; i < requestPointsNextProcesses[1]; i++) {
                points_for_next_process[i] = malloc(sizeof(double) * n_dims);
                for (int j = 0; j < n_dims; j++) {
                    points_for_next_process[i][j] = points_received[ortho_points_to_recv[end - requestPointsNextProcesses[1] + i].point_id][j];
                }
                MPI_Send(points_for_next_process[i], n_dims, MPI_DOUBLE, me + 1, REQUEST_POINTS_TAG, communicator);  
            } 
        
            // Merge all to points
            if (requestPoints[0] == 1) {
                //fprintf(stderr, "I AM PROCESS %d and I am merging my points - I SENT AND SENT\n", me);
                for (int i = 0; i < end; i++) {
                    for (int j = 0; j < n_dims; j++) {
                        points[i][j] = points_received[ortho_points_to_recv[i + requestPoints[1]].point_id][j];
                    }
                }

            } else if (requestPoints[0] == 2) {
                //fprintf(stderr, "I AM PROCESS %d and I am merging my points - I RECEIVED AND SENT\n", me);
                //fprintf(stderr, "Process %d ||| size = %d\n", me, requestPoints[1]);
                for (int i = 0; i < requestPoints[1]; i++) {
                    //fprintf(stderr, "Process %d ||| alive\n", me);
                    for (int j = 0; j < n_dims; j++) {
                        points[i][j] = points_from_previous_process[i][j];
                    }
                    //fprintf(stderr, "Process %d ||| dims ok\n", me);
                }
                for (int i = requestPoints[1]; i < end; i++) {
                    for (int j = 0; j < n_dims; j++) {
                        points[i][j] = points_received[ortho_points_to_recv[i - requestPoints[1]].point_id][j];
                    }
                }

            } else {
                //fprintf(stderr, "I AM PROCESS %d and I am merging my points - I X AND SENT\n", me);
                for (int i = 0; i < end; i++) {
                    for (int j = 0; j < n_dims; j++) {
                        points[i][j] = points_received[ortho_points_to_recv[i].point_id][j];
                    }
                }
            }

        } else {
            requestPointsNextProcesses[0] = 3;
            requestPointsNextProcesses[1] = 0;
            MPI_Send(requestPointsNextProcesses, 2, MPI_INT, me + 1, REQUEST_POINTS_TAG, communicator);

            // Merge all to points
            //fprintf(stderr, "I AM PROCESS %d and I am merging my points\n", me);
            if (requestPoints[0] == 1) {
                //fprintf(stderr, "I AM PROCESS %d and I am merging my points - I SENT AND X\n", me);
                for (int i = 0; i < end; i++) {
                    for (int j = 0; j < n_dims; j++) {
                        points[i][j] = points_received[ortho_points_to_recv[i + requestPoints[1]].point_id][j];
                    }
                }

            } else if (requestPoints[0] == 2) {
                //fprintf(stderr, "I AM PROCESS %d and I am merging my points - I RECEIVED AND X\n", me);
                for (int i = 0; i < requestPoints[1]; i++) {
                    for (int j = 0; j < n_dims; j++) {
                        points[i][j] = points_from_previous_process[i][j];
                    }
                }

                for (int i = requestPoints[1]; i < end; i++) {
                    for (int j = 0; j < n_dims; j++) {
                        points[i][j] = points_received[ortho_points_to_recv[i - requestPoints[1]].point_id][j];
                    }
                }

            } else {
                //fprintf(stderr, "I AM PROCESS %d and I am merging my points - I X AND X\n", me);
                for (int i = 0; i < end; i++) {
                    for (int j = 0; j < n_dims; j++) {
                        points[i][j] = points_received[ortho_points_to_recv[i].point_id][j];
                    }
                }
            }
        }
    } else {
        // Merge all to points
        if (requestPoints[0] == 1) {
            //fprintf(stderr, "I AM PROCESS %d and I am merging my points - I SENT AND X\n", me);
            for (int i = 0; i < end; i++) {
                for (int j = 0; j < n_dims; j++) {
                    points[i][j] = points_received[ortho_points_to_recv[i + requestPoints[1]].point_id][j];
                }
            }
        } else if (requestPoints[0] == 2) {
            //fprintf(stderr, "I AM PROCESS %d and I am merging my points - I RECEIVED AND X\n", me);
            for (int i = 0; i < requestPoints[1]; i++) {
                for (int j = 0; j < n_dims; j++) {
                    points[i][j] = points_from_previous_process[i][j];
                }
            }

            for (int i = requestPoints[1]; i < end; i++) {
                for (int j = 0; j < n_dims; j++) {
                    points[i][j] = points_received[ortho_points_to_recv[i - requestPoints[1]].point_id][j];
                }
            }
        } else {
            //fprintf(stderr, "I AM PROCESS %d and I am merging my points - I X AND X\n", me);
            for (int i = 0; i < end; i++) {
                for (int j = 0; j < n_dims; j++) {
                    points[i][j] = points_received[ortho_points_to_recv[i].point_id][j];
                }
            }
        }
    }

    // Sorting and exchange finished
    //fprintf(stderr, "I AM PROCESS %d and I finished my merge\n", me);
    for (int i = 0; i < end; i++) {
        ortho_points[i].point_id = i;
    }

    /*for (int i = 0; i < end; i++) {
        fprintf(stderr, "[PROCESS %d]: points after Alltoall: ", me);
        for (int j = 0; j < n_dims; j++) {
            fprintf(stderr, " %f", points[i][j]);
        }
        fprintf(stderr, "\n");
    }*/

    free(points_to_send);
    free(points_to_recv);
    free(ortho_points_to_recv);
    free(pivots);
    free(send_counts);
    free(recv_counts);
    free(send_displs);
    free(recv_displs);

    double *p1 = malloc(sizeof(double) * n_dims);
    double *p2 = malloc(sizeof(double) * n_dims);
    //fprintf(stderr, "I am process %d and i am sending my median values\n", me);
    if (initial_size % 2 == 0) {
        if (max_processes > 2) {
            if ((max_processes / 2 - 1) == me) MPI_Send(points[end - 1], n_dims, MPI_DOUBLE, 0, MEDIAN_POINT_TAG, communicator);
            else if ((max_processes / 2) == me) MPI_Send(points[0], n_dims, MPI_DOUBLE, 0, MEDIAN_POINT_TAG, communicator);
        } 
        else MPI_Send(points[0], n_dims, MPI_DOUBLE, 0, MEDIAN_POINT_TAG, communicator);
    } else {
        // 1 process
    }

    /*
     * Separate L and R sets
     */
    double *median_point = malloc(sizeof(double) * n_dims);
    //fprintf(stderr, "I am process %d and i am receiving my median point\n", me);
    MPI_Recv(median_point, n_dims, MPI_DOUBLE, 0, MEDIAN_POINT_TAG, communicator, MPI_STATUS_IGNORE);

    double furthest_distance = get_furthest_distance_parallel(median_point, start, end, threads);
    //fprintf(stderr, "I am process %d and i am sending my furthest distance = %f\n", me, furthest_distance);
    MPI_Send(&furthest_distance, 1, MPI_DOUBLE, 0, FURTHEST_DISTANCE_TAG, communicator);
    
    /*
     * Split processes in smaller groups or start serial version
     */
    /*for (int i = 0; i < end; i++) {
        fprintf(stderr, "Process %d point %d =", rank_on_world, i);
        for (int j = 0; j < n_dims; j++) {
            fprintf(stderr, " %f", points[i][j]);
        }
        fprintf(stderr, "\n");
    }*/

    node_t *tree = malloc(sizeof(node_t));
    if (max_processes > 2) {
        int half_max_processes = max_processes / 2;
        int ranks[half_max_processes];
        if (me < half_max_processes) for (int i = 0; i < half_max_processes; i++) ranks[i] = i;
        else                         for (int i = half_max_processes; i < max_processes; i++) ranks[i - half_max_processes] = i;

        MPI_Group world_group;
        MPI_Comm_group(WORLD, &world_group);

        MPI_Group new_group;
        MPI_Group_incl(world_group, half_max_processes, ranks, &new_group);

        // Create a new communicator based on the new group
        MPI_Comm new_communicator;
        MPI_Comm_create_group(MPI_COMM_WORLD, new_group, 0, &new_communicator);

        MPI_Comm_rank(new_communicator, &me);
    
        if (me == 0) {
            tree = build_tree_parallel_mpi(start, end, me, half_max_processes, new_communicator, threads);
        } else {
            wait_mpi(start, end, me, half_max_processes, new_communicator, threads);
        }
    } else {
        if (threads > 2) {
            tree = build_tree_parallel_omp(start, end, threads, me);
        } else {
            tree = build_tree(start, end);
        }
    }

    // Send end confirmation
    int confirmation = 1;
    MPI_Send(&confirmation, 1, MPI_INT, 0, CONFIRMATION_TAG, WORLD);
    //fprintf(stderr, "//////////////%d sent confirmation\n", me);

    MPI_Barrier(WORLD);
    MPI_Comm_rank(WORLD, &me);

    // Send node count
    long n_count = 0;
    count_nodes(tree, &n_count);
    //fprintf(stderr, "Process %d sending node count\n", me);
    MPI_Send(&n_count, 1, MPI_LONG, 0, COUNT_TAG, WORLD);

    MPI_Barrier(WORLD);

    // Receive node id and print type 
    long node_id[3];
    //fprintf(stderr, "Process %d waiting for print order\n", me);
    MPI_Recv(node_id, 3, MPI_LONG, MPI_ANY_SOURCE, NODE_TAG, WORLD, MPI_STATUS_IGNORE);

    long print_result;
    if (node_id[1]) {
        current_print_proc = me + 1;
        print_result = aux_print_tree_main_proc(tree, n_dims, points, n_count, node_id[0], me);
    } else {  
        print_result = aux_print_tree(tree, n_dims, points, n_count, node_id[0]);
    }

    MPI_Send(&print_result, 1, MPI_LONG, node_id[2], PRINT_TAG, WORLD);

    if (node_id[1]) {
        MPI_Send(&current_print_proc, 1, MPI_INT, node_id[2], PRINT_2_TAG, WORLD);
    }

    MPI_Finalize();
    exit(0);
}

int ballAlg_large_mpi(int argc, char *argv[], long n_samples) {
    double exec_time;
    double *first_point;
    long start, end;

    initial_size = n_samples;
    n_dims = atoi(argv[1]);

    omp_set_nested(1);
    omp_set_dynamic(1);
    max_threads = omp_get_max_threads();

    exec_time = -omp_get_wtime();

    int me;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    initial_procs = nprocs;

    MPI_Group world_group;
    MPI_Comm_group(WORLD, &world_group);
    //printf("Create world group\n");

    start = (n_samples / nprocs) * me;
    end = start + (n_samples / nprocs);

    first_point = malloc(sizeof(double) * n_dims);
    points = get_points_large(argc, argv, n_dims, n_samples, start, end, first_point);
    n_samples = end - start;

    /*fprintf(stderr, "[PROCESS %d] First points", me);
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
    }*/

    create_furthest_point_mpi(&furthest_point_mpi);

    /*
     * Get ortho projection of points in line ab
     */
    ortho_points = malloc(sizeof(node_t) * n_samples);

    #pragma omp parallel for 
    for (long i = 0; i < n_samples; i++) {
        ortho_points[i].center = malloc(sizeof(double) * n_dims);
        ortho_points[i].point_id = i;
    }
    
    if (me != 0)
        wait_mpi(0, n_samples, me, nprocs, WORLD, max_threads);

    node_t *tree = build_tree_parallel_mpi(0, n_samples, 0, nprocs, WORLD, max_threads);
    
    // Receive confirmation
    int confirmation;
    //fprintf(stderr, "I am process %d waiting for other processes to send me end confirmation\n", me);
    for (int i = 1; i < nprocs; i++) {
        MPI_Recv(&confirmation, 1, MPI_LONG, MPI_ANY_SOURCE, CONFIRMATION_TAG, WORLD, MPI_STATUS_IGNORE);
    }

    exec_time += omp_get_wtime();
    fprintf(stderr, "%lf\n", exec_time);
    fflush(stderr);
    MPI_Barrier(WORLD);

    // Receive node count
    int node_count = 0;
    int count;
    for (int i = 1; i < nprocs; i++){
        count = 0;
        MPI_Recv(&count, 1, MPI_LONG, MPI_ANY_SOURCE, COUNT_TAG, WORLD, MPI_STATUS_IGNORE);
        node_count += count;
    }

    //fprintf(stderr, "Process %d with node count %d\n", me, node_count);
    print_tree_mpi(tree, n_dims, points, node_count, me);
    
    //fprintf(stderr, "I AM ALMOSTT ENDIIIIIIIIIIIIIIIIIIING\n");
    free(ortho_points);
    free_node(tree);
    
    free(points[0]);
    free(points);

    //fprintf(stderr, "Process %d is ENDING\n", me);
    MPI_Finalize();
    exit(0);
}
