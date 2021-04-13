// OpenMP version

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "gen_points.h"
#include "quickSelect.h"

#define ELEM_SWAP(a,b) { register node_t t = (a); (a) = (b); (b) = t; }

typedef struct {
    double max_distance;
    long max;
} furthest_point;

double **points;
node_t *ortho_points;
furthest_point* furthest_points;
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

double distance_sqrd_med(double *pt1, double *pt2, double *median) {
  double dist = 0.0;

  for (int d = 0; d < n_dims; d++) {
    median[d] = pt2[d];
    dist += (pt1[d] - pt2[d]) * (pt1[d] - pt2[d]);
  }

  return dist;
}

double distance_sqrd_med2(double *pt1, double *pt2, double *pt3, double *median) {
  double dist = 0.0;

  for (int d = 0; d < n_dims; d++) {
    median[d] = (pt2[d] + pt3[d]) / 2;
    dist += (pt1[d] - median[d]) * (pt1[d] - median[d]);
  }
  return dist;
}

double distance(double *pt1, double *pt2, double *median) {
  return sqrt(distance_sqrd_med(pt1, pt2, median));
}

double distance2(double *pt1, double *pt2, double *pt3, double *median) {
  return sqrt(distance_sqrd_med2(pt1, pt2, pt3, median));
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

int get_furthest_point_parallel(long point, long start, long end, int threads) {
  double *point_point = points[ortho_points[point].point_id];

  for (int i = 0; i < threads; i++) {
    furthest_points[i].max = point;
    furthest_points[i].max_distance = 0.0;
  }

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
      // FIXME Isto dá mal porque está sempre a ir buscar o valor inicial de furthest_points, mas fica bem mais rápido.
      // É uma questão de consgeuir meter mais eficiente
      // Verificar o quickSelect - mudei os ifs de retorno porque nos estava a dar ao contrário
      /*int num_thread = omp_get_thread_num();
      furthest_point fp = furthest_points[num_thread];

      if (fp.max_distance < distance) {
        furthest_points[num_thread].max = i;
        furthest_points[num_thread].max_distance = distance;
      }*/
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
    return create_node(points[ortho_points[0].point_id], ortho_points[0].point_id, 0);

  } else if (start == end - 2) {  // 2 points
    double *median_point = malloc(sizeof(double) * n_dims);
    double *point1 = points[ortho_points[start].point_id];
    double *point2 = points[ortho_points[end - 1].point_id];
    double dist = distance2(point1, point1, point2, median_point);

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

  /*fprintf(stderr, "median = %f %f\n", point_median_1.center[0], point_median_2.center[0]);
  for (int i = start; i < end; i++) {
    fprintf(stderr, "value at %d = %f\n", i, ortho_points[i].center[0]);
  }*/
  /*long high = end - 1;
  ELEM_SWAP(ortho_points[median_ids.second], ortho_points[high]);

  // pivot (Element to be placed at right position)
  double pivot = point_median_2.center[0];  

  long i = start - 1;   // Index of smaller element and indicates the 
                        // right position of pivot found so far

  for (long j = start; j <= high - 1; j++) {
    // If current element is smaller than the pivot
    if ((ortho_points[j].center[0] - pivot) < 0)
    {
      i++;    // increment index of smaller element
      ELEM_SWAP(ortho_points[i], ortho_points[j]);
    }
  }

  ELEM_SWAP(ortho_points[i + 1], ortho_points[high]);*/

  #pragma omp task
  {
    /* Calc ortho projection of median points */
    double *p1 = points[point_median_1.point_id];
    double *p2 = points[point_median_2.point_id];
    calc_ortho_projection(point_a, point_b, p1, p2, ortho_points, median_ids.first, median_ids.second);

    /*
    * Get the radius of the ball (largest distance)
    */
    double distances[2] = {0.0, 0.0};
    if ((end - start) % 2 != 0) {
      distances[0] = distance(point_a, point_median_1.center, median_point);
      distances[1] = distance(point_b, point_median_1.center, median_point);
      
    } else {
      distances[0] = distance2(point_a, point_median_1.center, point_median_2.center, median_point);
      distances[1] = distance2(point_b, point_median_1.center, point_median_2.center, median_point);
    }

    tree->radius = ((distances[0] - distances[1]) > 0) ? distances[0] : distances[1];
  }

  //#pragma omp task
  //{
    
    if (start == median_ids.second) { 
      tree->L = NULL; 
    } else {
      #pragma omp task
        tree->L = build_tree(start, median_ids.second);
    }

    if (median_ids.second == end) {
      tree->R = NULL;
    } else {
      #pragma omp task
        tree->R = build_tree(median_ids.second, end);
    }
  //}

  //#pragma omp taskwait

  return tree;
}

// not inclusive
node_t *build_tree_parallel(long start, long end, int threads) {  
  if (start == end - 1) {  // 1 point
  return create_node(points[ortho_points[0].point_id], ortho_points[0].point_id, 0);

  } else if (start == end - 2) {  // 2 points
    double *median_point = malloc(sizeof(double) * n_dims);
    double *point1 = points[ortho_points[start].point_id];
    double *point2 = points[ortho_points[end - 1].point_id];
    double dist = distance2(point1, point1, point2, median_point);

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

  /*fprintf(stderr, "median = %f %f\n", point_median_1.center[0], point_median_2.center[0]);
  for (int i = start; i < end; i++) {
    fprintf(stderr, "value at %d = %f\n", i, ortho_points[i].center[0]);
  }*/

  /*long high = end - 1;
  ELEM_SWAP(ortho_points[median_ids.second], ortho_points[high]);

  // pivot (Element to be placed at right position)
  double pivot = ortho_points[high].center[0];  

  long i = start - 1;   // Index of smaller element and indicates the 
                        // right position of pivot found so far

  for (long j = start; j <= high - 1; j++) {
    // If current element is smaller than the pivot
    if ((ortho_points[j].center[0] - pivot) < 0)
    {
      i++;    // increment index of smaller element
      ELEM_SWAP(ortho_points[i], ortho_points[j]);
    }
  }

  ELEM_SWAP(ortho_points[i + 1], ortho_points[high]);*/

  //depth++;
  //fprintf(stderr, "%d - %ld => %d\n", threads, end - start, depth);
  #pragma omp parallel
  {
    #pragma omp single 
    {
      #pragma omp task
      {
        /* Calc ortho projection of median points */
        double *p1 = points[point_median_1.point_id];
        double *p2 = points[point_median_2.point_id];
        calc_ortho_projection(point_a, point_b, p1, p2, ortho_points, median_ids.first, median_ids.second);

        /*
        * Get the radius of the ball (largest distance)
        */
        double distances[2] = {0.0, 0.0};
        if ((end - start) % 2 != 0) {
          distances[0] = distance(point_a, point_median_1.center, median_point);
          distances[1] = distance(point_b, point_median_1.center, median_point);
          
        } else {
          distances[0] = distance2(point_a, point_median_1.center, point_median_2.center, median_point);
          distances[1] = distance2(point_b, point_median_1.center, point_median_2.center, median_point);
        }

        tree->radius = ((distances[0] - distances[1]) > 0) ? distances[0] : distances[1];
      }

      //#pragma omp task
      //{
        
      if (start == median_ids.second) { 
        tree->L = NULL; 
      //} else if ((threads / depth) > 1) {
        //fprintf(stderr, "Going left - %d\n", depth);
        //#pragma omp task
          //tree->L = build_tree_parallel(start, median_ids.second, depth, (threads + 1) / 2);
      } else {
        #pragma omp task
          tree->L = build_tree(start, median_ids.second);
      }

      if (median_ids.second == end) {
        tree->R = NULL;
      //} else if ((threads / depth) > 1) {
        //fprintf(stderr, "Going right - %d\n", depth);
        //#pragma omp task
          //tree->R = build_tree_parallel(median_ids.second, end, depth, threads / 2);
      } else {
        #pragma omp task
          tree->R = build_tree(median_ids.second, end);
      }
      //}      
    }
    #pragma omp taskwait
  }
  return tree;
}

char **aux_print_tree(node_t *tree, int n_dims, double **points,
                      long *n_count) {
  n_count[0]++;
  long my_id, left_id = 0, right_id = 0;
  long left_count = 0, right_count = 0;
  char **left = NULL, **right = NULL;

  my_id = *n_count;
  if (tree->L) {
    left_id = *n_count + 1;
    left = aux_print_tree(tree->L, n_dims, points, n_count);
    left_count = *n_count - left_id + 1;
  }
  if (tree->R) {
    right_id = *n_count + 1;
    right = aux_print_tree(tree->R, n_dims, points, n_count);
    right_count = *n_count - right_id + 1;
  }

  // 128 is MAX_SIZE for now
  char *my_line = malloc(sizeof(char) * 2048);
  sprintf(my_line, "%ld %ld %ld %f", my_id - 1, left_id - 1, right_id - 1,
          tree->radius);
  char aux_str[128];
  for (int i = 0; i < n_dims; i++) {
    sprintf(aux_str, " %f", tree->center[i]);
    strcat(my_line, aux_str);
  }
  char **result = malloc(sizeof(char *) * (1 + left_count + right_count));
  result[0] = my_line;
  for (int i = 0; i < left_count; i++) {
    result[i + 1] = left[i];
  }
  if (left) {
    free(left);
  }

  for (int i = 0; i < right_count; i++) {
    result[i + left_count + 1] = right[i];
  }
  if (right) {
    free(right);
  }

  return result;
}

void print_tree(node_t *tree, int n_dims, double **points) {
  long n_count = 0;
  char **result = aux_print_tree(tree, n_dims, points, &n_count);
  printf("%d %ld\n", n_dims, n_count);
  for (int i = 0; i < n_count; i++) {
    printf("%s\n", result[i]);
    free(result[i]);
  }
    free(result);
}

int main(int argc, char *argv[]) {
  double exec_time;
  long n_samples;
  int max_threads = omp_get_max_threads();

  omp_set_nested(1);
  omp_set_dynamic(1);

  exec_time = -omp_get_wtime();
  points = get_points(argc, argv, &n_dims, &n_samples);

  /*
   * Get ortho projection of points in line ab
   */
  ortho_points = malloc(sizeof(node_t) * n_samples);
  furthest_points = malloc((sizeof(furthest_point) + 2048) * max_threads);
  node_t * tree;
 
  #pragma omp parallel for 
  for (long i = 0; i < n_samples; i++) {
    ortho_points[i].center = malloc(sizeof(double) * n_dims);
    ortho_points[i].point_id = i;
  }

  for (int i = 0; i < max_threads; i++) {
    furthest_points[i].max = -1;
    furthest_points[i].max_distance = -1;
  }
  
  tree = build_tree_parallel(0, n_samples, max_threads);

  exec_time += omp_get_wtime();
  fprintf(stderr, "%lf\n", exec_time);

  if (n_samples < 1000)
  print_tree(tree, n_dims, points);

  free(ortho_points);
  free_node(tree);
  free(points[0]);
  free(points);
}