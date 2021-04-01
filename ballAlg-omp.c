// OpenMP version

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "gen_points.h"
#include "quickSelect.h"

double **points;
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

double distance_sqrd(double *pt1, double *pt2, node_t *ortho_point) {
  double dist = 0.0;

  for (int d = 0; d < n_dims; d++) {
    double tmp = pt1[d] - pt2[d];
    dist += tmp * tmp;
    ortho_point->center[d] = tmp;
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

int get_furthest_point(long point, long n_set, node_t* ortho_points) {
  long max = point;
  double *point_point = points[ortho_points[point].point_id];
  double distance, max_distance = 0.0;

  //#pragma omp parallel for
  for (long i = 0; i < n_set; i++) {
    distance = distance_sqrd(points[ortho_points[i].point_id], point_point, &ortho_points[i]);
    if (max_distance < distance) {
      max = i;
      max_distance = distance;
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

node_t *build_tree(long n_set, node_t* ortho_points) {
  if (n_set == 0) {
    return NULL;

  } else if (n_set == 1) {
    return create_node(points[ortho_points[0].point_id], ortho_points[0].point_id, 0);

  } else if (n_set == 2) {
    double *median_point = malloc(sizeof(double) * n_dims);
    double *point1 = points[ortho_points[0].point_id];
    double *point2 = points[ortho_points[1].point_id];
    double dist = distance2(point1, point1, point2, median_point);

    node_t *tree = create_node(median_point, -1, dist);
    
    if ((point1[0] - point2[0]) < 0) {
      tree->L = create_node(point1, ortho_points[0].point_id, 0);
      tree->R = create_node(point2, ortho_points[1].point_id, 0);
    } else {
      tree->L = create_node(point2, ortho_points[1].point_id, 0);
      tree->R = create_node(point1, ortho_points[0].point_id, 0);
    }
    
    return tree;
  }

  /*
   * Get furthest points
   */
  long a = get_furthest_point(0, n_set, ortho_points);
  long b = get_furthest_point(a, n_set, ortho_points);

  double *point_a = points[ortho_points[a].point_id];
  double *point_b = points[ortho_points[b].point_id];


  /*
   * Get projections to allow median calc
   */
  int d;
  //#pragma omp parallel for private(d) 
  for (int i = 0; i < n_set; i++) {
    double projection = 0.0;
    for (d = 0; d < n_dims; d++) {
      projection += ortho_points[i].center[d] * (point_b[d] - point_a[d]);
    }
    ortho_points[i].center[0] = projection * (point_b[0] - point_a[0]);
  }

  /*
   * Get median point which will be the center of the ball
   */
  medianValues median_ids = quickSelect(ortho_points, n_set);

  /*
   * Separate L and R sets
   */
  long left_set_count = n_set / 2;
  long right_set_count = n_set / 2;
  double *median_point = malloc(sizeof(double) * n_dims);
  node_t *right_set;
  node_t *left_set = malloc(sizeof(node_t) * left_set_count);
  node_t *tree = create_node(median_point, -1, -1);

  if (n_set % 2 != 0) {
    right_set_count += 1;
    right_set = malloc(sizeof(node_t) * right_set_count);

    double median_before_projection = ortho_points[median_ids.first].center[0];
    int j = 0, k = 0;
    for (int i = 0; i < n_set; i++) {
      if ((ortho_points[i].center[0] - median_before_projection) < 0)
        left_set[j++] = ortho_points[i];
      else 
        right_set[k++] = ortho_points[i];
    }
  } else {
    right_set = malloc(sizeof(node_t) * right_set_count);

    double median_before_projection = ortho_points[median_ids.first].center[0];
    int j = 0, k = 0;
    for (int i = 0; i < n_set; i++) {
      if ((ortho_points[i].center[0] - median_before_projection) <= 0) 
        left_set[j++] = ortho_points[i];
      else 
        right_set[k++] = ortho_points[i];
    }
  }

  #pragma omp parallel
  #pragma omp single
  {
    #pragma omp task
    {
      /* Calc ortho projection of median points */
      double *p1 = points[ortho_points[median_ids.first].point_id];
      double *p2 = points[ortho_points[median_ids.second].point_id];
      calc_ortho_projection(point_a, point_b, p1, p2, ortho_points, median_ids.first, median_ids.second);

      /*
      * Get the radius of the ball (largest distance)
      */
      double distances[2] = {0.0, 0.0};
      if (n_set % 2 != 0) {
        distances[0] = distance(point_a, ortho_points[median_ids.first].center, median_point);
        distances[1] = distance(point_b, ortho_points[median_ids.first].center, median_point);
        
      } else {
        distances[0] = distance2(point_a, ortho_points[median_ids.first].center, ortho_points[median_ids.second].center, median_point);
        distances[1] = distance2(point_b, ortho_points[median_ids.first].center, ortho_points[median_ids.second].center, median_point);
      }

      tree->radius = ((distances[0] - distances[1]) > 0) ? distances[0] : distances[1];
    }

    #pragma omp task depend(in: left_set)
    {
      tree->L = build_tree(left_set_count, left_set);
      free(left_set);
    }

    #pragma omp task depend(in: right_set, right_set_count)
    {
      tree->R = build_tree(right_set_count, right_set);
      free(right_set);
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

  exec_time = -omp_get_wtime();
  points = get_points(argc, argv, &n_dims, &n_samples);

  /*
   * Get ortho projection of points in line ab
   */
  node_t *ortho_points = malloc(sizeof(node_t) * n_samples);
  //#pragma omp parallel for
  for (long i = 0; i < n_samples; i++) {
    ortho_points[i].center = malloc(sizeof(double) * n_dims);
    ortho_points[i].point_id = i;
  }
  
  node_t *tree = build_tree(n_samples, ortho_points);

  exec_time += omp_get_wtime();
  fprintf(stderr, "%lf\n", exec_time);

  if (n_samples < 1000)
    print_tree(tree, n_dims, points);

  free(ortho_points);
  free_node(tree);
  free(points[0]);
  free(points);
}