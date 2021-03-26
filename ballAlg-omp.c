// OpenMP version

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "gen_points.h"

typedef struct _node {
  long point_id;
  double *center;
  double radius;
  struct _node *L;
  struct _node *R;
} node_t;

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

double distance_sqrd(int n_dims, double *pt1, double *pt2, node_t *ortho_point) {
  double dist = 0.0;

  #pragma omp parallel for
  for (int d = 0; d < n_dims; d++) {
    double tmp = pt1[d] - pt2[d];
    dist += tmp * tmp;
    ortho_point->center[d] = tmp;
  }

  return dist;
}

double distance_sqrd_med(int n_dims, double *pt1, double *pt2, double *median) {
  double dist = 0.0;

  #pragma omp parallel for
  for (int d = 0; d < n_dims; d++) {
    median[d] = pt2[d];
    dist += (pt1[d] - pt2[d]) * (pt1[d] - pt2[d]);
  }

  return dist;
}

double distance_sqrd_med2(int n_dims, double *pt1, double *pt2, double *pt3, double *median) {
  double dist = 0.0;

  #pragma omp parallel for
  for (int d = 0; d < n_dims; d++) {
    median[d] = (pt2[d] + pt3[d]) / 2;
    dist += (pt1[d] - median[d]) * (pt1[d] - median[d]);
  }

  return dist;
}

double distance(int n_dims, double *pt1, double *pt2, double *median) {
  return sqrt(distance_sqrd_med(n_dims, pt1, pt2, median));
}

double distance2(int n_dims, double *pt1, double *pt2, double *pt3, double *median) {
  return sqrt(distance_sqrd_med2(n_dims, pt1, pt2, pt3, median));
}

int get_furthest_point(double **points, long point, int n_dims, long n_set,
                       node_t* ortho_points) {
  long max = point;
  double *point_point = points[ortho_points[point].point_id];
  double distance, max_distance = 0.0;
  for (long i = 0; i < n_set; i++) {
    distance = distance_sqrd(n_dims, points[ortho_points[i].point_id], point_point, &ortho_points[i]);
    if (max_distance < distance) {
      max = i;
      max_distance = distance;
    }
  }
  return max;
}

void calc_ortho_projection(double **points, int n_set, int n_dims, long index_a, long index_b, long index_p1, long index_p2, node_t *ortho_points) {
  double top_inner_product1 = 0;
  double top_inner_product2 = 0;
  double bot_inner_product = 0;

  for (int i = 0; i < n_dims; i++) {
    double *value_a = &points[index_a][i];
    double *value_b = &points[index_b][i];
    double *value_p1 = &points[index_p1][i];
    double *value_p2 = &points[index_p2][i];
    double b_minus_a = *value_b - *value_a;
    top_inner_product1 += (*value_p1 - *value_a) * b_minus_a;
    top_inner_product2 += (*value_p2 - *value_a) * b_minus_a;
    bot_inner_product += b_minus_a * b_minus_a;
  }

  double inner_product1 = top_inner_product1 / bot_inner_product;
  double inner_product2 = top_inner_product2 / bot_inner_product;

  #pragma omp parallel for
  for (int i = 0; i < n_dims; i++) {
    double *value_a = &points[index_a][i];
    double *value_b = &points[index_b][i];
    ortho_points[n_set / 2].center[i] = inner_product1 * (*value_b - *value_a) + *value_a;
    ortho_points[n_set / 2 - 1].center[i] = inner_product2 * (*value_b - *value_a) + *value_a;
  }
}

int cmpfunc (const void * pa, const void * pb) {
  const node_t *a = (const node_t *)pa;
  const node_t *b = (const node_t *)pb;

  if (a->center[0] > b->center[0]) {
    return 1;
  } else if (a->center[0] < b->center[0]) {
    return -1;
  } else {
    return 0;
  }
}

node_t *build_tree(double **points, int n_dims, long n_set, node_t* ortho_points) {
  if (n_set == 0) {
    return NULL;

  } else if (n_set == 1) {
    return create_node(points[ortho_points[0].point_id], ortho_points[0].point_id, 0);

  } else if (n_set == 2) {
    double *median_point = malloc(sizeof(double) * n_dims);
    double *point1 = points[ortho_points[0].point_id];
    double *point2 = points[ortho_points[1].point_id];

    double dist = distance2(n_dims, point1, point1, point2, median_point);

    node_t *tree = create_node(median_point, -1, dist);
    
    if (point1[0] < point2[0]) {
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
  long a = get_furthest_point(points, 0, n_dims, n_set, ortho_points);
  long b = get_furthest_point(points, a, n_dims, n_set, ortho_points);

  double *point_a = points[ortho_points[a].point_id];
  double *point_b = points[ortho_points[b].point_id];

  for (int i = 0; i < n_set; i++) {
    double projection = 0.0;
    for (int d = 0; d < n_dims; d++) {
      projection += ortho_points[i].center[d] * (point_b[d] - point_a[d]);
    }
    ortho_points[i].center[0] = projection * (point_b[0] - point_a[0]);
  }

  /*
  * Sort ortho projection points
  */
  qsort(ortho_points, n_set, sizeof(node_t), cmpfunc);

  /*
   * Get median point which will be the center of the ball
   */
  double *median_point = malloc(sizeof(double) * n_dims);

  /*
   * Get the radius of the ball (largest distance)
   */
  double distances[2] = {0.0, 0.0};
  long left_set_count;
  long right_set_count;
  
  /* Calc ortho projection of median points */
  calc_ortho_projection(points, n_set, n_dims, ortho_points[0].point_id, ortho_points[n_set - 1].point_id,
      ortho_points[n_set / 2].point_id, ortho_points[n_set / 2 - 1].point_id, ortho_points);

  if (n_set % 2 != 0) {
    distances[0] = distance(n_dims, points[ortho_points[0].point_id], ortho_points[n_set / 2].center, median_point);
    distances[1] = distance(n_dims, points[ortho_points[n_set - 1].point_id], ortho_points[n_set / 2].center, median_point);
    left_set_count = n_set / 2;
    right_set_count = n_set / 2 + 1;

  } else {
    distances[0] = distance2(n_dims, points[ortho_points[0].point_id], ortho_points[n_set / 2].center, ortho_points[n_set / 2 - 1].center, median_point);
    distances[1] = distance2(n_dims, points[ortho_points[n_set - 1].point_id], ortho_points[n_set / 2].center, ortho_points[n_set / 2 - 1].center, median_point);
    left_set_count = n_set / 2;
    right_set_count = n_set / 2;
  }

  double radius = (distances[0] > distances[1]) ? distances[0] : distances[1];

  /*
   * Separate L and R sets
   */
  node_t *left_set = malloc(sizeof(node_t) * left_set_count);
  node_t *right_set = malloc(sizeof(node_t) * right_set_count);
  memcpy(left_set, ortho_points, sizeof(node_t) * left_set_count); 
  memcpy(right_set, &ortho_points[left_set_count], sizeof(node_t) * right_set_count);
  
  node_t *tree = create_node(median_point, -1, radius);
  
  tree->L = build_tree(points, n_dims, left_set_count, left_set);
  tree->R = build_tree(points, n_dims, right_set_count, right_set);

  free(left_set);
  free(right_set);

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

  exec_time = -omp_get_wtime();
  int n_dims;
  long n_samples;
  double **points = get_points(argc, argv, &n_dims, &n_samples);

  /*
   * Get ortho projection of points in line ab
   */
  node_t *ortho_points = malloc(sizeof(node_t) * n_samples);
 
  #pragma omp parallel for
  for (long i = 0; i < n_samples; i++) {
    ortho_points[i].center = malloc(sizeof(double) * n_dims);
    ortho_points[i].point_id = i;
  }
  
  node_t *tree = build_tree(points, n_dims, n_samples, ortho_points);

  exec_time += omp_get_wtime();
  fprintf(stderr, "%lf\n", exec_time);

  if (n_samples < 1000)
    print_tree(tree, n_dims, points);

  free(ortho_points);
  free_node(tree);
  free(points[0]);
  free(points);
}