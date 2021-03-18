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
  if (node->L) {
    free_node(node->L);
  }
  if (node->R) {
    free_node(node->R);
  }
  if (node->radius != 0) {
    free(node->center);
  }
  free(node);
}

double distance_sqrd(int n_dims, double *pt1, double *pt2, double *median) {
  double dist = 0.0;

  if (median != NULL) {
    for (int d = 0; d < n_dims; d++) {
      median[d] = pt2[d];
      dist += (pt1[d] - pt2[d]) * (pt1[d] - pt2[d]);
    }
  } else {
    for (int d = 0; d < n_dims; d++) 
      dist += (pt1[d] - pt2[d]) * (pt1[d] - pt2[d]);
  }

  return dist;
}

double distance_sqrd2(int n_dims, double *pt1, double *pt2, double *pt3, double *median) {
  double dist = 0.0;

  for (int d = 0; d < n_dims; d++) {
    median[d] = (pt2[d] + pt3[d]) / 2;
    dist += (pt1[d] - median[d]) * (pt1[d] - median[d]);
  }
  return dist;
}

double distance(int n_dims, double *pt1, double *pt2, double *median) {
  return sqrt(distance_sqrd(n_dims, pt1, pt2, median));
}

double distance2(int n_dims, double *pt1, double *pt2, double *pt3, double *median) {
  return sqrt(distance_sqrd2(n_dims, pt1, pt2, pt3, median));
}

int get_furthest_point(double **points, long point, int n_dims, long n_set,
                       long *set) {
  long max = point;
  double distance, max_distance = 0;
  for (long i = 0; i < n_set; i++) {
    long current_point = set[i];
    if (i != point) {
      distance = distance_sqrd(n_dims, points[set[point]], points[current_point], NULL);
      if (max_distance < distance) {
        max = i;
        max_distance = distance;
      }
    }
  }
  return max;
}

void calc_ortho_projection(double **points, int n_dims, int a, double* b_minus_a_sqr_set, double b_minus_a_sqr_sum, int p,
                           long *set, node_t *out_points) {
  long index_a = set[a];
  long index_p = set[p];

  for (int i = 0; i < n_dims; i++) {
    double *value_a = &points[index_a][i];
    double *value_p = &points[index_p][i];
    out_points[p].center[i] += (*value_p - *value_a) / b_minus_a_sqr_sum * b_minus_a_sqr_set[i] / n_dims + a;
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

node_t *build_tree(double **points, int n_dims, long n_set, long *set, node_t* ortho_points) {
  if (n_set < 1) {
    return NULL;
  } else if (n_set == 1) {
    return create_node(points[set[0]], set[0], 0);
  }

  /*
   * Get furthest points
   */
  long a = 0;
  long b = 1;
  if (n_set != 2) {
    a = get_furthest_point(points, 0, n_dims, n_set, set);
    b = get_furthest_point(points, a, n_dims, n_set, set);
  }

  double *b_minus_a_sqr_set = malloc(sizeof(double) * n_dims);
  double b_minus_a_sqr_sum = 0.0;
  for (int i = 0; i < n_dims; i++) {
    double *value_a = &points[a][i];
    double *value_b = &points[b][i];
    double b_minus_a = *value_b - *value_a;
    double b_minus_a_sqr = b_minus_a * b_minus_a;
    
    b_minus_a_sqr_set[i] = b_minus_a_sqr;
    b_minus_a_sqr_sum += b_minus_a_sqr;
  }

  for (long i = 0; i < n_set; i++) {
    ortho_points[i].center = malloc(sizeof(double) * n_dims);
    ortho_points[i].point_id = set[i];
    calc_ortho_projection(points, n_dims, a, b_minus_a_sqr_set, b_minus_a_sqr_sum, i, set, ortho_points);
  }
  free(b_minus_a_sqr_set);

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
  if (n_set % 2 == 0) {
    distances[0] = distance2(n_dims, points[ortho_points[0].point_id], ortho_points[n_set / 2].center, ortho_points[n_set / 2 - 1].center, median_point);
    distances[1] = distance2(n_dims, points[ortho_points[n_set - 1].point_id], ortho_points[n_set / 2].center, ortho_points[n_set / 2 - 1].center, median_point);
    left_set_count = n_set / 2;
    right_set_count = n_set / 2;
  } else {
    distances[0] = distance(n_dims, points[ortho_points[0].point_id], ortho_points[n_set / 2].center, median_point);
    distances[1] = distance(n_dims, points[ortho_points[n_set - 1].point_id], ortho_points[n_set / 2].center, median_point);
    left_set_count = n_set / 2;
    right_set_count = n_set / 2 + 1;
  }

  double radius = (distances[0] > distances[1]) ? distances[0] : distances[1];

  /*
   * Separate L and R sets
   */
  long *left_set = malloc(sizeof(long) * left_set_count);
  long *right_set = malloc(sizeof(long) * right_set_count);
  memcpy(left_set, set, sizeof(long) * left_set_count); 
  memcpy(right_set, &set[left_set_count], sizeof(long) * right_set_count);
  //free(set);

  node_t *tree = create_node(median_point, -1, radius);

  tree->L = build_tree(points, n_dims, left_set_count, left_set, ortho_points);
  tree->R = build_tree(points, n_dims, right_set_count, right_set, ortho_points);

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
  
  long n_set = n_samples;
  long *set = malloc(sizeof(long) * n_set);
  for (long i = 0; i < n_set; i++) {
    set[i] = i;
  }

  /*
   * Get ortho projection of points in line ab
   */
  node_t *n_ortho_points = malloc(sizeof(node_t) * n_samples);

  node_t *tree = build_tree(points, n_dims, n_set, set, n_ortho_points);

  exec_time += omp_get_wtime();
  fprintf(stderr, "%lf\n", exec_time);

  print_tree(tree, n_dims, points);

  free(n_ortho_points);
  free_node(tree);
  free(set);
  free(points[0]);
  free(points);
}
