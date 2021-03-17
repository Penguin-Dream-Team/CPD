#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "gen_points.h"
#include "quick_sort.h"

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

double distance_sqrd(int n_dims, double *pt1, double *pt2) {
  double dist = 0.0;

  for (int d = 0; d < n_dims; d++)
    dist += (pt1[d] - pt2[d]) * (pt1[d] - pt2[d]);
  return dist;
}

double distance(int n_dims, double *pt1, double *pt2) {
  return sqrt(distance_sqrd(n_dims, pt1, pt2));
}

int get_furthest_point(double **points, long point, int n_dims, long n_set,
                       long *set) {
  double *distances = malloc(sizeof(double) * n_set);

  long max = point;
  double max_distance = 0;
  for (long i = 0; i < n_set; i++) {
    long current_point = set[i];
    if (i != point) {
      //printf("i %ld n_dims %d set_point %ld current_point %ld point %ld\n", i, n_dims, set[point], current_point, point);
      points[set[point]];
      points[current_point];
      distances[i];
      distance(n_dims, points[set[point]], points[current_point]);
      distances[i] = distance(n_dims, points[set[point]], points[current_point]);
      if (max_distance < distances[i]) {
        max = i;
        max_distance = distances[i];
      }
    }
  }
  //printf("\n");
  free(distances);
  return max;
}

void calc_ortho_projection(double **points, int n_dims, int a, int b, int p,
                           long *set, node_t *out_points) {
  long index_a = set[a];
  long index_b = set[b];
  long index_p = set[p];
  double top_inner_product = 0;
  double bot_inner_product = 0;
  for (int i = 0; i < n_dims; i++) {
    double *value_a = &points[index_a][i];
    double *value_b = &points[index_b][i];
    double *value_p = &points[index_p][i];
    top_inner_product += (*value_p - *value_a) * (*value_b - *value_a);
    bot_inner_product += (*value_b - *value_a) * (*value_b - *value_a);
  }
  double inner_product = top_inner_product / bot_inner_product;

  for (int i = 0; i < n_dims; i++) {
    double *value_a = &points[index_a][i];
    double *value_b = &points[index_b][i];
    out_points[p].center[i] = inner_product * (*value_b - *value_a);
    out_points[p].center[i] += *value_a;
  }
}

int cmpfunc (const void * pa, const void * pb) {
  const node_t *a = (const node_t *)pa;
  const node_t *b = (const node_t *)pb;

  return (int) (a->center[0] - b->center[0]);
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

  for (long i = 0; i < n_set; i++) {
    ortho_points[i].center = malloc(sizeof(double) * n_dims);
    ortho_points[i].point_id = set[i];
    calc_ortho_projection(points, n_dims, a, b, i, set, ortho_points);
  }

  /*
   * Sort ortho projection points
   */
  //quick_sort(ortho_points, 0, n_set - 1, 0, n_dims, set);
  qsort(ortho_points, n_set, sizeof(node_t), cmpfunc);
  //qsort(points, n_set, sizeof(points[0]), cmpfunc);
  /*
   * Get median point which will be the center of the ball
   */
  double *median_point = malloc(sizeof(double) * n_dims);
  if (n_set % 2 == 0) {
    int mid_point_i = n_set / 2 - 1;
    for (int i = 0; i < n_dims; i++) {
      median_point[i] =
          (ortho_points[mid_point_i].center[i] + ortho_points[mid_point_i + 1].center[i]) / 2;
    }
  } else {
    int mid_point_i = n_set / 2;
    for (int i = 0; i < n_dims; i++) {
      median_point[i] = ortho_points[mid_point_i].center[i];
    }
  }

  /*
   * Get the radius of the ball (largest distance)
   */
  double distances[2] = {0.0, 0.0};

  distances[0] = distance(n_dims, points[ortho_points[0].point_id], median_point);
  distances[1] = distance(n_dims, points[ortho_points[n_set - 1].point_id], median_point);

  double radius = (distances[0] > distances[1]) ? distances[0] : distances[1];

  /*
   * Separate L and R sets
   */
  long *left_set = malloc(sizeof(long) * n_set);
  long *right_set = malloc(sizeof(long) * n_set);
  long left_set_count = 0;
  long right_set_count = 0;
  for (long i = 0; i < n_set; i++) {
    long point_index = ortho_points[i].point_id;
    if (ortho_points[i].center[0] < median_point[0]) {
      left_set[left_set_count++] = point_index;
    } else {
      right_set[right_set_count++] = point_index;
    }
  }

  node_t *tree = create_node(median_point, -1, radius);
  tree->L = build_tree(points, n_dims, left_set_count, left_set, ortho_points);
  tree->R = build_tree(points, n_dims, right_set_count, right_set, ortho_points);

//free(left_set);
//free(right_set);

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
