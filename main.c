#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "gen_points.h"
#include "quick_sort.h"

double distance(int n_dims, double *pt1, double *pt2) {
  double dist = 0.0;

  for (int d = 0; d < n_dims; d++)
    dist += (pt1[d] - pt2[d]) * (pt1[d] - pt2[d]);
  return sqrt(dist);
}

int get_furthest_point(double **points, int first_point, int *n_dims,
                       long *n_samples) {
  double distances[*n_samples];

  int max = first_point;
  double max_distance = 0;
  for (int i = 0; i < *n_samples; i++) {
    if (i != first_point) {
      distances[i] = distance(*n_dims, points[first_point], points[i]);
      if (max_distance < distances[i]) {
        max = i;
        max_distance = distances[i];
      }
#ifdef D_DISTANCE
      printf("Distance %d %d = %f\n", first_point, i, distances[i]);
#endif
    }
  }

  return max;
}

void calc_ortho_projection(double **points, int n_dims, int a, int b, int p,
                           double **out_points) {
  double top_inner_product = 0;
  double bot_inner_product = 0;
  for (int i = 0; i < n_dims; i++) {
    double *value_a = &points[a][i];
    double *value_b = &points[b][i];
    double *value_p = &points[p][i];
    top_inner_product += (*value_p - *value_a) * (*value_b - *value_a);
    bot_inner_product += (*value_b - *value_a) * (*value_b - *value_a);
  }
  double inner_product = top_inner_product / bot_inner_product;

  for (int i = 0; i < n_dims; i++) {
    double *value_a = &points[a][i];
    double *value_b = &points[b][i];
    out_points[p][i] = inner_product * (*value_b - *value_a);
    out_points[p][i] += *value_a;
  }
}

int main(int argc, char *argv[]) {
  int n_dims;
  long n_samples;

  double **points = get_points(argc, argv, &n_dims, &n_samples);

  printf("INITIAL POINTS:\n");
  print_points(points, n_dims, n_samples);

  /*
   * Get furthest points
   */
  int a = get_furthest_point(points, 0, &n_dims, &n_samples);
  int b = get_furthest_point(points, a, &n_dims, &n_samples);

#ifdef D_DISTANCE
  printf("A: %d | B: %d\n", a, b);
#endif
  /*
   * Get ortho projection of points in line ab
   */
  int n_set = n_samples;
  double **ortho_points = (double **)create_array_pts(n_dims, n_set);

  for (int i = 0; i < n_set; i++) {
    calc_ortho_projection(points, n_dims, a, b, i, ortho_points);
  }

  printf("ORTHO POINTS:\n");
  print_points(ortho_points, n_dims, n_set);

  /*
   * Sort ortho projection points
   */
  for (int i = 0; i < n_dims; i++) {
    quick_sort(ortho_points, 0, n_set - 1, i, n_dims);
  }
  printf("SORTED ORTHO POINTS:\n");
  print_points(ortho_points, n_dims, n_set);

  /*
   * Get median point which will be the center of the ball
   */
  double *median_point = malloc(sizeof(double) * n_dims);
  if (n_set % 2 == 0) {
    int mid_point_i = n_set / 2 - 1;
    for (int i = 0; i < n_dims; i++) {
      median_point[i] =
          (ortho_points[mid_point_i][i] + ortho_points[mid_point_i + 1][i]) / 2;
    }
  } else {
    int mid_point_i = n_set / 2;
    for (int i = 0; i < n_dims; i++) {
      median_point[i] = ortho_points[mid_point_i][i];
    }
  }

  printf("Median point:\n");
  print_point(median_point, n_dims);

  /*
   * Get the radius of the ball (largest distance)
   */
  double distances[2] = {0.0, 0.0};

  distances[0] = distance(n_dims, points[a], median_point);
  distances[1] = distance(n_dims, points[b], median_point);

  printf("Distance 1: %f\n", distances[0]);
  printf("Distance 2: %f\n", distances[1]);

  double radius = (distances[0] > distances[1]) ? distances[0] : distances[1];
  printf("RADIUS: %f\n", radius);

  free(median_point);

  free(ortho_points[0]);
  free(ortho_points);

  free(points[0]);
  free(points);
}
