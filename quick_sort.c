#include "quick_sort.h"

void swap(double **arr, int n_dims, long a, long b, long *set) {
  for (long i = 0; i < n_dims; i++) {
    double t = arr[a][i];
    arr[a][i] = arr[b][i];
    arr[b][i] = t;
  }
  long t = set[a];
  set[a] = set[b];
  set[b] = t;
}

long partition(double **arr, long low, long high, int dim, int n_dims,
              long *set) {
  double pivot = arr[high][dim];
  long i = (low - 1);

  for (long j = low; j <= high - 1; j++) {
    if (arr[j][dim] < pivot) {
      i++;
      swap(arr, n_dims, i, j, set);
    }
  }

  swap(arr, n_dims, i + 1, high, set);
  return (i + 1);
}

void quick_sort(double **arr, long low, long high, int dim, int n_dims,
                long *set) {
  if (low < high) {
    long pi = partition(arr, low, high, dim, n_dims, set);
    quick_sort(arr, low, pi - 1, dim, n_dims, set);
    quick_sort(arr, pi + 1, high, dim, n_dims, set);
  }
}
