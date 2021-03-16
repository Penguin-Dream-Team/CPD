#include "quick_sort.h"

void swap(double **arr, int n_dims, int a, int b) {
  for (int i = 0; i < n_dims; i++) {
    double t = arr[a][i];
    arr[a][i] = arr[b][i];
    arr[b][i] = t;
  }
}

int partition(double **arr, int low, int high, int dim, int n_dims) {
  double pivot = arr[high][dim];
  int i = (low - 1);

  for (int j = low; j <= high - 1; j++) {
    if (arr[j][dim] < pivot) {
      i++;
      swap(arr, n_dims, i, j);
    }
  }

  swap(arr, n_dims, i + 1, high);
  return (i + 1);
}

void quick_sort(double **arr, int low, int high, int dim, int n_dims) {
  if (low < high) {
    int pi = partition(arr, low, high, dim, n_dims);
    quick_sort(arr, low, pi - 1, dim, n_dims);
    quick_sort(arr, pi + 1, high, dim, n_dims);
  }
}
