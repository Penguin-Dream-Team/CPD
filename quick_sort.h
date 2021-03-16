#ifndef __H_QUICK_SORT__
#define __H_QUICK_SORT__

void swap(double **points, int n_dims, int a, int b);
int partition(double **arr, int low, int high, int dim, int n_dims);
void quick_sort(double **arr, int low, int high, int dim, int n_dims);

#endif
