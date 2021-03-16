#ifndef __H_QUICK_SORT__
#define __H_QUICK_SORT__

void swap(double **points, int n_dims, long a, long b, long* set);
long partition(double **arr, long low, long high, int dim, int n_dims, long* set);
void quick_sort(double **arr, long low, long high, int dim, int n_dims, long* set);

#endif
