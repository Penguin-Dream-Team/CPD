#define RANGE 10

extern void print_point(double *, int);
extern void print_point_indent(double *, int, int);
double **create_array_pts(int n_dims, long np);
double **get_points(int argc, char *argv[], int *n_dims, long *np);
void print_points(double** points, int np, long n_dims);