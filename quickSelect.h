typedef struct _node {
  long point_id;
  double *center;
  double radius;
  struct _node *L;
  struct _node *R;
} node_t;

typedef struct {
  long first, second;
} medianValues;

medianValues quickSelect(node_t *ortho_points, long n);