#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct _node {
    double radius;
    long id;
    long L;
    long R;
} node_t;

typedef struct _hash {
    long id;
    long index;
    struct _hash *next;
} hash_t;

int n_dims;
long n_nodes;
double *point;

node_t *tree;
double **center;
hash_t **hash;

long currBest;
double minDist = 1000000.0;


void allocate_hash()
{
    hash = (hash_t **) calloc(n_nodes, sizeof(hash_t *));
    if(hash == NULL){
        printf("Error allocating hash, exiting.\n");
        exit(20);
    }
}

void hash_insert(long tree_id, long array_idx)
{
    hash_t *new;
    long i;

    i = tree_id % n_nodes;

    new = (hash_t *) malloc(sizeof(hash_t));
    new->id = tree_id;
    new->index = array_idx;
    new->next = hash[i];
    hash[i] = new;
}

long hash_get_index(long id)
{
    long i;
    hash_t *h;

    i = id % n_nodes;

    for(h = hash[i]; h != NULL && h->id != id; h = h->next);

    if(h == NULL){
	printf("Id %d not found?!\n", id);
	exit(30);
    }

    return h->index;
}
    
void allocate_tree()
{
    double *_p_center;

    tree = (node_t *) malloc(n_nodes * sizeof(node_t));
    _p_center = (double *) malloc(n_nodes * n_dims * sizeof(double));
    center = (double **) malloc(n_nodes * sizeof(double *));
    if((_p_center == NULL) || (center == NULL) || (tree == NULL)){
        printf("Error allocating tree, exiting.\n");
        exit(10);
    }
    for(long i = 0; i < n_nodes; i++)
        center[i] = &(_p_center[i * n_dims]);
}


double distance(double *pt1, double *pt2)
{
    double dist = 0.0;

    for(int d = 0; d < n_dims; d++)
        dist += (pt1[d] - pt2[d]) * (pt1[d] - pt2[d]);
    return sqrt(dist);
}


void search_tree(long idx)
{
    double dist;
    long idxl, idxr;

    if(tree[idx].radius == 0.0){   // found leave
        dist = distance(center[idx], point);
        if(dist < minDist){
            minDist = dist;
            currBest = idx;
        }
        return;
    }
    
    idxl = hash_get_index(tree[idx].L);
    if(distance(center[idxl], point) - tree[idx].radius < minDist)
        search_tree(idxl);
    idxr = hash_get_index(tree[idx].R);
    if(distance(center[idxr], point) - tree[idx].radius < minDist)
        search_tree(idxr);
}

int main(int argc, char *argv[])
{
    FILE *fp;
    node_t *node;
    long i,
	 node_idx;
    int d;

    if(argc < 3){
        printf("Usage: %s <ball-tree-file> <point>\n", argv[0]);
        exit(1);
    }

    fp = fopen(argv[1], "r");
    if(fp == NULL){
        printf("Cannot open input file '%s'.\n", argv[1]);
        exit(2);
    }

    fscanf(fp, "%d %ld", &n_dims, &n_nodes);
    if(n_dims < 2){
        printf("Illegal number of dimensions (%d), must be above 1.\n", n_dims);
        exit(3);
    }
    if(n_nodes < 2){
        printf("Illegal number of nodes (%ld), must be above 1.\n", n_nodes);
        exit(2);
    }

    if(argc != n_dims + 2){
        printf("Wrong number of coordinates for <point>\n");
        exit(3);
    }

    point = (double *) malloc(n_dims * sizeof(double));
    if(point == NULL){
        printf("Error allocating point, exiting.\n");
        exit(4);
    }
    for(d = 0; d < n_dims; d++)
        point[d] = atof(argv[d + 2]);

    allocate_hash();
    allocate_tree();
    for(i = 0; i < n_nodes; i++){
        fscanf(fp, "%ld", &node_idx);
	hash_insert(node_idx, i);
        node = &(tree[i]);
        fscanf(fp, "%ld %ld %lf", &(node->L), &(node->R), &(node->radius));
        for(d = 0; d < n_dims; d++)
            fscanf(fp, "%lf", &(center[i][d]));
    }

    // tree and point are global, index 0 is root; currBest has result
    search_tree(hash_get_index(0));
    
    // print closest sample
    for(d = 0; d < n_dims; d++)
        printf("%lf ", center[currBest][d]);
    printf("\n");
}
