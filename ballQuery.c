#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct _node {
    double radius;
    long L;
    long R;
} node_t;

int n_dims;
node_t *tree;
double **center;
double *point;
double minDist = 1000000.0;
long currBest;


void allocate_tree(long n_nodes)
{
    double *_p_center;

    tree = (node_t *) malloc(n_nodes * sizeof(node_t));
    _p_center = (double *) malloc(n_nodes * n_dims * sizeof(double));
    center = (double **) malloc(n_nodes * sizeof(double *));
    if ((_p_center == NULL) || (center == NULL) || (tree == NULL)){
        printf("Error allocating tree, exiting.\n");
        exit(10);
    }
    for (long i = 0; i < n_nodes; i++)
        center[i] = &(_p_center[i * n_dims]);
}


double distance(double *pt1, double *pt2)
{
    double dist = 0.0;

    for (int d = 0; d < n_dims; d++)
        dist += (pt1[d] - pt2[d]) * (pt1[d] - pt2[d]);
    return sqrt(dist);
}


void search_tree(long idx)
{
    double dist;

    if (tree[idx].radius == 0.0){   // found leave
        dist = distance(center[idx], point);
        if(dist < minDist){
            minDist = dist;
            currBest = idx;
        }
        return;
    }

    if (distance(center[tree[idx].L], point) - tree[idx].radius < minDist)
        search_tree(tree[idx].L);
    if (distance(center[tree[idx].R], point) - tree[idx].radius < minDist)
        search_tree(tree[idx].R);
}

int main(int argc, char *argv[])
{
    FILE *fp;
    node_t *node;
    long i,
         node_idx,
         n_nodes;
    int d;

    if (argc < 3){
        printf("Usage: %s <ball-tree-file> <point>\n", argv[0]);
        exit(1);
    }

    fp = fopen(argv[1], "r");
    if (fp == NULL){
        printf("Cannot open input file '%s'.\n", argv[1]);
        exit(2);
    }

    fscanf(fp, "%d %ld", &n_dims, &n_nodes);
    if (n_dims < 2){
        printf("Illegal number of dimensions (%d), must be above 1.\n", n_dims);
        exit(3);
    }
    if (n_nodes < 2){
        printf("Illegal number of nodes (%ld), must be above 1.\n", n_nodes);
        exit(2);
    }

    if (argc != n_dims + 2){
        printf("Wrong number of coordinates for <point>\n");
        exit(3);
    }

    point = (double *) malloc(n_dims * sizeof(double));
    if (point == NULL){
        printf("Error allocating point, exiting.\n");
        exit(4);
    }
    for (d = 0; d < n_dims; d++)
        point[d] = atof(argv[d + 2]);

    allocate_tree(n_nodes);
    for (i = 0; i < n_nodes; i++){        
        fscanf(fp, "%ld", &node_idx);
        node = &(tree[node_idx]);
        fscanf(fp, "%ld %ld %lf", &(node->L), &(node->R), &(node->radius));
        for (d = 0; d < n_dims; d++)
            fscanf(fp, "%lf", &(center[node_idx][d]));
    }

    // tree and point are global, index 0 is root; currBest has result
    search_tree(0);
    
    // print point
    for (d = 0; d < n_dims; d++)
        printf("%lf ", center[currBest][d]);
    printf("\n");
}