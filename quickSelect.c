#include "quickSelect.h"

#define ELEM_SWAP(a,b) { register node_t t = (a); (a) = (b); (b) = t; }

medianValues quickSelect(node_t *ortho_points, long n) 
{
    medianValues result;
    long low, high ;
    long median;
    long middle, ll, hh;

    low = 0 ; high = n-1 ; median = (low + high) / 2;
    for (;;) {
        if (high <= low) { /* One element only */
            result.first = median;
            result.second = median;
            return result;
        }
        if (high == low + 1) {  /* Two elements only */
            result.first = low;
            result.second = high;
            return result;
        }

    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if ((ortho_points[middle].center[0] - ortho_points[high].center[0]) > 0)    ELEM_SWAP(ortho_points[middle], ortho_points[high]) ;
    if ((ortho_points[low].center[0] - ortho_points[high].center[0]) > 0)       ELEM_SWAP(ortho_points[low], ortho_points[high]) ;
    if ((ortho_points[middle].center[0] - ortho_points[low].center[0]) > 0)     ELEM_SWAP(ortho_points[middle], ortho_points[low]) ;

    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAP(ortho_points[middle], ortho_points[low+1]) ;

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
        do ll++; while ((ortho_points[low].center[0] - ortho_points[ll].center[0]) > 0) ;
        do hh--; while ((ortho_points[hh].center[0]  - ortho_points[low].center[0]) > 0) ;

        if (hh < ll)
        break;

        ELEM_SWAP(ortho_points[ll], ortho_points[hh]) ;
    }

    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP(ortho_points[low], ortho_points[hh]) ;

    /* Re-set active partition */
    if (hh <= median)
        low = ll;
        if (hh >= median)
        high = hh - 1;
    }
}