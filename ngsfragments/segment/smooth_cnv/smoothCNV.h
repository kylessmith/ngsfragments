#include <stdio.h>
#include <stdlib.h>

#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))


static int cmpfunc(const void * a, const void * b);
void smoothRL(int n, double *gdat, int nchr, long *cfrq, double *sgdat, int k, double oSD, double sSD);