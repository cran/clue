#include <R.h>
#include "assignment.h"
#include "clue.h"

void
solve_LSAP(double *c, Sint *n, Sint *p)
{
    AP *ap;
    double **data;

    data = clue_vector_to_square_matrix(c, *n);
    ap = ap_create_problem(data, *n);
    ap_hungarian(ap);
    ap_assignment(ap, p);
    ap_free(ap);
}

/*
void
solve_LSAP(double *c, Sint *n, Sint *p)
{
    AP *ap;
    double **data;
    int i, j;

    data = (double **) R_alloc(*n, sizeof(double));
    for(i = 0; i < *n; i++) {
	data[i] = (double *) R_alloc(*n, sizeof(double));
	for(j = 0; j < *n; j++)
	    data[i][j] = c[i + *n * j];
    }
    ap = ap_create_problem(data, *n);
    ap_hungarian(ap);
    ap_assignment(ap, p);
    ap_free(ap);
}
*/
