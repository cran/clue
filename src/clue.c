#include <R.h>
#include <Rmath.h>

#include "clue.h"

double **clue_vector_to_square_matrix(double *x, Sint n)
{
    double **data, *val;
    Sint i, j;
    data = (double **) R_alloc(n, sizeof(double));
    for(i = 0; i < n; i++) {
        data[i] = (double *) R_alloc(n, sizeof(double));
        val = x + i;
        for(j = 0; j < n; j++, val += n)
            data[i][j] = *val;
    }
    return(data);
}

void
clue_ultrametric_penalty_function(double *x, Sint *n, double *v)
{
    double **data, p, d;
    Sint i, j, k;

    data = clue_vector_to_square_matrix(x, *n);

    p = 0;
    for(i = 0; i < *n; i++)
	for(j = i + 1; j < *n; j++)
	    for(k = 0; k < *n; k++) {
		if(data[i][j] <= fmin2(data[i][k], data[j][k])) {
		    d = data[i][k] - data[j][k];
		    p += d * d;
		}
	    }

    *v = 2 * p;
}

void
clue_ultrametric_penalty_gradient(double *x, Sint *n, double *out)
{
    
    double **data, v;
    double d_i_j, d_i_k, d_j_k;
    Sint i, j, k;

    data = clue_vector_to_square_matrix(x, *n);    

    for(i = 0; i < *n; i++)
	for(j = 0; j < *n; j++) {
	    d_i_j = data[i][j];
	    v = 0;
	    for(k = 0; k < *n; k++) {
		d_i_k = data[i][k];
		d_j_k = data[j][k];
		if(d_i_k <= fmin2(d_i_j, d_j_k))
		    v += d_i_j - d_j_k;
		if(d_j_k <= fmin2(d_i_j, d_i_k))
		    v += d_i_j - d_i_k;
	    }
	    out[*n * i + j] = 4 * v;
	}

}

static int clue_sign(double x)
{
    if(x == 0) return(0);
    return((x > 0) ? 1 : -1);
}

void
clue_dissimilarity_count_inversions(double *x, double *y, Sint *n,
				    double *count)
{
    Sint i, j;
    for(i = 0; i < *n; i++)
	for(j = 0; j < *n; j++)
	    if((clue_sign(x[i] - x[j]) * clue_sign(y[i] - y[j])) < 0)
		(*count)++;
}
