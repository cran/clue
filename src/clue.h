#ifndef _CLUE_H
#define _CLUE_H

#include <R.h>

void solve_LSAP(double *c, Sint *n, Sint *p);

double **clue_vector_to_square_matrix(double *x, Sint n);

void clue_ultrametric_penalty_function(double *x, Sint *n, double *v);
void clue_ultrametric_penalty_gradient(double *x, Sint *n, double *out);

void clue_dissimilarity_count_inversions(double *x, double *y, Sint *n,
					 double *count);

#endif
