#include <stdio.h>
#include <stdlib.h>
#include <limits.h>		/* INT_MAX */
#include <float.h>		/* DBL_MAX */
#include <time.h>		/* time_t */

#define ROW             0
#define COLUMN          1
#define NOTYPE         -1

typedef struct {
    int n;		/* order of problem */
    double **t;		/* cost matrix */
    double **tc;	/* copy of cost matrix */
    int ztype;		/* ROW | COLUMN with minimum number of zeroes */
    int zr;		/* row index */
    int zc;		/* column index */
    int *ir, *ic;	/* cut row and column indicator */
    int *s;		/* solution */
    int runs;		/* number of iterations */
    double c;		/* minimum cost */
    time_t rtime;	/* time needed to solve */
} AP;

AP    *ap_create_problem(double **t, int n);
int    ap_assignment(AP *p, int *res);
double ap_mincost(AP *p);
void   ap_hungarian(AP *p);
int    ap_make_cuts(AP *p);
void   ap_preproc(AP *p);
void   ap_print_solution(AP *p);
AP    *ap_read_problem(char *file);
void   ap_reduce_cuts(AP *p);
void   ap_show_data(AP *p);

int    ap_iterations(AP *p);
int    ap_size(AP *p);
int    ap_time(AP *p);
int    ap_costmatrix(AP *p, double **m);
int    ap_datamatrix(AP *p, double **m);
int    ap_generate_solution(AP *p);

int    ap_free(AP *p);
