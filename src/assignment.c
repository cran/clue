#include "assignment.h"

void ap_print_solution(AP *p)
{
    int i;
  
    printf("%d itertations, %d secs.\n", p->runs, (int)p->rtime);
    printf("Min Cost: %10.4f\n", p->c);
  
    for(i = 0; i < p->n; i++)
	printf("%4d" ,p->s[i]);
    printf("\n");
}

void ap_reduce_cuts(AP *p)
{
    int i,j;
    double min;

    min = DBL_MAX;

    for(i = 0; i < p->n; i++)
	for(j = 0; j < p->n; j++)
	    if(p->ic[j] == 0 && p->ir[i] == 0 && p->t[i][j] < min)
		min = p->t[i][j];

    for(i = 0; i < p->n; i++)
	for(j = 0; j < p->n; j++){
	    if(p->ic[j] == 0 && p->ir[i] == 0)
		p->t[i][j]-= min;
	    if(p->ic[j] == 1 && p->ir[i] == 1)
		p->t[i][j]+= min;
	}
}

int ap_make_cuts(AP *p)
{
    int i, j;
    int ncut = 0;
    int **zeroes;
    int min;

    /* reset row and column indicators */
    for(i = 0; i < p->n; i++){
	p->ic[i] = 0;
	p->ir[i] = 0;
    }

    /* count zeroes */
    zeroes = (int **) malloc(p->n * sizeof(int *));
    for(i = 0; i < p->n; i++)
	zeroes[i] = (int *) calloc(2, sizeof(int));
    
    for(i = 0; i < p->n; i++)
	for(j = 0; j < p->n; j++)
	    if(p->ir[i] == 0 && p->ic[j] == 0 && p->t[i][j] == 0) {
		++zeroes[i][0];
		++zeroes[j][1];
	    }

    while(1){
	min = INT_MAX;
	p->zr = -1;
	p->zc = -1;
	p->ztype = NOTYPE;

	for(i = 0; i < p->n; i++) {
	    if(zeroes[i][0] > 0 && zeroes[i][0] < min && p->ir[i] == 0) {
		min = zeroes[i][0];
		p->zr = i;
		p->ztype = ROW;
	    }
	    if(zeroes[i][1] > 0 && zeroes[i][1] < min && p->ic[i] == 0) {
		min = zeroes[i][1];
		p->zc = i;
		p->ztype = COLUMN;
	    }
	}
    
	if(p->ztype == NOTYPE)
	    break;

	if(p->ztype == ROW) {
	    /* zeroes[p->zr][0] -= 1; */
	    for(i = 0; i < p->n; i++)
		if(p->t[p->zr][i] == 0 && p->ic[i] == 0) {
		    p->zc = i;
		    break;
		}
	} else {
	    /* zeroes[p->zc][1] -= 1; */
	    for(i = 0; i < p->n; i++)
		if(p->t[i][p->zc] == 0 && p->ir[i] == 0) {
		    p->zr = i;
		    break;
		}
	}

	/* apply cut */
	++ncut;
	if(p->ztype == ROW) {
	    /* printf("cut column %d\n",p->zc); */
	    p->ic[p->zc] = 1;
	    for(i = 0; i < p->n; i++)
		if(p->t[i][p->zc] == 0)
		    --zeroes[i][0];
	} else {
	    /* printf("cut row %d\n",p->zr); */
	    p->ir[p->zr] = 1;
	    for(i = 0; i < p->n; i++)
		if(p->t[p->zr][i] == 0)
		    --zeroes[i][1];
	}
    }

    /* free zeroes */
    for(i = 0; i < p->n; i++)
	free(zeroes[i]);
    free(zeroes);
    return ncut;
}

void ap_show_data(AP *p)
{
    int i,j;
    
    printf("      ");
    for(i = 0; i < p->n; i++)
	printf("%5d", p->ic[i]);
    printf("\n");
    for(i = 0; i < p->n; i++) {
	printf("%5d:", p->ir[i]);
	for(j = 0; j < p->n; j++)
	    printf("%5.2f", p->t[i][j]);
	printf("\n");
    }
}

void ap_hungarian(AP *p)
{
    int c;
    time_t start, end;

    start = time(0);
    ap_preproc(p);

    while((c = ap_make_cuts(p)) < p->n) {
	++p->runs;
	ap_reduce_cuts(p);
    }

    /* construct solution */
    ap_generate_solution(p);
    end = time(0);
    p->rtime = end - start;
    return;
}

void ap_preproc(AP *p)
{
    int i,j;
    double min;
    
    /* first reduce row minima */
    for(i = 0; i < p->n; i++){
	min = p->t[i][0];
	for(j = 1; j < p->n; j++)
	    min = (p->t[i][j] < min) ? p->t[i][j] : min;
	if(min > 0)
	    for(j = 0; j < p->n; j++)
		p->t[i][j]-= min;
    }
    
    /* do same with columns */
    for(i = 0; i < p->n; i++){
	min = p->t[0][i];
	for(j = 1; j < p->n; j++)
	    min = (p->t[j][i] < min) ? p->t[j][i] : min;
	if(min > 0)
	    for(j = 0; j < p->n; j++)
		p->t[j][i]-= min;
    }
}

AP *ap_read_problem(char *file)
{
    FILE *f;
    int i,j,c;
    int m,n;
    double x;
    double **t;
    int nrow,ncol;
    AP *p;

    f = fopen(file,"r");
    if(f == NULL)
	return NULL;

    t = (double **)malloc(sizeof(double*));

    m = 0; 
    n = 0;  

    nrow = 0;
    ncol = 0;
    
    while(EOF != (i = fscanf(f, "%lf", &x))) {
	if(i == 1){
	    if(n == 0){
		t = (double **) realloc(t,(m + 1) * sizeof(double *));
		t[m] = (double *) malloc(sizeof(double));
	    } else
		t[m] = (double *) realloc(t[m], (n + 1) * sizeof(double));

	    t[m][n++] = x;

	    ncol = (ncol < n) ? n : ncol;
	    c = fgetc(f);
	    if(c == '\n'){
		n = 0;
		++m;
		nrow = (nrow < m) ? m : nrow;
	    }
	}
    }
    fclose(f);

    /* prepare data */

    if(nrow != ncol) {
	printf("error: problem not quadratic\nrows = %d, cols = %d\n",
	       nrow, ncol);
	exit(1);
    }

    p = (AP*) malloc(sizeof(AP)); 
    p->n = ncol;
    p->runs = 0;
    
    p->t  = (double **) malloc((nrow)*sizeof(double *));
    p->tc = (double **) malloc((nrow)*sizeof(double *));
    if(p->t == NULL || p->tc == NULL)
	return NULL;

    for(i = 0; i < nrow; i++){
	p->t[i] = (double *) calloc(ncol, sizeof(double));
	p->tc[i] = (double *) calloc(ncol, sizeof(double));
	if(p->t[i] == NULL || p->tc[i] == NULL)
	    return NULL;
    }

    for(i = 0; i < nrow; i++)
	for( j = 0; j < ncol; j++){
	    p->t[i][j] = t[i][j];
	    p->tc[i][j] = t[i][j];
	}
    for(i = 0; i < nrow; i++)
	free(t[i]);
    free(t);
    
    /* prepare ir,ic: */
    p->ic = (int *) calloc(nrow,sizeof(int));
    p->ir = (int *) calloc(nrow,sizeof(int));
    if(p->ic == NULL || p->ir == NULL)
	return NULL;
    
    /* set cost to maximum to indicate that the problem has not been
       solved */
    p->c = DBL_MAX;
    p->s = NULL;
    
    p->rtime = 0;
    return p;
}

AP *ap_create_problem(double **t, int n)
{
    int i,j;
    AP *p;

    p = (AP*) malloc(sizeof(AP)); 
    if(p == NULL)
	return NULL;
    
    p->n = n;
    p->runs = 0;
    
    p->t  = (double **) malloc(n * sizeof(double *));
    p->tc = (double **) malloc(n * sizeof(double *));
    if(p->t == NULL || p->tc == NULL)
	return NULL;
    
    for(i = 0; i < n; i++){
	p->t[i] = (double *) calloc(n, sizeof(double));
	p->tc[i] = (double *) calloc(n, sizeof(double));
	if(p->t[i] == NULL || p->tc[i] == NULL)
	    return NULL;
    }

    for(i = 0; i < n; i++)
	for( j = 0; j < n; j++){
	    p->t[i][j] = t[i][j];
	    p->tc[i][j] = t[i][j];
	}  if(t == NULL)
	    return NULL;

  
    /* prepare ir,ic: */
    p->ic = (int *) calloc(n, sizeof(int));
    p->ir = (int *) calloc(n, sizeof(int));
    
    if(p->ic == NULL || p->ir == NULL)
	return NULL;
    
    /* set cost to maximum to indicate that the problem has not been
       solved */
    p->c = DBL_MAX;
    p->s = NULL;
    
    p->rtime = 0;
    return p;
}

int ap_assignment(AP *p, int *res)
{
    int i;
    
    if(p->s == NULL)
	ap_hungarian(p);
    
    for(i = 0; i < p->n; i++)
	res[i] = p->s[i];
    
    return p->n;
}

double ap_mincost(AP *p)
{
    if(p->s == NULL)
	ap_hungarian(p);
    
    return p->c;
}

int ap_size(AP *p)
{
    return p->n;
}

int ap_time(AP *p)
{
    return (int) p->rtime; 
}

int ap_iterations(AP *p)
{
    return p->runs; 
}

int ap_costmatrix(AP *p, double **m)
{
  int i,j;
  
  for(i = 0; i < p->n; i++)
    for(j = 0; j < p->n; j++)
      m[i][j] = p->tc[i][j];

  return p->n;
}

int ap_datamatrix(AP *p, double **m)
{
    int i,j;
  
    for(i = 0; i < p->n; i++)
	for(j = 0; j < p->n; j++)
	    m[i][j] = p->t[i][j];
    
    return p->n;
}

int ap_generate_solution(AP *p)
{
    int i, j;
    int ncut = 0;
    int **zeroes;
    int min;
    
    p->s = (int *) calloc(p->n, sizeof(int));
    
    /* reset row and column indicators */
    for(i = 0; i < p->n; i++) {
	p->ic[i] = 0;
	p->ir[i] = 0;
    }

    /* count zeroes */
    zeroes = (int **) malloc(p->n * sizeof(int *));
    for(i = 0; i < p->n; i++)
	zeroes[i] = (int *) calloc(2,sizeof(int));

    for(i = 0; i < p->n; i++)
	for(j = 0; j < p->n; j++)
	    if(p->ir[i] == 0 && p->ic[j] == 0 && p->t[i][j] == 0) {
		++zeroes[i][0];
		++zeroes[j][1];
	    }

    p->c = 0;
    while(1) {
	min = INT_MAX;
	p->zr = -1;
	p->zc = -1;
	p->ztype = NOTYPE;

	for(i = 0; i < p->n; i++){
	    if(zeroes[i][0] > 0 && zeroes[i][0] < min && p->ir[i] == 0) {
		min = zeroes[i][0];
		p->zr = i;
		p->ztype = ROW;
	    }
	    if(zeroes[i][1] > 0 && zeroes[i][1] < min && p->ic[i] == 0) {
		min = zeroes[i][1];
		p->zc = i;
		p->ztype = COLUMN;
	    }
	}
	
	if(p->ztype == NOTYPE)
	    break;

	if(p->ztype == ROW){
	    for(i = 0; i < p->n; i++)
		if(p->t[p->zr][i] == 0 && p->ic[i] == 0){
		    p->zc = i;
		    break;
		}
	} else {
	    for(i = 0; i < p->n; i++)
		if(p->t[i][p->zc] == 0 && p->ir[i] == 0){
		    p->zr = i;
		    break;
		}
	}
	
	/*apply cut */
	++ncut;
	p->ic[p->zc] = 1;
	p->ir[p->zr] = 1;
	for(i = 0; i < p->n; i++)
	    if(p->t[i][p->zc] == 0)
		--zeroes[i][0];
	for(i = 0; i < p->n; i++)
	    if(p->t[p->zr][i] == 0)
		--zeroes[i][1];
	p->s[p->zr] = p->zc;
	p->c+= p->tc[p->zr][p->zc];
    }

    /* free zeroes */
    for(i = 0; i < p->n; i++)
	free(zeroes[i]);
    free(zeroes);
    return ncut;
}

int ap_free(AP *p)
{
    int i;
    
    for(i = 0; i < p->n; i++){
	free(p->t[i]);
	free(p->tc[i]);
    }
    free(p->t);
    free(p->tc);
    
    free(p->ir);
    free(p->ic);
    free(p->s);
    free(p);
    return 0;
}
