cl_prototypes <-
function(x)
    UseMethod("cl_prototypes")

## No default method.

## Package stats: kmeans() (R 2.1.0 or better).
cl_prototypes.kmeans <-
function(x)
    x$centers

## Package cluster: clara() always gives prototypes.
cl_prototypes.clara <-
function(x)
    x$medoids
## Package cluster: fanny() never gives prototypes.
## Package cluster: pam() does not give prototypes if given a
## dissimilarity matrix.
cl_prototypes.pam <-
function(x)
{
    p <- x$medoids
    if(!is.matrix(p))
        stop("Cannot determine prototypes.")
    p
}

## Package cclust: cclust().
cl_prototypes.cclust <- cl_prototypes.kmeans

## Package e1071: cmeans() gives objects of class "fclust".
cl_prototypes.fclust <- cl_prototypes.kmeans
## Package e1071: cshell().
cl_prototypes.cshell <- cl_prototypes.kmeans
## Package e1071: bclust().
cl_prototypes.bclust <- cl_prototypes.kmeans

## Package mclust: Mclust().
cl_prototypes.Mclust <-
function(x)
{
    p <- x$mu
    ## For multidimensional models, we get a matrix whose columns are
    ## the means of each group in the best model, and hence needs to be
    ## transposed.
    if(is.matrix(p))
        p <- t(p)
    p
}

## Package clue: cl_pclust().
cl_prototypes.cl_pclust <-
function(x)
    x$prototypes
