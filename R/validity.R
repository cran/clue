## A slightly polymorphic function, similar to cluster::silhouette() and
## its methods.

cl_validity <-
function(x, ...)
    UseMethod("cl_validity")

cl_validity.default <-
function(x, d, ...)
{
    ## Note that providing methods for classes "cl_partition" and
    ## "cl_hierarchy" is not good enough ...
    out <- list()
    if(inherits(x, "cl_partition")) {
        v <- .cl_validity_partition_d_a_f(cl_membership(x),
                                          as.matrix(d))
        out <- list("Dissimilarity accounted for" = v)
    }
    else if(inherits(x, "cl_hierarchy")) {
        ## Currently nothing.
        ## Consider adding e.g. the Agglomerative Coefficient or
        ## Divisive Coeffcient for more than cluster::agnes() and
        ## cluster::diana(), respectively.
    }
    class(out) <- "cl_validity"
    out
}

## Package cluster: agnes().
cl_validity.agnes <-
function(x, ...)
{
    out <- list("Agglomerative coefficient" = x$ac)
    class(out) <- "cl_validity"
    out
}
## Package cluster: diana().
cl_validity.diana <-
function(x, ...)
{
    out <- list("Divisive coefficient" = x$dc)
    class(out) <- "cl_validity"
    out
}

cl_validity.cl_pclust <-
function(x, ...)
    x$validity

print.cl_validity <-
function(x, ...)
{
    for(nm in names(x))
        cat(nm, ": ", x[[nm]], "\n", sep = "")
    invisible(x)
}

.cl_validity_partition_d_a_f <-
function(m, d)
{
    ## "Dissimilarity accounted for".
    ## Internal function for computing
    ##   \sum_{i,j} \sum_k m_{ik}m_{jk} d(i,j) / \sum_{i,j} d(i,j)
    ## where m is the membership matrix and d a *symmetric* matrix of
    ## dissimilarities.
    within <- sum(sapply(seq(length = ncol(m)),
                         function(k) {
                             z <- m[, k]
                             sum(outer(z, z, "*") * d)
                         }))
    1 - within / sum(d)
}
