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
    if(is.cl_partition(x)) {
        v <- .cl_validity_partition_d_a_f(cl_membership(x),
                                          as.matrix(d))
        out <- list("Dissimilarity accounted for" = v)
    }
    else if(is.cl_hierarchy(x)) {
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
    ## Internal function for computing 1 - a / mean(d), where the
    ## "average within dissimilarity" a is given by
    ##   \frac{\sum_{i,j} \sum_k m_{ik}m_{jk} d(i,j)}
    ##        {\sum_{i,j} \sum_k m_{ik}m_{jk}}
    ## where m is the membership matrix and d a *symmetric* matrix of
    ## dissimilarities.
    within_sums <-
        rowSums(sapply(seq(length = ncol(m)),
                       function(k) {
                           z <- m[, k]
                           w <- outer(z, z, "*")
                           c(sum(w * d), sum(w))
                       }))
    average_within_d <- within_sums[1] / within_sums[2]
    1 - average_within_d / mean(d)
}
