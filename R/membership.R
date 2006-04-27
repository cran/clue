### * cl_membership

## Get the class membership matrix from a partition.

## <NOTE>
## We could use sparse matrices for the memberships of hard partitions.
## Not sure if this is really that important, though, as we typically
## use memberships in a context where dense matrices (memberships of
## soft partitions) occur.
## </NOTE>

## <NOTE>
## Currently, the number of classes to be used for the memberships must
## not be less than the number of classes in the partition.  We might
## eventually change this so that "optimal" collapsing of classes is
## performed (but note that optimality needs to be relative to some
## dissimilarity measure) ...
## However, from the discussion of the second method in Gordon and Vichi
## (2001) we note that whereas optimal assignment is "simple", optimal
## collapsing (equivalent to partitioning into an arbitrary number of
## partitions) is of course very hard.
## </NOTE>

cl_membership <-
function(x, k = n_of_classes(x))
{
    if(k < n_of_classes(x))
        stop("k cannot be less than the number of classes in x.")
    UseMethod("cl_membership")
}

## Default method.
cl_membership.default <-
function(x, k = n_of_classes(x))
    .cl_membership_from_class_ids(cl_class_ids(x), k)

## Package stats: kmeans() (R 2.1.0 or better).
cl_membership.kmeans <- cl_membership.default

## Package cluster: clara(), fanny(), and pam() give objects of the
## respective class inheriting from class "partition".
cl_membership.fanny <-
function(x, k = n_of_classes(x))
    .cl_membership_from_memberships(x$membership, k)
cl_membership.partition <- cl_membership.default

## Package cclust: cclust().
cl_membership.cclust <- cl_membership.default

## Package e1071: cmeans() gives objects of class "fclust".
cl_membership.fclust <- cl_membership.fanny
## Package e1071: cshell().
cl_membership.cshell <- cl_membership.fanny
## Package e1071: bclust().
cl_membership.bclust <- cl_membership.default

## Package flexmix: class "flexmix".
cl_membership.flexmix <-
function(x, k = n_of_classes(x))
    .cl_membership_from_memberships(flexmix::posterior(x), k)

## Package mclust: Mclust().
cl_membership.Mclust <-
function(x, k = n_of_classes(x))
    .cl_membership_from_memberships(x$z, k)

## Package clue: Memberships.
cl_membership.cl_membership <-
function(x, k = n_of_classes(x))
    .cl_membership_from_memberships(x, k)
## (Note: we cannot simply return x in case k equals n_of_classes(x),
## because ncol(x) might be different.)

## Package clue: cl_pclust().
cl_membership.cl_pclust <- cl_membership.fanny
## Package clue: (virtual) class "cl_partition".
cl_membership.cl_partition <-
function(x, k = n_of_classes(x))
    cl_membership(.get_representation(x), k)

### * .cl_membership_from_class_ids

.cl_membership_from_class_ids <-
function(x, k = NULL)
{
    x <- factor(x)
    n_of_objects <- length(x)
    n_of_classes <- nlevels(x)
    if(is.null(k))
        k <- n_of_classes
    else if(k < n_of_classes)
        stop("k cannot be less than the number of classes in x.")
    ## <TODO>
    ## Should really use a sparse encoding of this ...
    M <- matrix(0, n_of_objects, k)
    ## (Could also use .one_entry_per_column(M, as.numeric(x)) <- 1 for
    ## the time being.)
    M[cbind(seq(length = n_of_objects), as.numeric(x))] <- 1
    ## </TODO>
    if(nlevels(x) == k)
        colnames(M) <- levels(x)
    if(!is.null(nm <- names(x)))
        rownames(M) <- nm
    attr(M, "n_of_classes") <- n_of_classes
    attr(M, "is_cl_hard_partition") <- TRUE
    class(M) <- "cl_membership"
    M
}

### * .cl_membership_from_memberships

.cl_membership_from_memberships <-
function(x, k = NULL)
{
    n_of_objects <- nrow(x)
    x <- x[ , colSums(x) > 0, drop = FALSE]
    n_of_classes <- ncol(x)
    if(!is.null(k)) {
        if(k < n_of_classes)
            stop("k cannot be less than the number of classes in x.")
        if(k > n_of_classes) {
            ## Fill up with zero columns.
            x <- cbind(x, matrix(0, nrow(x), k - n_of_classes))
        }
    }
    attr(x, "n_of_classes") <- n_of_classes
    attr(x, "is_cl_hard_partition") <- all(rowSums(x == 1))
    class(x) <- "cl_membership"
    x
}

### * as.cl_membership

as.cl_membership <-
function(x)
    UseMethod("as.cl_membership")
as.cl_membership.default <-
function(x)
{
    if(inherits(x, "cl_membership"))
        x
    else if(is.atomic(x))
        .cl_membership_from_class_ids(x)
    else
        cl_membership(x)
}
as.cl_membership.matrix <-
function(x)
    .cl_membership_from_memberships(x)

### * .memberships_from_cross_dissimilarities

.memberships_from_cross_dissimilarities <-
function(d, power = 2)
{
    ## For a given matrix of cross-dissimilarities [d_{bj}], return a
    ## matrix [u_{bj}] such that \sum_{b,j} u_{bj}^p d_{bj}^q => min!
    ## under the constraint that u is a stochastic matrix.
    ## If only one power is given, it is taken as p, with q as 1.
    ## <NOTE>
    ## This returns a plain matrix of membership values and not a
    ## cl_membership object (so that it does not deal with possibly
    ## dropping or re-introducing unused classes).
    ## </NOTE>
    exponent <- if(length(power) == 1)
        1 / (power - 1)
    else
        power[2] / (power[1] - 1)
    u <- matrix(0, nrow(d), ncol(d))
    FUN <- function(s, t) (s / t) ^ exponent
    zero_incidences <- (d == 0)
    n_of_zeroes <- rowSums(zero_incidences)
    if(any(ind <- (n_of_zeroes > 0)))
        u[ind, ] <- zero_incidences[ind, ] / n_of_zeroes[ind]
    if(any(!ind)) {
        u[!ind, ] <-
            t(apply(d[!ind, ], 1,
                    function(s) 1 / rowSums(outer(s, s, FUN))))
    }
    u
}

### * print.cl_membership

print.cl_membership <-
function(x, ...)
{
    writeLines("Memberships:")
    print(matrix(as.vector(x), nr = nrow(x), dimnames = dimnames(x)),
          ...)
    invisible(x)
}

### .has_object_memberships

## Be nice to users when computing proximities: all measures for
## "partitions" we currently consider really only assume that we can
## compute memberships and/or class ids.

## Note that the cl_membership() default method works for cl_class_ids.

.has_object_memberships <-
function(x)
    (is.cl_partition(x)
     || inherits(x, "cl_membership")
     || inherits(x, "cl_class_ids"))


### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
