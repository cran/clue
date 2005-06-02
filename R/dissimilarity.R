### * cl_dissimilarity

cl_dissimilarity <-
function(x, y = NULL, method = "euclidean")
{
    x <- as.cl_ensemble(x)
    is_partition_ensemble <- inherits(x, "cl_partition_ensemble")

    ## Be nice.
    if(is.character(y) || is.function(y)) {
        method <- y
        y <- NULL
    }
    
    if(!is.function(method)) {
        builtin_methods <- if(is_partition_ensemble)
            c("euclidean", "comemberships")
        else
            c("euclidean", "cophenetic", "gamma")
        builtin_method_names <- if(is_partition_ensemble)
            c("minimal euclidean membership distances",
              "euclidean comembership distances")
        else
            c("euclidean ultrametric distances",
              "cophenetic correlations",
              "rate of inversions")
        if(is.na(ind <- pmatch(tolower(method),
                               tolower(builtin_methods))))

            stop(gettextf("Value '%s' is not a valid abbreviation for a dissimilarity method.",
                          method),
                 domain = NA)
        method <- paste(".cl_dissimilarity",
                        if(is_partition_ensemble)
                        "partition"                        
                        else
                        "hierarchy",
                        builtin_methods[ind],
                        sep = "_")
        method_name <- builtin_method_names[ind]
    }
    else {
        method_name <- "user-defined method"
    }
              
    if(!is.null(y)) {
        y <- as.cl_ensemble(y)
        if(inherits(y, "cl_partition_ensemble")
           != is_partition_ensemble)
            stop("Cannot mix partitions and hierarchies.")
        if(n_of_objects(x) != n_of_objects(y))
            stop("All clusterings must have the same number of objects.")
        ## Build a cross-proximity object of cross-dissimilarities.
        d <- matrix(0, length(x), length(y))
        for(j in seq(along = y))
            d[, j] <- sapply(x, method, y[[j]])
        dimnames(d) <-
            list(if(is.null(names(x))) seq(along = x) else names(x),
                 if(is.null(names(y))) seq(along = y) else names(y))
        description <- paste("Dissimilarities using", method_name)
        return(cl_cross_proximity(d, description,
                                  class = "cl_cross_dissimilarity"))
    }
    
    ## Otherwise, build a proximity object of dissimilarities.
    n <- length(x)
    d <- vector("list", length = n - 1)
    ind <- seq(length = n)
    while(length(ind) > 1) {
        j <- ind[1]
        ind <- ind[-1]
        d[[j]] <- sapply(x[ind], method, x[[j]])
    }
    
    cl_proximity(unlist(d),
                 paste("Dissimilarities using", method_name),
                 labels = {
                     if(is.null(names(x)))
                         seq(along = x)
                     else
                         names(x)},
                 size = n,
                 class = c("cl_dissimilarity", "cl_proximity", "dist"))
}

### ** .cl_dissimilarity_partition_euclidean

.cl_dissimilarity_partition_euclidean <-
function(x, y)
{
    k <- max(n_of_classes(x), n_of_classes(y))
    M_x <- cl_membership(x, k)
    M_y <- cl_membership(y, k)
    ## Match classes from conforming memberships.
    ind <- solve_LSAP(crossprod(M_x, M_y), max = TRUE)
    sum((M_x - M_y[, ind]) ^ 2) / (2 * n_of_objects(x))
}

### ** .cl_dissimilarity_partition_comemberships

.cl_dissimilarity_partition_comemberships <-
function(x, y)
{
    ## We used to have the straightforward
    ##   C_x <- crossprod(t(cl_membership(x))) # M_x M_x'
    ##   C_y <- crossprod(t(cl_membership(y))) # M_y M_y'
    ##   sum((C_x - C_y) ^ 2) / n_of_objects(x) ^ 2
    ## But note that
    ##   \| AA' - BB' \|^2
    ##      = tr((AA' - BB')'(AA' - BB')
    ##      = tr(A'A A'A) - 2 tr(A'B B'A) + tr(B'B B'B)
    ##      = \| A'A \|^2 - 2 \| A'B \|^2 + \| B'B \|^2
    ## which can be computed much more efficiently as all involved cross
    ## product matrices are "small" ...
    k <- max(n_of_classes(x), n_of_classes(y))
    M_x <- cl_membership(x, k)
    M_y <- cl_membership(y, k)
    (sum(crossprod(M_x) ^ 2)
     - 2 * sum(crossprod(M_x, M_y) ^ 2)
     + sum(crossprod(M_y) ^ 2)) / (n_of_objects(x) ^ 2)
}

### ** .cl_dissimilarity_hierarchy_euclidean

.cl_dissimilarity_hierarchy_euclidean <-
function(x, y)
{
    mean((cl_ultrametric(x) - cl_ultrametric(y)) ^ 2)
}

### ** .cl_dissimilarity_hierarchy_cophenetic

.cl_dissimilarity_hierarchy_cophenetic <-
function(x, y)
    1 - cor(cl_ultrametric(x), cl_ultrametric(y)) ^ 2

### ** .cl_dissimilarity_hierarchy_gamma

.cl_dissimilarity_hierarchy_gamma <-
function(x, y)
{
    ## <NOTE>
    ## This is a dissimilarity measure that works for arbitrary
    ## dissimilarities, see e.g. Bock.
    u <- cl_ultrametric(x)
    v <- cl_ultrametric(y)
    n <- length(u)
    .C("clue_dissimilarity_count_inversions",
       u, v, n, count = double(1),
       PACKAGE = "clue") $ count / choose(n, 2)
}

### * as.dist.cl_dissimilarity

as.dist.cl_dissimilarity <-
function(m, diag = FALSE, upper = FALSE)
{
    y <- c(m)
    ## Fill non-inherited attributes with default values.
    attributes(y) <- c(attributes(m)[c("Size", "Labels")],
                       Diag = diag, Upper = upper, call = match.call())
    ## (Note that as.dist.default() does not automatically add
    ## 'method'.)
    class(y) <- "dist"
    y
}

### * [.cl_dissimilarity

"[.cl_dissimilarity" <-
function(x, i, j)
{
    y <- NextMethod("[")
    if(!inherits(y, "cl_dissimilarity")) {
        description <- attr(x, "description")
        return(cl_cross_proximity(y, description = description,
                                  class = "cl_cross_dissimilarity"))
    }
    y
}

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
