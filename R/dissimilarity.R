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

    if(is.function(method))
        method_name <- "user-defined method"
    else {
        if(!inherits(method, "cl_dissimilarity_method")) {
            ## Get the method definition and description from the
            ## registry. 
            type <- ifelse(is_partition_ensemble,
                           "partition", "hierarchy")
            method <- get_cl_dissimilarity_method(method, type)
        }
        method_name <- method$description
        method <- method$definition
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
        dimnames(d) <- list(names(x), names(y))
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
                 labels = names(x),
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
    sqrt(sum((M_x - M_y[, ind]) ^ 2))
}

### ### ** .cl_dissimilarity_partition_manhattan

.cl_dissimilarity_partition_manhattan <-
function(x, y)
{
    k <- max(n_of_classes(x), n_of_classes(y))
    M_x <- cl_membership(x, k)
    M_y <- cl_membership(y, k)
    C <- .cxdist(M_x, M_y, "manhattan")
    ind <- solve_LSAP(C)
    sum(C[cbind(seq(along = ind), ind)])
}    

### ** .cl_dissimilarity_partition_comemberships

.cl_dissimilarity_partition_comemberships <-
function(x, y)
{
    ## We used to have the straightforward
    ##   C_x <- tcrossprod(cl_membership(x)) # M_x M_x'
    ##   C_y <- tcrossprod(cl_membership(y)) # M_y M_y'
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
    sqrt(sum(crossprod(M_x) ^ 2)
         - 2 * sum(crossprod(M_x, M_y) ^ 2)
         + sum(crossprod(M_y) ^ 2))
}

### ** .cl_dissimilarity_partition_symdiff

.cl_dissimilarity_partition_symdiff <-
function(x, y)
{
    ## Cardinality of the symmetric difference of the partitions
    ## regarded as binary equivalence relations, i.e., the number of
    ## discordant pairs.

    ## Handle soft partitions using the corresponding hard ones.
    ## (At least, for the time being.)

    n <- n_of_objects(x)
    .cl_dissimilarity_partition_Rand(x, y) * choose(n, 2)
}
    
### ** .cl_dissimilarity_partition_Rand

.cl_dissimilarity_partition_Rand <-
function(x, y)    
{
    ## Handle soft partitions using the corresponding hard ones.
    ## (At least, for the time being.)

    1 - .cl_agreement_partition_Rand(x, y)
}

### ** .cl_dissimilarity_partition_GV1

.cl_dissimilarity_partition_GV1 <-
function(x, y)
{    
    k_x <- n_of_classes(x)
    k_y <- n_of_classes(y)
    M_x <- cl_membership(x, k_x)
    M_y <- cl_membership(y, k_y)
    C <- outer(colSums(M_x ^ 2), colSums(M_y ^ 2), "+") -
        2 * crossprod(M_x, M_y)
    if(k_x < k_y)
        C <- rbind(C, matrix(0, nr = k_y - k_x, nc = k_y))
    else if(k_x > k_y)
        C <- cbind(C, matrix(0, nr = k_x, nc = k_x - k_y))
    ind <- solve_LSAP(C)
    sqrt(sum(C[cbind(seq(along = ind), ind)]))
    ## (Note that this sum really only includes matched non-dummy
    ## classes.)
}

### ** .cl_dissimilarity_partition_BA_A

.cl_dissimilarity_partition_BA_A <-
function(x, y)
{
    .cl_dissimilarity_partition_manhattan(as.cl_hard_partition(x),
                                          as.cl_hard_partition(y)) / 2
    ## Could to this more efficiently, of course ...
}

### ** .cl_dissimilarity_partition_BA_C

.cl_dissimilarity_partition_BA_C <-
function(x, y)
{
    n_of_classes(x) + n_of_classes(y) - 2 * n_of_classes(cl_join(x, y))
}

### ** .cl_dissimilarity_partition_BA_D

.cl_dissimilarity_partition_BA_D <-
    .cl_dissimilarity_partition_Rand

### ** .cl_dissimilarity_partition_BA_E

.cl_dissimilarity_partition_BA_E <-
function(x, y)
{
    z <- table(cl_class_ids(x), cl_class_ids(y))
    z <- z / sum(z)

    ## Average mutual information between the partitions.
    y <- outer(rowSums(z), colSums(z))
    i <- which((z > 0) & (y > 0))
    I <- sum(z[i] * log(z[i] / y[i]))

    ## Entropy of meet(x, y).
    i <- which(z > 0)
    H <- - sum(z[i] * log(z[i]))

    1 - I / H
}

### ** .cl_dissimilarity_hierarchy_euclidean

.cl_dissimilarity_hierarchy_euclidean <-
function(x, y)
{
    if(!.has_object_dissimilarities(x) ||
       !.has_object_dissimilarities(y))
        return(NA)
    u <- cl_object_dissimilarities(x)
    v <- cl_object_dissimilarities(y)
    sqrt(sum((u - v) ^ 2))
}

### ** .cl_dissimilarity_hierarchy_manhattan

.cl_dissimilarity_hierarchy_manhattan <-
function(x, y)
{
    if(!.has_object_dissimilarities(x) ||
       !.has_object_dissimilarities(y))
        return(NA)
    u <- cl_object_dissimilarities(x)
    v <- cl_object_dissimilarities(y)
    sum(abs(u - v))
}

### ** .cl_dissimilarity_hierarchy_cophenetic

.cl_dissimilarity_hierarchy_cophenetic <-
function(x, y)
{
    if(!.has_object_dissimilarities(x) ||
       !.has_object_dissimilarities(y))
        return(NA)
    u <- cl_object_dissimilarities(x)
    v <- cl_object_dissimilarities(y)
    1 - cor(u, v) ^ 2
}

### ** .cl_dissimilarity_hierarchy_gamma

.cl_dissimilarity_hierarchy_gamma <-
function(x, y)
{
    ## <NOTE>
    ## This is a dissimilarity measure that works for arbitrary
    ## dissimilarities, see e.g. Bock.
    ## (And the current implementation finally respects this ...)
    ## </NOTE>
    if(!.has_object_dissimilarities(x) ||
       !.has_object_dissimilarities(y))
        return(NA)
    u <- cl_object_dissimilarities(x)
    v <- cl_object_dissimilarities(y)
    n <- length(u)
    .C("clue_dissimilarity_count_inversions",
       u, v, n, count = double(1),
       PACKAGE = "clue") $ count / choose(n, 2)
}

### ** .cl_dissimilarity_hierarchy_symdiff

.cl_dissimilarity_hierarchy_symdiff <-
function(x, y)
{
    ## Cardinality of the symmetric difference of the n-trees when
    ## regarded as sets of subsets (classes) of the set of objects.

    n <- n_of_objects(x)
    x <- cl_classes(x)
    y <- cl_classes(y)
    lx <- sapply(x, length)
    ly <- sapply(y, length)
    s <- 0
    for(i in seq(length = n)) {
        sx <- x[lx == i]
        sy <- y[lx == i]
        s <- s + sum(is.na(match(sx, sy))) + sum(is.na(match(sy, sx)))
    }
    s
}

### ** .cl_dissimilarity_hierarchy_Chebyshev

.cl_dissimilarity_hierarchy_Chebyshev <-
function(x, y)
{
    if(!.has_object_dissimilarities(x) ||
       !.has_object_dissimilarities(y))
        return(NA)
    u <- cl_object_dissimilarities(x)
    v <- cl_object_dissimilarities(y)
    max(abs(u - v))
}

### ** .cl_dissimilarity_hierarchy_Lyapunov

.cl_dissimilarity_hierarchy_Lyapunov <-
function(x, y)
{
    if(!.has_object_dissimilarities(x) ||
       !.has_object_dissimilarities(y))
        return(NA)
    q <- cl_object_dissimilarities(x) / cl_object_dissimilarities(y)
    if(is.matrix(q)) q <- q[lower.tri(q)]
    log(max(q) / min(q))
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

### .cxdist

.cxdist <-
function(A, B, method = "manhattan")
{
    ## Return the column cross distance matrix of A and B.
    ## I.e., the matrix C = [c_{j,k}] with
    ##   c_{j,k} = distance(A[, j], B[, k])
    ## Currently, only Manhattan (L1) distances are provided.
    ## Extensions to Minkowski or even more distances (a la dist())
    ## could be added eventually.

    ## <NOTE>
    ## Possible implementations include
    ##
    ## foo_a <- function(A, B)
    ##   apply(B, 2, function(u) colSums(abs(A - u)))
    ## foo_d <- function(A, B) {
    ##   out <- as.matrix(dist(rbind(t(A), t(B)), "manhattan"))
    ##   dimnames(out) <- NULL
    ##   nc_B <- NCOL(B)
    ##   out[seq(from = NCOL(A) + 1, length = nc_B), seq(length = nc_B)]
    ## }
    ## foo_f <- function(A, B) {
    ##   out <- matrix(0, NCOL(A), NCOL(B))
    ##   for(j in seq(length = NCOL(A)))
    ##     for(k in seq(length = NCOL(B)))
    ##       out[j, k] = sum(abs(A[, j] - B[, k]))
    ##   out
    ## }
    ##
    ## The one actually used seems to be the best performer, with the
    ## "for" version a close second (note that "typically", A and B have
    ## much fewer columns than rows).
    ## only few columns 
    
    out <- matrix(0, NCOL(A), NCOL(B))
    for(k in seq(length = NCOL(B)))
        out[, k] <- colSums(abs(A - B[, k]))
    out
}

### .rxdist

.rxdist <-
function(A, B, method = c("euclidean", "manhattan"))
{
    ## Return the row cross distance matrix of A and B.
    ## I.e., the matrix C = [c_{j,k}] with
    ##   c_{j,k} = distance(A[j, ], B[k, ])
    
    ## <NOTE>
    ## Could also do something like
    ##   ind <- seq(length = NROW(B))
    ##   as.matrix(dist(rbind(B, A)))[-ind, ind]
    ## but that is *very* inefficient for the "usual" data by prototype
    ## case (where NROW(B) << NROW(A)).
    ## </NOTE>
    
    ## No fancy pmatching for methods for the time being.
    method <- match.arg(method)
    
    ## Workhorse: Full A, single row of b.
    FOO <- if(method == "euclidean")
        function(A, b) sqrt(rowSums(sweep(A, 2, b) ^ 2))
    else
        function(A, b) rowSums(abs(sweep(A, 2, b)))
        
    out <- matrix(0, NROW(A), NROW(B))
    for(k in seq(length = NROW(B)))
        out[, k] <- FOO(A, B[k, ])
    out
}

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
