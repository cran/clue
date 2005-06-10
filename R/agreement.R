### * cl_agreement

cl_agreement <-
function(x, y = NULL, method = "euclidean")
{
    ## <NOTE>
    ## This code is repeated from cl_dissimilarity(), mutatis mutandis.
    ## Not really a big surprise ...
    ## </NOTE>
    
    x <- as.cl_ensemble(x)
    is_partition_ensemble <- inherits(x, "cl_partition_ensemble")

    ## Be nice.
    if(is.character(y) || is.function(y)) {
        method <- y
        y <- NULL
    }

    if(!is.function(method)) {
        builtin_methods <- if(is_partition_ensemble)
            c("euclidean", "manhattan", "Rand", "cRand", "NMI", "KP",
              "angle", "diag")
        else
            c("euclidean", "manhattan", "cophenetic", "angle", "gamma")
        builtin_method_names <- if(is_partition_ensemble)
            c("minimal euclidean membership distances",
              "minimal manhattan membership distances",
              "Rand index",
              "corrected Rand index",
              "normalized mutual information",
              "Katz-Powell index",
              "maximal angle between memberships",
              "maximal co-classification rate")
        else
            c("euclidean ultrametric distances",
              "manhattan ultrametric distances",              
              "cophenetic correlations",
              "angle between ultrametrics",
              "rate of inversions")
        if(is.na(ind <- pmatch(tolower(method),
                               tolower(builtin_methods))))
            stop(gettextf("Value '%s' is not a valid abbreviation for an agreement method.",
                          method),
                 domain = NA)
        method <- paste(".cl_agreement",
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
        ## Build a cross-proximity object of cross-agreements.
        d <- matrix(0, length(x), length(y))
        for(j in seq(along = y))
            d[, j] <- sapply(x, method, y[[j]])
        dimnames(d) <-
            list(if(is.null(names(x))) seq(along = x) else names(x),
                 if(is.null(names(y))) seq(along = y) else names(y))
        description <- paste("Agreements using", method_name)
        return(cl_cross_proximity(d, description,
                                  class = "cl_cross_agreement"))
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

    ## <NOTE>
    ## We assume that self-agreements are always one ...
    ## </NOTE>
    cl_proximity(unlist(d),
                 paste("Agreements using", method_name),
                 labels = {
                     if(is.null(names(x)))
                         seq(along = x)
                     else
                         names(x)},
                 self = rep.int(1, length(x)),
                 size = n, class = "cl_agreement")
}

### ** .cl_agreement_partition_euclidean

.cl_agreement_partition_euclidean <-
function(x, y)
{
    ## <NOTE>
    ## Upper bound for maximal dissimilarity, maybe improve eventually.
    d_max <- sqrt(2 * n_of_objects(x))
    ## </NOTE>
    1 - .cl_dissimilarity_partition_euclidean(x, y) / d_max
}

### ** .cl_agreement_partition_manhattan

.cl_agreement_partition_manhattan <-
function(x, y)
{
    ## <NOTE>
    ## Upper bound for maximal dissimilarity, maybe improve eventually.
    d_max <- 2 * n_of_objects(x)
    ## </NOTE>
    1 - .cl_dissimilarity_partition_manhattan(x, y) / d_max
}
    
### ** .cl_agreement_partition_Rand

.cl_agreement_partition_Rand <-
function(x, y)
{
    n <- n_of_objects(x)
    
    ## Handle soft partitions using the corresponding hard ones.
    ## (At least, for the time being.)
    x <- table(cl_class_ids(x), cl_class_ids(y))

    ## <NOTE>
    ## The number A of concordant pairs is given by
    ##   A = choose(n,2) + \sum_{i,j} x_{ij}^2
    ##       - (1/2) * (\sum_i x_{i.}^2 + \sum_j x_{.j}^2)
    ##     = choose(n,2) + 2 \sum_{i,j} choose(x_{ij},2)
    ##       - (\sum_i choose(x_{i.},2) + \sum_j choose(x_{.j},2)
    ## with the first version certainly much faster to compute.
    ## </NOTE>
    1 + (sum(x^2) -
         (sum(rowSums(x)^2) + sum(colSums(x)^2)) / 2) / choose(n, 2)
}

### ** .cl_agreement_partition_cRand

.cl_agreement_partition_cRand <-
function(x, y)
{
    if(!is.cl_hard_partition(x) || !is.cl_hard_partition(y))
        stop("Can only handle hard partitions.")

    n <- n_of_objects(x)    
    x <- table(cl_class_ids(x), cl_class_ids(y))

    ## <NOTE>
    ## The basic formula is
    ##   (Sxy - E) / ((Sx. + S.y) / 2 - E)
    ## where
    ##   Sxy = \sum_{i,j} choose(x_{ij}, 2)
    ##   Sx. = \sum_i     choose(x_{i.}, 2)
    ##   S.y = \sum_j     choose(x_{.j}, 2)
    ## and
    ##   E = Sx. * S.y / choose(n, 2)
    ## We replace the bincoefs by the corresponding sums of squares,
    ## getting
    ##   (Txy - E) / ((Tx. + T.y) / 2 - E)
    ## where
    ##   Txy = \sum_{i,j} x_{ij}^2 - n
    ##   Tx. = \sum_i     x_{i.}^2 - n
    ##   T.y = \sum_j     x_{.j}^2 - n
    ## and
    ##   E = Tx. * T.y / (n^2 - n)
    ## </NOTE>
    Tx. <- sum(rowSums(x)^2)
    T.y <- sum(colSums(x)^2)
    E <- Tx. * T.y / (n^2 - n)
    (sum(x^2) - E) / ((Tx. + T.y) / 2 - E)
}


### ** .cl_agreement_partition_NMI

.cl_agreement_partition_NMI <-
function(x, y)
{
    if(!is.cl_hard_partition(x) || !is.cl_hard_partition(y))
        stop("Can only handle hard partitions.")
    x <- table(cl_class_ids(x), cl_class_ids(y))
    x <- x / sum(x)
    
    m_x <- rowSums(x)
    m_y <- colSums(x)
    y <- outer(m_x, m_y)

    i <- which((x > 0) & (y > 0))
    out <- sum(x[i] * log(x[i] / y[i]))
    e_x <- sum(m_x * log(ifelse(m_x > 0, m_x, 1)))
    e_y <- sum(m_y * log(ifelse(m_y > 0, m_y, 1)))

    out / sqrt(e_x * e_y)
}

### ** .cl_agreement_partition_KP

.cl_agreement_partition_KP <-
function(x, y)
{
    ## Agreement measure due to Katz & Powell (1953, Psychometrika), see
    ## also Messatfa (1992, Journal of Classification).

    n <- n_of_objects(x)    

    ## Handle soft partitions using the corresponding hard ones.
    ## (At least, for the time being.)
    x <- table(cl_class_ids(x), cl_class_ids(y))

    A_xy <- sum(x ^ 2)
    A_x. <- sum(rowSums(x) ^ 2)
    A_.y <- sum(colSums(x) ^ 2)

    (n^2 * A_xy - A_x. * A_.y) /
        sqrt(A_x. * (n^2 - A_x.) * A_.y * (n^2 - A_.y))
}

### ** .cl_agreement_partition_angle

.cl_agreement_partition_angle <-
function(x, y)
{
    ## Maximal angle between the matched memberships.

    k <- max(n_of_classes(x), n_of_classes(y))
    M_x <- cl_membership(x, k)
    M_y <- cl_membership(y, k)
    ## Match classes from conforming memberships.
    ind <- solve_LSAP(crossprod(M_x, M_y), max = TRUE)
    sum(M_x * M_y[, ind]) / sqrt(sum(M_x ^ 2) * sum(M_y ^ 2))
}

### ** .cl_agreement_partition_diag

.cl_agreement_partition_diag <-
function(x, y)
{
    ## Maximal co-classification rate.

    k <- max(n_of_classes(x), n_of_classes(y))
    M_x <- cl_membership(x, k)
    M_y <- cl_membership(y, k)
    ## Match classes from conforming memberships.
    ind <- solve_LSAP(crossprod(M_x, M_y), max = TRUE)
    sum(M_x * M_y[, ind]) / n_of_objects(x)
}

### ** .cl_agreement_hierarchy_euclidean

.cl_agreement_hierarchy_euclidean <-
function(x, y)
    1 / (1 + .cl_dissimilarity_hierarchy_euclidean(x, y))

### ** .cl_agreement_hierarchy_manhattan

.cl_agreement_hierarchy_manhattan <-
function(x, y)
    1 / (1 + .cl_dissimilarity_hierarchy_manhattan(x, y))

### ** .cl_agreement_hierarchy_cophenetic

.cl_agreement_hierarchy_cophenetic <-
function(x, y)
{
    ## Cophenetic correlation.
    cor(cl_ultrametric(x), cl_ultrametric(y))
}

### ** .cl_agreement_hierarchy_angle

.cl_agreement_hierarchy_angle <-
function(x, y)
{
    ## Angle between ultrametrics.
    u_x <- cl_ultrametric(x)
    u_y <- cl_ultrametric(y)
    sum(u_x * u_y) / sqrt(sum(u_x ^ 2) * sum(u_y ^ 2))
}

### ** .cl_agreement_hierarchy_gamma

.cl_agreement_hierarchy_gamma <-
function(x, y)
    1 - .cl_dissimilarity_hierarchy_gamma(x, y)
    
### * [.cl_agreement

"[.cl_agreement" <-
function(x, i, j)
{
    y <- NextMethod("[")
    if(!inherits(y, "cl_agreement")) {
        description <- attr(x, "description")
        return(cl_cross_proximity(y, description = description,
                                  class = "cl_cross_agreement"))
    }
    y
}


### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
