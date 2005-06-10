### * cl_consensus

cl_consensus <-
function(x, method = NULL, weights = 1, control = list())
{
    ## <NOTE>
    ## Interfaces are a matter of taste.
    ## E.g., one might want to have a 'type' argument indication whether
    ## hard or soft partitions are sought.  One could then do
    ##   cl_consensus(x, method = "euclidean", type = "hard")
    ## to look for an optimal median (or least squares) hard partition
    ## (for euclidean dissimilarity).
    ## For us, "method" really indicates a certain algorithm, with its
    ## bells and whistles accessed via the 'control' argument.
    ## </NOTE>

    clusterings <- as.cl_ensemble(x)

    if(!length(clusterings))
        stop("Cannot compute consensus of empty ensemble.")

    weights <- rep(weights, length = length(clusterings))
    if(any(weights < 0))
        stop("Argument 'weights' has negative elements.")
    if(!any(weights > 0))
        stop("Argument 'weights' has no positive elements.")

    is_partition_ensemble <-
        inherits(clusterings, "cl_partition_ensemble")

    if(!is.function(method)) {
        builtin_methods <-
            if(is_partition_ensemble)
                c("DWH", "GV1", "GV3")
            else
                c("cophenetic")
        if(is.null(method))
            ind <- 1
        else if(is.na(ind <- pmatch(tolower(method),
                                    tolower(builtin_methods))))
            stop(paste("Value", sQuote(method),
                       "is not a valid abbreviation",
                       "for a consensus method."))
        method <- get(paste(".cl_consensus",
                            if(is_partition_ensemble)
                            "partition"                        
                            else
                            "hierarchy",
                            builtin_methods[ind],
                            sep = "_"))
    }

    method(clusterings, weights, control)
}

### ** .cl_consensus_partition_DWH

.cl_consensus_partition_DWH <-
function(clusterings, weights, control)
{
    ## <TODO>
    ## Could make things more efficient by subscripting on positive
    ## weights.
    ## (Note that this means control$order has to be subscripted as
    ## well.)
    ## </TODO>
    
    max_n_of_classes <- max(sapply(clusterings, n_of_classes))

    ## Control parameters.
    k <- control$k
    if(is.null(k))
        k <- max_n_of_classes
    order <- control$order
    if(is.null(order))
        order <- sample(seq(along = clusterings))

    clusterings <- clusterings[order]
    weights <- weights[order]

    k_max <- max(k, max_n_of_classes)    
    s <- weights / cumsum(weights)
    
    M <- cl_membership(clusterings[[1]], k_max)
    for(b in seq(along = clusterings)[-1]) {
        mem <- cl_membership(clusterings[[b]], k_max)
        ## Match classes from conforming memberships.
        ind <- solve_LSAP(crossprod(M, mem), max = TRUE)
        M <- (1 - s[b]) * M + s[b] * mem[, ind]
        if(k < k_max)
            M <- .project_to_leading_columns(M, k)
    }

    cl_membership(as.cl_membership(M[, 1 : k]), k)    
}

### ** .cl_consensus_partition_GV1

.cl_consensus_partition_GV1 <-
function(clusterings, weights, control)
{
    ## <TODO>
    ## Could make things more efficient by subscripting on positive
    ## weights.
    ## </TODO>
    
    match_memberships <- function(M, N) {
        ## Return the M[, ind] column permutation of M optimally
        ## matching N.
        M[, solve_LSAP(crossprod(N, M), max = TRUE)]
    }

    value <- function(memberships, M, w) {
        sum(w * sapply(memberships, function(u) sum((u - M) ^ 2)))
    }

    B <- length(clusterings)    
    n <- n_of_objects(clusterings)
    max_n_of_classes <- max(sapply(clusterings, n_of_classes))

    ## Control parameters.
    k <- control$k
    if(is.null(k))
        k <- max_n_of_classes
    maxiter <- control$maxiter
    if(is.null(maxiter))
        maxiter <- 100
    reltol <- control$reltol
    if(is.null(reltol))
        reltol <- sqrt(.Machine$double.eps)
    start <- control$start
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")

    if(is.null(start)) {
        ## Random starting value.
        M <- matrix(runif(n * k), n, k)
        M <- M / rowSums(M)
    }
    else
        M <- start

    w <- weights / sum(weights)
    k_max <- max(k, max_n_of_classes)
    if(k_max > k)
        M <- cbind(M, matrix(0, nrow(M), k_max - k))
    
    ## This algorithm works by iteratively finding permutations P_i so
    ## that M_1 P_1 \approx M, ..., M_B P_B \approx M, and then updating
    ## M <- \sum_b w_b M_b P_b.

    ## As successive column permutations are column permutations of the
    ## original matrices, we work with the *current* M_b P_b.  We could
    ## use a 3-dimensional array so that M_b <=> memberships[ , , b],
    ## but prefer to work with a list of memberships.

    memberships <- lapply(clusterings, cl_membership, k_max)
    memberships <- lapply(memberships, match_memberships, M)
    old_value <- value(memberships, M, w)
    iter <- 1

    while(iter <= maxiter) {
        ## Update M as \sum w_b M_b P_b.
        M <- matrix(rowSums(mapply("*", memberships, w)), n)
        ## If k_max > k, "project" as indicated in Gordon & Vichi
        ## (2001), p. 238.
        if(k_max > k)
            M <- .project_to_leading_columns(M, k)
        ## Match the M_b P_b to M.
        memberships <- lapply(memberships, match_memberships, M)
        new_value <- value(memberships, M, w)
        if(verbose)
            cat("Iteration:", iter,
                "Old value:", old_value,
                "New value:", new_value,
                "\n")
        if(abs(old_value - new_value)
           < reltol * (abs(old_value) + reltol))
            break
        old_value <- new_value
        iter <- iter + 1
    }

    attr(M, "converged") <- (iter <= maxiter)
    rownames(M) <- rownames(memberships[[1]])
    cl_membership(as.cl_membership(M[, 1 : k]), k)
}

### ** .cl_consensus_partition_GV3

.cl_consensus_partition_GV3 <-
function(clusterings, weights, control)
{
    ## Use a SUMT similar to the one in ls_fit_ultrametric() to solve
    ##   \| Y - M M' \|_F^2 => min
    ## where M is a membership matrix and Y = \sum_b w_b M_b M_b'.

    B <- length(clusterings)    
    n <- n_of_objects(clusterings)
    max_n_of_classes <- max(sapply(clusterings, n_of_classes))

    ## Control parameters.
    eps <- control$eps
    if(is.null(eps))
        eps <- .Machine$double.eps
    k <- control$k
    if(is.null(k))
        k <- max_n_of_classes
    method <- control$method
    ## We use optim(method = "CG") as default, as nlm() may be
    ## computationally infeasible.
    if(is.null(method))
        method <- "CG"
    q <- control$q
    if(is.null(q))
        q <- 10
    start <- control$start
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")
    ## Do this at last ...
    control <- as.list(control$control)

    w <- weights / sum(weights)

    comemberships <-
        lapply(clusterings, function(x) {
            ## No need to force a common k here.
            M <- cl_membership(x)
            M %*% t(M)
        })

    Y <- matrix(rowSums(mapply("*", comemberships, w)), n)
    M <- if(is.null(start)) {
        e <- eigen(Y, symmetric = TRUE)
        ## Use M <- U_k \lambda_k^{1/2}.
        e$vectors[, 1:k] * rep(sqrt(e$values[1:k]), each = n)
    }
    else
        start
    y <- c(Y)

    L <- function(m) sum((y - crossprod(t(matrix(m, n)))) ^ 2)
    P <- function(m) {
        sum(pmin(m, 0) ^ 2) + sum((rowSums(matrix(m, n)) - 1) ^ 2)
    }
    Phi <- function(rho, m) L(m) + rho * P(m)
    grad_Phi <- function(rho, m) {
        M <- matrix(m, n)
        4 * c((crossprod(t(M)) - Y) %*% M) +
            2 * rho * (pmin(m, 0) + rep.int(rowSums(M) - 1, k))
    }
    
    make_Phi <- if(method == "nlm") {
        function(rho) {
            function(m) {
                y <- Phi(rho, m)
                attr(y, "gradient") <- grad_Phi(rho, m)
                y
            }
        }
    }
    else
        function(rho) {
            function(m)
                Phi(rho, m)
        }
    make_grad_Phi <- function(rho) { function(m) grad_Phi(rho, m) }

    optimize_with_penalty <- if(method == "nlm")
        function(rho, m)
            nlm(make_Phi(rho), m, check.analyticals = FALSE) $ estimate
    else
        function(rho, m)
            optim(m, make_Phi(rho), gr = make_grad_Phi(rho),
                  method = method, control = control) $ par

    m <- c(M)
    ## <TODO>
    ## Better upper/lower bounds for rho?
    rho <- max(L(m), 0.00001) / max(P(m), 0.00001)
    ## </TODO>
    iter <- 1
    repeat {
        if(verbose)
            cat("Iteration:", iter,
                "Rho:", rho,
                "P:", P(m),
                "\n")
        m_old <- m
        m <- optimize_with_penalty(rho, m)
        if(sum((m_old - m) ^ 2) < eps)
            break
        iter <- iter + 1        
        rho <- q * rho
    }

    ## Ensure that a stochastic matrix is returned.
    M <- matrix(pmax(m, 0), n)
    M <- M / rowSums(M)
    rownames(M) <- rownames(cl_membership(clusterings[[1]]))
    cl_membership(as.cl_membership(M), k)
}

### ** .cl_consensus_hierarchy_cophenetic

.cl_consensus_hierarchy_cophenetic <-
function(clusterings, weights, control)
{
    w <- weights / sum(weights)
    B <- length(clusterings)
    ultrametrics <- lapply(clusterings, cl_ultrametric)
    w <- rep.int(w, rep.int(length(ultrametrics[[1]]), B))
    dissimilarities <- rowSums(matrix(w * unlist(ultrametrics), nc = B))
    ## Cannot use as.cl_ultrametric() as we only have a dissimilarity
    ## (and are trying to fit an ultrametric to it) ...
    ## <FIXME 2.1.0>
    ## as.dist() is generic in R >= 2.1.0, maybe turn .dist_from_vector
    ## into an S3 method for it.
    ## Note that we need to ensure that labels get preserved ...
    ## (which is also why we use lapply() rather than sapply() above).
    labels <- attr(ultrametrics[[1]], "Labels")
    d <- clue:::.dist_from_vector(dissimilarities, labels = labels)
    ## </FIXME>
    ls_fit_ultrametric(d, control)
}

### * .project_to_leading_columns

.project_to_leading_columns <-
function(x, k)
{
    ## For a given matrix stochastic matrix x, return the stochastic
    ## matrix y which has columns from k+1 on all zero which is closest
    ## to x in the Frobenius distance.
    y <- x[, 1 : k]
    y <- cbind(pmax(y + (1 - rowSums(y)) / k, 0),
               matrix(0, nrow(y), ncol(x) - k))
    ## (Use the pmax to ensure that entries remain nonnegative.)
}

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
