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
        if(!inherits(method, "cl_consensus_method")) {
            ## Get the method definition from the registry.
            type <- ifelse(is_partition_ensemble,
                           "partition", "hierarchy")
            if(is.null(method)) {
                ## Defaults: change to SBH (Soft Boehm-Hornik)
                ## eventually. 
                method <- ifelse(is_partition_ensemble,
                                 "DWH", "cophenetic")
            }
            method <- get_cl_consensus_method(method, type)
        }
        method <- method$definition
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
    s[is.na(s)] <- 0                    # Division by zero ...
    
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

### ** .cl_consensus_partition_AO

.cl_consensus_partition_AO <-
function(clusterings, weights, control, type = c("soft", "hard"))
{
    ## The start of a general purpose optimizer for determining
    ## consensus partitions by minimizing
    ##   \sum_b d(M, M_b) ^ p
    ##     = \sum_b \min_{P_b} f(M, M_b P_b) ^ p
    ## for dissimilarity measures d involving explicitly matching
    ## labels via permutation matrices P_b.  The AO ("alternative
    ## optimization") proceeds by alternatively by matching the M_b to
    ## M by minimizing f(M, M_b P_b) over P_b, and fitting M by
    ## minimizing \sum_b f(M, M_b P_b) ^ p for fixed matchings.
    ##
    ## Such a procedure requires three ingredients: a function for
    ## matching M_b to M (in fact simply replacing M_b by the matched
    ## M_b P_b); a function for fitting M to the \{M_b P_b\}, and a
    ## function for computing the value of the criterion function
    ## corresponding to this fit (so that one can stop if the relative
    ## improvement is small enough).
    ##
    ## For the time being, we only use this to determine soft and hard
    ## Euclidean least squares consensus partitions (soft and hard
    ## Euclidean means), so the interface does not yet reflect the
    ## generality of the approach (which would either pass the three
    ## functions, or even set up family objects encapsulating the three
    ## functions).

    ## <TODO>
    ## Could make things more efficient by subscripting on positive
    ## weights.
    ## </TODO>

    ## For the time being ...
    type <- match.arg(type)
    
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

    ## <NOTE>
    ## As successive column permutations are column permutations of the
    ## original matrices, we work with the *current* M_b P_b.  We could
    ## use a 3-dimensional array so that M_b <=> memberships[ , , b],
    ## but prefer to work with a list of memberships.
    ## </NOTE>

    ## Function for fitting M to (fixed) memberships M_b P_b.
    ## <NOTE>
    ## In general, this clearly needs to know how many classes/columns
    ## to use, and hence needs an explicit argument 'k'.  Both 'n' and
    ## 'k_max' can always be determined from the memberships ensemble,
    ## but could also be passed on explicitly for efficiency reasons.
    ## </NOTE>
    if(type == "soft")
        fit_memberships <- function(memberships) {
            ## Update M as \sum w_b M_b P_b.
            M <- matrix(rowSums(mapply("*", memberships, w)), n)
            ## If k_max > k, "project" as indicated in Gordon & Vichi
            ## (2001), p. 238.
            if(k_max > k)
                M <- .project_to_leading_columns(M, k)
            M
        }
    else {
        ## For hard partitions, we currently cannot handle the case
        ## where k < k_max.
        if(k < k_max)
            stop("Currently cannot compute soft means for reduced 'k'.")
        fit_memberships <- function(memberships) {
            ## Compute M as \sum w_b M_b P_b.
            M <- matrix(rowSums(mapply("*", memberships, w)), n)
            ## And compute a closest hard partition H(M) from that.
            cl_membership(as.cl_membership(max.col(M)), k)
        }
    }
            
    memberships <- lapply(clusterings, cl_membership, k_max)
    memberships <- lapply(memberships, match_memberships, M)
    old_value <- value(memberships, M, w)
    iter <- 1

    while(iter <= maxiter) {
        ## Fit M to the M_b P_b.
        M <- fit_memberships(memberships)
        ## Match the M_b P_b to M.
        memberships <- lapply(memberships, match_memberships, M)
        ## Update value.
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
    

### ** .cl_consensus_partition_GV1

.cl_consensus_partition_GV1 <-
function(clusterings, weights, control)
    .cl_consensus_partition_AO(clusterings, weights, control, "soft")

### ** .cl_consensus_partition_HBH

.cl_consensus_partition_HBH <-
function(clusterings, weights, control)
    .cl_consensus_partition_AO(clusterings, weights, control, "hard")

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
    d <- .dist_from_vector(dissimilarities, labels = labels)
    ## </FIXME>
    ls_fit_ultrametric(d, control)
}

### ** .cl_consensus_hierarchy_majority

.cl_consensus_hierarchy_majority <-
function(clusterings, weights, control)
{
    ## Have no use for control arguments.

    w <- weights / sum(weights)

    classes <- lapply(clusterings, .get_classes_in_hierarchy)
    all_classes <- unique(unlist(classes, recursive = FALSE))
    gamma <- double(length = length(all_classes))
    for(i in seq(along = classes))
        gamma <- gamma + w[i] * !is.na(match(all_classes, classes[[i]]))

    maj_classes <- all_classes[gamma > 1 / 2]
    attr(maj_classes, "labels") <- attr(classes[[1]], "labels")

    .cl_ultrametric_from_classes(maj_classes)
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
