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

### * .cl_consensus_partition_DWH

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

### * .cl_consensus_partition_AOS

.cl_consensus_partition_AOS <-
function(clusterings, weights, control, type = c("SE", "HE"))
{
    ## The start of a general purpose optimizer for determining
    ## consensus partitions by minimizing
    ##   \sum_b w_b d(M, M_b) ^ e
    ##     = \sum_b \min_{P_b} w_b f(M, M_b P_b) ^ e
    ## for the special case where the criterion function is based on
    ## M and M_b P_b (i.e., column permutations of M_b), as opposed to
    ## the general case where d(M, M_b) = \min_{P_b} f(M, P_b, M_b)
    ## handled by .cl_consensus_partition_AOG().
    ##
    ## The AO ("alternative optimization") proceeds by alternatively
    ## matching the M_b to M by minimizing f(M, M_b P_b) over P_b, and
    ## fitting M by minimizing \sum_b w_b f(M, M_b P_b) ^ e for fixed
    ## matchings.
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
    ##
    ## This special case is provided for efficiency and convenience.
    ## Using the special form of the criterion function, we can simply
    ## always work memberships with the same maximal number of columns,
    ## and with the permuted \{ M_b P_b \}.

    ## For the time being ...
    type <- match.arg(type)

    w <- weights / sum(weights)    
    n <- n_of_objects(clusterings)
    k_max <- max(sapply(clusterings, n_of_classes))

    ## Control parameters.
    k <- control$k
    if(is.null(k))
        k <- k_max
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

    ## The maximal (possible) number of classes in M and the \{ M_b \}.
    k_all <- max(k, k_max)

    if(k < k_all)
        M <- cbind(M, matrix(0, nrow(M), k_all - k))

    ## Currently, only Euclidean dissimilarities ...
    value <- function(M, memberships, w) {
        sum(w * sapply(memberships, function(u) sum((u - M) ^ 2)))
    }
    match_memberships <- function(M, N) {
        ## Return the M[, ind] column permutation of M optimally
        ## matching N.
        M[, solve_LSAP(crossprod(N, M), max = TRUE)]
    }
    ## Function for fitting M to (fixed) memberships \{ M_b P_b \}.
    ## As we use a common number of columns for all membership matrices
    ## involved, we need to pass the desired 'k' ...
    if(type == "SE")
        fit_M <- function(memberships, w, k) {
            ## Update M as \sum w_b M_b P_b.
            M <- matrix(rowSums(mapply("*", memberships, w)), nrow(M))
            ## If k < k_all, "project" as indicated in Gordon & Vichi
            ## (2001), p. 238.
            if(k < ncol(M))
                M <- .project_to_leading_columns(M, k)
            M
        }
    else if(type == "HE")
        fit_M <- function(memberships, w, k) {
            ## Compute M as \sum w_b M_b P_b.
            M <- matrix(rowSums(mapply("*", memberships, w)), nrow(M))
            ## And compute a closest hard partition H(M) from that,
            ## using the first k columns of M.
            cl_membership(as.cl_membership(max.col(M[ , 1 : k])),
                          ncol(M))
            ## <FIXME>
            ## Perhaps more efficiently:
            ##   .cl_membership_from_class_ids(max.col(M[ , 1 : k]),
            ##                                 ncol(M))
            ## </FIXME>
        }

    memberships <- lapply(clusterings, cl_membership, k_all)
    memberships <- lapply(memberships, match_memberships, M)
    old_value <- value(M, memberships, w)
    iter <- 1

    while(iter <= maxiter) {
        ## Fit M to the M_b P_b.
        M <- fit_M(memberships, w, k)
        ## Match the \{ M_b P_b \} to M.
        memberships <- lapply(memberships, match_memberships, M)
        ## Update value.
        new_value <- value(M, memberships, w)
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

    rownames(M) <- rownames(memberships[[1]])    
    M <- cl_membership(as.cl_membership(M[, 1 : k]), k)

    ## Add these attributes here, as the above would not preserve them.
    attr(M, "converged") <- (iter <= maxiter)
    attr(M, "value") <- new_value

    M
}

### ** .cl_consensus_partition_SE

.cl_consensus_partition_SE <-
function(clusterings, weights, control)
    .cl_consensus_partition_AOS(clusterings, weights, control, "SE")

### ** .cl_consensus_partition_HE

.cl_consensus_partition_HE <-
function(clusterings, weights, control)
    .cl_consensus_partition_AOS(clusterings, weights, control, "HE")

### * .cl_consensus_partition_AOG

.cl_consensus_partition_AOG <-
function(clusterings, weights, control, type = c("GV1"))
{
    ## The start of a general purpose optimizer for determining
    ## consensus partitions by minimizing
    ##   \sum_b w_b d(M, M_b) ^ p
    ##     = \sum_b \min_{P_b} w_b f(M, M_b, P_b) ^ e
    ## for general dissimilarity matrices which involve class matching
    ## via permutation matrices P_b.
    ##
    ## The AO ("Alternative Optimization") proceeds by alternating
    ## between determining the optimal permutations P_b by minimizing
    ##   f(M, M_b, P_b)
    ## for fixed M, and fitting M by minimizing
    ##   \sum_b w_b f(M, M_b, P_b) ^ e
    ## for fixed \{ P_b \}.
    ##
    ## We encapsulate this into functions fit_P() and fit_M() (and a
    ## value() function for the criterion function to be minimized with
    ## respect to both M and \{ P_b \}, even though the current
    ## interface does not yet reflect the generality of the approach.
    ##
    ## Note that rather than passing on information about the numbers of
    ## classes (e.g., needed for GV1) and representing all involved
    ## membership matrices with the same maximal number of columns, we
    ## use "minimal" representations with no dummy classes (strictly
    ## speaking, with the possible exception of M, for which the given k
    ## is used).

    ## For the time being ...
    type <- match.arg(type)

    w <- weights / sum(weights)    
    n <- n_of_objects(clusterings)
    k_max <- max(sapply(clusterings, n_of_classes))
    
    ## Control parameters.
    k <- control$k
    if(is.null(k))
        k <- k_max
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

    ## <NOTE>
    ## For the given memberships, we can simply use ncol() in the
    ## computations (rather than n_of_classes(), because we used
    ## cl_membership() to create them.  For M, the number of classes
    ## could be smaller than the given k "target".
    ## </NOTE>
    
    value <- function(M, permutations, memberships, w) {

        k <- .n_of_nonzero_columns(M)
        d <- function(u, p) {
            ## Compute the squared GV1 dissimilarity between M and u
            ## based on the M->u class matching p.
            nc_u <- ncol(u)
            if(nc_u == k) {
                ## Simple case: all classes are matched.
                sum((u[, p] - M) ^ 2)
            }
            else {
                ## Only include the matched non-dummy classes of M ..
                ind <- 1 : k
                ## ... which are matched to non-dummy classes of u.
                ind <- ind[p[ind] <= nc_u]
                sum((u[, p[ind]] - M[, ind]) ^ 2)
            }
        }
        
        sum(w * mapply(d, memberships, permutations))
    }

    fit_P <- function(u, M) {
        
        ## Return a permutation representing a GV1 optimal matching of
        ## the columns of M to the columns of u (note the order of the
        ## arguments), using a minimal number of dummy classes (i.e., p
        ## has max(.n_of_nonzero_columns(M), n_of_classes(u)) entries).

        ## See also .cl_dissimilarity_partition_GV1().
        
        C <- outer(colSums(M ^ 2), colSums(u ^ 2), "+") -
            2 * crossprod(M, u)
        nc_M <- .n_of_nonzero_columns(M)
        nc_u <- ncol(u)
        ## (See above for ncol() vs n_of_classes().)
        if(nc_M < nc_u)
            C <- rbind(C, matrix(0, nr = nc_u - nc_M, nc = nc_u))
        else if(nc_M > nc_u)
            C <- cbind(C, matrix(0, nr = nc_M, nc = nc_M - nc_u))
        
        solve_LSAP(C)
    }

    fit_M <- function(permutations, memberships, w) {
        
        ## Here comes the trickiest part ...
        ##
        ## In general, M = [m_{iq}] is determined as follows.
        ## Write value(M, permutations, memberships, w) as
        ##   \sum_b \sum_i \sum_{p=1}^{k_b} \sum_{q=1}^k
        ##      w_b (u_{ip}(b) - m_{iq})^2 x_{pq}(b)
        ## where U(b) and X(b) are the b-th membership matrix and the
        ## permutation matrix representing the M->U(b) non-dummy class
        ## matching (as always, note the order of the arguments).
        ##
        ## Let
        ##   \beta_{iq} = \sum_b \sum_{p=1}^{k_b} w_b u_{ip}(b) x_{pq}(b)
        ##   \alpha_q   = \sum_b \sum_{p=1}^{k_b} w_b x_{pq}(b)
        ## and
        ##   \bar{m}_{iq} =
        ##     \cases{\beta_{iq}/\alpha_q, & $\alpha_q > 0$ \cr
        ##            0                    & otherwise}.
        ## Then, as the cross-product terms cancel out, the value
        ## function rewrites as
        ##   \sum_b \sum_i \sum_{p=1}^{k_b} \sum_{q=1}^k
        ##      w_b (u_{ip}(b) - \bar{m}_{iq})^2 x_{pq}(b)
        ##   + \sum_i \sum_q \alpha_q (\bar{m}_{iq} - m_{iq}) ^ 2,
        ## where the first term is a constant, and the minimum is found
        ## by solving
        ##   \sum_q \alpha_q (\bar{m}_{iq} - m_{iq}) ^ 2 => min!
        ## s.t.
        ##   m_{i1}, ..., m_{ik} >= 0, \sum_{iq} m_{iq} = 1.
        ##
        ## We can distinguish three cases.
        ## A. If S_i = \sum_q \bar{m}_{iq} = 1, things are trivial.
        ## B. If S_i = \sum_q \bar{m}_{iq} < 1.
        ##    B1. If some \alpha_q are zero, then we can choose
        ##        m_{iq} = \bar{m}_{iq} for those q with \alpha_q = 0;
        ##        m_{iq} = 1 / number of zero \alpha's, otherwise.
        ##    B2. If all \alpha_q are positive, we can simply
        ##        equidistribute 1 - S_i over all classes as written
        ##        in G&V.
        ## C. If S_i > 1, things are not so clear (as equidistributing
        ##    will typically result in violations of the non-negativity
        ##    constraint).  We currently revert to using solve.QP() from
        ##    package quadprog, as constrOptim() already failed in very
        ##    simple test cases.
        ##
        ## Now consider \sum_{p=1}^{k_b} x_{pq}(b).  If k <= k_b for all
        ## b, all M classes from 1 to k are matched to one of the k_b
        ## classes in U(b), hence the sum and also \alpha_q are one.
        ## But then
        ##    \sum_q \bar{m}_{iq}
        ##      =  \sum_b \sum_{p=1}^{k_b} w_b u_{ip}(b) x_{pq}(b)
        ##      <= \sum_b \sum_{p=1}^{k_b} w_b u_{ip}(b)
        ##      =  1
        ## with equality if k = k_b for all b.  I.e.,
        ## * If k = \min_b k_b = \max k_b, we are in case A.
        ## * If k <= \min_b k_b, we are in case B2.
        ## And it makes sense to handle these cases explicitly for
        ## efficiency reasons.

        ## And now for something completely different ... the code.

        k <- .n_of_nonzero_columns(M)
        nr_M <- nrow(M)
        nc_M <- ncol(M)
        nc_memberships <- sapply(memberships, ncol)

        if(k <= min(nc_memberships)) {
            ## Compute the weighted means \bar{M}.
            M <- matrix(rowSums(mapply("*",
                                       mapply(function(u, p)
                                              u[ , p[1 : k]],
                                              memberships,
                                              permutations,
                                              SIMPLIFY = FALSE),
                                       w)),
                        nr_M)
            ## And add dummy classes if necessary.
            if(k < nc_M)
                M <- cbind(M, matrix(0, nr_M, nc_M - k))
            ## If we always got the same number of classes, we are
            ## done.  Otherwise, equidistribute ...
            if(k < max(nc_memberships))
                M <- pmax(M + (1 - rowSums(M)) / nc_M, 0)
            return(M)
        }

        ## Here comes the general case.

        ## First, compute the \alpha and \beta.
        alpha <- rowSums(rep(w, each = k) *
                         mapply(function(p, n) p[1 : k] <= n,
                                permutations, nc_memberships))
        ## Alternatively (more literally):
        ##   X <- lapply(permutations, .make_X_from_p)
        ##   alpha1 <- double(length = k)
        ##   for(b in seq(along = permutations)) {
        ##       alpha1 <- alpha1 +
        ##           w[b] * colSums(X[[b]][1 : nc_memberships[b], ])
        ##   }
        
        ## A helper function giving suitably permuted memberships.
        pmem <- function(u, p) {
            ## Only matched classes, similar to the one used in value(),
            ## maybe merge eventually ...
            v <- matrix(0, nr_M, k)
            ind <- 1 : k
            ind <- ind[p[ind] <= ncol(u)]
            if(any(ind))
                v[ , ind] <- u[ , p[ind]]
            v
        }
        beta <- matrix(rowSums(mapply("*",
                                       mapply(pmem,
                                              memberships,
                                              permutations,
                                              SIMPLIFY = FALSE),
                                       w)),
                       nr_M)
        ## Alternatively (more literally):
        ##   beta1 <- matrix(0, nr_M, nc_M)
        ##   for(b in seq(along = permutations)) {
        ##     ind <- seq(length = nc_memberships[b])
        ##     beta1 <- beta1 +
        ##       w[b] * memberships[[b]][, ind] %*% X[[b]][ind, ]
        ##   }
        
        ## Compute the weighted means \bar{M}.
        M <- .cscale(beta, ifelse(alpha > 0, 1 / alpha, 0))
        ## Alternatively (see comments for .cscale()):
        ##   M1 <- beta %*% diag(ifelse(alpha > 0, 1 / alpha, 0))
        
        ## And add dummy classes if necessary.
        if(k < nc_M)
            M <- cbind(M, matrix(0, nr_M, nc_M - k))

        S <- rowSums(M)
        ## Take care of those rows with row sums < 1.
        ind <- (S < 1)
        if(any(ind)) {
            i_0 <- alpha == 0
            if(any(i_0))
                M[ind, i_0] <- 1 / sum(i_0)
            else
                M[ind, ] <- pmax(M[ind, ] + (1 - S[ind]) / nc_M, 0)
        }
        ## Take care of those rows with row sums > 1.
        ind <- (S > 1)
        if(any(ind)) {
            require("quadprog")
            ## Argh.  Call solve.QP() for each such i.  Alternatively,
            ## could set up on very large QP, but is this any better?
            Dmat <- diag(alpha, nc_M)
            Amat <- t(rbind(rep(-1, nc_M), diag(1, nc_M)))
            bvec <- c(-1, rep(0, nc_M))
            for(i in which(ind))
                M[i, ] <- quadprog::solve.QP(Dmat, alpha * M[i, ],
                                             Amat, bvec)$solution
        }

        M
    }

    memberships <- lapply(clusterings, cl_membership)
    permutations <- lapply(memberships, fit_P, M)
    old_value <- value(M, permutations, memberships, w)
    iter <- 1

    while(iter <= maxiter) {
        ## Fit M.
        M <- fit_M(permutations, memberships, w)
        ## Fit \{ P_b \}.
        permutations <- lapply(memberships, fit_P, M)
        ## Update value.
        new_value <- value(M, permutations, memberships, w)
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

    rownames(M) <- rownames(memberships[[1]])    
    M <- cl_membership(as.cl_membership(M[, 1 : k]), k)

    ## Add these attributes here, as the above would not preserve them.
    attr(M, "converged") <- (iter <= maxiter)
    attr(M, "value") <- new_value

    M
}

### ** .cl_consensus_partition_GV1

.cl_consensus_partition_GV1 <-
function(clusterings, weights, control)
    .cl_consensus_partition_AOG(clusterings, weights, control, "GV1")


### * .cl_consensus_partition_GV3

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

### * .cl_consensus_hierarchy_cophenetic

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

### * .cl_consensus_hierarchy_majority

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

### * Utilities

### ** .project_to_leading_columns

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

### ** .make_X_from_p

.make_X_from_p <-
function(p) {
    ## X matrix corresponding to permutation p as needed for the AO
    ## algorithms.  I.e., x_{ij} = 1 iff j->p(j)=i.
    X <- matrix(0, length(p), length(p))
    i <- seq(along = p)
    X[cbind(p[i], i)] <- 1
    X
}

### ** .n_of_nonzero_columns

## <NOTE>
## Could turn this into n_of_classes.matrix().
.n_of_nonzero_columns <-
function(x)
    sum(colSums(x) > 0)
## </NOTE>

### ** .cscale

## <FIXME>
## Move to utilities eventually ...
.cscale <-
function(A, x)
{
    ## Scale the columns of matrix A by the elements of vector x.
    ## Formally, A %*% diag(x), but faster.
    ## Could also use sweep(A, 2, x, "*")
    rep(x, each = nrow(A)) * A
}

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
