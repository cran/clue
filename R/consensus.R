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

    as.cl_partition(cl_membership(as.cl_membership(M[, 1 : k]), k))
}

### * .cl_consensus_partition_AOS

.cl_consensus_partition_AOS <-
function(clusterings, weights, control,
         type = c("SE", "HE", "SM", "HM"))
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

    value <-
        switch(type,
               SE = , HE = function(M, memberships, w) {
                   sum(w * sapply(memberships,
                                  function(u) sum((u - M) ^ 2)))
               },
               SM = , HM = function(M, memberships, w) {
                   sum(w * sapply(memberships,
                                  function(u) sum(abs(u - M))))
               })
    ## Return the M[, ind] column permutation of M optimally matching N.
    match_memberships <-
        switch(type,
               SE = , HE = function(M, N) {
                   M[, solve_LSAP(crossprod(N, M), max = TRUE)]
               },
               SM = , HM = function(M, N) {
                   M[, solve_LSAP(.cxdist(N, M, "manhattan"))]
               })
    ## Function for fitting M to (fixed) memberships \{ M_b P_b \}.
    ## As we use a common number of columns for all membership matrices
    ## involved, we need to pass the desired 'k' ...
    fit_M <-
        switch(type,
               SE = function(memberships, w, k) {
                   ## Update M as \sum w_b M_b P_b.
                   M <- .weighted_sum_of_matrices(memberships, w, nrow(M))
                   ## If k < k_all, "project" as indicated in Gordon &
                   ## Vichi (2001), p. 238.
                   if(k < ncol(M))
                       M <- .project_to_leading_columns(M, k)
                   M
               },
               HE = , HM = function(memberships, w, k) {
                   ## Compute M as \sum w_b M_b P_b.
                   M <- .weighted_sum_of_matrices(memberships, w, nrow(M))
                   ## And compute a closest hard partition H(M) from
                   ## that, using the first k columns of M.
                   cl_membership(as.cl_membership(max.col(M[ , 1 : k])),
                                 ncol(M))
                   ## <FIXME>
                   ## Perhaps more efficiently:
                   ##   .cl_membership_from_class_ids(max.col(M[ , 1 : k]),
                   ##                                 ncol(M))
                   ## </FIXME>
               },
               SM = .l1_fit_M)

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

    as.cl_partition(M)
}

.l1_fit_M <-
function(memberships, w, k)
{
    ## Determine stochastic matrix M with at most k leading nonzero
    ## columns such that
    ##
    ##    \sum_b w_b \sum_{i,j} | m_{ij}(b) - m_{ij} | => min
    ##
    ## where the sum over j goes from 1 to k.
    ##
    ## Clearly, this can be done separately for each row, where we need
    ## to minimize
    ##
    ##    \sum_b w_b \sum_j | y_j(b) - x_j | => min
    ##
    ## over all probability vectors x.  Such problems can e.g. be solved
    ## via the following linear program:
    ##
    ##    \sum_b \sum_j w_b e'(u(b) + v(b)) => min
    ##
    ## subject to
    ##
    ##    u(1), v(1), ..., u(B), v(B), x >= 0
    ##                   x + u(b) - v(b)  = y(b),    b = 1, ..., B
    ##                               e'x  = 1
    ##
    ## (where e = [1, ..., 1]).
    ##
    ## So we have one long vector z of "variables":
    ##
    ##    z = [u(1)', v(1)', ..., u(B)', v(B)', x']'
    ##
    ## of length (2B + 1) k, with x the object of interest.
    
    ## Rather than providing a separate function for weighted L1 fitting
    ## of probability vectors we prefer doing "everything" at once, in
    ## order to avoid recomputing the coefficients and constraints of
    ## the associated linear program.

    B <- length(memberships)
    L <- (2 * B + 1) * k

    ## Set up associated linear program.

    ## Coefficients in the objective function.
    objective_in <- c(rep(w, each = 2 * k), rep.int(0, k))
    
    ## Constraints.
    constr_mat <- rbind(diag(1, L),
                        cbind(kronecker(diag(1, B),
                                        cbind(diag(1, k),
                                              diag(-1, k))),
                              kronecker(rep.int(1, B),
                                        diag(1, k))),
                        c(rep.int(0, 2 * B * k), rep.int(1, k)))
    constr_dir <- c(rep.int(">=", L), rep.int("==", B * k + 1))

    ind <- seq(from = 2 * B * k + 1, length = k)
    nr <- NROW(memberships[[1]])
    nc <- NCOL(memberships[[1]])
    M <- matrix(0, nr = nr, nc = k)

    ## Put the memberships into one big array so that we can get their
    ## rows more conveniently (and efficiently):

    memberships <- array(unlist(memberships), c(nr, nc, B))

    require("lpSolve")

    for(i in seq(length = nr)) {
        out <- lpSolve::lp("min",
                           objective_in,
                           constr_mat,
                           constr_dir,
                           c(rep.int(0, L), memberships[i, 1 : k, ], 1))
        M[i, ] <- out$solution[ind]
    }

    ## Add zero columns if necessary.
    if(k < nc)
        M <- cbind(M, matrix(0, nr, nc - k))

    M
}

### ** .cl_consensus_partition_soft_euclidean

.cl_consensus_partition_soft_euclidean <-
function(clusterings, weights, control)
    .cl_consensus_partition_AOS(clusterings, weights, control, "SE")

### ** .cl_consensus_partition_hard_euclidean

.cl_consensus_partition_hard_euclidean <-
function(clusterings, weights, control)
    .cl_consensus_partition_AOS(clusterings, weights, control, "HE")

### ** .cl_consensus_partition_soft_manhattan

.cl_consensus_partition_soft_manhattan <-
function(clusterings, weights, control)
    .cl_consensus_partition_AOS(clusterings, weights, control, "SM")

### ** .cl_consensus_partition_hard_manhattan

.cl_consensus_partition_hard_manhattan <-
function(clusterings, weights, control)
    .cl_consensus_partition_AOS(clusterings, weights, control, "HM")

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
            M <- .weighted_sum_of_matrices(mapply(function(u, p)
                                                  u[ , p[1 : k]],
                                                  memberships,
                                                  permutations,
                                                  SIMPLIFY = FALSE),
                                           w, nr_M)
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
        beta <- .weighted_sum_of_matrices(mapply(pmem,
                                                 memberships,
                                                 permutations,
                                                 SIMPLIFY = FALSE),
                                          w, nr_M)
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

    as.cl_partition(M)
}

### ** .cl_consensus_partition_GV1

.cl_consensus_partition_GV1 <-
function(clusterings, weights, control)
    .cl_consensus_partition_AOG(clusterings, weights, control, "GV1")


### * .cl_consensus_partition_GV3

.cl_consensus_partition_GV3 <-
function(clusterings, weights, control)
{
    ## Use a SUMT to solve
    ##   \| Y - M M' \|_F^2 => min
    ## where M is a membership matrix and Y = \sum_b w_b M_b M_b'.

    B <- length(clusterings)    
    n <- n_of_objects(clusterings)
    max_n_of_classes <- max(sapply(clusterings, n_of_classes))

    ## Control parameters:
    ## k,
    k <- control$k
    if(is.null(k))
        k <- max_n_of_classes
    ## nruns,
    nruns <- control$nruns
    ## start.
    start <- control$start

    w <- weights / sum(weights)

    comemberships <-
        lapply(clusterings, function(x) {
            ## No need to force a common k here.
            tcrossprod(cl_membership(x))
        })

    Y <- .weighted_sum_of_matrices(comemberships, w, n)

    ## Handle start values and number of runs.
    if(!is.null(start)) {
        if(!is.list(start)) {
            ## Be nice to users.
            start <- list(start)
        }
    } else {
        if(is.null(nruns)) {
            ## Use nruns only if start is not given.
            nruns <- 1
        }
        e <- eigen(Y, symmetric = TRUE)
        ## Use M <- U_k \lambda_k^{1/2}, or random perturbations
        ## thereof.
        M <- e$vectors[, 1:k] * rep(sqrt(e$values[1:k]), each = n)
        m <- c(M)
        start <- c(list(m),
                   replicate(nruns - 1,
                             m + rnorm(length(m), sd = sd(m) / sqrt(3)),
                             simplify = FALSE))
    }

    y <- c(Y)

    L <- function(m) sum((y - tcrossprod(matrix(m, n))) ^ 2)
    P <- .make_penalty_function_membership(n, k)
    grad_L <- function(m) {
        M <- matrix(m, n)
        4 * c((tcrossprod(M) - Y) %*% M)
    }
    grad_P <- .make_penalty_gradient_membership(n, k)

    m <- SUMT(start, L, P, grad_L, grad_P,
              method = control$method, eps = control$eps,
              q = control$q, verbose = control$verbose,
              control = as.list(control$control))

    ## Ensure that a stochastic matrix is returned.
    M <- matrix(pmax(m, 0), n)
    M <- M / rowSums(M)
    rownames(M) <- rownames(cl_membership(clusterings[[1]]))
    as.cl_partition(cl_membership(as.cl_membership(M), k))
}

### * .cl_consensus_partition_soft_symdiff

.cl_consensus_partition_soft_symdiff <-
function(clusterings, weights, control)
{
    ## Use a SUMT to solve
    ##   \sum_b w_b \sum_{ij} | c_{ij}(b) - c_{ij} | => min
    ## where C(b) = comembership(M(b)) and C = comembership(M) and M is
    ## a membership matrix.

    ## Control parameters:
    ## gradient,
    gradient <- control$gradient
    if(is.null(gradient))
        gradient <- TRUE
    ## k,
    k <- control$k
    ## nruns,
    nruns <- control$nruns
    ## start.
    start <- control$start

    ## Handle start values and number of runs.
    if(!is.null(start)) {
        if(!is.list(start)) {
            ## Be nice to users.
            start <- list(start)
        }
    } else if(is.null(nruns)) {
        ## Use nruns only if start is not given.
        nruns <- 1
    }
    
    max_n_of_classes <- max(sapply(clusterings, n_of_classes))
    if(is.null(k))
        k <- max_n_of_classes

    B <- length(clusterings)    
    n <- n_of_objects(clusterings)
    
    w <- weights / sum(weights)

    comemberships <-
        lapply(clusterings, function(x) {
            ## No need to force a common k here.
            tcrossprod(cl_membership(x))
        })

    ## Handle start values and number of runs.
    if(!is.null(start)) {
        if(!is.list(start)) {
            ## Be nice to users.
            start <- list(start)
        }
    } else {
        if(is.null(nruns)) {
            ## Use nruns only if start is not given.
            nruns <- 1
        }
        ## Try using a rank k "root" of the weighted median of the
        ## comemberships as starting value.
        Y <- apply(array(unlist(comemberships), c(n, n, B)), c(1, 2),
                   weighted_median, w)
        e <- eigen(Y, symmetric = TRUE)
        ## Use M <- U_k \lambda_k^{1/2}, or random perturbations
        ## thereof.
        M <- e$vectors[, 1:k] * rep(sqrt(e$values[1:k]), each = n)
        m <- c(M)
        start <- c(list(m),
                   replicate(nruns - 1,
                             m + rnorm(length(m), sd = sd(m) / sqrt(3)),
                             simplify = FALSE))
    }

    L <- function(m) {
        M <- matrix(m, n)
        C_M <- tcrossprod(M)
        sum(w * sapply(comemberships,
                       function(C) sum(abs(C_M - C))))
    }
    P <- .make_penalty_function_membership(n, k)
    if(gradient) {
        grad_L <- function(m) {
            M <- matrix(m, n)
            C_M <- tcrossprod(M)
            .weighted_sum_of_matrices(lapply(comemberships,
                                             function(C)
                                             2 * sign(C_M - C) %*% M),
                                      w, n)
        }
        grad_P <- .make_penalty_gradient_membership(n, k)
    }
    else
        grad_L <- grad_P <- NULL

    m <- SUMT(start, L, P, grad_L, grad_P,
              method = control$method, eps = control$eps,
              q = control$q, verbose = control$verbose,
              control = as.list(control$control))

    ## Ensure that a stochastic matrix is returned.
    M <- matrix(pmax(m, 0), n)
    M <- M / rowSums(M)
    rownames(M) <- rownames(cl_membership(clusterings[[1]]))
    as.cl_partition(cl_membership(as.cl_membership(M), k))
}

### * .cl_consensus_hierarchy_cophenetic

.cl_consensus_hierarchy_cophenetic <-
function(clusterings, weights, control)
{
    w <- weights / sum(weights)
    B <- length(clusterings)
    ultrametrics <- lapply(clusterings, cl_ultrametric)
    dissimilarities <- .weighted_sum_of_vectors(ultrametrics, w)
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
    as.cl_dendrogram(ls_fit_ultrametric(d, control = control))
}

### * .cl_consensus_hierarchy_manhattan

.cl_consensus_hierarchy_manhattan <-
function(clusterings, weights, control)
{
    ## Control parameters:
    ## gradient,
    gradient <- control$gradient
    if(is.null(gradient))
        gradient <- TRUE
    ## nruns,
    nruns <- control$nruns
    ## start.
    start <- control$start

    ## Handle start values and number of runs.
    if(!is.null(start)) {
        if(!is.list(start)) {
            ## Be nice to users.
            start <- list(start)
        }
    } else if(is.null(nruns)) {
        ## Use nruns only if start is not given.
        nruns <- 1
    }
    
    w <- weights / sum(weights)
    B <- length(clusterings)
    ultrametrics <- lapply(clusterings, cl_ultrametric)
    
    if(B == 1)
        return(as.cl_dendrogram(ultrametrics[[1]]))

    n <- n_of_objects(ultrametrics[[1]])
    labels <- cl_object_names(ultrametrics[[1]])

    ## We need to do
    ##
    ##   \sum_b w_b \sum_{i,j} | u_{ij}(b) - u_{ij} | => min
    ##
    ## over all ultrametrics u.  Let's use a SUMT (for which "gradients"
    ## can optionally be switched off) ...
    
    L <- function(d) {
        sum(w * sapply(ultrametrics, function(u) sum(abs(u - d))))
        ## Could also do something like
        ##   sum(w * sapply(ultrametrics, cl_dissimilarity, d,
        ##                  "manhattan"))
    }
    P <- .make_penalty_function_ultrametric(n)
    if(gradient) {
        grad_L <- function(d) {
            ## "Gradient" is \sum_b w_b sign(d - u(b)).
            .weighted_sum_of_vectors(lapply(ultrametrics,
                                            function(u) sign(d - u)),
                                     w)
        }
        grad_P <- .make_penalty_gradient_ultrametric(n)
    } else
        grad_L <- grad_P <- NULL
            
    if(is.null(start)) {
        ## Initialize by "random shaking" of the weighted median of the
        ## ultrametrics.  Any better ideas?
        ## <FIXME>
        ## Using var(x) / 3 is really L2 ...
        ## </FIXME>
        x <- apply(matrix(unlist(ultrametrics), nc = B), 1, 
                   weighted_median, w)
        start <- replicate(nruns,  
                           x + rnorm(length(x), sd = sd(x) / sqrt(3)),
                           simplify = FALSE)
    }

    d <- SUMT(start, L, P, grad_L, grad_P,
              method = control$method, eps = control$eps,
              q = control$q, verbose = control$verbose,
              control = as.list(control$control))

    ## Round to enforce ultrametricity, and hope for the best ...
    d <- .cl_ultrametric_from_ultrametric_approximation(d, size = n,
                                                        labels = labels)

    as.cl_dendrogram(d)
}

### * .cl_consensus_hierarchy_majority

.cl_consensus_hierarchy_majority <-
function(clusterings, weights, control)
{
    w <- weights / sum(weights)

    p <- control$p
    if(is.null(p))
        p <- 1 / 2
    ## <FIXME>
    ## Add bounds checking for [1/2, 1] eventually.
    ## </FIXME>

    classes <- lapply(clusterings, cl_classes)
    all_classes <- unique(unlist(classes, recursive = FALSE))
    gamma <- double(length = length(all_classes))
    for(i in seq(along = classes))
        gamma <- gamma + w[i] * !is.na(match(all_classes, classes[[i]]))

    maj_classes <- if(p == 1) {
        ## Strict consensus tree.
        all_classes[gamma == 1]
    }
    else
        all_classes[gamma > p]
    attr(maj_classes, "labels") <- attr(classes[[1]], "labels")

    ## <FIXME>
    ## Stop auto-coercing that to dendrograms once we have suitable ways
    ## of representing n-trees.
    as.cl_hierarchy(.cl_ultrametric_from_classes(maj_classes))
    ## </FIXME>
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
## </FIXME>

## .make_penalty_function_membership

.make_penalty_function_membership <-
function(nr, nc)
    function(m) {
        sum(pmin(m, 0) ^ 2) + sum((rowSums(matrix(m, nr)) - 1) ^ 2)
    }

## .make_penalty_gradient_membership

.make_penalty_gradient_membership <-
function(nr, nc)
    function(m) {
        2 * (pmin(m, 0) + rep.int(rowSums(matrix(m, nr)) - 1, nc))
    }


### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
