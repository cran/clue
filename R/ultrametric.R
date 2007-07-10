### * cl_ultrametric

cl_ultrametric <-
function(x, size = NULL, labels = NULL) 
{
    if(inherits(x, "cl_hierarchy")) {
        ## <NOTE>
        ## Strictly, not every hierarchy corresponds to an ultrametric.
        ## </NOTE>
        return(cl_ultrametric(.get_representation(x),
                              size = size, labels = labels))
    }
    else if(!inherits(x, "cl_ultrametric")) {
        ## Try using cophenetic().
        ## This starts by coercing to hclust, which has methods for all
        ## currently supported hierarchical classification methods.
        ## To support others, either provide as.hclust methods for
        ## these, or make cl_ultrametric() generic and add methods.
        ## Or use the fact that in R >= 2.1.0, stats::cophenetic() is
        ## generic.
        out <- cophenetic(x)
    }
    else {
        out <- x
        if(is.null(labels))
            labels <- attr(x, "Labels")
    }
    .cl_ultrametric_from_veclh(out, labels = labels, size = size)
}

.cl_ultrametric_from_veclh <-
function(x, size = NULL, labels = NULL)
{
    cl_proximity(x, "Ultrametric distances",
                 labels = labels, size = size,
                 class = c("cl_ultrametric", "cl_dissimilarity",
                 "cl_proximity", "dist"))
}

### * as.cl_ultrametric

as.cl_ultrametric <-
function(x)
    UseMethod("as.cl_ultrametric")
as.cl_ultrametric.default <-
function(x)
{
    if(inherits(x, "cl_ultrametric"))
        x
    else if(is.atomic(x))
        .cl_ultrametric_from_veclh(x)
    else
        cl_ultrametric(x)
}
as.cl_ultrametric.matrix <-
function(x)
    .cl_ultrametric_from_veclh(x[row(x) > col(x)],
                               labels = rownames(x))

### * as.dendrogram.cl_ultrametric

as.dendrogram.cl_ultrametric <-
function(object, ...)
    as.dendrogram(as.hclust(object), ...)

### * as.hclust.cl_ultrametric

as.hclust.cl_ultrametric <-
function(x, ...)
{
    ## Hierarchical clustering with single linkage gives the minimal
    ## ultrametric dominated by a dissimilarity, see e.g. Bock (1974,
    ## Theorem 39.2).  Hence, hclust(method = "single") on an
    ## ultrametric gives the hclust representation of the associated
    ## dendrogram.
    hclust(x, "single")
}

### * cophenetic.cl_ultrametric

cophenetic.cl_ultrametric <-
function(x)
    as.dist(x)

### * plot.cl_ultrametric

plot.cl_ultrametric <-
function(x, ...)
    plot(as.dendrogram(x), ...)

### * ls_fit_ultrametric

ls_fit_ultrametric <-
function(x, method = c("SUMT", "IP", "IR"), weights = 1,
         control = list())
{
    if(inherits(x, "cl_ultrametric"))
        return(x)

    ## Handle weights.
    ## This is somewhat tricky ...
    if(is.matrix(weights)) {
        weights <- as.dist(weights)
        if(length(weights) != length(x))
            stop("Arguments 'weights' must be compatible with 'x'.")
    }
    else
        weights <- rep(weights, length.out = length(x))
    if(any(weights < 0))
        stop("Argument 'weights' has negative elements.")
    if(!any(weights > 0))
        stop("Argument 'weights' has no positive elements.")

    method <- match.arg(method)
    switch(method,
           SUMT = .ls_fit_ultrametric_by_SUMT(x, weights, control),
           IP = {
               .ls_fit_ultrametric_by_iterative_projection(x, weights,
                                                           control)
           },
           IR = {
               .ls_fit_ultrametric_by_iterative_reduction(x, weights,
                                                          control)
           })
}

### ** .ls_fit_ultrametric_by_SUMT
           
.ls_fit_ultrametric_by_SUMT <-
function(x, weights = 1, control = list())
{
    ## Fit an ultrametric to a dissimilarity by minimizing euclidean
    ## dissimilarity subject to the ultrametric constraint, using the
    ## sequential algorithm of de Soete (1984) with a slight change: we
    ## try to ensure that what we obtain satisfies the constraints
    ## "exactly" rather than approximately.  We (currently?) do that via
    ## rounding ...

    ## <NOTE>
    ## This fits and hence returns an ultrametric, *not* the hierarchy
    ## corresponding to the ultrametric.
    ## </NOTE>

    w <- weights / sum(weights)

    ## Control parameters:
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

    ## If x is an ultrametric, or satisfies the ultrametricity
    ## constraints, return it.
    if(inherits(x, "cl_ultrametric")
       || (.non_ultrametricity(x, max = TRUE) == 0))
        return(as.cl_ultrametric(x))

    ## For the time being, use a simple minimizer.

    n <- attr(x, "Size")
    if(n <= 2) return(as.cl_ultrametric(x))
    labels <- attr(x, "Labels")

    ## Handle missing values in x along the lines of de Soete (1984):
    ## set the corresponding weights to 0, and impute by the weighted
    ## mean.
    ind <- which(is.na(x))
    if(any(ind)) {
        w[ind] <- 0
        x[ind] <- weighted.mean(x, w, na.rm = TRUE)
    }
    
    ## We follow de Soete's notation, and use the veclh's (vector of
    ## lower half, in S the same as x[lower.tri(x)]) of the respective
    ## proximity objects.

    L <- function(d) sum(w * (x - d) ^ 2)
    P <- .make_penalty_function_ultrametric(n)
    grad_L <- function(d) 2 * w * (d - x)
    grad_P <- .make_penalty_gradient_ultrametric(n)

    if(is.null(start)) {
        ## Initialize by "random shaking".  Use sd() for simplicity.
        start <- replicate(nruns,  
                           x + rnorm(length(x), sd = sd(x) / sqrt(3)),
                           simplify = FALSE)
    }
    
    ## And now ...
    d <- sumt(start, L, P, grad_L, grad_P,
              method = control$method, eps = control$eps,
              q = control$q, verbose = control$verbose,
              control = as.list(control$control))
    
    ## Round to enforce ultrametricity, and hope for the best ...
    ## Alternatively, we could try running one more optimization step
    ## with just the penalty function.
    .cl_ultrametric_from_ultrametric_approximation(d, size = n,
                                                   labels = labels)
}

.make_penalty_function_ultrametric <-
function(n)
    function(d) {
        ## Smooth penalty function measuring the extent of violation of
        ## the ultrametricity constraint.  Also ensure nonnegativity ...
        (.non_ultrametricity(.symmetric_matrix_from_veclh(d, n))
         + sum(pmin(d, 0) ^ 2))
    }

.make_penalty_gradient_ultrametric <-
function(n)
    function(d) {
        gr <- matrix(.C("deviation_from_ultrametricity_gradient",
                        as.double(.symmetric_matrix_from_veclh(d, n)),
                        n,
                        gr = double(n * n),
                        PACKAGE = "clue")$gr,
                     n, n)
        gr[row(gr) > col(gr)] + 2 * sum(pmin(d, 0))
    }

### ** .ls_fit_ultrametric_by_iterative_projection

## <NOTE>
## Functions
##   .ls_fit_ultrametric_by_iterative_projection()
##   .ls_fit_ultrametric_by_iterative_reduction()
## are really identical apart from the name of the C routine they call.
## (But will this necessarily always be the case in the future?)
## Merge maybe ...
## </NOTE>

.ls_fit_ultrametric_by_iterative_projection <-
function(x, weights = 1, control = list())
{
    if(any(diff(weights)))
        warning("Non-identical weights currently not supported.")

    if(attr(x, "Size") <= 2) return(as.cl_ultrametric(x))

    labels <- attr(x, "Labels")
    x <- as.matrix(x)
    n <- nrow(x)

    ## Control parameters:
    ## maxiter,
    maxiter <- control$maxiter
    if(is.null(maxiter))
        maxiter <- 10000
    ## nruns,
    nruns <- control$nruns
    ## order,
    order <- control$order
    ## tol,
    tol <- control$tol
    if(is.null(tol))
        tol <- 1e-8
    ## verbose.
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")

    ## Handle order and nruns.
    if(!is.null(order)) {
        if(!is.list(order))
            order <- as.list(order)
        if(!all(sapply(order,
                       function(o) all(sort(o) == seq_len(n)))))
            stop("All given orders must be valid permutations.")
    }
    else {
        if(is.null(nruns))
            nruns <- 1
        order <- replicate(nruns, sample(n), simplify = FALSE)
    }

    L <- function(d) sum(weights * (x - d) ^ 2)

    d_opt <- NULL
    v_opt <- Inf
    for(run in seq_along(order)) {
        if(verbose)
            cat("Iterative projection run:", run, "\n")
        d <- .C("ls_fit_ultrametric_by_iterative_projection",
                as.double(as.matrix(x)),
                as.integer(n),
                as.integer(order[[run]] - 1),
                as.integer(maxiter),
                iter = integer(1),
                as.double(tol),
                as.logical(verbose))[[1]]
        v <- L(d)
        if(v < v_opt) {
            v_opt <- v
            d_opt <- d
        }
    }
        
    d <- as.dist(matrix(d_opt, n))
    .cl_ultrametric_from_ultrametric_approximation(d, n, labels)
}

### ** .ls_fit_ultrametric_by_iterative_reduction

.ls_fit_ultrametric_by_iterative_reduction <-
function(x, weights = 1, control = list())
{
    if(any(diff(weights)))
        warning("Non-identical weights currently not supported.")

    if(attr(x, "Size") <= 2) return(as.cl_ultrametric(x))

    labels <- attr(x, "Labels")
    x <- as.matrix(x)
    n <- nrow(x)

    ## Control parameters:
    ## maxiter,
    maxiter <- control$maxiter
    if(is.null(maxiter))
        maxiter <- 10000
    ## nruns,
    nruns <- control$nruns
    ## order,
    order <- control$order
    ## tol,
    tol <- control$tol
    if(is.null(tol))
        tol <- 1e-8
    ## verbose.
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")

    ## Handle order and nruns.
    if(!is.null(order)) {
        if(!is.list(order))
            order <- as.list(order)
        if(!all(sapply(order,
                       function(o) all(sort(o) == seq_len(n)))))
            stop("All given orders must be valid permutations.")
    }
    else {
        if(is.null(nruns))
            nruns <- 1
        order <- replicate(nruns, sample(n), simplify = FALSE)
    }

    L <- function(d) sum(weights * (x - d) ^ 2)

    d_opt <- NULL
    v_opt <- Inf
    for(run in seq_along(order)) {
        if(verbose)
            cat("Iterative reduction run:", run, "\n")
        d <- .C("ls_fit_ultrametric_by_iterative_reduction",
                as.double(x),
                as.integer(n),
                as.integer(order[[run]] - 1),
                as.integer(maxiter),
                iter = integer(1),
                as.double(tol),
                as.logical(verbose),
                PACKAGE = "clue")[[1]]
        v <- L(d)
        if(v < v_opt) {
            v_opt <- v
            d_opt <- d
        }
    }
        
    d <- as.dist(matrix(d_opt, n))
    ## <NOTE>
    ## If we want to add an attribute with the convergence info, we need
    ## to do this after using as.cl_ultrametric().
    ## </NOTE>
    .cl_ultrametric_from_ultrametric_approximation(d, n, labels)
}

### * Ultrametric Target Fitters.

### ** ls_fit_ultrametric_target

ls_fit_ultrametric_target <-
function(x, y, weights = 1)
{
    fitter <- if(identical(weights, 1)) # Default.
        function(x, w) mean(x)
    else
        function(x, w) weighted.mean(x, w)
    .fit_ultrametric_target(x, y, weights, fitter)
}

### ** l1_fit_ultrametric_target

l1_fit_ultrametric_target <-
function(x, y, weights = 1)
{
    fitter <- if(identical(weights, 1)) # Default.
        function(x, w) median(x)
    else
        function(x, w) weighted_median(x, w)
    .fit_ultrametric_target(x, y, weights, fitter)
}

### ** .fit_ultrametric_target    

.fit_ultrametric_target <-
function(x, y, w, fitter)
{
    w <- .handle_weights_for_ultrametric_target_fitters(w, x)
    x <- as.matrix(x)
    y <- as.hclust(y)
    n <- length(y$order)
    ilist <- vector("list", n)
    out <- matrix(0, n, n)
    for(i in seq_len(n - 1)) {
        inds <- y$merge[i, ]
        ids1 <- if(inds[1] < 0) -inds[1] else ilist[[inds[1]]]
        ids2 <- if(inds[2] < 0) -inds[2] else ilist[[inds[2]]]
        ilist[[i]] <- c(ids1, ids2)
        out[ids1, ids2] <- fitter(x[ids1, ids2], w[ids1, ids2])
    }
    rownames(out) <- y$labels
    as.cl_ultrametric(out + t(out))
}

### ** .handle_weights_for_ultrametric_target_fitters

.handle_weights_for_ultrametric_target_fitters <-
function(weights, x)
{
    ## Handle weights for the ultrametric target fitters.
    ## This is somewhat tricky ...
    if(is.matrix(weights)) {
        if(any(dim(weights) != attr(x, "Size")))
            stop("Arguments 'weights' must be compatible with 'x'.")
    }
    else
        weights <- as.matrix(.dist_from_vector(rep(weights, length =
                                                   length(x)))) 
    if(any(weights < 0))
        stop("Argument 'weights' has negative elements.")
    if(!any(weights > 0))
        stop("Argument 'weights' has no positive elements.")
    weights
}

### l1_fit_ultrametric

l1_fit_ultrametric <-
function(x, method = c("SUMT", "IRIP"), weights = 1,
         control = list())
{
    if(inherits(x, "cl_ultrametric"))
        return(x)

    ## Handle weights.
    ## This is somewhat tricky ...
    if(is.matrix(weights)) {
        weights <- as.dist(weights)
        if(length(weights) != length(x))
            stop("Arguments 'weights' must be compatible with 'x'.")
    }
    else
        weights <- rep(weights, length.out = length(x))
    if(any(weights < 0))
        stop("Argument 'weights' has negative elements.")
    if(!any(weights > 0))
        stop("Argument 'weights' has no positive elements.")

    method <- match.arg(method)
    switch(method,
           SUMT = .l1_fit_ultrametric_by_SUMT(x, weights, control),
           IRIP = .l1_fit_ultrametric_by_IRIP(x, weights, control))
}

### ** .l1_fit_ultrametric_by_SUMT

.l1_fit_ultrametric_by_SUMT <-
function(x, weights = 1, control = list())
{
    ## Try a SUMT with "pseudo-gradients".

    w <- weights / sum(weights)

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

    ## If x is an ultrametric, or satisfies the ultrametricity
    ## constraints, return it.
    if(inherits(x, "cl_ultrametric")
       || (.non_ultrametricity(x, max = TRUE) == 0))
        return(as.cl_ultrametric(x))

    ## For the time being, use a simple minimizer.

    n <- attr(x, "Size")
    if(n <= 2) return(as.cl_ultrametric(x))
    labels <- attr(x, "Labels")

    L <- function(d) sum(w * abs(d - x))
    P <- .make_penalty_function_ultrametric(n)
    if(gradient) {
        grad_L <- function(d) w * sign(d - x)
        grad_P <- .make_penalty_gradient_ultrametric(n)
    } else
        grad_L <- grad_P <- NULL

    if(is.null(start)) {
        ## Initialize by "random shaking".  Use sd() for simplicity.
        start <- replicate(nruns,  
                           x + rnorm(length(x), sd = sd(x) / sqrt(3)),
                           simplify = FALSE)
    }

    ## And now ...
    d <- sumt(start, L, P, grad_L, grad_P,
              method = control$method, eps = control$eps,
              q = control$q, verbose = control$verbose,
              control = as.list(control$control))
    
    as.cl_ultrametric(d)
}

### ** .l1_fit_ultrametric_by_IRIP

.l1_fit_ultrametric_by_IRIP <-
function(x, weights = 1, control = list())
{
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")
    
    eps <- 1e-6                         # Make settable later.

    ## Initialize by "random shaking" as for the L2 SUMT, but perhaps we
    ## should not do this?  [Or do it differently?]
    u <- x + rnorm(length(x), sd = sd(x) / 3)
    ## u <- x

    iter <- 1

    repeat {
        if(verbose)
            cat("Outer iteration:", iter)
        u_old <- u
        w <- weights / pmax(abs(u - x), 0.000001)
        ## u <- ls_fit_ultrametric(x, weights = w)
        u <- .ls_fit_ultrametric_by_SUMT(x, weights = w,
                                         control = as.list(control$control))
        ## Use some control arguments lateron ...
        D <- sum(abs(u - u_old))
        if(verbose)
            cat(", Change:", D, "\n")
        if(D < eps) break
        iter <- iter + 1
    }

    ## <FIXME>
    ## Use
    ##   .cl_ultrametric_from_ultrametric_approximation(u)
    ## eventually ...
    as.cl_ultrametric(u)
    ## </FIXME>
}

## * ls_fit_sum_of_ultrametrics

ls_fit_sum_of_ultrametrics <-
function(x, nterms = 1, weights = 1, control = list())
{
    ## Control parameters:
    ## eps,
    eps <- control$eps
    if(is.null(eps))
        eps <- 1e-8
    ## method,
    method <- control$method
    if(is.null(method))
        method <- "SUMT"
    ## verbose.
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")
    ## Do this at last.
    control <- as.list(control$control)
    ## And be nice ...
    if(identical(method, "SUMT") && is.null(control$nruns))
        control$nruns <- 10
    
    ## Init.
    u <- rep.int(list(as.cl_ultrametric(0 * x)), nterms)

    ## Loop.
    iter <- 1
    repeat {
        if(verbose)
            cat("Iteration:", iter, "\n")
        delta <- 0
        for(i in seq_len(nterms)) {
            if(verbose)
                cat("Term:", i)
            u_old <- u[[i]]
            ## Compute residual r = x - \sum_{j: j \ne i} u(j)
            r <- x - rowSums(matrix(unlist(u[-i]), ncol = nterms - 1))
            ## Fit residual.
            u[[i]] <- ls_fit_ultrametric(r, method, weights, control)
            ## Accumulate change.
            change <- sum((u[[i]] - u_old) ^ 2)
            if(verbose)
                cat(" Change:", change, "\n")
            delta <- delta + change
        }
        if(verbose)
            cat("Total change:", delta, "\n\n")
        if(delta < eps)
            break
        iter <- iter + 1
    }

    u
}

### * .non_ultrametricity

.non_ultrametricity <-
function(x, max = FALSE)
{
    if(!is.matrix(x))
        x <- .symmetric_matrix_from_veclh(x)
    .C("deviation_from_ultrametricity",
       as.double(x),
       nrow(x),
       fn = double(1),
       as.logical(max),
       PACKAGE = "clue")$fn
}

### * .cl_ultrametric_from_ultrametric_approximation

.cl_ultrametric_from_ultrametric_approximation <-
function(x, size = NULL, labels = NULL)
{
    ## Turn x into an ultrametric after possibly rounding to
    ## non-ultrametric significance.
    mnum <- .non_ultrametricity(x, max = TRUE)
    cl_ultrametric(as.cl_ultrametric(round(x, floor(abs(log10(mnum))))), 
                   size = size, labels = labels)
}

### * .cl_ultrametric_from_classes

.cl_ultrametric_from_classes <-
function(x)
{

    ## Compute an ultrametric from a hierarchy of classes (i.e., an
    ## n-tree).
    
    ## .get_classes_in_hierarchy() orders according to cardinality, but
    ## a consensus method may forget to ...
    x[] <- x[order(sapply(x, length))]

    ## Get the objects (unique codes in the classes).
    objects <- sort(unique(unlist(x)))
    ## (Could also look at the classes of length 1.)

    ## Recursively compute the heights of the classes.
    heights <- double(length = length(x))
    for(i in which(sapply(x, length) > 1)) {
        ## Find the relevant classes.
        j <- sapply(x[seq_len(i - 1)],
                    function(s) all(s %in% x[[i]]))
        heights[i] <- max(heights[j]) + 1
    }

    ## Next, create an incidence matrix (objects by classes).
    incidences <- sapply(x, function(s) objects %in% s)

    ## Now that we have the heights and incidences, we can compute
    ## distances, using the idea that
    ##     distance(i, j) = min(height(A): A contains i and j)
    n <- length(objects)
    d <- matrix(0, n, n)
    for(i in objects)
        d[i, ] <- heights[apply((rep(incidences[i, ], each = n)
                                 & incidences),
                                1, which.max)]
    dimnames(d) <- rep(list(attr(x, "labels")), 2)

    as.cl_ultrametric(d)
}


### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
