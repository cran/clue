### * cl_ultrametric

cl_ultrametric <-
function(x, size = NULL, labels = NULL) 
{
    if(!inherits(x, "cl_ultrametric")) {
        ## Try using cophenetic().
        ## This starts by coercing to hclust, which has methods for all
        ## currently supported hierarchical classification methods.
        ## To support others, either provide as.hclust methods for
        ## these, or make cl_ultrametric() generic and add methods.
        ## Or use the fact that in R >= 2.1.0, stats::cophenetic() is
        ## generic.
        out <- cophenetic(x)
        ## <FIXME 2.1.0>
        ## In R <= 2.0.1, cophenetic() does not preserve labels.
        ## Remove eventually.
        if(is.null(labels))
            labels <- as.hclust(x)$labels
        ## </FIXME>
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
        weights <- rep(weights, length = length(x))
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

    ## Control parameters.
    eps <- control$eps
    ## <TODO>
    ## Maybe divide by length(x) to force maximal non-ultrametricity
    ## down even more.
    if(is.null(eps))
        eps <- .Machine$double.eps
    ## </TODO>
    method <- control$method
    ## We use optim(method = "CG") as default, as nlm() may be
    ## computationally infeasible.
    if(is.null(method))
        method <- "CG"
    q <- control$q
    if(is.null(q))
        q <- 10
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")
    ## Do this at last ...
    control <- as.list(control$control)

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
    P <- function(d) {
        ## Smooth penalty function measuring the extent of violation of
        ## the ultrametricity constraint.
        .non_ultrametricity(.symmetric_matrix_from_veclh(d, n))
    }
    grad_P <- function(d) {
        d <- .symmetric_matrix_from_veclh(d, n)
        n <- nrow(d)
        gr <- matrix(.C("deviation_from_ultrametricity_gradient",
                        as.double(d), n, gr = double(n * n),
                        PACKAGE = "clue")$gr, n, n)
        gr[row(gr) > col(gr)]
    }
    Phi <- function(rho, d) L(d) + rho * P(d)
    grad_Phi <- function(rho, d) {
        2 * w * (d - x) + rho * grad_P(d)
    }

    make_Phi <- if(method == "nlm") {
        function(rho) {
            function(d) {
                y <- Phi(rho, d)
                attr(y, "gradient") <- grad_Phi(rho, d)
                y
            }
        }
    }
    else
        function(rho) {
            function(d)
                Phi(rho, d)
        }
    make_grad_Phi <- function(rho) { function(d) grad_Phi(rho, d) }

    ## <NOTE>
    ## For the penalized minimization, the Newton-type nlm() may be
    ## computationally infeasible (although it works much faster on the
    ## Phonemes data).
    ## De Soete recommends using Conjugate Gradients.
    ## We provide a simple choice: by default, optim(method = "CG") is
    ## used.  If control$method is non-null and not "nlm", we use
    ## optim() with this method.  In both cases, control$control gives
    ## the control parameters for optim().
    ## If control$method is "nlm", nlm() is used, in which case
    ## control$control is ignored.  Note that we *must* call nlm() which
    ## checking analyticals turned off, as in fact the penalty function
    ## is not even continuous ...
    optimize_with_penalty <- if(method == "nlm")
        function(rho, d)
            nlm(make_Phi(rho), d, check.analyticals = FALSE) $ estimate
    else {
        function(rho, d)
            optim(d, make_Phi(rho), gr = make_grad_Phi(rho),
                  method = method, control = control) $ par
    }
    ## </NOTE>
    
    ## Initialize by "random shaking".  Use sd() for simplicity.
    d_old <- x
    d <- x + rnorm(length(x), sd = sd(x) / 3)
    ## <TODO>
    ## Better upper/lower bounds for rho?
    rho <- L(d) / max(P(d), 0.00001)
    ## </TODO>

    iter <- 1
    ## <TODO>
    ## Shouldnt't we also have maxiter, just in case ...?
    ## </TODO>
    while(sum((d_old - d) ^ 2) >= eps) {
        if(verbose)
            cat("Iteration:", iter,
                "Rho:", rho,
                "P:", P(d),
                "\n")
        d_old <- d
        d <- optimize_with_penalty(rho, d)
        iter <- iter + 1        
        rho <- q * rho
    }

    ## Round to enforce ultrametricity, and hope for the best ...
    ## Alternatively, we could try running one more optimization step
    ## with just the penalty function.
    .cl_ultrametric_from_ultrametric_approximation(d, size = n,
                                                   labels = labels)
}

### * .ls_fit_ultrametric_by_iterative_projection

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
    ## order,
    order <- control$order
    if(is.null(order))
        order <- 1 : n
    else {
        if(!all(sort(order) == 1 : n))
            stop("Given order is not a valid permutation.")
    }
    ## tol,
    tol <- control$tol
    if(is.null(tol))
        tol <- 1e-8
    ## verbose.
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")
    
    out <- .C("ls_fit_ultrametric_by_iterative_projection",
              as.double(as.matrix(x)),
              as.integer(n),
              as.integer(order - 1),
              as.integer(maxiter),
              iter = integer(1),
              as.double(tol),
              as.logical(verbose))
    d <- as.dist(matrix(out[[1]], n))
    .cl_ultrametric_from_ultrametric_approximation(d, n, labels)
}

### * .ls_fit_ultrametric_by_iterative_reduction

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
    ## order,
    order <- control$order
    if(is.null(order))
        order <- 1 : n
    else {
        if(!all(sort(order) == 1 : n))
            stop("Given order is not a valid permutation.")
    }
    ## tol,
    tol <- control$tol
    if(is.null(tol))
        tol <- 1e-8
    ## verbose.
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")

    out <- .C("ls_fit_ultrametric_by_iterative_reduction",
              as.double(x),
              as.integer(n),
              as.integer(order - 1),
              as.integer(maxiter),
              iter = integer(1),
              as.double(tol),
              as.logical(verbose),
              PACKAGE = "clue")
    d <- as.dist(matrix(out[[1]], n))
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
    if(identical(weights, 1))           # Default.
        fitter <- function(x, w) mean(x)
    else
        fitter <- function(x, w) weighted.mean(x, w)
    .fit_ultrametric_target(x, y, weights, fitter)
}

### ** l1_fit_ultrametric_target

l1_fit_ultrametric_target <-
function(x, y, weights = 1)
{
    if(identical(weights, 1))           # Default.
        fitter <- function(x, w) median(x)
    else
        fitter <- function(x, w) weighted_median(x, w)
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
    for(i in 1 : (n - 1)) {
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

### * .dist_from_vector
        
.dist_from_vector <-
function(x, n = NULL, labels = NULL)
{
    ## This might be useful as as.dist.vector, perhaps without the extra
    ## argument n then which we only have for minimal performance gains.
    if(is.null(n))
        n <- as.integer((sqrt(1 + 8 * length(x)) + 1) / 2)
    attr(x, "Size") <- n
    if(!is.null(labels))
        attr(x, "Labels") <- labels
    class(x) <- "dist"
    x
}

### * .symmetric_matrix_from_veclh

.symmetric_matrix_from_veclh <-
function(x, n = NULL)
{
    ## In essence the same as as.matrix.dist, but without handling the
    ## additional attributes that dist objects might have.
    if(is.null(n))
        n <- as.integer((sqrt(1 + 8 * length(x)) + 1) / 2)
    M <- matrix(0, n, n)
    M[row(M) > col(M)] <- x
    M + t(M)
}

### .non_ultrametricity

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

### .cl_ultrametric_from_ultrametric_approximation

.cl_ultrametric_from_ultrametric_approximation <-
function(x, size = NULL, labels = NULL)
{
    ## Turn x into an ultrametric after possibly rounding to
    ## non-ultrametric significance.
    mnum <- .non_ultrametricity(x, max = TRUE)
    cl_ultrametric(as.cl_ultrametric(round(x, floor(abs(log10(mnum))))), 
                   size = size, labels = labels)
}

### .cl_ultrametric_from_classes

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
        j <- sapply(x[seq(length = i - 1)],
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
