ls_fit_addtree <-
function(x, method = c("IP", "IR"), weights = 1, control = list())
{
    ## Handle argument 'weights'.
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
           IP = {
               .ls_fit_addtree_by_iterative_projection(x, weights,
                                                       control)
           },
           IR = {
               .ls_fit_addtree_by_iterative_reduction(x, weights,
                                                      control)
           })
}

## <NOTE>
## Functions
##   .ls_fit_addtree_by_iterative_projection()
##   .ls_fit_addtree_by_iterative_reduction()
## are really identical apart from the name of the C routine they call.
## (But will this necessarily always be the case in the future?)
## Merge maybe ...
## </NOTE>

.ls_fit_addtree_by_iterative_projection <-
function(x, weights = 1, control = list())
{
    if(any(diff(weights)))
        warning("Non-identical weights currently not supported.")

    if(attr(x, "Size") <= 3)
        return(x)

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
    
    out <- .C("ls_fit_addtree_by_iterative_reduction",
              as.double(x),
              as.integer(n),
              as.integer(order - 1),
              as.integer(maxiter),
              iter = integer(1),
              as.double(tol),
              as.logical(verbose),
              PACKAGE = "clue")
    ## Nothing fancy for the time being, as we have no data structures
    ## for representing "additive trees".
    ## Also need some rounding a la
    ## .cl_ultrametric_from_ultrametric_approximation().
    d <- matrix(out[[1]], n)
    dimnames(d) <- list(labels, labels)
    as.dist(d)
}

.ls_fit_addtree_by_iterative_reduction <-
function(x, weights = 1, control = list())
{
    if(any(diff(weights)))
        warning("Non-identical weights currently not supported.")

    if(attr(x, "Size") <= 3)
        return(x)

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

    out <- .C("ls_fit_addtree_by_iterative_reduction",
              as.double(x),
              as.integer(n),
              as.integer(order - 1),
              as.integer(maxiter),
              iter = integer(1),
              as.double(tol),
              as.logical(verbose),
              PACKAGE = "clue")
    ## Nothing fancy for the time being, as we have no data structures
    ## for representing "additive trees".
    ## Also need some rounding a la
    ## .cl_ultrametric_from_ultrametric_approximation().
    d <- matrix(out[[1]], n)
    dimnames(d) <- list(labels, labels)
    as.dist(d)
}

.non_additivity <-
function(x, max = FALSE)
{
    if(!is.matrix(x))
        x <- .symmetric_matrix_from_veclh(x)
    .C("deviation_from_additivity",
       as.double(x),
       nrow(x),
       fn = double(1),
       as.logical(max),
       PACKAGE = "clue")$fn
}


### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
