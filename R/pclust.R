cl_pclust <-
function(x, k, m = 1, control = list())
{
    ## Partition a cluster ensemble x into (at most) k classes by
    ## minimizing
    ##   \sum_b \sum_j u_{bj}^m d(x_b, p_j) ^ e
    ## for "suitable" prototypes p_1, ..., p_k, where 1 <= m < \infty,
    ## with 1 corresponding to hard (secondary) partitions, and d a
    ## dissimilarity measure (such as Euclidean dissimilarity of
    ## partitions or hierarchies).
    ##
    ## The algorithm works whenever there is a consensus method for
    ## solving
    ##   \sum_b u_{bj}^m d(x_b, p) ^ e => \min_p
    ##
    ## As we refer to consensus methods by their *name* (e.g., 'HBH'),
    ## we rely on the registration mechanism (set_cl_consensus_method())
    ## to provide the required information about d and e.

    maxiter <- control$maxiter
    if(is.null(maxiter))
        maxiter <- 100
    reltol <- control$reltol
    if(is.null(reltol))
        reltol <- sqrt(.Machine$double.eps)
    ## To be passed on to cl_consensus().
    method <- control$method
    ## Do this at last ...
    control <- as.list(control$control)

    clusterings <- as.cl_ensemble(x)
    B <- length(clusterings)

    is_partition_ensemble <-
        inherits(clusterings, "cl_partition_ensemble")

    if(is_partition_ensemble) {
        ## Canonicalize by turning into an ensemble of membership
        ## matrices with the same (minimal) number of columns.
        memberships <- lapply(clusterings, cl_membership,
                              max(sapply(clusterings, n_of_classes)))
        clusterings <- cl_ensemble(list = memberships)
    }

    if(is.null(method)) {
        ## Default: use squared Euclidean dissimilarities.  Note that
        ## there is a small catch, as the default consensus methods for
        ## partitions and hierarchies are *not* named "euclidean".
        ## Hence, fill in d and e ourselves ...
        e <- 2
        d <- function(x, y = NULL)
            cl_dissimilarity(x, y, method = "euclidean")
    }
    else {
        ## Get required information on d and e from the registry.        
        if(!inherits(method, "cl_consensus_method")) {
            type <- if(is_partition_ensemble)
                "partition"
            else
                "hierarchy" 
            method <- get_cl_consensus_method(method, type)
            ## Note that this avoids registry lookup in subsequent
            ## calls to cl_consensus().
        }
        if(is.null(method$exponent))
            stop("No information on exponent in consensus method used.")
        e <- method$exponent
        if(is.null(method$dissimilarity))
            stop("No information on dissimilarity in consensus method used.")
        d <- function(x, y = NULL)
            cl_dissimilarity(x, y, method = method$dissimilarity)
    }
    
    ## Take random clusterings as prototypes.
    ## For partitions, it may be better to use random soft partitions.
    prototypes <- clusterings[sample(1 : B, k)]
    dissimilarities <- d(clusterings, prototypes) ^ e

    if(m == 1) {
        ## Hard secondary partitions.
        class_ids <- max.col( - dissimilarities )
        old_value <-
            sum(.one_entry_per_column(dissimilarities, class_ids))
        iter <- 1
        while(iter <= maxiter) {
            ## Compute new prototypes.
            ## <NOTE>
            ## Splitting lists is broken in R versions up to 2.0.1, so
            ## we use a loop here.  Something based on
            ##   lapply(split(clusterings, class_ids), cl_consensus, ...)
            ## would be nicer.
            for(j in unique(class_ids))
                prototypes[[j]] <-
                    cl_consensus(clusterings[class_ids %in% j],
                                 method = method, control = control)
            ## </NOTE>
            ## Update the class ids.
            dissimilarities <- d(clusterings, prototypes) ^ e
            class_ids <- max.col( - dissimilarities )
            new_value <-
                sum(.one_entry_per_column(dissimilarities, class_ids))
            if(abs(old_value - new_value)
               < reltol * (abs(old_value) + reltol))
                break
            old_value <- new_value
            iter <- iter + 1
        }
        u <- matrix(0, B, k)
        u[cbind(seq(length = B), class_ids)] <- 1
    }
    else {
        ## Soft secondary partitions.
        value <- function(u, dissimilarities)
            sum(u ^ m * dissimilarities)
        u <- .memberships_from_cross_dissimilarities(dissimilarities, m)
        old_value <- value(u, dissimilarities)
        iter <- 1        
        while(iter <= maxiter) {
            ## Update the prototypes.
            ## This amounts to solving, for each j:
            ##   \sum_b u_{bj}^m d(x_b, p) ^ e => \min_p
            ## I.e., p_j is the *weighted* consensus clustering of the
            ## x_b with corresponding weights u_{bj}^m.
            for(j in 1 : k) {
                prototypes[[j]] <-
                    cl_consensus(clusterings, weights = u[, j] ^ m,
                                 method = method, control = control)
            }
            ## Update u.
            dissimilarities <- d(clusterings, prototypes) ^ e
            u <- .memberships_from_cross_dissimilarities(dissimilarities, m)
            new_value <- value(u, dissimilarities)
            if(abs(old_value - new_value)
               < reltol * (abs(old_value) + reltol))
                break
            old_value <- new_value
            iter <- iter + 1
        }
        class_ids <- max.col(u)
    }

    dissimilarities <- as.matrix(d(clusterings) ^ e)
    ## Note that our dissimilarities inherit from "cl_proximity" but not
    ## "dist", and as.dist() is not a generic function.
    u <- cl_membership(as.cl_membership(u), k)
                                          
    out <- list(prototypes = prototypes,
                membership = u,
                cluster = class_ids,
                silhouette = silhouette(class_ids,
                                        dmatrix = dissimilarities),
                validity = cl_validity(u, dissimilarities),
                m = m, d = d, e = e,
                call = match.call())
    attr(out, "converged") <- (iter <= maxiter)
    class(out) <- "cl_pclust"
    out
}

print.cl_pclust <-
function(x, ...)
{
    is_hard <- x$m == 1
    class_ids <- cl_class_ids(x)
    ## <TODO>
    ## Note that for m > 1 we could also get hard partitions ...
    ## Can we really?
    txt <- paste("A",
                 if(is_hard) "hard" else "soft",
                 "partition",
                 if(!is_hard) paste("(degree m = ", x$m, ")", sep = ""),
                 "of a cluster ensemble with",
                 length(class_ids), 
                 "elements into",
                 n_of_classes(x),
                 "classes.")
    ## </TODO>
    writeLines(strwrap(txt))
    if(is_hard) {
        writeLines("Class ids:")
        print(class_ids, ...)
    }
    else {
        writeLines("Class memberships:")
        print(cl_membership(x), ...)
        writeLines("Class ids of closest hard partition:")
        print(class_ids, ...)
    }
    print(cl_validity(x), ...)
    invisible(x)
}

silhouette.cl_pclust <-
function(x, ...)
    x$silhouette
