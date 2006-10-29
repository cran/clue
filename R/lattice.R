cl_meet <-
function(x, y)
{
    ## General case.
    ## x either an ensemble, or x and y two clusterings with the same
    ## number of objects.
    if(!inherits(x, "cl_ensemble")) {
        ## Be nice about error messages.
        if(n_of_objects(x) != n_of_objects(y))
            stop("Arguments 'x' and 'y' must have the same number of objects.")
        x <- cl_ensemble(x, y)
    }

    if(inherits(x, "cl_partition_ensemble"))
        .cl_meet_partition(x)
    else if(inherits(x, "cl_dendrogram_ensemble"))
        .cl_meet_dendrogram(x)
    else
        stop("Cannot compute meet of given clusterings.")
}

.cl_meet_partition <-
function(x)    
{
    x <- unique(x)
    if(length(x) == 1)
        return(cl_partition_by_class_ids(cl_class_ids(x[[1]])))

    ids <- seq_len(n_of_objects(x[[1]]))
    ## Cross-classify the objects.
    z <- split(ids, lapply(x, cl_class_ids))
    ## Subscript on the non-empty cells to get adjacent class ids.
    lens <- sapply(z, length)
    pos <- which(lens > 0)
    ids[unlist(z, use.names = FALSE)] <-
        rep(seq_along(z[pos]), lens[pos])
    cl_partition_by_class_ids(ids)
}

.cl_meet_dendrogram <-
function(x)
{
    ## Meet of an ensemble of dendrograms.
    ## We need the maximal ultrametric dominated by the given ones,
    ## which can be obtained by hierarchical clustering with single
    ## linkage on the pointwise minima of the ultrametrics.
    as.cl_dendrogram(hclust(as.dist(do.call(pmin,
                                            lapply(x, cl_ultrametric))),
                            "single"))
}

cl_join <-
function(x, y)
{
    ## General case.
    ## x either an ensemble, or x and y two clusterings with the same
    ## number of objects.
    if(!inherits(x, "cl_ensemble")) {
        ## Be nice about error messages.
        if(n_of_objects(x) != n_of_objects(y))
            stop("Arguments 'x' and 'y' must have the same number of objects.")
        x <- cl_ensemble(x, y)
    }

    if(inherits(x, "cl_partition_ensemble"))
        .cl_join_partition(x)
    else if(inherits(x, "cl_dendrogram_ensemble"))
        .cl_join_dendrogram(x)
    else
        stop("Cannot compute join of given clusterings.")
}

.cl_join_partition <-
function(x)
{
    x <- unique(x)
    if(length(x) == 1)
        return(cl_partition_by_class_ids(cl_class_ids(x[[1]])))

    ## Canonicalize: ensure that class ids are always the integers from
    ## one to the number of classes.
    n <- sapply(x, n_of_classes)
    ids <- mapply(function(p, ncp) match(cl_class_ids(p),
                                         seq_len(ncp)),
                  x, n, SIMPLIFY = FALSE)
    ## Order according to the number of classes.
    ids <- ids[order(n)]

    ## And now incrementally build the join.
    jcids <- ids[[1]]                   # Class ids of the current join.
    jnc <- length(unique(jcids))        # Number of classes of this.
    for(b in seq(from = 2, to = length(x))) {
        z <- table(jcids, ids[[b]])
        ## It is faster to work on the smaller partition, but this
        ## should be ensured by the reordering ...
        C_new <- C_old <- C <- (tcrossprod(z) > 0)
        repeat {
            C_new <- (C_old %*% C) > 0
            if(all(C_new == C_old)) break
        }
        ## This should now have the connected components.
        ## Next, compute the map of the join class ids to the ids of
        ## these components. 
        cnt <- 0
        map <- remaining_ids <- seq_len(jnc)
        while(length(remaining_ids)) {
            cnt <- cnt + 1
            pos <- which(C[remaining_ids[1], remaining_ids] > 0)
            map[remaining_ids[pos]] <- cnt
            remaining_ids <- remaining_ids[-pos]
        }
        ## And update the join:
        jcids <- map[jcids]
        jnc <- cnt
    }

    cl_partition_by_class_ids(jcids)    
}

.cl_join_dendrogram <-
function(x)
{
    ## Join of an ensemble of dendrograms.
    as.cl_dendrogram(do.call(pmax, lapply(x, cl_ultrametric)))
}
