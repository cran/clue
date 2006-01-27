cl_meet <-
function(x, y)
{
    ## General case.
    ## x either a partition ensemble, or x and y two partitions with the
    ## same number of objects.
    if(!inherits(x, "cl_partition_ensemble")) {
        ## Be nice about error messages.
        if(!all(c(is.cl_partition(x), is.cl_partition(y))))
            stop("'x' and 'y' must both be partitions.")
        if(n_of_objects(x) != n_of_objects(y))
            stop("'x' and 'y' must have the same number of objects")
        x <- cl_ensemble(x, y)
    }

    x <- unique(x)
    if(length(x) == 1)
        return(cl_membership(x[[1]]))

    ids <- seq(length = n_of_objects(x[[1]]))
    ## Cross-classify the objects.
    z <- split(ids, lapply(x, cl_class_ids))
    ## Subscript on the non-empty cells to get adjacent class ids.
    lens <- sapply(z, length)
    pos <- which(lens > 0)
    ids[unlist(z, use.names = FALSE)] <-
        rep(seq(along = z[pos]), lens[pos])
    .cl_membership_from_class_ids(ids)
}

cl_join <-
function(x, y)
{
    ## General case.
    ## x either a partition ensemble, or x and y two partitions with the
    ## same number of objects.
    if(!inherits(x, "cl_partition_ensemble")) {
        ## Be nice about error messages.
        if(!all(c(is.cl_partition(x), is.cl_partition(y))))
            stop("'x' and 'y' must both be partitions.")
        if(n_of_objects(x) != n_of_objects(y))
            stop("'x' and 'y' must have the same number of objects")
        x <- cl_ensemble(x, y)
    }

    x <- unique(x)
    if(length(x) == 1)
        return(cl_membership(x[[1]]))

    ## Canonicalize: ensure that class ids are always the integers from
    ## one to the number of classes.
    n <- sapply(x, n_of_classes)
    ids <- mapply(function(p, ncp) match(cl_class_ids(p),
                                         seq(length = ncp)),
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
        C_new <- C_old <- C <- (z %*% t(z)) > 0
        repeat {
            C_new <- (C_old %*% C) > 0
            if(all(C_new == C_old)) break
        }
        ## This should now have the connected components.
        ## Next, compute the map of the join class ids to the ids of
        ## these components. 
        cnt <- 0
        map <- remaining_ids <- seq(length = jnc)
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

    .cl_membership_from_class_ids(jcids)    
}

## <FIXME>
## These methods "work", but really do too much.
## The current design is really based on the fact that memberships and
## ultrametrics are numeric objects with certain attributes, and that
## one can thus "directly" perform numeric computations on them.
## In fact, the methods should be fine for basically all cl_partition
## objects which are *not* cl_membership objects ...
## </FIXME>

## Ops.cl_membership <-
## function(e1, e2)
## {
##     if(nargs() == 1)
##         stop(gettextf("Unary '%s' not defined for \"%s\" objects.",
##                       .Generic, .Class))
##
##     ## Only comparisons are supprorted.
##     if(!(as.character(.Generic) %in% c("<", "<=", ">", ">=",
##                                        "==", "!=")))
##         stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
##                       .Generic, .Class))
##
##     ci1 <- cl_class_ids(e1)
##     ci2 <- cl_class_ids(e2)
##     if(length(ci1) != length(ci2))
##         stop("Partitions must have the same number of objects.")
##     z <- table(ci1, ci2) > 0
##     switch(.Generic,
##            "<=" = all(rowSums(z) == 1),
##            "<"  = all(rowSums(z) == 1) && any(colSums(z) > 1),
##            ">=" = all(colSums(z) == 1),
##            ">"  = all(colSums(z) == 1) && any(rowSums(z) > 1),
##            "==" = all(rowSums(z) == 1) && all(colSums(z) == 1),
##            "!=" = any(rowSums(z) > 1) || any(colSums(z) > 1))
##
## }

## Summary.cl_membership <-
## function(x, ...)
## {
##     ok <- switch(.Generic, max = , min = , range = TRUE, FALSE)
##     if(!ok)
##         stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
##                       .Generic, .Class))
##     args <- list(x, ...)
##     ## Remove 'na.rm' component:
##     args$na.rm <- NULL
##     switch(.Generic,
##            "min" = cl_meet(cl_ensemble(list = args)),
##            "max" = cl_join(cl_ensemble(list = args)),
##            "range" = {
##                cl_ensemble(min = cl_meet(cl_ensemble(list = args)),
##                            max = cl_join(cl_ensemble(list = args)))
##            })
## }
