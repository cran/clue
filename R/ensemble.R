cl_ensemble <-
function(..., list = NULL)
{
    clusterings <- c(list(...), list)

    if(!length(clusterings)) {
        ## Return an empty cl_ensemble.
        ## In this case, we cannot additionally know whether it contains
        ## partitions or hierarchies ...
        attr(clusterings, "n_of_objects") <- as.integer(NA)
        class(clusterings) <- "cl_ensemble"
        return(clusterings)
    }

    ## Otherwise, the elements of the ensemble must either all be
    ## partitions, or all be hierarchies.
    if(all(sapply(clusterings, is.cl_partition)))
        class(clusterings) <- c("cl_partition_ensemble", "cl_ensemble")
    else if(all(sapply(clusterings, is.cl_hierarchy)))
        class(clusterings) <- c("cl_hierarchy_ensemble", "cl_ensemble")
    else
        stop("Ensemble elements must be all partitions or all hierarchies.")

    n <- sapply(clusterings, n_of_objects)
    if(any(diff(n)))
        stop("All clusterings must have the same number of objects.")
    attr(clusterings, "n_of_objects") <- n[1]

    clusterings
}

is.cl_ensemble <-
function(x)
    inherits(x, "cl_ensemble")
as.cl_ensemble <-
function(x)
    if(is.cl_ensemble(x)) x else cl_ensemble(x)

c.cl_ensemble <-
function(..., recursive = FALSE)
{
    clusterings <- unlist(lapply(list(...), as.cl_ensemble),
                          recursive = FALSE)
    cl_ensemble(list = clusterings)
}

"[.cl_ensemble" <-
function(x, i)
{
    clusterings <- unclass(x)[i]
    ## <NOTE>
    ## This is not very efficient.
    ## But note that we need to special-case zero length clusterings.
    cl_ensemble(list = clusterings)
    ## </NOTE>
}

rep.cl_ensemble <-
function(x, times, ...)
    cl_ensemble(list = NextMethod("rep"))

print.cl_partition_ensemble <-
function(x, ...)
{
    writeLines(strwrap(paste("An ensemble of", length(x),
                             if(length(x) == 1)
                             "partition of" else "partitions of",
                             n_of_objects(x),
                             "objects")))
    invisible(x)
}

print.cl_hierarchy_ensemble <-
function(x, ...)
{
    writeLines(strwrap(paste("An ensemble of", length(x),
                             if(length(x) == 1)
                             "hierarchy of" else "hierarchies of",
                             n_of_objects(x),
                             "objects")))
    invisible(x)
}
