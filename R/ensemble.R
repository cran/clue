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

    ## Previously, we used to require that the elements of the ensemble
    ## either all be partitions, or all be hierarchies.  We no longer do
    ## this, as it makes sense to also allow e.g. object dissimilarities
    ## (raw "dist" objects or additive distances) as elements (e.g.,
    ## when computing proximities), and it is rather cumbersome to
    ## decide in advance which combinations of elements might be useful
    ## and hence should be allowed.  All we enforce is that all elements
    ## correspond to the same number of objects (as we typically cannot
    ## verify that they relate to the *same* objects).  For "pure"
    ## ensembles of partitions or hierarchies we add additional class
    ## information.

    if(all(sapply(clusterings, is.cl_partition)))
        class(clusterings) <- c("cl_partition_ensemble", "cl_ensemble")
    else if(all(sapply(clusterings, is.cl_dendrogram)))
        class(clusterings) <- c("cl_dendrogram_ensemble",
                                "cl_hierarchy_ensemble", "cl_ensemble")
    else if(all(sapply(clusterings, is.cl_hierarchy)))
        class(clusterings) <- c("cl_hierarchy_ensemble", "cl_ensemble")
    else
        class(clusterings) <- "cl_ensemble"

    n <- sapply(clusterings, n_of_objects)
    if(any(diff(n)))
        stop("All elements must have the same number of objects.")
    attr(clusterings, "n_of_objects") <- as.integer(n[1])

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
    msg <- sprintf(ngettext(length(x),
                            "An ensemble of %d partition of %d objects.",
                            "An ensemble of %d partitions of %d objects."),
                   length(x), n_of_objects(x))
    writeLines(strwrap(msg))
    invisible(x)
}

Summary.cl_partition_ensemble <-
function(..., na.rm = FALSE)
{
    ok <- switch(.Generic, max = , min = , range = TRUE, FALSE)
    if(!ok)
        stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                      .Generic, .Class))
    args <- list(...)
    ## Combine the given partition ensembles.
    x <- do.call(c, args)
    switch(.Generic,
           "min" = cl_meet(x),
           "max" = cl_join(x),
           "range" = cl_ensemble(min = cl_meet(x), max = cl_join(x)))
}

print.cl_dendrogram_ensemble <-
function(x, ...)
{
    msg <- sprintf(ngettext(length(x),
                            "An ensemble of %d dendrogram of %d objects.",
                            "An ensemble of %d dendrograms of %d objects."),
                   length(x), n_of_objects(x))
    writeLines(strwrap(msg))
    invisible(x)
}

print.cl_hierarchy_ensemble <-
function(x, ...)
{
    msg <- sprintf(ngettext(length(x),
                            "An ensemble of %d hierarchy of %d objects.",
                            "An ensemble of %d hierarchies of %d objects."),
                   length(x), n_of_objects(x))
    writeLines(strwrap(msg))
    invisible(x)
}

print.cl_ensemble <-
function(x, ...)
{
    writeLines(sprintf(ngettext(length(x),
                                "An ensemble with %d element.",
                                "An ensemble with %d elements."),
                       length(x)))
    invisible(x)
}
                        
unique.cl_ensemble <-
function (x, incomparables = FALSE, ...)
    cl_ensemble(list = NextMethod("unique"))
