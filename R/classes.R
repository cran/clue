cl_classes <-
function(x)
    UseMethod("cl_classes")

cl_classes.default <-
function(x)
{
    ## Be nice to users ...
    if(is.cl_partition(x))
        cl_classes(as.cl_partition(x))
    else if(is.cl_dendrogram(x))
        cl_classes(as.cl_dendrogram(x))
    else
        stop("Can only determine classes of partitions or hierarchies.")
}

cl_classes.cl_partition <-
function(x)
{
    n <- n_of_objects(x)
    out <- split(seq(length = n), cl_class_ids(x))
    class(out) <- c("cl_classes_of_partition_of_objects",
                    "cl_classes_of_objects")
    attr(out, "n_of_objects") <- n
    attr(out, "labels") <- cl_object_labels(x)
    out
}

cl_classes.cl_dendrogram <-
function(x)
{
    x <- as.hclust(x)
    n <- n_of_objects(x)
    labels <- seq(length = n)
    groups <- cutree(x, labels)
    ## Give a list with the (unique) sets of numbers of the objects.
    out <- unique(unlist(sapply(split(groups, col(groups)),
                                function(k) split(labels, k)),
                         recursive = FALSE,
                         use.names = FALSE))
    ## Preserve labels if possible, and re-order according to
    ## cardinality.
    out <- out[order(sapply(out, length))]
    class(out) <- c("cl_classes_of_hierarchy_of_objects",
                    "cl_classes_of_objects")
    attr(out, "n_of_objects") <- n
    attr(out, "labels") <- x$labels
    out
}

## Be nice to users of ultrametric fitters ... which should really fit
## dendrograms.
cl_classes.cl_ultrametric <- cl_classes.cl_dendrogram

## <FIXME>
## Add a decent print method for class "cl_classes_of_objects"
## eventually, giving the usual
##   { {A}, {B}, {C}, {A, B}, {A, B, C} }
## output.
## </FIXME>
