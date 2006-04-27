### * is.cl_hierarchy

## Determine whether an object is a hierarchy.
## Note that strictly speaking we only deal with indexed hierarchies so
## that there is always a corresponding ultrametric (and dendrogram).

is.cl_hierarchy <-
function(x)
    UseMethod("is.cl_hierarchy")

## Default method.
is.cl_hierarchy.default <- .false

## Package stats: hclust().
is.cl_hierarchy.hclust <- .true

## Package cluster: agnes() and diana() give objects inheriting from
## class "twins".
is.cl_hierarchy.twins <- .true
## Package cluster: mona().
is.cl_hierarchy.mona <- .true

## Package clue: (virtual) class "cl_hierarchy".
## Note that "raw" cl_ultrametric objects are *not* hierarchies, as
## these are meant for numeric computations.
is.cl_hierarchy.cl_hierarchy <- .true

### * as.cl_hierarchy

## Note that cl_partition conceptually is a virtual class, so there are
## no prototypes and no cl_hierarchy() creator.

.cl_hierarchy_classes <- "cl_hierarchy"

as.cl_hierarchy <-
function(x)
{
    if(is.cl_hierarchy(x)) {
        if(!inherits(x, "cl_hierarchy"))
            .make_container(x, .cl_hierarchy_classes)
        else
            x
    }
    else
        .make_container(as.cl_ultrametric(x),
                        .cl_hierarchy_classes)
}

### * print.cl_hierarchy

print.cl_hierarchy <-
function(x, ...)
    .print_container(x, "cl_hierarchy", ...)

### * Complex.cl_hierarchy

## No Complex() for any kind of hierarchy.
Complex.cl_hierarchy <-
function(z)
    stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                  .Generic, .Class))

### * Math.cl_hierarchy

## No Math() for any kind of hierarchy.
Math.cl_hierarchy <-
function(x, ...)
    stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                  .Generic, .Class))

### * is.cl_dendrogram

## <FIXME>
## Once we have cl_dendrogram testing, we can simplify cl_hierarchy
## testing.  E.g.,
##   is.cl_hierachy.default <- is.cl_dendrogram
## should be ok, and we can add cl_hierarchy predicates for hierarchies
## which are not dendrograms on top of that.
## </FIXME>

is.cl_dendrogram <-
function(x)
    UseMethod("is.cl_dendrogram")
## Default method.
is.cl_dendrogram.default <- .false
## Package stats: hclust().
is.cl_dendrogram.hclust <- .true
## Package cluster: agnes() and diana() give objects inheriting from
## class "twins".
is.cl_dendrogram.twins <- .true
## Package cluster: mona().
is.cl_dendrogram.mona <- .true
## Package clue: (virtual) class "cl_dendrogram".
is.cl_dendrogram.cl_dendrogram <- .true

### * as.cl_dendrogram

.cl_dendrogram_classes <- c("cl_dendrogram", "cl_hierarchy")

as.cl_dendrogram <-
function(x)
{
    if(is.cl_dendrogram(x)) {
        if(!inherits(x, "cl_dendrogram"))
            .make_container(x, .cl_dendrogram_classes)
        else
            x
    }
    else
        .make_container(as.cl_ultrametric(x),
                        .cl_dendrogram_classes)
}

### * print.cl_dendrogram

print.cl_dendrogram <-
function(x, ...)
    .print_container(x, "cl_dendrogram", ...)

### * plot.cl_dendrogram

plot.cl_dendrogram <-
function(x, ...)
    plot(cl_ultrametric(.get_representation(x)), ...)

### * Group methods for cl_dendrogram objects.

Ops.cl_dendrogram <-
function(e1, e2)
{
    if(nargs() == 1)
        stop(gettextf("Unary '%s' not defined for \"%s\" objects.",
                      .Generic, .Class))

    ## Only comparisons are supprorted.
    if(!(as.character(.Generic) %in% c("<", "<=", ">", ">=",
                                       "==", "!=")))
        stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                      .Generic, .Class))

    u1 <- cl_ultrametric(e1)
    u2 <- cl_ultrametric(e2)
    if(length(u1) != length(u2))
        stop("Dendrograms must have the same number of objects.")
    switch(.Generic,
           "<=" = all(u1 <= u2),
           "<"  = all(u1 <= u2) && any(u1 < u2),
           ">=" = all(u1 >= u2),
           ">"  = all(u1 >= u2) && any(u1 > u2),
           "==" = all(u1 == u2),
           "!=" = any(u1 != u2))
}

### * Summary.cl_dendrogram

## <NOTE>
## This is really the same as Summary.cl_partition().
## Of course, we could make things slighly more efficient by calling
## .cl_meet_dendrogram() and .cl_join_dendrogram() directly instead of
## calling cl_meet() and cl_join().
## </NOTE>

Summary.cl_dendrogram <-
function(x, ...)
{
    ok <- switch(.Generic, max = , min = , range = TRUE, FALSE)
    if(!ok)
        stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                      .Generic, .Class))
    args <- list(x, ...)
    ## Remove 'na.rm' component ...
    args$na.rm <- NULL
    switch(.Generic,
           "min" = cl_meet(cl_ensemble(list = args)),
           "max" = cl_join(cl_ensemble(list = args)),
           "range" = {
               cl_ensemble(min = cl_meet(cl_ensemble(list = args)),
                           max = cl_join(cl_ensemble(list = args)))
           })
}

### * as.hclust.cl_dendrogram

as.hclust.cl_dendrogram <-
function(x, ...)
    as.hclust(.get_representation(x), ...)
  

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
