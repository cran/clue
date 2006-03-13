## * Matrix utilities

.one_entry_per_column <-
function(x, j)
{
    ## For a matrix x and a vector of column indices j_1, ..., j_n where
    ## n is the number of rows of x, get x[1,j_1], ..., x[n,j_n].
    if(!is.matrix(x))
        stop("Argument 'x' must be a matrix.")
    x[cbind(seq(length = nrow(x)), j)]
}

".one_entry_per_column<-" <-
function(x, j, value)
{
    if(!is.matrix(x))
        stop("Argument 'x' must be a matrix.")        
    x[cbind(seq(length = nrow(x)), j)] <- value
    x
}

## <FIXME>
## Transition function.
## Remove once R 2.3.0 is out and required.
.tcrossprod <- if(exists("tcrossprod", envir = baseenv())) {
    tcrossprod
} else {
    function(x, y = NULL) {
        if(is.null(y))
            x %*% t(x)
        else
            x %*% t(y)
    }
}
## </FIXME>

## * Containers

## Creator.
.make_container <-
function(x, classes, properties = NULL)
    structure(list(.Data = x, .Meta = properties),
              class = unique(classes))

## Getters.
.get_representation <-
function(x)
    x$.Data
.get_properties <-
function(x)
    x$.Meta
.get_property <-
function(x, which)
    x$.Meta[[which]]
.has_property <-
function(x, which)
    which %in% names(x$.Meta)
.get_property_from_object_or_representation <-
function(x, which, getter)
{
    if(.has_property(x, which))
        .get_property(x, which)
    else {
        if(missing(getter)) getter <- get(which)
        getter(.get_representation(x))
    }
}

## Methods (sort of).
.print_container <-
function(x, cls, ...)
{
    writeLines(gettextf("An object of virtual class '%s', with representation:\n",
                        cls))
    print(.get_representation(x), ...)
    invisible(x)
}
    
## * Others

weighted_median <-
function(x, w = 1, na.rm = FALSE)
{
    w <- rep(w, length = length(x))
    if(na.rm && any(ind <- is.na(x))) {
        x <- x[!ind]
        w <- w[!ind]
    }
    if(any(is.na(x)) || !length(x))
        return(NA)
    w <- w / sum(w)    
    ind <- order(x)
    x <- x[ind]
    w <- w[ind]
    x[which.min(x * (cumsum(w) - 0.5) - cumsum(w * x))]
}


### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
