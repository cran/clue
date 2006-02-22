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

.make_container <-
function(x, classes)
    structure(x, class = unique(c(classes, class(x))))

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
