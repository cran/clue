cl_predict <-
function(object, newdata = NULL, ...)
    UseMethod("cl_predict")

## Default method.
## Should also work for kcca() from package flexclust.
cl_predict.default <-
function(object, newdata = NULL, ...)
    as.cl_membership(predict(object, newdata, ...))

## Package stats: kmeans() (R 2.1.0 or better).
cl_predict.kmeans <-
function(object, newdata = NULL, ...)
{
    if(is.null(newdata))
        return(cl_membership(object))
    as.cl_membership(max.col(- .rxdist(newdata, object$centers)))
}

## Package cluster:
## * fanny() cannot make "new" predictions.
## * clara() gives medoids, and takes metric data using Euclidean or
##   Manhattan dissimilarities (and we can figure out which by looking
##   at the call and the default values).
## * pam() gives medoids, but might have been called with dissimilarity
##   data, so is tricky.  We can always find out which by looking at the
##   medoids: as in the dissimilarity input case this is a vector of
##   class labels, and a matrix with in each row the coordinates of one
##   medoid otherwise.  We then still need to figure out whether
##   Euclidean or Manhattan distances were used by looking at the call
##   and the default values.
## Both pam() and clara() show that the interfaces could be improved to
## accomodate modern needs, e.g., for bagging.

cl_predict.fanny <-
function(object, newdata = NULL, ...)
{
    if(is.null(newdata))
        return(cl_membership(object))
    stop("Cannot make new predictions.")
}    

cl_predict.clara <-
function(object, newdata = NULL, ...)
{
    if(is.null(newdata))
        return(cl_membership(object))
    ## <FIXME>
    ## Add support eventually ...
    if(identical(object$call$stand, TRUE))
        warning("Standardization is currently not supported.")
    ## </FIXME>
    method <- object$call$metric
    if(is.null(method)) {
        ## Not given in the call, hence use default value.
        method <- formals(cluster::clara)$metric
        ## (Or hard-wire the default value: "euclidean".)
    }
    d <- .rxdist(newdata, object$medoids, method)
    as.cl_membership(max.col(-d))
}

cl_predict.pam <- 
function(object, newdata = NULL, ...)
{
    if(is.null(newdata))
        return(cl_membership(object))
    prototypes <- object$medoids
    if(!is.matrix(prototypes))
        stop("Cannot make new predictions.")
    ## <FIXME>
    ## Add support eventually ...
    if(identical(object$call$stand, TRUE))
        warning("Standardization is currently not supported.")
    ## </FIXME>
    method <- object$call$metric
    if(is.null(method)) {
        ## Not given in the call, hence use default value.
        method <- formals(cluster::pam)$metric
        ## (Or hard-wire the default value: "euclidean".)
    }
    d <- .rxdist(newdata, object$medoids, method)
    as.cl_membership(max.col(-d))
}

## Package RWeka: clusterers return objects inheriting from
## "Weka_clusterer".
cl_predict.Weka_clusterer <-
function(object, newdata = NULL, ...)
{
    if(is.null(newdata))
        return(cl_membership(object))
    as.cl_membership(predict(object, newdata = newdata,
                             type = "memberships", ...))
}

## Package cba: ccfkms().
cl_predict.ccfkms <-
function(object, newdata = NULL, ...)
{
    if(is.null(newdata))
        return(cl_membership(object))
    as.cl_membership(as.vector(predict(object, newdata)$cl))
}
## Package cba: rockCluster() returns objects of class "rock".
## Currently, no predictions.
## Note that if x is a Rock object, fitted(x) and predict(x, newdata)
## can result in missing classifications, as
##   In the case a 'drop' value greater than zero is specified, all
##   clusters with size equal or less than this value are removed from
##   the classifier. Especially, 'fitted' uses a threshold of one
##   because for singleton clusters the neighborhood is empty.
## so we need to discuss with CeeBoo whether this is really what we
## want, or if e.g. cl_predict() should use drop = 0.

## Package cclust: cclust().
cl_predict.cclust <-
function(object, newdata = NULL, ...)
{
    ## Package cclust provides predict.cclust() which returns (again) an
    ## object of class "cclust", but does not give the labels of the
    ## original data in case no new data are given.
    if(is.null(newdata))
        return(cl_membership(object))
    cl_membership(predict(object, newdata))
}

## Package e1071: cmeans() gives objects of class "fclust".
cl_predict.fclust <-
function(object, newdata = NULL, ...)
{
    if(is.null(newdata))
        return(cl_membership(object))

    ## Note that the 'fclust' objects returned by cmeans() do not always
    ## directly contain the information on the fuzzification parameter m
    ## and the distance (Euclidean/Manhattan) employed, so we have to
    ## engineer this from the matched call and the default arguments.
    nms <- names(object$call)
    ## Note that we cannot directly use object$call$m, as this could
    ## give the 'method' argument if 'm' was not given.
    m <- if("m" %in% nms)
        object$call$m
    else {
        ## Not given in the call, hence use default value.
        formals(e1071::cmeans)$m
        ## (Or hard-wire the default value: 2.)
    }
    method <- if("dist" %in% nms)
        object$call$dist
    else {
        ## Not given in the call, hence use default value.        
        formals(e1071::cmeans)$dist
        ## (Or hard-wire the default value: "euclidean".)
    }
    
    d <- .rxdist(newdata, object$centers, method)
    power <- c(m, if(method == "euclidean") 2 else 1)
    as.cl_membership(.memberships_from_cross_dissimilarities(d, power))
}

## Package e1071: cshell().
cl_predict.cshell <-
function(object, newdata = NULL, ...)
{
    if(is.null(newdata))
        return(cl_membership(object))

    ## Not surprisingly, this is rather similar to what we do for fclust
    ## objects.  Only dissimiliraties (and exponents) need to be
    ## computed differently ...
    nms <- names(object$call)
    m <- if("m" %in% nms)
        object$call$m
    else {
        ## Not given in the call, hence use default value.
        formals(e1071::cshell)$m
        ## (Or hard-wire the default value: 2.)
    }
    method <- if("dist" %in% nms)
        object$call$dist
    else {
        ## Not given in the call, hence use default value.        
        formals(e1071::cshell)$dist
        ## (Or hard-wire the default value: "euclidean".)
    }

    d <- .rxdist(newdata, object$centers, method)
    d <- sweep(d, 2, object$radius) ^ 2
    as.cl_membership(.memberships_from_cross_dissimilarities(d, m))
}

## Package e1071: bclust().
## <NOTE>
## One might argue that it would be better to use the 'dist.method'
## employed for the hierarchical clustering, but it seems that class
## labels ("clusters") are always assigned using Euclidean distances.
cl_predict.bclust <- cl_predict.kmeans
## </NOTE>

# Package flexclust: kcca() returns objects of S4 class "kcca" which
## extends S4 class "flexclust".
cl_predict.kcca <- cl_predict.default

## Package mclust: Ron Wehrens will add a predict method for Mclust
## objects eventually (as suggested by us).

## Package clue: cl_pclust().
cl_predict.cl_pclust <-
function(object, newdata = NULL, ...)
{
    if(is.null(newdata))
        return(cl_membership(object))

    d <- object$d(newdata, object$prototypes)
    power <- c(object$m, object$e)
    as.cl_membership(.memberships_from_cross_dissimilarities(d, power))
}

## Package clue: (virtual) class "cl_partition".
cl_predict.cl_partition <-
function(object, newdata = NULL, ...)
    cl_predict(.get_representation(object), newdata = newdata, ...)
