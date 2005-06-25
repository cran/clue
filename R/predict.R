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

## Package e1071: cmeans(), cshell(), and ufcl() give objects of class
## "fclust".
cl_predict.fclust <-
function(object, newdata = NULL, ...)
{
    if(is.null(newdata))
        return(cl_membership(object))

    ## ARGH, the 'fclust' objects returned by cmeans() do not directly
    ## contain the information on the fuzzification parameter m and the
    ## distance (euclidean/manhattan) employed, so we have to engineer
    ## this from the matched call and the default arguments.
    nms <- names(object$call)
    ## Note that we cannot directly use object$call$m, as this could
    ## the 'method' argument if 'm' was not given.
    m <- if("m" %in% nms)
        object$call$m
    else
        formals(cmeans)$m
    method <- if("dist" %in% nms)
        object$call$dist
    else
        formals(cmeans)$dist

    d <- .rxdist(newdata, object$centers, method)    
    power <- c(m, if(method == "euclidean") 2 else 1)

    as.cl_membership(.memberships_from_cross_dissimilarities(d, power))
}

## <FIXME>
## Currently, there is no support for partitioning algorithms from
## packages cluster and Mclust.
##
## Re cluster:
## * fanny() cannot make predictions.
## * clara() gives medoids, and takes metric data using Euclidean or
##   Manhattan dissimilarities, so should be straightforward.
## * pam() gives medoids, but might have been called with dissimilarity
##   data, so is tricky.  We can always find out which by looking at the
##   medoids: as in the dissimilarity input case this is a vector of
##   class labels, and a matrix with in each row the coordinates of one
##   medoid otherwise.  We then still need to figure out whether
##   Euclidean or Manhattan distances were used by looking at the call
##   and the default values.
## Both pam() and clara() show that the interfaces could be improved to
## accomodate modern needs, e.g. for bagging.
##
## Re mclust: not clear why there are no predict() methods for Mclust
## objects ...?
## </FIXME>
