cl_boot <-
function(x, B, k = NULL, 
         algorithm = if(is.null(k)) "hclust" else "kmeans",
         parameters = list(), resample = FALSE)
{
    if(resample)
        stop("Resampling is currently not supported.")
    x <- rep.int(list(x), B)
    clusterings <-
        eval(as.call(c(list(as.name("lapply"), x, algorithm),
                       if(!is.null(k)) list(k),
                       parameters)))
    cl_ensemble(list = clusterings)
}
