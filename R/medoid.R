cl_medoid <-
function(x, method = "euclidean")
{
    ## <NOTE>
    ## In principle we can get the same using pam(k = 1)$medoids.
    ## </NOTE>

    clusterings <- as.cl_ensemble(x)
    if(!length(clusterings))
        stop("Cannot compute medoid of empty ensemble.")
    dissimilarities <-
        as.matrix(cl_dissimilarity(clusterings, method = method))
    clusterings[[which.min(rowSums(dissimilarities))]]
}
