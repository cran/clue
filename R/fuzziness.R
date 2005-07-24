cl_fuzziness <-
function(x, method = NULL, normalize = TRUE)
{
    x <- as.cl_ensemble(x)
    if(inherits(x, "cl_hierarchy_ensemble")) {
        ## Currently, no fuzzy hierarchies ... hence, always 0.
        out <- double(length(x))
        attr(out, "description") <- "Fuzziness"
        class(out) <- "cl_fuzziness"
        return(out)
    }

    if(!is.function(method)) {
        builtin_methods <- c("PC", "PE")
        builtin_method_names <-
            c("partition coefficient", "partition entropy")
        if(is.null(method))
            ind <- 1
        else if(is.na(ind <- pmatch(tolower(method),
                                    tolower(builtin_methods))))
            stop(gettextf("Value '%s' is not a valid abbreviation for a fuzziness method.",
                          method),
                 domain = NA)
        method <- paste(".cl_fuzziness_partition", builtin_methods[ind],
                        sep = "_")
        method_name <- builtin_method_names[ind]
        if(normalize)
            method_name <- paste("normalized", method_name)
    }
    else
        method_name <- "user-defined method"
    out <- as.numeric(sapply(x, method, normalize))
    attr(out, "description") <- paste("Fuzziness using", method_name)
    class(out) <- "cl_fuzziness"
    out 
}

.cl_fuzziness_partition_PC <-
function(x, normalize = TRUE)
{
    ## Dunn's Partition Coefficient, see also ?fanny.
    ## Note that we normalize differently ...
    if(!.maybe_is_proper_soft_partition(x) && is.cl_hard_partition(x))
        return(1 - normalize)
    pc <- sum(cl_membership(x) ^ 2) / n_of_objects(x)
    if(normalize)
        pc <- (1 - pc) / (1 - 1 / n_of_classes(x))
    pc
}

.cl_fuzziness_partition_PE <-
function(x, normalize = TRUE)
{
    ## Bezdek's Partition Entropy.
    ## Note that we normalize differently ...
    if(!.maybe_is_proper_soft_partition(x) && is.cl_hard_partition(x))
        return(0)
    M <- cl_membership(x)
    pe <- - sum(ifelse(M > 0, M * log(M), 0)) / n_of_objects(x)
    if(normalize)
        pe <- pe / log(n_of_classes(x))
    pe
}

print.cl_fuzziness <-
function(x, ...)
{
    cat(attr(x, "description"), ":\n", sep = "")
    print(as.vector(x), ...)
    invisible(x)
}
