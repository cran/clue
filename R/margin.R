## <FIXME>
## Improve this.
## For hard partitions, all margins are one, and we can do this much
## more efficiently.
## </FIXME>

cl_margin <- function(x) {
    x <- cl_membership(x)
    i <- seq(length = nrow(x))
    j <- cbind(i, max.col(x))
    out <- x[j]
    x[j] <- 0
    out <- out - x[cbind(i, max.col(x))]
    if(!is.null(rn <- rownames(x)))
        names(out) <- rn
    out
}
