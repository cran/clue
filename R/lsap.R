solve_LSAP <-
function(x, maximum = FALSE)
{
    if(!is.matrix(x)
       || ((n <- nrow(x)) != ncol(x))
       || any(x < 0))
        stop("x must be a square matrix with nonnegative entries.")
    if(maximum) x <- max(x) - x
    storage.mode(x) <- "double"
    out <- .C("solve_LSAP", x, n, p = integer(n),
              PACKAGE = "clue")$p + 1
    class(out) <- "solve_LSAP"
    out
}

print.solve_LSAP <-
function(x, ...)
{
    writeLines(c("Optimal assignment:",
                 gsub("x", " ",
                      strwrap(paste(seq(length = length(x)), x,
                                    sep = "x=>x", collapse = ", ")))))
    invisible(x)
}
