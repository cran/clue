\name{solve-LSAP}
\alias{solve_LSAP}
\title{Solve Linear Sum Assignment Problem}
\description{
  Solve the linear sum assignment problem using the Hungarian method.
}
\usage{
solve_LSAP(x, maximum = FALSE)
}
\arguments{
  \item{x}{a square matrix with nonnegative entries.}
  \item{maximum}{a logical indicating whether to minimize of maximize
    the sum of assigned costs.}
}
\details{
  If \eqn{n} is the number of rows and columns of \code{x},
  \code{solve_LSAP} finds a permutation \code{p} of the numbers from 1
  to \eqn{n} such that \eqn{\sum_{i=1}^n x[i, p[i]]} is minimized or
  maximized.

  This permutation can be found using a linear program (and package
  \pkg{lpSolve} provides a function \code{lp.assign} for doing so), but
  typically more efficiently and provably in polynomial time
  \eqn{O(n^3)} using primal-dual methods such as the so-called Hungarian
  method (see the references).
}
\value{
  An object of class \code{"solve_LSAP"} with the optimal assignment of
  rows to columns.
}
\references{
  C. Papadimitriou and K. Steiglitz (1982)
  \emph{Combinatorial Optimization: Algorithms and Complexity}.
  Englewood Cliffs: Prentice Hall.
}
\author{
  Walter Boehm \email{Walter.Boehm@wu-wien.ac.at} kindly provided C code
  implementing the Hungarian method.
}
\examples{
x <- matrix(c(5, 1, 4, 3, 5, 2, 2, 4, 4), nr = 3)
solve_LSAP(x)
solve_LSAP(x, max = TRUE)
}
\keyword{optimize}
