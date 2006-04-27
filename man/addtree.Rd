\name{addtree}
\encoding{UTF-8}
\alias{ls_fit_addtree}
\alias{ls_fit_centroid}
\title{Least Squares Fit of Additive Tree Distances to Dissimilarities}
\description{
  Find the additive tree distance or centroid distance minimizing least
  squares distance (Euclidean dissimilarity) to a given dissimilarity
  object.
}
\usage{
ls_fit_addtree(x, method = c("SUMT", "IP", "IR"), weights = 1,
               control = list())
ls_fit_centroid(x)
}
\arguments{
  \item{x}{a dissimilarity object inheriting from class
    \code{"\link{dist}"}.}
  \item{method}{a character string indicating the fitting method to be
    employed.  Must be one of \code{"SUMT"} (default), \code{"IP"}, or
    \code{"IR"}, or a unique abbreviation thereof.}
  \item{weights}{a numeric vector or matrix with non-negative weights
    for obtaining a weighted least squares fit.  If a matrix, its
    numbers of rows and columns must be the same as the number of
    objects in \code{x}, and the lower diagonal part is used.
    Otherwise, it is recycled to the number of elements in \code{x}.}
  \item{control}{a list of control parameters.  See \bold{Details}.}
}
\value{
  An object of class \code{"cl_addtree"} containing the optimal additive
  tree distances.
}
\details{
  Additive tree distances are object dissimilarities \eqn{d} satisfying
  the so-called \emph{additive tree conditions}, also known as
  \emph{four-point conditions} \eqn{d_{ij} + d_{kl} \le \max(d_{ik} +
    d_{jl}, d_{il} + d_{jk})} for all quadruples \eqn{i, j, k, l}.
  Equivalently, for each such quadruple, the largest two values of the
  sums \eqn{d_{ij} + d_{kl}}, \eqn{d_{ik} + d_{jl}}, and \eqn{d_{il} +
    d_{jk}} must be equal.
  Centroid distances are additive tree distances where the inequalities
  in the four-point conditions are strengthened to equalities (such that
  all three sums are equal), and can be represented as \eqn{d_{ij} = g_i
    + g_j}, i.e., as sums of distances from a \dQuote{centroid}.
  See, e.g., Barthélémy and Guénoche (1991) for more details on additive
  tree distances.

  With \eqn{L(d) = \sum w_{ij} (x_{ij} - d_{ij})^2}, the problem to be
  solved by \code{ls_fit_addtree} is minimizing \eqn{L} over all
  additive tree distances \eqn{d}.  This problem is known to be NP
  hard.

  We provide three heuristics for solving this problem.

  Method \code{"SUMT"} implements the SUMT (Sequential Unconstrained
  Minimization Technique, Fiacco and McCormick, 1968) approach of de
  Soete (1983).  Incomplete dissimilarities are currently not supported.

  Methods \code{"IP"} and \code{"IR"} implement the Iterative
  Projection and Iterative Reduction approaches of Hubert and Arabie
  (1995) and Roux (1988), respectively.  Non-identical weights and
  incomplete dissimilarities are currently not supported.

  See \code{\link{ls_fit_ultrametric}} for details on these methods and
  available control parameters.

  It should be noted that all methods are heuristics which can not be
  guaranteed to find the global minimum.  Standard practice would
  recommend to use the best solution found in \dQuote{sufficiently many}
  replications of the base algorithm.

  \code{ls_fit_centroid} finds the centroid distance \eqn{d} minimizing
  \eqn{L(d)} (currently, only for the case of identical weights).  This
  optimization problem has a closed-form solution.

  There is a \code{\link{plot}} method for the fitted additive tree
  distance.
}
\references{
  J.-P. Barthélémy and A. Guénoche (1991).
  \emph{Trees and proximity representations}.
  Chichester: John Wiley \& Sons.
  ISBN 0-471-92263-3.

  A. V. Fiacco and G. P. McCormick (1968).
  \emph{Nonlinear programming: Sequential unconstrained minimization
    techniques}.
  New York: John Wiley & Sons.

  L. Hubert and P. Arabie (1995).
  Iterative projection strategies for the least squares fitting of tree
  structures to proximity data.
  \emph{British Journal of Mathematical and Statistical Psychology},
  \bold{48}, 281--317.

  M. Roux (1988).
  Techniques of approximation for building two tree structures.
  In C. Hayashi and E. Diday and M. Jambu and N. Ohsumi (Eds.),
  \emph{Recent Developments in Clustering and Data Analysis}, pages
  151--170.
  New York: Academic Press.
  
  G. de Soete (1983).
  A least squares algorithm for fitting additive trees to proximity
  data.
  \emph{Psychometrika}, \bold{48}, 621--626.
}
\keyword{cluster}
\keyword{optimize}
