\name{ls_fit_ultrametric}
\encoding{UTF-8}
\alias{ls_fit_ultrametric}
\title{Least Squares Fit of Ultrametrics to Dissimilarities}
\description{
  Find the ultrametric minimizing least squares distance (Euclidean
  dissimilarity) to a given dissimilarity object.
}
\usage{
ls_fit_ultrametric(x, method = c("SUMT", "IP", "IR"), weights = 1,
                   control = list())
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
  An object of class \code{"\link{cl_ultrametric}"} containing the
  optimal ultrametric distances.
}
\details{
  With \eqn{L(u) = \sum w_{ij} (x_{ij} - u_{ij})^2}, the problem to be
  solved is minimizing \eqn{L} over all \eqn{u} satisfying the
  ultrametric constraints (i.e., for all \eqn{i, j, k}, \eqn{u_{ij} \le
    \max(u_{ik}, u_{jk})}).  This problem is known to be NP hard
  (Krivanek and Moravek, 1986).

  We provide three heuristics for solving this problem.

  Method \code{"SUMT"} implements the SUMT (Sequential Unconstrained
  Minimization Technique, Fiacco and McCormick, 1968) approach of de
  Soete (1986) which in turn simplifies the suggestions in Carroll and
  Pruzansky (1980).  (See \code{\link{sumt}} for more information on the
  SUMT approach.)  We then use a final rounding step to ensure that the
  returned object exactly satisfies the ultrametric constraints.  The
  starting value \eqn{u_0} is obtained by \dQuote{random shaking} of the
  given dissimilarity object (if not given).  If there are missing
  values in \code{x}, i.e., the given dissimilarities are
  \emph{incomplete}, we follow a suggestion of de Soete (1984), imputing
  the missing values by the weighted mean of the non-missing ones, and
  setting the corresponding weights to zero.

  Availabe control parameters are \code{method}, \code{control},
  \code{eps}, \code{q}, and \code{verbose}, which have the same roles as
  for \code{\link{sumt}}, and the following.
  
  \describe{
    \item{\code{nruns}}{an integer giving the number of runs to be
      performed.  Defaults to 1.}
    \item{\code{start}}{a single dissimilarity, or a list of
      dissimilarities to be employed as starting values.}
  }

  The default optimization using conjugate gradients should work
  reasonably well for medium to large size problems.  For \dQuote{small}
  ones, using \code{nlm} is usually faster.  Note that the number of
  ultrametric constraints is of the order \eqn{n^3}, where \eqn{n} is
  the number of objects in the dissimilarity object, suggesting to use
  the SUMT approach in favor of \code{\link[stats]{constrOptim}}.

  If starting values for the SUMT are provided via \code{start}, the
  number of starting values gives the number of runs to be performed,
  and control option \code{nruns} is ignored.  Otherwise, \code{nruns}
  starting values are obtained by random shaking of the dissimilarity to
  be fitted.  In the case of multiple SUMT runs, the (first) best
  solution found is returned.

  Method \code{"IP"} implements the Iterative Projection approach of
  Hubert and Arabie (1995).  This iteratively projects the current
  dissimilarities to the closed convex set given by the ultrametric
  constraints (3-point conditions) for a single index triple \eqn{(i, j,
    k)}, in fact replacing the two largest values among \eqn{d_{ij},
    d_{ik}, d_{jk}} by their mean.  The following control parameters can
  be provided via the \code{control} argument.

  \describe{
    \item{\code{nruns}}{an integer giving the number of runs to be
      performed.  Defaults to 1.}
    \item{\code{order}}{a permutation of the numbers from 1 to the
      number of objects in \code{x}, specifying the order in which the
      ultrametric constraints are considered, or a list of such
      permutations.}
    \item{\code{maxiter}}{an integer giving the maximal number of
      iterations to be employed.}
    \item{\code{tol}}{a double indicating the maximal convergence
      tolerance.  The algorithm stops if the total absolute change in
      the dissimilarities in an iteration is less than \code{tol}.}
    \item{\code{verbose}}{a logical indicating whether to provide some
      output on minimization progress.  Defaults to
      \code{getOption("verbose")}.}
  }

  If permutations are provided via \code{order}, the number of these
  gives the number of runs to be performed, and control option
  \code{nruns} is ignored.  Otherwise, \code{nruns} randomly generated
  orders are tried.  In the case of multiple runs, the (first) best
  solution found is returned.

  Non-identical weights and incomplete dissimilarities are currently not
  supported.

  Method \code{"IR"} implements the Iterative Reduction approach
  suggested by Roux (1988), see also Barthélémy and Guénoche (1991).
  This is similar to the Iterative Projection method, but modifies the
  dissimilarities between objects proportionally to the aggregated
  change incurred from the ultrametric projections.  Available control
  parameters are identical to those of method \code{"IP"}.

  Non-identical weights and incomplete dissimilarities are currently not
  supported.

  It should be noted that all methods are heuristics which can not be
  guaranteed to find the global minimum.  Standard practice would
  recommend to use the best solution found in \dQuote{sufficiently many}
  replications of the base algorithm.
}
\references{
  J.-P. Barthélémy and A. Guénoche (1991).
  \emph{Trees and proximity representations}.
  Chichester: John Wiley \& Sons.
  ISBN 0-471-92263-3.
  
  J. D. Carroll and S. Pruzansky (1980).
  Discrete and hybrid scaling models.
  In E. D. Lantermann and H. Feger (eds.), \emph{Similarity and Choice}.
  Bern (Switzerland): Huber.

  L. Hubert and P. Arabie (1995).
  Iterative projection strategies for the least squares fitting of tree
  structures to proximity data.
  \emph{British Journal of Mathematical and Statistical Psychology},
  \bold{48}, 281--317.
  
  M. Krivanek and J. Moravek (1986).
  NP-hard problems in hierarchical tree clustering.
  \emph{Acta Informatica}, \bold{23}, 311--323.
  
  M. Roux (1988).
  Techniques of approximation for building two tree structures.
  In C. Hayashi and E. Diday and M. Jambu and N. Ohsumi (Eds.),
  \emph{Recent Developments in Clustering and Data Analysis}, pages
  151--170.
  New York: Academic Press.

  G. de Soete (1984).
  Ultrametric tree representations of incomplete dissimilarity data.
  \emph{Journal of Classification}, \bold{1}, 235--242.

  G. de Soete (1986).
  A least squares algorithm for fitting an ultrametric tree to a
  dissimilarity matrix.
  \emph{Pattern Recognition Letters}, \bold{2}, 133--137.
}
\seealso{
  \code{\link{cl_consensus}} for computing least squares consensus
  hierarchies by least squares fitting of average ultrametric
  distances.
}
\examples{
## Least squares fit of an ultrametric to the Miller-Nicely consonant
## phoneme confusion data.
data("Phonemes")
## Note that the Phonemes data set has the consonant misclassification
## probabilities, i.e., the similarities between the phonemes.
d <- as.dist(1 - Phonemes)
u <- ls_fit_ultrametric(d, control = list(verbose = TRUE))
## Cophenetic correlation:
cor(d, u)
## Plot:
plot(u)
## ("Basically" the same as Figure 1 in de Soete (1986).)
}
\keyword{cluster}
\keyword{optimize}
