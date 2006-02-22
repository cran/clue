\name{pclust}
\alias{cl_pclust}
\title{Prototype-Based Partitions of Clusterings}
\description{
  Compute prototype-based partitions of a cluster ensemble by minimizing
  \eqn{\sum u_{bj}^m d(x_b, p_j)^e}, the sum of the membership-weighted
  \eqn{e}-th powers of the dissimilarities between the elements
  \eqn{x_b} of the ensemble and the prototypes \eqn{p_j}, for suitable
  dissimilarities \eqn{d} and exponents \eqn{e}.
}
\usage{
cl_pclust(x, k, m = 1, control = list())
}
\arguments{
  \item{x}{an ensemble of partitions or hierarchies, or something
    coercible to that (see \code{\link{cl_ensemble}}).}
  \item{k}{an integer giving the number of classes to be used in the
    partition.}
  \item{m}{a number not less than 1 controlling the softness of the
    partition (as the \dQuote{fuzzification parameter} of the fuzzy
    \eqn{c}-means algorithm).  The default value of 1 corresponds to
    hard partitions obtained from a generalized \eqn{k}-means problem;
    values greater than one give partitions of increasing softness
    obtained from a generalized fuzzy \eqn{c}-means problem.}
  \item{control}{a list of control parameters.  See \bold{Details}.}
}
\value{
  An object of class \code{"cl_pclust"} representing the obtained
  \dQuote{secondary} partition, which is a list with the following
  components.
  \item{prototypes}{a cluster ensemble with the \eqn{k} prototypes.}
  \item{membership}{an object of class \code{"\link{cl_membership}"}
    with the membership values \eqn{u_{bj}}.}
  \item{cluster}{the class ids of the nearest hard partition.}
  \item{silhouette}{Silhouette information for the partition, see
    \code{\link[cluster]{silhouette}}.} 
  \item{validity}{precomputed validity measures for the partition.}
  \item{m}{the softness control argument.}
  \item{call}{the matched call.}
  \item{d}{the dissimilarity function \eqn{d = d(x, p)} employed.}
  \item{e}{the exponent \eqn{e} employed.}
}
\details{
  For \eqn{m = 1}, a generalization of the Lloyd-Forgy variant of the
  \eqn{k}-means algorithm is used, which iterates between reclassifying
  objects to their closest prototypes, and computing new prototypes as
  consensus clusterings for the classes.  This may result in degenerate
  solutions, and will be replaced by a Hartigan-Wong style algorithm
  eventually.

  For \eqn{m > 1}, a generalization of the fuzzy \eqn{c}-means recipe
  (e.g., Bezdek (1981)) is used, which alternates between computing
  optimal memberships for fixed prototypes, and computing new prototypes
  as the consensus clusterings for the classes.

  This procedure is repeated until convergence occurs, or the maximal
  number of iterations is reached.

  Consensus clusterings are computed using \code{\link{cl_consensus}}.

  Available control parameters are as follows.
  \describe{
    \item{\code{maxiter}}{an integer giving the maximal number of
      iterations to be performed.  Defaults to 100.}
    \item{\code{reltol}}{the relative convergence tolerance.  Defaults
      to \code{sqrt(.Machine$double.eps)}.}
    \item{\code{method}}{the method to be used in
      \code{\link{cl_consensus}}.}
    \item{\code{control}}{control parameters to be used in
      \code{\link{cl_consensus}}.}
  }

  The dissimilarities \eqn{d} and exponent \eqn{e} are implied by the
  consensus method employed, and inferred via a registration mechanism
  currently only made available to built-in consensus methods.  The
  default methods compute Least Squares Euclidean consensus clusterings,
  i.e., use Euclidean dissimilarity \eqn{d} and \eqn{e = 2}.

  The fixed point approach employed is a heuristic which cannot be
  guaranteed to find the global minimum (as this is already true for the
  computation of consensus clusterings).  Standard practice would
  recommend to use the best solution found in \dQuote{sufficiently many}
  replications of the base algorithm.
}
\references{
  J. C. Bezdek (1981).
  \emph{Pattern recognition with fuzzy objective function algorithms}.
  New York: Plenum.
}

\seealso{
  \code{\link[stats]{kmeans}},
  \code{\link[e1071]{cmeans}}.
}  
\examples{
## Use a precomputed ensemble of 50 k-means partitions of the
## Cassini data.
data("CKME")
CKME <- CKME[1 : 30]		# for saving precious time ...
diss <- cl_dissimilarity(CKME)
hc <- hclust(diss)
plot(hc)
## This suggests using a partition with three classes, which can be
## obtained using cutree(hc, 3).  Could use cl_consensus() to compute
## prototypes as the least squares consensus clusterings of the classes,
## or alternatively:
set.seed(123)
x1 <- cl_pclust(CKME, 3, m = 1)
x2 <- cl_pclust(CKME, 3, m = 2)
## Agreement of solutions.
cl_dissimilarity(x1, x2)
table(cl_class_ids(x1), cl_class_ids(x2))
}
\keyword{cluster}
