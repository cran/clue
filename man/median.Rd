\name{median}
\alias{cl_median}
\title{Median Partitions and Hierarchies}
\description{
  Compute the median of an ensemble of partitions or hierarchies.  The
  median minimizes the sum of dissimilarities between itself and the
  elements of the ensemble over a suitable class of partitions or
  hierarchies.
}
\usage{
cl_median(x, method = NULL, weights = 1, control = list())
}
\arguments{
  \item{x}{an ensemble of partitions or hierarchies, or something
    coercible to that (see \code{\link{cl_ensemble}}).}
  \item{method}{a character string specifying one of the built-in
    methods for computing medians, or a function to be taken as a
    user-defined method, or \code{NULL} (default value).  If a character
    string, its lower-cased version is matched against the lower-cased
    names of the available built-in methods using \code{\link{pmatch}}.
    See \bold{Details} for available built-in methods and defaults.}
  \item{weights}{a numeric vector with non-negative case weights.
    Recycled to the number of elements in the ensemble given by \code{x}
    if necessary.}
  \item{control}{a list of control parameters.  See \bold{Details}.}
}
\value{
  The median partition or hierarchy.
}
\details{
  Median clusterings are special cases of \dQuote{consensus} clusterings
  characterized as the solutions of an optimization problem.  See Gordon
  (2001) for more information.
  
  If all elements of the ensemble are partitions, the built-in methods
  for obtaining medians proceed by minimizing \eqn{L(m) = \sum w_b
    d(x_b, m)} for a suitable dissimilarity measure \eqn{d} (see
  \code{\link{cl_dissimilarity}}) over all soft partitions with \eqn{k}
  classes, where \eqn{w_b} is the case weight given to element \eqn{x_b}
  of the ensemble.

  Available methods are as follows.
  \describe{
    \item{\code{"DWH"}}{an extension of the greedy algorithm in
      Dimitriadou, Weingessel and Hornik (2002) for approximately
      minimizing \eqn{L} with \eqn{d} being Euclidean dissimilarity.
      The reference provides some structure theory relating finding
      the median to an instance of the multiple assignment problem,
      which is known to be NP-hard, and suggests a simple heuristic
      based on successively matching an individual partitions \eqn{x_b}
      to the current approximation to the median, and compute the
      memberships of the next approximation as a weighted average of
      those of the current one and of \eqn{x_b} after permuting its
      columns for the optimal matching of class ids.

      The following control parameters are available for this method.
      \describe{
	\item{\code{k}}{an integer giving he number of classes to be
	  used for the median partition.  By default, the maximal number
	  of classes in the ensemble is used.}
	\item{\code{order}}{a permutation of the integers from 1 to the
	  size of the ensemble, specifying the order in which the
	  partitions in the ensemble should be aggregated.  Defaults to
	  using a random permutation (unlike the reference, which does
	  not permute at all).}
      }
    }

    \item{\code{"GV1"}}{the fixed-point algorithm for the \dQuote{first
	model} in Gordon and Vichi (2001) for minimizing \eqn{L} with
      \eqn{d} again being Euclidean dissimilarity.  This iterates
      between individually matching all partitions to the current
      approximation to the median, and computing the next approximation
      as a weighted average of the memberships of all partitions after
      permuting their columns for the optimal matchings of class ids.

      The following control parameters are available for this method.
      \describe{
	\item{\code{k}}{an integer giving he number of classes to be
	  used for the median partition.  By default, the maximal number
	  of classes in the ensemble is used.}
	\item{\code{maxiter}}{an integer giving the maximal number of
	  iterations to be performed.  Defaults to 100.}
	\item{\code{reltol}}{the relative convergence tolerance.
	  Defaults to \code{sqrt(.Machine$double.eps)}.}
	\item{\code{start}}{a matrix with number of rows equal to the
	  size of the cluster ensemble, and \eqn{k} columns, to be used
	  as a starting value.  By default, suitable random membership
	  matrices are used.}
	\item{\code{verbose}}{a logical indicating whether to provide
	  some output on minimization progress.  Defaults to
	  \code{getOption("verbose")}.}
      }

    }

    \item{\code{"GV3"}}{a SUMT algorithm for the \dQuote{third model} in
      Gordon and Vichi (2001) for minimizing \eqn{L} with \eqn{d} being
      co-membership dissimilarity.  See \code{\link{ls_fit_ultrametric}}
      for more information on the SUMT approach.  This optimization
      problem is equivalent to finding the membership matrix \eqn{m} for
      which the sum of the squared differences between \eqn{C(m) = m m'}
      and the weighted average co-membership matrix \eqn{\sum_b w_b
	C(m_b)} of the partitions is minimal.

      Availabe control parameters are \code{method}, \code{control},
      \code{eps}, \code{q}, and \code{verbose}, which have the same
      roles as for \code{\link{ls_fit_ultrametric}}, and the following.
      \describe{
	\item{\code{k}}{an integer giving he number of classes to be
	  used for the median partition.  By default, the maximal number
	  of classes in the ensemble is used.}
	\item{\code{start}}{a matrix with number of rows equal to the
	  size of the cluster ensemble, and \eqn{k} columns, to be used
	  as a starting value.  By default, a membership based on a rank
	  \eqn{k} approximation to the weighted average co-membership
	  matrix is used.}
      }

    }

  }
  
  By default, method \code{"DWH"} is used.
  
  If all elements of the ensemble are hierarchies, the built-in method
  (named \code{"cophenetic"}) for computing medians is based on
  minimizing \eqn{L(u) = \sum w_b d(x_b, u)} over all ultrametrics,
  where \eqn{d} is Euclidean dissimilarity.  This is equivalent to
  finding the best least squares ultrametric approximation of the
  weighted average \eqn{d = \sum w_b u_b} of the ultrametrics \eqn{u_b}
  of the hierarchies \eqn{x_b}, which is attempted by calling
  \code{\link{ls_fit_ultrametric}} on \eqn{d} with appropriate control
  parameters.

  If a user-defined agreement method is to be employed, it must be a
  function taking the cluster ensemble, the case weights, and a list of
  control parameters as its arguments.

  All built-in methods use heuristics for solving hard optimization
  problems, and cannot be guaranteed to find a global minimum.  Standard
  practice would recommend to use the best solution found in
  \dQuote{sufficiently many} replications of the methods.
}
\references{
  E. Dimitriadou and A. Weingessel and K. Hornik (2002).
  A combination scheme for fuzzy clustering.
  \emph{International Journal of Pattern Recognition and Artificial
    Intelligence}, \bold{16}, 901--912.

  A. D. Gordon and M. Vichi (2001).
  Fuzzy partition models for fitting a set of partitions.
  \emph{Psychometrika}, \bold{66}, 229--248.

  A. D. Gordon (1999).
  \emph{Classification} (2nd edition).
  Boca Raton, FL: Chapman \& Hall/CRC.
}
\seealso{
  \code{\link{cl_medoid}}
}
\examples{
## Median partition for the Rosenberg-Kim kinship terms partition
## data based on co-membership dissimilarities.
data("Kinship82")
m1 <- cl_median(Kinship82, method = "GV3",
                control = list(k = 3, verbose = TRUE))
## (Note that one should really use several replicates of this.)
## Total co-membership dissimilarity:
sum(cl_dissimilarity(Kinship82, m1, "comem"))
## Compare to the consensus solution given in Gordon & Vichi (2001).
data("Kinship82_Consensus")
m2 <- Kinship82_Consensus[["JMF"]]
sum(cl_dissimilarity(Kinship82, m2, "comem"))
## Seems we get a better solution ...
## How dissimilar are these solutions?
cl_dissimilarity(m1, m2, "comem")
## How "fuzzy" are they?
cl_fuzziness(cl_ensemble(m1, m2))
## Do the "nearest" hard partitions fully agree?
cl_dissimilarity(as.cl_hard_partition(m1),
                 as.cl_hard_partition(m2))
## Hmm ...

## Median partition for the Gordon and Vichi (2001) macroeconomic
## partition data based on Euclidean dissimilarities.
data("Macro")
set.seed(1)
m1 <- cl_median(Macro, method = "GV1",
                control = list(k = 2, verbose = TRUE))
## (Note that one should really use several replicates of this.)
## Total Euclidean dissimilarity:
sum(cl_dissimilarity(Macro, m1))
## Compare to the consensus solution given in Gordon & Vichi (2001).
data("Macro_Consensus")
m2 <- Macro_Consensus[["MF1"]]
sum(cl_dissimilarity(Macro, m2))
## Seems we get a better solution ...
## And in fact, it is qualitatively different:
table(as.cl_hard_partition(m1),
      as.cl_hard_partition(m2))
## Hmm ...
}
\keyword{cluster}
