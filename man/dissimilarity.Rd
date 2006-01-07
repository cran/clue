\name{dissimilarity}
\alias{cl_dissimilarity}
\title{Dissimilarity Between Partitions or Hierarchies}
\description{Compute the dissimilarity between (ensembles) of partitions
  or hierarchies.}
\usage{
cl_dissimilarity(x, y = NULL, method = "euclidean")
}
\arguments{
  \item{x}{an ensemble of partitions or hierarchies, or something
    coercible to that (see \code{\link{cl_ensemble}}).}
  \item{y}{\code{NULL} (default), or as for \code{x}.}
  \item{method}{a character string specifying one of the built-in
    methods for computing dissimilarity, or a function to be taken as
    a user-defined method.  If a character string, its lower-cased
    version is matched against the lower-cased names of the available
    built-in methods using \code{\link{pmatch}}.  See \bold{Details} for
    available built-in methods.}
}
\value{
  If \code{y} is \code{NULL}, an object of class
  \code{"cl_dissimilarity"} containing the dissimilarities between all
  pairs of components of \code{x}.  Otherwise, an object of class
  \code{"cl_cross_dissimilarity"} with the dissimilarities between the
  components of \code{x} and the components of \code{y}.
}
\details{
  If \code{y} is given, its components must be of the same kind as those
  of \code{x} (i.e., components must either all be partitions, or all be
  hierarchies).

  If all components are partitions, the following built-in methods for
  measuring dissimilarity between two partitions with respective
  membership matrices \eqn{u} and \eqn{v} (brought to a common number of
  columns) are available:

  \describe{
    \item{\code{"euclidean"}}{the Euclidean dissimilarity of the
      memberships, i.e., the square root of the minimal sum of the
      squared differences of \eqn{u} and all column permutations of
      \eqn{v}.  See Dimitriadou, Weingessel and Hornik (2002).}
    \item{\code{"manhattan"}}{the Manhattan dissimilarity of the
      memberships, i.e., the minimal sum of the absolute differences of
      \eqn{u} and all column permutations of \eqn{v}.}
    \item{\code{"comemberships"}}{the Euclidean dissimilarity of the
      elements of the co-membership matrices \eqn{C(u) = u u'} and
      \eqn{C(v)}, i.e., the square root of the sum of the squared
      differences of \eqn{C(u)} and \eqn{C(v)}.}
    \item{\code{"symdiff"}}{the cardinality of the symmetric set
      difference of the sets of co-classified pairs of distinct objects
      in the partitions.  I.e., the number of distinct pairs of objects
      in the same class in exactly one of the partitions.
      (Alternatively, the cardinality of the symmetric set difference
      between the (binary) equivalence relations corresponding to the
      partitions.)  For soft partitions, (currently) the symmetric set
      difference of the corresponding \dQuote{nearest} hard partitions
      is used.}
    \item{\code{"Rand"}}{the Rand distance, i.e., the rate of distinct
      pairs of objects in the same class in exactly one of the
      partitions.  For soft partitions, (currectly) the Rand distance of
      the of the corresponding \dQuote{nearest} hard partitions is used.}
    \item{\code{"GV1"}}{the square root of the dissimilarity
      \eqn{\Delta_1}{Delta_1} used for the first model in Gordon and
      Vichi (2001), i.e., the square root of the minimal sum of the
      squared differences of the \emph{matched} non-zero columns of
      \eqn{u} and \eqn{v}.}
  }

  If all components are hierarchies, available built-in methods for
  measuring agreement between two hierarchies with respective
  ultrametrics \eqn{u} and \eqn{v} are as follows.

  \describe{
    \item{\code{"euclidean"}}{the Euclidean dissimilarity of the
      ultrametrics (i.e., the square root of the sum of the squared
      differences of \eqn{u} and \eqn{v}).}
    \item{\code{"manhattan"}}{the Manhattan dissimilarity of the
      ultrametrics (i.e., the sum of the absolute differences of \eqn{u}
      and \eqn{v}).}
    \item{\code{"cophenetic"}}{\eqn{1 - c^2}, where \eqn{c} is the
      cophenetic correlation coefficient (i.e., the product-moment
      correlation of the ultrametrics).}
    \item{\code{"gamma"}}{the rate of inversions between the
      ultrametrics (i.e., the rate of pairs \eqn{(i,j)} and \eqn{(k,l)}
      for which \eqn{u_{ij} < u_{kl}} and \eqn{v_{ij} > v_{kl}}).}
    \item{\code{"symdiff"}}{the cardinality of the symmetric set
      difference of the sets of classes (hierarchies in the strict
      sense) induced by the dendrograms.  I.e., the number of sets of
      objects obtained by a split in exactly one of the hierarchies.}
  }

  If a user-defined agreement method is to be employed, it must be a
  function taking two clusterings as its arguments.

  Symmetric dissimilarity objects of class \code{"cl_dissimilarity"} are
  implemented as symmetric proximity objects with self-proximities
  identical to zero, and inherit from class \code{"cl_proximity"}.  They
  can be coerced to dense square matrices using \code{as.matrix}.  It
  is possible to use 2-index matrix-style subscripting for such objects;
  unless this uses identical row and column indices, this results in a
  (non-symmetric dissimilarity) object of class
  \code{"cl_cross_dissimilarity"}.

  Symmetric dissimilarity objects also inherit from class
  \code{"\link{dist}"} (although they currently do not \dQuote{strictly}
  extend this class), thus making it possible to use them directly for
  clustering algorithms based on dissimilarity matrices of this class,
  see the examples.
}
\references{
  E. Dimitriadou and A. Weingessel and K. Hornik (2002).
  A combination scheme for fuzzy clustering.
  \emph{International Journal of Pattern Recognition and Artificial
    Intelligence}, \bold{16}, 901--912.

  A. D. Gordon and M. Vichi (2001).
  Fuzzy partition models for fitting a set of partitions.
  \emph{Psychometrika}, \bold{66}, 229--248.
}  
\seealso{
  \code{\link{cl_agreement}}
}
\examples{
## An ensemble of partitions.
data("CKME")
pens <- CKME[1 : 30]
diss <- cl_dissimilarity(pens)
summary(c(diss))
cl_dissimilarity(pens[1:5], pens[6:7])
## Equivalently, using subscripting.
diss[1:5, 6:7]
## Can use the dissimilarities for "secondary" clustering
## (e.g. obtaining hierarchies of partitions):
hc <- hclust(diss)
plot(hc)
}
\keyword{cluster}
