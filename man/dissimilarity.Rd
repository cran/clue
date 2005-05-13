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
      memberships, i.e., the minimal sum of the squared differences of
      \eqn{u} and all column permutations of \eqn{v}, divided by twice
      the number of objects in the partitions.  See Dimitriadou,
      Weingessel and Hornik (2002).}
    \item{\code{"comemberships"}}{the Euclidean dissimilarity of the
      elements of the co-membership matrices \eqn{C(u) = u u'} and
      \eqn{C(v)}, i.e., the sum of the squared differences of \eqn{C(u)}
      and \eqn{C(v)}, divided by the squared number of objects in the
      partitions.}
  }

  If all components are hierarchies, available built-in methods for
  measuring agreement between two hierarchies with respective
  ultrametrics \eqn{u} and \eqn{v} are as follows.

  \describe{
    \item{\code{"euclidean"}}{the Euclidean dissimilarity of the
      ultrametrics (i.e., the sum of the squared differences of \eqn{u}
      and \eqn{v}).}
    \item{\code{"cophenetic"}}{\eqn{1 - c^2}, where \eqn{c} is the
      cophenetic correlation coefficient (i.e., the product-moment
      correlation of the ultrametrics).}
    \item{\code{"gamma"}}{the rate of inversions between the
      ultrametrics (i.e., the rate of pairs \eqn{(i,j)} and \eqn{(k,l)}
      for which \eqn{u_{ij} < u_{kl}} and \eqn{v_{ij} > v_{kl}}).}
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
