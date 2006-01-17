\name{clustering}
\alias{is.cl_partition}
\alias{is.cl_hard_partition}
\alias{is.cl_soft_partition}
\alias{is.cl_hierarchy}
\alias{as.cl_partition}
\alias{as.cl_hard_partition}
\alias{as.cl_hierarchy}
\title{Partitions and Hierarchies}
\description{Determine whether an R object represents a partition or
  hierarchy of objects, or coerce to an R object representing such.}
\usage{
is.cl_partition(x)
is.cl_hard_partition(x)
is.cl_soft_partition(x)
is.cl_hierarchy(x)

as.cl_partition(x)
as.cl_hard_partition(x)
as.cl_hierarchy(x)
}
\arguments{
  \item{x}{an R object.}
}
\value{
  For the testing functions, a logical indicating whether the given
  object represents a clustering of objects of the respective kind.

  For the coercion functions, the object itself if it already represents
  a clustering of the respective type.  Otherwise, an object inheriting
  from \code{"\link{cl_membership}"} or \code{"\link{cl_ultrametric}"}
  when coercing to partitions or hierarchies, respectively.
}
\details{
  Apart from \code{is.cl_soft_partition}, the testing functions are
  generic functions.

  The methods provided in package \pkg{clue} handle the partitions and
  hierarchies obtained from clustering functions in the base R
  distribution, as well as packages \pkg{cba}, \pkg{cclust},
  \pkg{cluster}, \pkg{e1071}, \pkg{flexclust}, \pkg{flexmix},
  \pkg{kernlab}, and \pkg{mclust} (and of course, \pkg{clue} itself).

  \code{is.cl_soft_partition} gives true iff \code{is.cl_partition} is
  true and \code{is.cl_hard_partition} is false.

  For \code{as.cl_partition} and \code{as.cl_hierarchy}, the given
  object is returned if it already represents a partition or hierarchy
  (i.e., the corresponding test returns true).  Otherwise,
  \code{\link{as.cl_membership}} or \code{\link{as.cl_ultrametric}} are
  called, creating suitable membership or ultrametric objects if
  possible.

  \code{as.cl_hard_partition(x)} returns \code{x} if this represents a
  hard partition (i.e., \code{is.cl_hard_partition(x)} is true).
  Otherwise, it returns an object of class \code{"cl_membership"} with
  the memberships of a hard partition with classes either obtained
  directly from \code{x} if this is an atomic vector of raw class ids,
  or, if \code{x} represents a soft partition or is a raw matrix of
  membership values, as the class ids of the \emph{closest hard
    partition}, defined by taking the class ids of the (first) maximal
  membership values.

  Conceptually, (hard) partitions and hierarchies are \emph{virtual}
  classes.
}
\examples{
data("Cassini")
pcl <- kmeans(Cassini$x, 3)
is.cl_partition(pcl)
is.cl_hard_partition(pcl)
is.cl_soft_partition(pcl)

hcl <- hclust(dist(USArrests))
is.cl_partition(hcl)
is.cl_hierarchy(hcl)
}
\keyword{cluster}
