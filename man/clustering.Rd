\name{clustering}
\alias{is.cl_partition}
\alias{is.cl_hard_partition}
\alias{is.cl_soft_partition}
\alias{is.cl_hierarchy}
\title{Partitions and Hierarchies}
\description{Determine whether an R object represents a partition or
  hierarchy of objects.}
\usage{
is.cl_partition(x)
is.cl_hard_partition(x)
is.cl_soft_partition(x)
is.cl_hierarchy(x)
}
\arguments{
  \item{x}{an R object.}
}
\value{
  A logical indicating whether the given object represents a clustering
  of objects of the respective kind.
}
\details{
  Apart from \code{is.cl_soft_partition}, these functions are generic
  functions.

  The methods provided in package \pkg{clue} handle the partitions and
  hierarchies obtained from clustering functions in the base R
  distribution, as well as packages \pkg{cclust}, \pkg{cluster},
  \pkg{e1071}, and \pkg{mclust} (and of course, \pkg{clue} itself).

  \code{is.cl_soft_partition} gives true iff \code{is.cl_partition} is
  true and \code{is.cl_hard_partition} is false.
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
