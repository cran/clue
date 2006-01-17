\name{n_of_objects}
\alias{n_of_objects}
\title{Number of Objects in a Partition or Hierarchy}
\description{Determine the number of objects from which a partition or
  hierarchy was obtained.}
\usage{
n_of_objects(x)
}
\arguments{
  \item{x}{an object representing a (hard of soft) partition or a
    hierarchy.}
}
\value{
  An integer giving the number of objects in the clustering.
}
\details{
  This is a generic function.

  The methods provided in package \pkg{clue} handle the partitions and
  hierarchies obtained from clustering functions in the base R
  distribution, as well as packages \pkg{cba}, \pkg{cclust},
  \pkg{cluster}, \pkg{e1071}, \pkg{flexclust}, \pkg{flexmix},
  \pkg{kernlab}, and \pkg{mclust} (and of course, \pkg{clue} itself).
}
\seealso{
  \code{\link{is.cl_partition}},
  \code{\link{is.cl_hierarchy}}
}
\examples{
data("Cassini")
pcl <- kmeans(Cassini$x, 3)
n_of_objects(pcl)
hcl <- hclust(dist(USArrests))
n_of_objects(hcl)
}
\keyword{cluster}
