\name{ensemble}
\alias{cl_ensemble}
\alias{as.cl_ensemble}
\alias{is.cl_ensemble}
\title{Cluster Ensembles}
\description{Creation and manipulation of cluster ensembles.}
\usage{
cl_ensemble(..., list = NULL)
as.cl_ensemble(x)
is.cl_ensemble(x)
}
\arguments{
  \item{\dots}{R objects representing clusterings of the same kind
    (either all partitions or all hierarchies) on the same objects.}
  \item{list}{a list of R objects as in \code{\dots}.}
  \item{x}{for \code{as.cl_ensemble}, an R object as in \code{\dots};
    for \code{is.cl_ensemble}, an arbitrary R object.}
}
\details{
  \code{cl_ensemble} creates cluster ensembles, which are realized as
  lists of clusterings with additional class information.   All
  clusterings in an ensemble must be of the same kind, and have the same
  number of objects.  If all clusterings are partitions, this list has
  class \code{"cl_partition_ensemble"} and inherits from
  \code{"cl_ensemble"}; if all clusterings are hierarchies, it has class
  \code{"cl_hierarchy_ensemble"} and inherits from
  \code{"cl_ensemble"}.   Empty ensembles
  cannot be categorized according to the kind of clusterings they
  contain, and hence only have class \code{"cl_ensemble"}.

  The list representation makes it possible to use \code{lapply} for
  computations on the individual clusterings in (i.e., the components
  of) a cluster ensemble.

  Available methods for cluster ensembles include those for
  subscripting, \code{c}, \code{rep}, and \code{print}.

  Note that (currently), \code{as.cl_ensemble} assumes that unclassed
  lists represent \emph{single} clusterings, as this in particular holds
  true for \code{\link{kmeans}} in versions of R prior to 2.1.0.  This
  may change eventually.
}
\value{
  \code{cl_ensemble} returns a list of the given clusterings, with
  additional class information (see \bold{Details}).
}
\examples{
d <- dist(USArrests)
hclust_methods <- c("ward", "single", "complete", "average",
                    "mcquitty", "median", "centroid")
hclust_results <- lapply(hclust_methods, function(m) hclust(d, m))
## Now create an ensemble from the results.
hens <- cl_ensemble(list = hclust_results)
names(hens) <- hclust_methods 
hens
## Subscripting.
hens[1 : 3]
## Replication.
rep(hens, 3)
## And continue to analyze the ensemble, e.g.
cl_dissimilarity(hens, method = "gamma")
}
\keyword{cluster}
