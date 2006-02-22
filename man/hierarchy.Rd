\name{hierarchy}
\alias{cl_hierarchy}		% class ...
\alias{is.cl_hierarchy}
\alias{as.cl_hierarchy}
\alias{cl_dendrogram}		% class ...
\alias{is.cl_dendrogram}
\alias{as.cl_dendrogram}
\title{Hierarchies}
\description{
  Determine whether an R object represents a hierarchy of objects, or
  coerce to an R object representing such.}
\usage{
is.cl_hierarchy(x)
is.cl_dendrogram(x)

as.cl_hierarchy(x)
as.cl_dendrogram(x)
}
\arguments{
  \item{x}{an R object.}
}
\value{
  For the testing functions, a logical indicating whether the given
  object represents a clustering of objects of the respective kind.

  For the coercion functions, a container object inheriting from
  \code{"cl_hierarchy"}, with a suitable representation of the hierarchy
  given by \code{x}.
}
\details{
  These functions are generic functions.

  The methods provided in package \pkg{clue} handle the partitions and
  hierarchies obtained from clustering functions in the base R
  distribution, as well as packages \pkg{RWeka}, \pkg{cba},
  \pkg{cclust}, \pkg{cluster}, \pkg{e1071}, \pkg{flexclust},
  \pkg{flexmix}, \pkg{kernlab}, and \pkg{mclust} (and of course,
  \pkg{clue} itself).

  The hierarchies considered by \pkg{clue} are \emph{\eqn{n}-trees}
  (hierarchies in the strict sense) and \emph{dendrograms} (also known
  as valued \eqn{n}-trees or total indexed hierarchies), which are
  represented by the virtual classes \code{"cl_hierarchy"} and
  \code{"cl_dendrogram"} (which inherits from the former),
  respectively.

  \eqn{n}-trees on a set \eqn{X} of objects correspond to collections
  \eqn{H} of subsets of \eqn{X}, usually called \emph{classes} of the
  hierarchy, which satisfy the following properties:
  \itemize{
    \item \eqn{H} contains all singletons with objects of \eqn{X},
    \eqn{X} itself, but not the empty set;
    \item The intersection of two sets \eqn{A} and \eqn{B} in \eqn{H} is
    either empty or one of the sets.
  }

  The classes of a hierarchy can be obtained by
  \code{\link{cl_classes}}.

  Dendrograms correspond to dissimilarities on \eqn{X} which satisfy the
  ultrametric (3-point) conditions \eqn{u_{ij} \le \max(u_{ik}, u_{jk})}
  for all triples \eqn{(i, j, k)} of objects.

  The ultrametric dissimilarities of a dendrogram can be obtained by
  \code{\link{cl_ultrametric}}.

  \code{as.cl_hierarchy} returns an object of class
  \code{"cl_hierarchy"} \dQuote{containing} the given object \code{x} if
  this already represents a hierarchy (i.e., \code{is.cl_hierarchy(x)}
  is true), or the ultrametric obtained from \code{x} via
  \code{\link{as.cl_ultrametric}}.

  \code{as.cl_dendrogram} returns an object which has class
  \code{"cl_dendrogram"} and inherits from \code{"cl_hierarchy"},
  and contains \code{x} if it represents a dendrogram (i.e.,
  \code{is.cl_dendrogram(x)} is true), or the ultrametric obtained from
  \code{x}.
  
  Conceptually, hierarchies and dendrograms are \emph{virtual} classes,
  allowing for a variety of representations.

  There are group methods for comparing dendrograms and computing their
  minimum, maximum, and range based on the meet and join operations, see
  \code{\link{cl_meet}}.
}
\examples{
hcl <- hclust(dist(USArrests))
is.cl_dendrogram(hcl)
is.cl_hierarchy(hcl)
}
\keyword{cluster}
