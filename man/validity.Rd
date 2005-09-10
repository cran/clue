\name{validity}
\alias{cl_validity}
\alias{cl_validity.default}
\title{Validity Measures for Partitions and Hierarchies}
\description{
  Compute validity measures for partitions and hierarchies, attempting
  to measure how well these clusterings capture the underlying structure
  in the data they were obtained from.
}
\usage{
cl_validity(x, ...)
\method{cl_validity}{default}(x, d, ...)
}
\arguments{
  \item{x}{an object representing a partition or hierarchy.}
  \item{d}{a dissimilarity object from which \code{x} was obtained.}
  \item{\dots}{arguments to be passed to or from methods.}
}
\value{
  A list of class \code{"cl_validity"} with the computed validity
  measures.
}
\details{
  \code{cl_validity} is a generic function.

  For partitions, its default method gives the \dQuote{dissimilarity
    accounted for}, defined as \eqn{1 - a_w / a_t}, where \eqn{a_t} is
  the average total dissimilarity, and the \dQuote{average within
    dissimilarity} \eqn{a_w} is given by
  \deqn{\sum_{i,j} \sum_k m_{ik}m_{jk} d_{ij} /
        \sum_{i,j} \sum_k m_{ik}m_{jk}} 
  where \eqn{d} and \eqn{m} are the dissimilarities and memberships,
  respectively, and the sums are over all pairs of objects and all
  classes.

  For hierarchies, currently no validity measures are computed by
  default.  For the results of using \code{\link[cluster]{agnes}} and
  \code{\link[cluster]{diana}}, the agglomerative and divisive
  coefficients are provided.
}
\seealso{
  \code{\link[fpc]{cluster.stats}} in package \pkg{fpc} for a variety of
  cluster validation statistics;
  \code{\link[e1071]{fclustIndex}} in package \pkg{e1071} for several
  fuzzy cluster indexes;
  \code{\link[cclust]{clustIndex}} in package \pkg{cclust};
  \code{\link[cluster]{silhouette}} in package \pkg{cluster}.
}
\keyword{cluster}
