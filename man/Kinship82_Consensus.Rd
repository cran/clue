\name{Kinship82_Consensus}
\alias{Kinship82_Consensus}
\title{Gordon-Vichi Kinship82 Consensus Partition Data}
\description{
  The soft (\dQuote{fuzzy}) consensus partitions for the Rosenberg-Kim
  kinship terms partition data given in Gordon and Vichi (2001).
}
\usage{data("Kinship82_Consensus")}
\format{
  A named cluster ensemble of three soft partitions of the 15 kinship
  terms into three classes.
}
\details{
  The elements of the ensemble are named \code{"MF1"}, \code{"MF2"}, and
  \code{"JMF"}, and correspond to the median partitions obtained by
  applying models 1, 2, and 3 in Gordon and Vichi (2001) to the kinship
  terms partition data in Rosenberg (1982), which are available as data
  set \code{\link{Kinship82}}.
}
\source{
  Table 6 in Gordon and Vichi (2001).
}
\references{
  A. D. Gordon and M. Vichi (2001).  
  Fuzzy partition models for fitting a set of partitions.
  \emph{Psychometrika}, \bold{66}, 229--248.

  S. Rosenberg (1982).
  The method of sorting in multivariate research with applications
  selected from cognitive psychology and person perception.
  In N. Hirschberg and L. G. Humphreys (eds.),
  \emph{Multivariate Applications in the Social Sciences}, 117--142.
  Hillsdale, NJ: Erlbaum.
}
\examples{
## Load the median consensus partitions.
data("Kinship82_Consensus")
## Fuzziness using the Partition Coefficient.
cl_fuzziness(Kinship82_Consensus)
## (Correponds to 1 - F in the source.)
## Dissimilarities:
cl_dissimilarity(Kinship82_Consensus)
cl_dissimilarity(Kinship82_Consensus, method = "comem")
}
\keyword{datasets}
