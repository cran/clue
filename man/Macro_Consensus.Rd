\name{Macro_Consensus}
\alias{Macro_Consensus}
\title{Gordon-Vichi Macroeconomic Consensus Partition Data}
\description{
  The soft (\dQuote{fuzzy}) consensus partitions for the macroeconomic
  partition data given in Gordon and Vichi (2001).
}
\usage{data("Macro_Consensus")}
\format{
  A named cluster ensemble of four soft partitions of 21 countries
  terms into two classes.
}
\details{
  The elements of the ensemble are named \code{"MF1"}, \code{"MF2"},
  \code{"JMF"}, and \code{"S&S"}, and correspond to the median
  partitions obtained by applying models 1, 2, and 3 in Gordon and Vichi
  (2001) and the approach in Sato and Sato (1994) to the macroeconomic
  partition data in Gordon and Vichi (2001), which are available as data
  set \code{\link{Macro}}.
}
\source{
  Table 4 in Gordon and Vichi (2001).
}
\references{
  A. D. Gordon and M. Vichi (2001).
  Fuzzy partition models for fitting a set of partitions.
  \emph{Psychometrika}, \bold{66}, 229--248.

  M. Sato and Y. Sato (1994).
  On a multicriteria fuzzy clustering method for 3-way data.
  \emph{International Journal of Uncertainty, Fuzziness and
    Knowledge-Based Systems}, \bold{2}, 127--142.
}
\examples{
## Load the median consensus partitions.
data("Macro_Consensus")
## Fuzziness using the Partition Coefficient.
cl_fuzziness(Macro_Consensus)
## (Correponds to 1 - F in the source.)
## Dissimilarities:
cl_dissimilarity(Macro_Consensus)
cl_dissimilarity(Macro_Consensus, method = "comem")
}
\keyword{datasets}
