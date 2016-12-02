\name{GVME_Consensus}
\alias{GVME_Consensus}
\title{Gordon-Vichi Macroeconomic Consensus Partition Data}
\description{
  The soft (\dQuote{fuzzy}) consensus partitions for the macroeconomic
  partition data given in Gordon and Vichi (2001).
}
\usage{data("GVME_Consensus")}
\format{
  A named cluster ensemble of eight soft partitions of 21 countries
  terms into two or three classes.
}
\details{
  The elements of the ensemble are consensus partitions for the
  macroeconomic partition data in Gordon and Vichi (2001), which are
  available as data set \code{\link{GVME}}.  Element names are of the
  form \code{"\var{m}/\var{k}"}, where \var{m} indicates the consensus
  method employed (one of \samp{MF1}, \samp{MF2}, \samp{JMF}, and
  \samp{S&S}, corresponding to the application of models 1, 2, and 3
  in Gordon and Vichi (2001) and the approach in Sato and Sato (1994),
  respectively), and \var{k} denotes the number classes (2 or 3).
}
\source{
  Tables 4 and 5 in Gordon and Vichi (2001).
}
\references{
  A. D. Gordon and M. Vichi (2001).
  Fuzzy partition models for fitting a set of partitions.
  \emph{Psychometrika}, \bold{66}, 229--248.
  \doi{10.1007/BF02294837}.

  M. Sato and Y. Sato (1994).
  On a multicriteria fuzzy clustering method for 3-way data.
  \emph{International Journal of Uncertainty, Fuzziness and
    Knowledge-Based Systems}, \bold{2}, 127--142. \cr
  \doi{10.1142/S0218488594000122}.
}
\examples{
## Load the consensus partitions.
data("GVME_Consensus")
## Pick the partitions into 2 classes.
GVME_Consensus_2 <- GVME_Consensus[1 : 4]
## Fuzziness using the Partition Coefficient.
cl_fuzziness(GVME_Consensus_2)
## (Corresponds to 1 - F in the source.)
## Dissimilarities:
cl_dissimilarity(GVME_Consensus_2)
cl_dissimilarity(GVME_Consensus_2, method = "comem")
}
\keyword{datasets}
