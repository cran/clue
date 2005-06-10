\name{ls_fit_ultrametric}
\alias{ls_fit_ultrametric}
\title{Least Squares Fit of Ultrametrics to Dissimilarities}
\description{
  Find the ultrametric minimizing least squares distance (Euclidean
  dissimilarity) to a given dissimilarity object.
}
\usage{
ls_fit_ultrametric(x, control = list())
}
\arguments{
  \item{x}{a dissimilarity object inheriting from class
    \code{"\link{dist}"}.}
  \item{control}{a list of control parameters.  See \bold{Details}.}
}
\value{
  An object of class \code{"\link{cl_ultrametric}"} containing the
  optimal ultrametric distances.
}
\details{
  With \eqn{L(u) = \sum (x - u)^2}, the problem to be solved is
  minimizing \eqn{L} over all \eqn{u} satisfying the ultrametric
  constraints (i.e., for all \eqn{i, j, k}, \eqn{u_{ij} \le \max(u_{ik},
    u_{jk})}).  This problem is known to be NP hard (Krivanek and
  Moravek, 1986).

  We follow de Soete (1986) to use a heuristic based on an SUMT
  (Sequential Unconstrained Minimization Technique) approach in turn
  simplifying the suggestions in Carroll and Pruzansky (1980).  One
  iteratively minimizes \eqn{L(u) + \rho_k P(u)}, where \eqn{P(u)} is a
  non-negative function penalizing violations of the ultrametric 
  constraints such that \eqn{P(u)} is zero iff \eqn{u} is an
  ultrametric.  The \eqn{\rho} values are increased according to the
  rule \eqn{\rho_{k+1} = q \rho_k} for some constant \eqn{q > 1}, until
  convergence is obtained in the sense that the Euclidean distance
  between successive solutions \eqn{u_k} and \eqn{u_{k+1}} is small
  enough.  We then use a final rounding step to ensure that the returned
  object exactly satisfies the ultrametric constraints.  The starting
  value \eqn{u_0} is obtained by \dQuote{random shaking} of the given
  dissimilarity object.

  The unconstrained minimizations are carried out using either
  \code{\link[stats]{optim}} or \code{\link[stats]{nlm}}, using the
  analytic gradients given in Carroll and Pruzansky (1980).  The
  following control parameters can be provided via the \code{control}
  argument.

  \describe{
    \item{\code{method}}{a character string, or \code{NULL}.  If not
      given, \code{"CG"} is used.  If equal to \code{"nlm"},
      minimization is carried out using \code{\link[stats]{nlm}}.
      Otherwise, \code{\link[stats]{optim}} is used with \code{method}
      as the given method.}
    \item{\code{control}}{a list of control parameters to be passed to
      the minimization routine in case \code{optim} is used.}
    \item{\code{eps}}{the absolute convergence tolerance.
      Defaults to \code{sqrt(.Machine$double.eps)}.}
    \item{\code{q}}{a double greater than one controlling the growth of
      the \eqn{\rho_k} as described above.  Defaults to 10.}
    \item{\code{verbose}}{a logical indicating whether to provide some
      output on minimization progress.  Defaults to
      \code{getOption("verbose")}.}
  }

  The default optimization using conjugate gradients should work
  reasonably well for medium to large size problems.  For \dQuote{small}
  ones, using \code{nlm} is usually faster.  Note that the number of
  ultrametric constraints is of the order \eqn{n^3}, where \eqn{n} is
  the number of objects in the dissimilarity object, suggesting to use
  the SUMT approach in favor of \code{\link[stats]{constrOptim}}.

  It should be noted that the SUMT approach is a heuristic which can not
  be guaranteed to find the global minimum.  Standard practice would
  recommend to use the best solution found in \dQuote{sufficiently many}
  replications of the base algorithm.
}
\references{
  J. D. Carroll and S. Pruzansky (1980).
  Discrete and hybrid scaling models.
  In E. D. Lantermann and H. Feger (eds.), \emph{Similarity and Choice}.
  Bern (Switzerland): Huber.
  
  M. Krivanek and J. Moravek (1986).
  NP-hard problems in hierarchical tree clustering.
  \emph{Acta Informatica}, \bold{23}, 311--323.

  G. de Soete (1986).
  A least squares algorithm for fitting an ultrametric tree to a
  dissimilarity matrix.
  \emph{Pattern Recognition Letters}, \bold{2}, 133--137.
}
\seealso{
  \code{\link{cl_consensus}} for computing least squares consensus
  hierarchies by least squares fitting of average ultrametric
  distances.
}
\examples{
## Least squares fit of an ultrametric to the Miller-Nicely consonant
## phoneme confusion data.
data("Phonemes")
## Note that the Phonemes data set has the consonant misclassification
## probabilities, i.e., the similarities between the phonemes.
d <- 1 - as.dist(Phonemes)
u <- ls_fit_ultrametric(d, control = list(verbose = TRUE))
## Cophenetic correlation:
cor(d, u)
## Dendrogram:
plot(u)
## ("Basically" the same as Figure 1 in de Soete (1986).)
}
\keyword{cluster}
\keyword{optimize}
