\name{predict}
\alias{cl_predict}
\title{Predict Memberships}
\description{
  Predict class memberships from R objects representing partitions.
}
\usage{
cl_predict(object, newdata = NULL, ...)
}
\arguments{
  \item{object}{an R object representing a partition of objects.}
  \item{newdata}{an optional data set giving the objects to make
    predictions for.  This must be of the same \dQuote{kind} as the data
    set employed for obtaining the partition.  If omitted, the original
    data are used.}
  \item{\dots}{arguments to be passed to and from methods.}
}
\value{
  An object of class \code{"cl_membership"} with the matrix of predicted
  membership values.
}
\details{
  Many algorithms resulting in partitions of a given set of objects can
  be taken to induce a partition of the underlying feature space for the
  measurements on the objects, so that class memberships for
  \dQuote{new} objects can be obtained from the induced partition.
  Examples include partitions based on assigning objects to their
  \dQuote{closest} prototypes, or providing mixture models for the
  distribution of objects in feature space.

  This is a generic function.  The methods provided in package
  \pkg{clue} handle the partitions obtained from clustering functions in
  the base R distribution, as well as packages \pkg{RWeka}, \pkg{cba},
  \pkg{cclust}, \pkg{cluster}, \pkg{e1071}, \pkg{flexclust},
  \pkg{flexmix}, \pkg{kernlab}, and \pkg{mclust} (and of course,
  \pkg{clue} itself).
}
\examples{
## Run kmeans on a random subset of the Cassini data, and predict the
## memberships for the "test" data set.
data("Cassini")
nr <- NROW(Cassini$x)
ind <- sample(nr, 0.9 * nr, replace = FALSE)
party <- kmeans(Cassini$x[ind, ], 3)
table(cl_class_ids(cl_predict(party, Cassini$x[-ind, ])),
      Cassini$classes[-ind])
}
\keyword{cluster}
