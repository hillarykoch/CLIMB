\name{get_prior_weights}
\alias{get_prior_weights}
\title{Compute priors on the mixing weights of each candidate latent class.}
\description{This is the third step of the CLIMB procedure. It is useful because one can eliminate latent classes whose prior weights are too small. This function is parallelizable.}
\usage{
get_prior_weights(reduced_classes, fits, parallel = FALSE, ncores = 20, delta = NULL)
}
\arguments{
  \item{reduced_classes}{Matrix output from \link{get_reduced_classes}.}
  \item{fits}{Fits output from \link{get_pairwise_fits}.}
  \item{parallel}{Boolean. Should this function be parallelized?}
  \item{ncores}{If parallel=TRUE, how many cores should this be parallelized over?}
  \item{delta}{Which value of integer tuning parameters delta should be checked? Defaults to 0:choose(D,2)}
}
\value{
List of vectors. Each element of the list is a vector corresponds to the computed prior weights for a given value delta.
}
\author{\email{hbk5086@psu.edu}}
\examples{
data(sim)
data(fits)

#------------------------------------------------------------
# Not run:
#
# D <- ncol(sim$data)  
# reduced_classes <- get_reduced_classes(fits, D, "lgf.txt")
# prior_weights <- get_prior_weights(reduced_classes, fits, parallel = FALSE)
# plot(
#    seq_along(prior_weights)-1,
#    sapply(prior_weights, function(X)
#    sum(X * n > D)),
#    ylab = "number of retained classes",
#    xlab =  expression(delta))
#------------------------------------------------------------
}