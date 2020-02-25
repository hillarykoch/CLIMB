\name{get_hyperparameters}
\alias{get_hyperparameters}
\title{Get prior hyperparameters for MCMC.}
\description{This is the fourth and final step of CLIMB. For a list of candidate latent classes, this function computes all necessary prior hyperparameters for the Bayesian analysis of the constrained normal mixture model.}
\usage{
get_hyperparameters(z, fits, reduced_classes, prior_weights)
}
\arguments{
  \item{z}{n by D matrix or data frame of appropriately pre-processed observations.}
  \item{fits}{Output from \link{get_pairwise_fits}}
  \item{reduced_classes}{Final list of candidate latent classes, after eliminating classes whose prior weights are too small.}
  \item{prior_weights}{Prior weights for a given delta, as output from \link{get_prior_weights}}
}
\value{
    \item{Psi0}{D by D by M array of prior cluster covariance matrices.}
    \item{mu0}{M by D matrix of prior cluster means.}
    \item{alpha}{M-vector of prior cluster mixing weights.}
    \item{kappa0}{alpha * n -- this controls the informativity of the prior, and can be changed by the user.}
}
\author{\email{hbk5086@psu.edu}}
\examples{
data(sim)
data(fits)
data(true_association)
p <- rep(1/nrow(true_assoc), times = nrow(true_assoc))
hyp <- get_hyperparameters(sim$data, fits, true_assoc, p)
}