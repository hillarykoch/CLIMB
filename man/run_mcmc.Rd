\name{run_mcmc}
\alias{run_mcmc}
\title{Fit constrained normal mixture model via MCMC.}
\description{This is the fourth and final step of CLIMB. It fits a constrained normal mixture model to the data, given a final list of candidate latent classes and prior hyperparameters.}
\usage{
run_mcmc(dat, hyp, nstep, retained_classes)
}
\arguments{
  \item{dat}{n by D matrix or data frame of appropriately pre-processed observations.}
  \item{hyp}{Hyperparameters output from \link{get_hyperparameters}.}
  \item{nstep}{Integer. Number of MCMC iterations.}
  \item{retained_classes}{Final list of candidate latent classes, after eliminating classes whose prior weights are too small.}
}
\details{
The proposals for each cluster in the MCMC are adaptively tuned such that the acceptance rates converge to about 0.3
}
\value{
\item{chain}{A Julia object. The estimated parameters over the course of nstep iterations}
\item{acceptance_rate_chain}{an M by nstep matrix of the acceptance rates for each cluster covariance.}
\item{tune_df_chain}{The tuning degrees of freedom across the chain, adjusted to yield optimal acceptance rates.}
}
\author{\email{hbk5086@psu.edu}}