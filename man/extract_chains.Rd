\name{extract_chains}
\alias{extract_chains}
\title{Extract results from MCMC output.}
\description{
\link{run_mcmc} returns a Julia Object. This function extracts the results to be usable within R.}
\usage{
extract_chains(mcmc)
}
\arguments{
  \item{mcmc}{Object output from \link{run_mcmc}.}
}
\value{

}
\references{}
\item{mu_chains}{List with M elements, each element a matrix of dimension iterations x D. mu_chains[[i]] is the MCMC samples for the mean vector of cluster i.}
 \item{Sigma_chains}{List with M elements, each element an array of dimension D x D x iterations. Sigma_chains[[i]] is the MCMC samples for the covariance matrix of cluster i.}
\item{prop_chain}{Matrix of dimension M x iterations, containing the MCMC samples for the mixing proportions of each class.}
\item{z_chain}{Matrix of dimension n x iterations, containing the MCMC samples for the class labels of each observation. These labels correspond to the row indices of the matrix of candidate latent classes.}
\author{\email{hbk5086@psu.edu}}