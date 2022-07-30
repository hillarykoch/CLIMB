\name{calculate_gelmanRubin}
\alias{calculate_gelmanRubin}
\title{Calculate the potential scale reduction factor from parallel MCMCs}
\description{The Gelman-Rubin convergence diagnostic calculates the "potential
    scale reduction factor" (PSRF) for a series of MCMC chains. If the collection
    of parallel MCMCs have converged for a given parameter, the PSRF will be
    close to 1.}
\usage{
calculate_gelmanRubin(chain_list, burnin_list)
}
\arguments{
  \item{chain_list}{A list of MCMC runs output from \link{extract_chains}.}
  \item{burnin_list}{A vector or list of vectors specifying the iterations to discard as burn-in for each chain.}
}
\details{
    NAs are returns for parameters associated with a 0 label, since these are not
    sampled during MCMC.
}
\value{
\item{PSRF_mu}{Calculated PSRF for cluster mean estimates.}
\item{PSRF_Sigma}{Calculated PSRF for cluster covariance estimates.}
\item{PSRF_prop}{Calculated PSRF for cluster proportion estimates.}
}
\author{\email{hbk5086@psu.edu}}
\references{Gelman, A., and D. B. Rubin. 1992. Inference from Iterative Simulation Using Multiple Sequences. Statistical Science 7: 457–511.}