\name{get_pairwise_fits}
\alias{get_pairwise_fits}
\title{Fit all D choose 2 pairwise models to D-dimensional data}
\description{This is the first step of the CLIMB procedure. It fits a penalized constrained normal mixture model to D-dimensional data, automatically selecting the number of clusters in each pairwise fit. This can be parallelized across up to D choose 2 cores.}
\usage{
get_pairwise_fits(z, nlambda = 10, parallel = TRUE, ncores = 10, bound = 0)
}
\arguments{
  \item{z}{n by D matrix or data frame of appropriately pre-processed observations.}
  \item{nlambda}{Integer. Number of penalty terms to try}
  \item{parallel}{Boolean. Should the analysis be executed in parallel?}
  \item{ncores}{Integer. If parallel=TRUE, number of cores to parallelize over.}
  \item{bound}{Numeric. Is there a minimum value for the mean of the non-null distribution?}
}
\details{CLIMB automatically determines an optimal grid of lambdas, given a number of lambdas. In all analyses in the original CLIMB paper, bound was always set to 0.}
\value{
CLIMB returns a list of pairwise fits, each with a name in the style of "d1_d2", implying it is the pairwise fit between dimensions d1 and d2. Each fit contains 12 elements.
\item{k}{Number of clusters estimated to be present in the pairwise fit.}
\item{prop}{Cluster mixing weights.}
\item{mu}{Non-null cluster mean.}
\item{sigma}{Non-null cluster standard deviation.}
\item{rho}{Correlation between replicates with association.}
\item{df}{Degrees of freedom for each cluster.}
\item{cluster}{Cluster labels for each of the n observations.}
\item{BIC}{BIC of optimal fit.}
\item{lambda}{Lambda of optimal fit.}
\item{ll}{Log likelihood of the data given the optimal fit.}
\item{post_prob}{Posterior probability of each observation arising from each cluster.}
\item{combos}{Used internally with fconstr_pGMCM. Can be ignored by the user.}
}
\author{\email{hbk5086@psu.edu}}
\examples{
data(sim)
#------------------------------------------------------------
# Not run:
#
# fits <- get_pairwise_fits(z = sim$data, parallel = FALSE)
# plot(sim$data[,1:2], col = fits[[1]]$cluster)
#------------------------------------------------------------
}