\name{fconstr_pGMM}
\alias{fconstr_pGMM}
\title{
Fit constrained penalized Gaussian mixture model
}
\description{
For a given maximum number of clusters in the data, find the optimal penalization parameter lambda and the optimal clustering. Lambda penalizes the mixing proportions. Model optimality is determined by best BIC. The constraints here are that, for any number of replicates, there are 3 components h = {-1,0,1}. h=0 implies no associations, while 1 and -1 imply positive and negative associations. The mean of the no association component is restricted to 0, while the positive and negative associations have mean mu and -mu. The standard deviation of the no association component is restricted to 1, while the positive and negative associations have standard deviation sigma. For two replicates where both have positive association or both have negative association, they have correlation rho. For two replicates, one with positive association and one with negative association, they have correlation -rho. Otherwise, correlation between replicates is restricted to 0.
}
\usage{
fconstr_pGMM(x, lambda = NULL, tol=1e-06,
                itermax = 300, penaltyType = c("SCAD", "LASSO"))
}
\arguments{
  \item{x}{An n by d numeric matrix of n observations with dimension d.}
  \item{lambda}{Penalty parameter. If unspecified, a grid is automatically generated.}
  \item{tol}{Tolerance for the stopping rule for the EM algorithm. A lower tolerance will require more iterations. Defaults to 1e-06.}
  \item{itermax}{Maximum number of iterations of the EM algorithm to perform. Defaults to 300.}
  \item{penaltyType}{Character string specifying which term to use to penalize the cluster mixing proportions. Defaults to SCAD.}
}
\details{The model is fit using the expectation-maximization algorithm. Clusters are initialized using kmeans with 3^d clusters. Initial values of mixing proportions, means, and variance (covariance) matrices for EM are computed from these clusters.}
\value{
    \item{k}{Optimal number of clusters.}
    \item{prop}{Mixing proportions in each cluster.}
    \item{mu}{Cluster means.}
    \item{sigma}{Cluster variances.}
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
\seealso{\link{fconstr_pGMCM}}
\examples{
library(mvtnorm)
set.seed(234)
pal <- sample(get_pals(4), 9, replace = FALSE)

n <- 3600
prop <- c(0.4,0.05,0.1,0.05,0.1,0.05,0.1,0.05,0.1)
mu <- 4
sigma <- 1.3
rho <- 0.8
sim <- rconstr_GMCM(n, prop, mu, sigma, rho, 2)

################################################################################
# Not run:
# fit <- fconstr_pGMM(sim$data, itermax = 100)
#
# par(mfrow = c(1,2))
# plot(sim$data, col = pal[sim$cluster], main = "observed")
# plot(sim$data, col = pal[fit$cluster], main = "GMM classification")
################################################################################
}
