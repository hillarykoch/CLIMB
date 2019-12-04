\name{fpGMM}
\alias{fpGMM}
\title{
Fit penalized Gaussian mixture model
}
\description{
For a given maximum number of clusters in the data, find the optimal penalization parameter lambda and the optimal clustering. Lambda penalized the mixing proportions. Model optimality is determined by best BIC.
}
\usage{
fpGMM(x, kmax, lambda = NULL, tol=1e-06, itermax = 200)
}
\arguments{
  \item{x}{An n by d numeric matrix of n observations with dimension d.}
  \item{kmax}{Maximum number of clusters to consider, and the number of clusters the algorithm is initialized with.}
  \item{lambda}{Penalty parameter. If unspecified, a grid is automatically generated.}
  \item{tol}{Tolerance for the stopping rule for the EM algorithm. A lower tolerance will require more iterations. Defaults to 1e-06.}
  \item{itermax}{Maximum number of iterations of the EM algorithm to perform. Defaults to 200.}

}
\details{The model is fit using the expectation-maximization algorithm. Clusters are initialized using kmeans with kmax clusters. Initial values of mixing proportions, means, and variance (covariance) matrices for EM are computed from these clusters.}
\value{
\item{k}{Optimal number of clusters.}
\item{prop}{Mixing proportions in each cluster.}
\item{mu}{Cluster means.}
\item{sigma}{Cluster variance (covariance) matrices}
\item{cluster}{Cluster labels for each of the n observations.}
\item{BIC}{BIC of optimal fit.}
\item{lambda}{Lambda of optimal fit.}
\item{ll}{Log likelihood of the data given the optimal fit.}
}
\author{\email{hbk5086@psu.edu}}
\examples{
library(mvtnorm)
set.seed(123)
pal <- get_pals(1)
sim <- rbind(rmvnorm(500, mean = c(-1,0),
                sigma = matrix(c(.65,0,0,2), nrow = 2, byrow = T)),
         rmvnorm(500, mean = c(1,-3),
                sigma = matrix(c(2.25,.1,.1,.65), nrow = 2, byrow = T)),
         rmvnorm(100, mean = c(4,2),
                sigma = matrix(c(1,.5,.5,1), nrow = 2, byrow = T)))

################################################################################
# Not run:
# fit <- fpGMM(sim, lambda = 1, kmax = 10)
#
# par(mfrow = c(1,2))
# plot(sim, col = pal[rep(1:3,times = c(500, 500, 100))], pch = 16, cex = 0.75,
#     xlab = expression(x[1]), ylab = expression(x[2]), main = "true clusters")
# plot(sim, col = pal[fit$cluster], pch = 16, cex = 0.75,
#     xlab = expression(x[1]), ylab = expression(x[2]), main = "fitted clusters")
################################################################################
}
