\name{GMM_kmeans}
\alias{GMM_kmeans}
\title{
Initialize parameters for modified GMM using kmeans clustering.
}
\description{
For a given number of clusters, find initial clusters and compute initial parameter estimates.
}
\usage{
GMM_kmeans(y, k, iter.max = 30)
}
\arguments{
  \item{y}{An n by d numeric matrix with n observations of dimension d.}
  \item{k}{Number of cluster centers, passed to \link[stats]{kmeans}.}
  \item{iter.max}{The maximum number of iterations allowed, passed to \link[stats]{kmeans}.}
}
\value{
    \item{prop}{Proportion of observations assigned to each cluster.}
    \item{mu}{Vector of cluster means.}
    \item{sigma}{List of cluster variance covariance matrices.}
    \item{cluster}{Vecor of cluster labels.}
}
\author{\email{hbk5086@psu.edu}}
\examples{
sim <- rGMCM(1500, prop = c(1/3,1/3,1/3),
            mu = list(c(-1,1), c(1,1), c(0,-sqrt(2))),
            sigma = list(rbind(c(0.65, 0.7794), c(0.7794, 1.55)),
                        rbind(c(0.65, -0.7794), c(-0.7794, 1.55)),
                        diag(c(2, 0.2))))
init <- GMM_kmeans(sim$data, k = 3)
}
