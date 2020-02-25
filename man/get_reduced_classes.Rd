\name{get_reduced_classes}
\alias{get_reduced_classes}
\title{Find D-dimensional candidate latent classes from pairwise latent classes.}
\description{This is the second step of the CLIMB procedure. For higher dimensions, this function is fast but requires a lot of memory.}
\usage{
get_reduced_classes(fits, d, filepath = "lgf.txt")
}
\arguments{
  \item{fits}{Output from \link{get_pairwise_fits}}
  \item{d}{Dimension of data}
  \item{filepath}{Path to write the LEMON graph format fit, to be used internally by CLIMB.}
}
\value{
M by D matrix of candidate latent classes. Each row is one candidate class.
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
#------------------------------------------------------------
}