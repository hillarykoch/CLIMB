\name{get_pals}
\alias{get_pals}
\title{Stored color palettes.}
\description{A collection of pleasant color palettes, good for color contrast or color gradients.}
\usage{
get_pals(pal)
}
\arguments{
  \item{pal}{An integer corresponding to a color palette.}
}
\value{Hex codes of a color palette.}
\author{\email{hbk5086@psu.edu}}
\examples{
pal <- get_pals(4)
n <- length(pal)
plot(seq(n), runif(n), col = pal, cex = 15)
}
