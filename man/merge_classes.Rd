\name{merge_classes}
\alias{merge_classes}
\title{
Merge CLIMB classes into parent groups.
}
\description{
Use a theoretical distance between estimated classes to merge them together. This can be useful for more parsimonious visualizations.
}
\usage{
merge_classes(n_groups, chain, burnin)
}
\arguments{
  \item{n_groups}{Integer. How many parent groups to merge the classes into?}
  \item{chain}{An object output from \link[CLIMB]{extract_chains} function.}
  \item{burnin}{Iterations to be used as burn-in.}
}
\value{
    \item{merged_z}{New group labels assigned to each observation}
    \item{merged_mu}{A matrix of merged class means with each column corresponding to a dimension and each row corresponding to a group.}
    \item{merged_prop}{Merged class mixing weights.}
    \item{An object output from \link[stats]{hclust}, showing hierarchical relationship among the non-empty classes.}
}
\author{\email{hbk5086@psu.edu}}