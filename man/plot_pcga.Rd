\name{plot_pcga}
\alias{plot_pcga}
\title{
Plot the loadings of PCGA, optionally colored according to some categorical sampling design.
}
\description{
Visualizes the most important features from PCGA, i.e. loadings and explained variance of the first two PCs.
}
\usage{
plot_pcga(x,col.vec = NULL,...)
}
\arguments{
  \item{x}{
output from \link{pcga}
}
  \item{col.vec}{
optinal: a character vector indicating the colors to be used for loadings according to some specific sampling design
}
  \item{\dots}{
additional parameters to change the appearance of the plot
}

}
\references{
Buras et al., 2016: Tuning the voices of a choir: detecting ecological gradients in time-series populations. PLoS ONE 11(7): e0158346. DOI:10.1371/journal.pone.0158346
}
\author{
Allan Buras
}
\note{
for further details, please note the dendroschool video-tutorial on how to compute PCGA on dendroschool.org
}
