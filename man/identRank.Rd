\name{identRank}
\alias{identRank}
\title{
Identify the rank of tree-ring series from the loadings of the PCGA plot.
}
\description{
This function allows for an interactive identification of single tree-ring series from the PCGA plot. 
Application of this function to a pcga-object will automatically draw the PCGA plot. The rank of the desired tree-ring series is displayed and returned
by clicking on the corresponding loading in the PCGA plot.
}
\usage{
identRank(pcga.obj)
}
\arguments{
  \item{pcga.obj}{
An object as returned from the function pcga.
}
}
\value{
Returns the rank of the identified tree-ring series within the Principal Component Gradient.
}
\references{
Buras et al., 2016: Tuning the voices of a choir: detecting ecological gradients in time-series populations. PLoS ONE 11(7): e0158346. DOI:10.1371/journal.pone.0158346
}
\author{
Allan Buras
}
\note{
Please note that the virtual learning platform at dendroschool.org contains a video-tutorial on how to compute and evaluate pcga.
}


