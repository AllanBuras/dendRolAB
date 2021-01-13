\name{SGC.default}
\alias{SGC.default}
\title{
Computes Standardized Growth Changes of tree-ring data.
}
\description{
Standardized growth changes can be considered a novel tool for pointer year detection, which was shown by Buras et al. (2020) to feature a higher sensitivity to extraordinary growth.
}
\usage{
SGC.default(rwl, p.thresh = 0.05, make.plot = T, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rwl}{
An rwl.data.frame as obtained using the function read.rwl of the dplR package.
}
  \item{p.thresh}{
The significance level for pointer year detection. Defaults to 0.05.
}
  \item{make.plot}{
Logical. Whether to draw a plot which depicts the identified pointer years.
}
  \item{\dots}{
Additional arguments to be passed to the plot.
}
}
\details{
As recommended by Buras et al. (2020) make sure to only use the common overlap period of tree-ring data.
}
\value{
\item{negative.PY}{P-values of identified negative pointer years.}
\item{positive.PY}{P-values of identified positive pointer years.}
}
\references{
Buras, Rammig, and Zang (2020): A novel approach for the detection of pointer years. Dendrochronologia 63, 125746. https://doi.org/10.1016/j.dendro.2020.125746.
}
\author{
Allan Buras
}
\note{
I am currently working on a refined version of SGC which is scheduled for publication in 2021. 
The refined version will also detect pointer years following gradual growth decline and allow for defining the length of the recovery period after an event.
}


