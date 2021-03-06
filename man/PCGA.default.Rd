\name{PCGA.default}
\alias{PCGA.default}
\title{
compute a Principal Component Gradient Analysis
}
\description{
pcga makes use of an ordinary PCA to compute the similarity of tree-ring series based on the loadings of the first two principal components. 
x-y-coordinates of PC-loadings are transformed into polar coordinates which are used to define a gradient of similarity along the loadings. 
}
\usage{
PCGA.default(DATA,plot.CI=T)
}
\arguments{
  \item{DATA}{
An rwl.data.frame such as the one generated by read.rwl in the dplR-package.
}
\item{plot.CI}{
Specify whether the common overlap period shall be plotted. Defaults to TRUE.
}
}
\details{
pcga automatically extracts the common overlap period of the rwl.data.frame. If this period is too short, results may be meaningless -->
consider removing too short series. Moreover, detrending has an important impact on the results obtained by pcga. Depending on research
questions you may consider using rigorous detrending techniques to emphasize high-frequency variability (e.g. polynomial spline or autoregressive models) 
or representations of absolute growth such as basal area increments (bai.in and bai.out in 'dplR').
}
\value{
\item{pca}{result of the ordinary PCA computed within pcga}
\item{imp}{importance of the first two principal components expressed as proportions - the higher, the more meaningful the results of pcga}
\item{rank}{ranking according to the loadings of the PCA, i.e. the more similar the rank the more similar the variability of tree-ring series}
\item{pol.coord}{polar coordinates of the loadings}
\item{pop}{input data used for the computation}
\item{period}{common overlap period of the input data}
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



