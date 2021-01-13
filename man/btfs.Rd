\name{btfs}
\alias{btfs}
\title{
Compute a Bootstrapped Transfer Function Stability test BTFS.
}
\description{
The function bootstraps linear regression parameters between a predictor and predictand for two different periods of equal length.
Bootstrapped regression parameter ratios are computed between the two periods and tested for equality (equalling a ratio of 1).
If the true regression parameter ratio is unlikely to be 1 as derived from empirical cumulative density functions, the regression parameter is deemed instable over the two periods.
In addition to BTFS the classic tests reduction of error and coefficient of efficiency are also bootstrapped and displayed. 
However, RE and CE were shown to be less sensitive to instability if compared to BTFS (see Buras et al., 2017).
}
\usage{
btfs(x, y, boot.n = 1000)
}
\arguments{
  \item{x}{
The predictor.
}
  \item{y}{
The predictand.
}
  \item{boot.n}{
Number of iterations for the bootstrapping.
}
}
\details{
Make sure x and y have equal lenght, otherwise an error will occur.
}
\value{
\item{stats}{Output matrix containing the regression parameter ratios as well as the corresponding p-values. If p-values undergo 0.05, the parameter can be considered instable over time
To allow for compatibility with classic tests such as reduction of error and coefficient of efficiency, these metrics are also bootstrapped and displayed.}
\item{boot.mat}{Matrix containing the regression parameter estimates for each iteration which can be visualized if plotting the output of the function.}
}
\references{
Buras, Zang, Menzel, 2017: Testing the stability of transfer functions. Dendrochronologia 42: 56-62. http://dx.doi.org/10.1016/j.dendro.2017.01.005
}
\author{
Allan Buras
}
\note{
Please check out dendroschool.org for video tutorials on how to apply BTFS.
}