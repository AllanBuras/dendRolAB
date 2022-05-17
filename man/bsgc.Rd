\name{bsgc}
\alias{bsgc}
\title{
Compute the bias-adjusted standardized growth change method for pointer year analysis.
}
\description{
bsgc represents a refined version of the sgc method (Buras et al., 2020). It now also allows for detecting pointer years following successive growth decline as well as an estimate for the duration of the growth relaxation period, i.e. the period until when growth values have reached expected rates.
}
\usage{
bsgc(rwl, p.thresh=0.05, make.plot=T, maxlag=20, rm.succ=TRUE, w1=0.93, w2=0.53,...)
}
\arguments{
  \item{rwl}{
An rwl data frame as created by the read.rwl function of dplR
}
  \item{p.thresh}{
The significance level. Defaults to 0.05 and must not be increased to conform with international statistical regulations.
}
  \item{make.plot}{
Logical indicating whether or not to draw a plot that indicates identified pointer years and related growth relaxation periods. Defaults to TRUE.
}
  \item{maxlag}{
Numeric indicating the maximum lag years to take into consideration when computing the duration of the growth relaxation period. Defaults to 20 and can at most be n/2 where n corresponds with the number of years in the rwl.data.frame.
}
  \item{rm.succ}{
Logical indicating whether or not to correct for possible biases in succeeding pointer years, i.e. pointer years identified right after a pointer year of opposite sign (for details see Buras et al., 2022).
}
  \item{w1}{
Weighing factor applied to global standardized growth changes (GSCG) in the conflation of p-values (for details see section 2.3 in Buras et al., 2022). Changing w1 is strongly discouraged.
}
  \item{w2}{
Weighing factor applied to standardized growth changes (SCG) in the conflation of p-values (for details see section 2.3 in Buras et al., 2022). Changing w2 is strongly discouraged.
}
  \item{\dots}{
Additional arguments passed to change the appearance of the plot in case make.plot equals TRUE.
}
}
\details{
P.thresh must not be increased in order to conform with statistical rules and keep type I error probability low. w1 and w2 should not be changed since this changes the sensitivity of the method.
}
\value{
\item{neg}{A vector containing the p-values of identified negative pointer years. The vector names refer to corresponding years.}
\item{pos}{A vector containing the p-values of identified positive pointer years. The vector names refer to corresponding years.}
\item{pvals.sgc}{A matrix containing the p-values related to SGC. Rows refer to years, columns to lag-time in years.}
\item{pvals.gsgc}{A matrix containing the p-values related to GSGC. Rows refer to years, columns to lag-time in years.}
\item{pvals.conf}{A matrix containing conflated p-values combining pvals.sgc and pvals.gsgc which is used to identify pointer years and relaxation periods. Rows refer to years, columns to lag-time in years.}
\item{def.neg}{A vector containing the duration of deflection periods for all identified negative pointer years. The vector names refer to corresponding years.}
\item{def.pos}{A vector containing the duration of deflection periods for all identified positive pointer years. The vector names refer to corresponding years.}
}
\references{
Buras, Rammig, and Zang (2020): A novel approach for the detection of pointer years. Dendrochronologia 63, 125746. https://doi.org/10.1016/j.dendro.2020.125746.
Buras, Ovenden, Rammig, and Zang (2022): Refining the standardized growth change method: bias-adjustment and definition of the relaxation period. Dendrochronologia 65, 125964. https://doi.org/10.1016/j.dendro.2022.125964
}
\author{
Allan Buras
}
\note{
bsgc replaces the former function sgc, which is outdated.
A video tutorial exemplifying the application of bsgc is provided at https://dendroschool.org/bias-adjusted-standardized-growth-change-method-bsgc/
}

