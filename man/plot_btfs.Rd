\name{plot_btfs}
\alias{plot_btfs}
\title{
Boxplot the bootstrapped regression parameter estimates.
}
\description{
This plot visualizes the stability of the regression parameter estimates over the two periods.
}
\usage{
plot_btfs(x,...)
}
\arguments{
  \item{x}{
An object as generated using the btfs function.
}
\item{\dots}{
Additional parameters to change the plot appearance.
}
}
\details{
In case of instable regression parameters, significance stars will be added denoting p < 0.05 (*), p < 0.01 (**), and p < 0.05 (***).
}
\author{
Allan Buras
}
