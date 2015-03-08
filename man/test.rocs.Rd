\name{test.rocs}
\alias{test.rocs}
\title{
Testing for significant difference between two ROCS
}
\description{
The null hypothesis being tested is that the two ROCS are generated from populations with equal class-separating power. This function uses a bootstrap-based test.
}
\usage{
test.rocs(x0, x1, z0, z1, B = 1000, do.plot = TRUE)
}
\arguments{
  \item{x0}{Vector; the raw data of the true-negative class in study 1.}
  \item{x1}{Vector; the raw data of the true-positive class in study 1.}
  \item{z0}{Vector; the raw data of the true-negative class in study 2.}
  \item{z1}{Vector; the raw data of the true-positive class in study 2.}
  \item{B}{The number of bootstrap samples to be used in order to estimate the spread of the distribution under the null hypothesis.}
  \item{do.plot}{Whether to plot the spread based on bootstrap samples.}
}
\value{The p-value of the test is returned.}
\references{Yu T (2012) ROCS: Receiver Operating Characteristic Surface for Class-Skewed High-Throughput Data. PLoS ONE 7(7): e40598. }
\author{Tianwei Yu. Email: tianwei.yu@emory.edu. }
\examples{
x0<-rnorm(1000, mean=0, sd=1.5)
x1<-rnorm(100, mean=3, sd=1)

z0<-rnorm(800, mean=0, sd=1.5)
z1<-rnorm(80, mean=3.25, sd=1)

test.rocs(x0, x1, z0,z1, B=1000)
}
\keyword{classif}
