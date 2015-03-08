\name{fcauc.x}
\alias{fcauc.x}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generating ROC plot with shading based on FDR
}
\description{
The function plots the ROC curve. It shades the area corresponding to acceptable FDR level. The FDR-controlled area under the curve (FCAUC) is calculated for the shaded area.
}
\usage{
fcauc.x(x0, x1, FDR.cut = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x0}{Vector; the raw data of the null class.}
  \item{x1}{Vector; the raw data of the non-null class.}
  \item{FDR.cut}{The FDR level at which to shade the AUC and calculate the FCAUC.}
}
\value{
The FCAUC value is returned.
}
\references{Yu T (2012) ROCS: Receiver Operating Characteristic Surface for Class-Skewed High-Throughput Data. PLoS ONE 7(7): e40598. }
\author{Tianwei Yu. Email: tianwei.yu@emory.edu. }
\seealso{
fcauc.fptp
}
\examples{
#perfect separation
x0<-runif(1000)
x1<-runif(100)+2

fcauc.x(x0, x1)

#partial separation
x0<-rnorm(1000, mean=0, sd=1.5)
x1<-rnorm(100, mean=3, sd=1)

fcauc.x(x0, x1)
}
\keyword{classif}
