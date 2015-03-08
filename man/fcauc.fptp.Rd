\name{fcauc.fptp}
\alias{fcauc.fptp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generating ROC plot with shading based on FDR using user-provided FPR, TPR, FDR values
}
\description{
The function plots the ROC curve using user-provided FPR, TPR, FDR values. It shades the area corresponding to acceptable FDR level. The FDR-controlled area under the curve (FCAUC) is calculated for the shaded area.
}
\usage{
fcauc.fptp(FP, TP, TDR, FDR.cut = 0.2, do.plot=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FP}{Vector of false positive rate.}
  \item{TP}{Vector of true positive rate.}
  \item{TDR}{Vector of true discovery rate.}
  \item{FDR.cut}{The cutoff value of FDR, for the purpose of FCAUC cualculation.}
  \item{do.plot}{Whether to generate the plot.}
}
\value{The FCAUC value is returned.}
\references{Yu T (2012) ROCS: Receiver Operating Characteristic Surface for Class-Skewed High-Throughput Data. PLoS ONE 7(7): e40598. }
\author{Tianwei Yu. Email: tianwei.yu@emory.edu. }
\seealso{
fcauc.x
}
\keyword{classif}
