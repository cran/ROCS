\name{fcauc.comp.fptp}
\alias{fcauc.comp.fptp}
\title{
Plotting two ROC curves, highlighting the segments of the curves corresponding to acceptable FDR levels.
}
\description{
The function draws two colored ROC curves in one plot. For each curve, the segment that corresponds to acceptable FDR (e.g. <= 0.2) is colored differently than the rest of the curve. 
}
\usage{
fcauc.comp.fptp(roc.1, roc.2, FDR.cut = 0.2, lwd = 3, 
colors = c("blue", "green", "cyan", "red"), 
exp.labels = c("experiment 1", "experiment 2"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{roc.1}{ A list object with three components: TP, FP and TDR, each being a vector. }
  \item{roc.2}{ A list object with three components: TP, FP and TDR, each being a vector. }
  \item{FDR.cut}{The cutoff value of FDR.}
  \item{lwd}{The line width parameter to be passed on to the function lines().}
  \item{colors}{A vector of at least four colors. It is used for the coloring of the two ROC curves.}
  \item{exp.labels}{Experiment labels in the plot.}
}
\value{No value is returned.}
\references{Yu T (2012) ROCS: Receiver Operating Characteristic Surface for Class-Skewed High-Throughput Data. PLoS ONE 7(7): e40598. }
\author{Tianwei Yu. Email: tianwei.yu@emory.edu. }
\seealso{fcauc.fptp()}
\keyword{classif}
