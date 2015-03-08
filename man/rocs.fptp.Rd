\name{rocs.fptp}
\alias{rocs.fptp}
\title{
Plotting Receiver Operating Characteristics Surface with user-provided FPR, TPR & TDR
}
\description{
The function plots the colored ROCS in 3 dimensions using the rgl utilities. It shades the FDR-controlled AUC in the FPR-TPR plain. The VUS of the surface and the FCAUC area are calculated.
}
\usage{
rocs.fptp(FP, TP, TDR, FDR.cut = 0.2)
}
\arguments{
  \item{FP}{Vector of false positive rate.}
  \item{TP}{Vector of true positive rate.}
  \item{TDR}{Vector of true discovery rate.}
  \item{FDR.cut}{The cutoff value of FDR, for the purpose of FCAUC cualculation.}
}
\details{
The ROCS is the surface spanned by the TPR-FPR-TDR and its projection to the TPR-TDR plain. The vulume is between the surface and its projection on the TPR-FDR plain (i.e. the AUC of ROC).
}
\value{
The volume under the surface (VUS) is returned.
}
\references{Yu T (2012) ROCS: Receiver Operating Characteristic Surface for Class-Skewed High-Throughput Data. PLoS ONE 7(7): e40598. }
\author{Tianwei Yu. Email: tianwei.yu@emory.edu. }
\seealso{
rocs.x
}
\keyword{classif}
