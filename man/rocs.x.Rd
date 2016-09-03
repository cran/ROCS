\name{rocs.x}
\alias{rocs.x}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculating Receiver Operating Characteristics Surface from one-dimensional data.
}
\description{
The function plots the colored ROCS in 3 dimensions using the rgl utilities. It shades the FDR-controlled AUC in the FPR-TPR plain. The VUS of the surface and the FCAUC area are calculated. The significance of the VUX is assessed by one-sided permutation test.
}
\usage{
rocs.x(x0, x1, s0=NULL, s1=NULL, n.perm = 1000, do.plot = TRUE, FDR.cut=0.2)
}
\arguments{
  \item{x0}{Vector; the raw data of the null class.}
  \item{x1}{Vector; the raw data of the non-null class.}
  \item{s0}{Vector; the confidence level of the class assignment of the null class observations. The length should be the same as x0. The default is NULL, in which case all observations are considered to be assigned without uncertainty. }
  \item{s1}{Vector; the confidence level of the class assignment of the non-null class observations. The length should be the same as x0. The default is NULL, in which case all observations are considered to be assigned without uncertainty. }
  \item{n.perm}{The number of permutations to assess the significance of the VUX.}
  \item{do.plot}{Whether to plot the 3D surface, or just return the VUX.}
  \item{FDR.cut}{The FDR level at which to shade the AUC and calculate the FCAUC.}
}
\details{
The ROCS is the surface spanned by the TPR-FPR-TDR and its projection to the TPR-TDR plain. The volume is between the surface and its projection on the TPR-FDR plain (i.e. the AUC of ROC).
}
\value{
The volume under the surface (VUS) is returned.
}
\references{Yu T (2012) ROCS: Receiver Operating Characteristic Surface for Class-Skewed High-Throughput Data. PLoS ONE 7(7): e40598. 
Yu T, Jones DP.  (2014) Improving peak detection in high-resolution LC/MS metabolomics data using preexisting knowledge and machine learning approach. Bioinformatics. 30(20):  2941-2948.}
\author{Tianwei Yu. Email: tianwei.yu@emory.edu. }
\seealso{
rocs.fptp
}
\examples{

x0<-rnorm(500, mean=0, sd=1.5)
x1<-rnorm(50, mean=4, sd=1)

rocs.x(x0, x1)$vus

}
\keyword{classif}
