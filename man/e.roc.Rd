\name{e.roc}
\alias{e.roc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculating expected Receiver Operating Characteristics Curve in the presence of imperfect reference standards.
}
\description{
This function builds an eROC curve and returns the expected values of TPR, FPR, and TDR. It also calculates the AUC of the eROC curve. The 95\% bootstrap percentile confidence interval for the AUC is provided.}
\usage{
e.roc(x, mu, method='RNA', bt.ci=TRUE, bt.nreps=100, do.plot=TRUE)
}
\arguments{
  \item{x}{Vector; the scores yielded by the classifier.}
  \item{mu}{Vector; the probabilistic confidence assigned by the imperfect reference standard.}
  \item{method}{The method used to compute the cumulative distribution function for the Poisson binomial distribution. ``DFT-CF'' for the DFT-CF method, ``RF'' for the recursive formula, ``RNA'' for the refined normal approximation, ``NA'' for the normal approximation. }
  \item{bt.ci}{Whether to compute the bootstrap confidence interval.}
  \item{bt.nreps}{The number of bootstrap replicates.}
  \item{do.plot}{Whether to plot the eROC curve.}
}
\details{
The eROC curve is a generalization of ROC curve given the class membership uncertainties. See the reference for the definition of the eROC curve.}
\value{
Returns the area under the eROC curve, the expected values of TPR, FPR, and TDR. 
}
\references{Liao P, Wu H, and Yu T (2016). ROC Curve Analysis in the Presence of Imperfect
Reference Standards. Stat Biosci doi:10.1007/s12561-016-9159-7.}
\author{Peizhou Liao. Email: pliao3@emory.edu. }


\examples{
## normal scores
sample.p <- rnorm(100, mean=1, sd=sqrt(2))  
sample.n <- rnorm(100, mean=-1, sd=sqrt(2))  

## probabilistic confidence
mu.p <- rbeta(100, shape1=5, shape2=1)
mu.n <- rbeta(100, shape1=1, shape2=5)

## combine the sample
x.sample <- c(sample.p, sample.n)
mu.sample <- c(mu.p, mu.n)

## build eROC curve
e.roc.fit <- e.roc(x=x.sample, mu=mu.sample)

}
\keyword{classif}

