\name{Segmentor3IsBack-package}
\alias{Segmentor3IsBack-package}
\alias{Segmentor3IsBack}
\docType{package}
\title{
Implementation of the Pruned Dynamic Programming Algorithm for the exact optimal segmentation of profiles
}
\description{
Exact change-point algorithm for the segmentation of profiles according to the log-likelihood criterion for 5 possible models: Poisson, Gaussian homoscedastic, negative binomial, Gaussian with constant mean and Exponential.
}
\details{
\tabular{ll}{
Package: \tab Segmentor3IsBack\cr
Type: \tab Package\cr
Version: \tab 1.5 \cr
Date: \tab 2013-03-25\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}

}
\author{
Alice Cleynen

Maintainer: Alice Cleynen <alice.cleynen@agroparistech.fr>
}

\references{
  PDPA: Rigaill, G. Pruned dynamic programming for optimal multiple change-point detection: \emph{Submitted} \url{http://arxiv.org/abs/1004.0887}

  PDPA: Cleynen, A. and Koskas, M. and Lebarbier, E. and Rigaill, G. and Robin, S. Segmentor3IsBack (2014): an R package for the fast and exact
segmentation of Seq-data \emph{Algorithms for Molecular Biology}

  overdispersion parameter: Johnson, N. and Kemps, A. and Kotz, S. (2005) Univariate Discrete Distributions: \emph{John Wiley & Sons} 
  
  variance parameter: Hall, P. and Kay, J. and Titterington, D. (1990): Asymptotically optimal difference-based estimation of variance in non-parametric regression \emph{Biometrika}
  
  Selection criterion for counts: Cleynen, A. and Lebarbier, E. (2014) Segmentation of the Poisson and negative binomial rate models: a penalized estimator: \emph{Esaim : Probability and Statistics} 
  
  Selection criterion for Gaussian distribution: Lebarbier, E. (2005) Detecting multiple change-points in the mean of Gaussian process by model selection: \emph{Signal Processing}
  
  Slope heuristic: Arlot, S. and Bach, F. (2009) Data-driven calibration of penalties for least-squares regression: \emph{Journal of Machine Learning Research}
  
  modified BIC: Zhang, N. and Siegmund, D. (2007) A modified Bayes information criterion with applications to the analysis of comparative genomic hybridization data: \emph{Biometrics}
}

\keyword{changepoint}
\keyword{segmentation}
\keyword{exact}
\keyword{counts}
\keyword{selection}
\examples{
N=1000 
x=rnbinom(5*N, size=1.3, prob=rep(c(0.7,0.2,0.01,0.2,0.8),each=N))
res=Segmentor(data=x, model=3, Kmax=20, keep=TRUE);  
# Finds the optimal segmentation in up to 20 segments with respect to 
#the negative binomial model, keeping cost matrix.
Cr<-SelectModel(res, penalty='oracle', keep=FALSE)
Cr
#chooses the number of segments in the segmentation of x, not keeping
# values of constants for slope heuristic.
Best<-BestSegmentation(res, K=Cr, t=2*N)
matplot(Best$bestCost, type='l', lty=2)
points(apply(Best$bestCost,2,which.min),apply(Best$bestCost,2,min),pch=20,col=1:(Cr-1))
apply(Best$bestCost, 2, which.min)
getBreaks(res)[Cr,1:(Cr-1)]
#computes and plots cost of best segmentation in Cr segments with 
#change-point t, and compares result with change-point estimates.
Best$bestSeg
#returns the optimal segmentation in Cr segments with t as a
#change-point
}
