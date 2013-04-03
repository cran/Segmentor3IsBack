\name{SelectModel}
\alias{SelectModel}
\title{SelectModel}
\description{
  This function is used to select the number of segments in the segmentation of the data according to several criterion including BIC, AIC, mBIC and oracle penalties. 
}

\usage{
SelectModel(x,penalty="oracle",seuil=n/log(n),keep=FALSE)
}

\arguments{
  \item{x}{An object of class Segmentor returned by function Segmentor}
  \item{penalty}{An object of class string stating which penalty to use in the model selection criteria. Accepted penalties are BIC, AIC, Zhang's modified BIC: mBIC and oracle penalties: 'oracle'. In the case of Normal segmentation, criterion developed by Lebarbier, in the case of Poisson and Negative Binomial segmentation, criterion developed by Cleynen and Lebarbier. Default is oracle}
  \item{seuil}{If penalty='oracle', an Integer for the threshold to use for the slope heuristic. Default value is n/log(n) }
  \item{keep}{A Boolean stating whether or not to keep the values of the criterion. Default is FALSE.}
}

\value{
  \item{K}{The number of segments selected.} 
  \item{crit}{If keep=TRUE, a vector of criterion value for each possible K.}
}



\references{
  PDPA: Rigaill, G. Pruned dynamic programming for optimal multiple change-point detection: \emph{Submitted} \url{http://arxiv.org/abs/1004.0887}

  PDPA: Cleynen, A. and Koskas, M. and Rigaill, G. A Generic Implementation of the Pruned Dynamic Programing Algorithm: \emph{Submitted} \url{http://arxiv.org/abs/1204.5564}

  overdispersion parameter: Johnson, N. and Kemps, A. and Kotz, S. (2005) Univariate Discrete Distributions: \emph{John Wiley & Sons} 
  
  Selection criterion for counts: Cleynen, A. and Lebarbier, E. Segmentation of the Poisson and negative binomial rate models: a penalized estimator: \emph{Submitted} \url{http://arxiv.org/abs/1301.2534}
  
  Selection criterion for Gaussian distribution: Lebarbier, E. (2005) Detecting multiple change-points in the mean of Gaussian process by model selection: \emph{Signal Processing}
  
  Slope heuristic: Arlot, S. and Bach, F. (2009) Data-driven calibration of penalties for least-squares regression: \emph{Journal of Machine Learning Research}
  
  modified BIC: Zhang, N. and Siegmund, D. (2007) A modified Bayes information criterion with applications to the analysis of comparative genomic hybridization data: \emph{Biometrics}
}

\details{
  \tabular{ll}{
    Package: \tab Segmentor3IsBack\cr
    Type: \tab Package\cr
    Version: \tab 1.4\cr
    Date: \tab 2013-03-25\cr
    License: \tab GPL (>= 2)\cr
  }
}

\author{
Alice Cleynen, Michel Koskas and Guillem Rigaill

Maintainer: Who to complain to <alice.cleynen@agroparistech.fr>
}

\keyword{ R-package }
\keyword{ PDPA }
\keyword{ Segmentation }


\examples{
require(Segmentor3IsBack);
N=2000 
x=c(rnbinom(N,size=1.3,prob=0.7),rnbinom(N,size=1.3,prob=0.2),
rnbinom(N,size=1.3,prob=0.01), rnbinom(N,size=1.3,prob=0.2),
rnbinom(N,size=1.3,prob=0.8));
res=Segmentor(data=x,model=3,Kmax=20);  
# Finds the optimal segmentation in up to 20 segments with respect to 
#the negative binomial model.
Cr<-SelectModel(res,penalty='oracle',keep=FALSE)
Cr
#chooses the number of segments in the segmentation of x
}


