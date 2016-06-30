\name{SelectModel}
\alias{SelectModel}
\title{SelectModel}
\description{
  This function is used to select the number of segments in the segmentation of the data according to several criterion including BIC, AIC, mBIC and oracle penalties. 
}

\usage{
SelectModel(x,penalty="oracle",seuil=n/log(n),keep=FALSE,greatjump=FALSE)
}

\arguments{
  \item{x}{An object of class Segmentor returned by function Segmentor}
  \item{penalty}{An object of class string stating which penalty to use in the model selection criteria. Accepted penalties are BIC, AIC, Zhang's modified BIC: mBIC and oracle penalties: 'oracle'. In the case of Normal segmentation, criterion developed by Lebarbier, in the case of Poisson and Negative Binomial segmentation, criterion developed by Cleynen and Lebarbier. Default is oracle}
  \item{seuil}{If penalty='oracle', an Integer for the threshold to use for the slope heuristic. Default value is n/log(n) }
  \item{keep}{A Boolean stating whether or not to keep the values of the criterion. Default is FALSE.}
  \item{greatjump}{If penalty is "oracle", a boolean stating whether to use the greatest jump (TRUE) or the threshold for the slope heuristic. Default is FALSE.}
}

\value{
  \item{K}{The number of segments selected.} 
  \item{crit}{If keep=TRUE, a vector of criterion value for each possible K.}
}



\references{
  PDPA: Rigaill, G. Pruned dynamic programming for optimal multiple change-point detection: \emph{Submitted} \url{http://arxiv.org/abs/1004.0887}

  PDPA: Cleynen, A. and Koskas, M. and Lebarbier, E. and Rigaill, G. and Robin, S. Segmentor3IsBack (2014): an R package for the fast and exact
segmentation of Seq-data \emph{Algorithms for Molecular Biology}

  overdispersion parameter: Johnson, N. and Kemps, A. and Kotz, S. (2005) Univariate Discrete Distributions: \emph{John Wiley & Sons} 
  
  Selection criterion for counts: Cleynen, A. and Lebarbier, E. (2014) Segmentation of the Poisson and negative binomial rate models: a penalized estimator: \emph{ESAIM: Probability and Statistics} 
  
  Selection criterion for Gaussian distribution: Lebarbier, E. (2005) Detecting multiple change-points in the mean of Gaussian process by model selection: \emph{Signal Processing}
  
  Slope heuristic: Arlot, S. and Bach, F. (2009) Data-driven calibration of penalties for least-squares regression: \emph{Journal of Machine Learning Research}
  
  modified BIC: Zhang, N. and Siegmund, D. (2007) A modified Bayes information criterion with applications to the analysis of comparative genomic hybridization data: \emph{Biometrics}
}

\details{
  \tabular{ll}{
    Package: \tab Segmentor3IsBack\cr
    Type: \tab Package\cr
    Version: \tab 1.5\cr
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
x=rnbinom(5*N, size=1.3, prob=rep(c(0.7,0.2,0.01,0.2,0.8),each=N))
res=Segmentor(data=x,model=3,Kmax=20);  
# Finds the optimal segmentation in up to 20 segments with respect to 
#the negative binomial model.
Cr<-SelectModel(res,penalty='oracle',keep=FALSE)
Cr
#chooses the number of segments in the segmentation of x using
# an oracle-inequality approach

N=250 
x=rpois(10*N, rep(c(8,1,5,3,16,33,2,12,7,1),each=N))
res=Segmentor(data=x,model=3,Kmax=40);  
# Finds the optimal segmentation in up to 40 segments with respect to 
#the poisson model.
Cr<-SelectModel(res,penalty='BIC',keep=FALSE)
Cr
#chooses the number of segments in the segmentation of x using
# the BIC approach
}


