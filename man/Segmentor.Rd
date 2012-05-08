\name{Segmentor}
\alias{Segmentor}
\alias{print.Segmentor}
\alias{Segmentor.default}
\title{Segmentor}
\description{
  The functions are used for change-point problems. Given a loss function (Poisson, Normal homoscedastic, Negative Binomial or Normal Heteroscedastic (with given constant mean)), the function Segmentor estimates the optimal segmentation with respect to the log-likelihood. The Segmentor gives estimates of the breakpoint locations as well as the loss function parameter of each segment.
}

\usage{
Segmentor(data=numeric(), model=1, Kmax = 15, theta = numeric(), m = numeric())
}

\arguments{
  \item{data}{A vector of observations to be segmented. Must have no missing values.}
  \item{model}{Integer between 1 and 4 giving the modelisation of the observed data, 1: poisson (default), 2: normal-homoscedastic, 3: negative binomial or 4: normal-heteroscedastic}
  \item{Kmax}{The maximum number of segments wanted for the data. The Segmentor will find all optimal segmentations in 1 to Kmax segments.}
  \item{theta}{Needed only for the Negative Binomial distribution: the value of the inverse of the overdispersion parameter. If the user does not enter a value, the package uses a modified version of Johnson and Kotz's estimator where the mean is replaced by the median. }
  \item{m}{Needed only for Normal Heteroscedastic distribution: the value of the constant mean. If not entered, the function uses the empirical mean of the data.}
}

\value{
  \item{model}{Emission distribution (Poisson, Normal Homoscedastic, Negative Binomial or Normal Heteroscedastic)} 
  \item{breaks}{Matrix of size Kmax*Kmax of estimated change-point locations for each optimal segmentation in 1 to Kmax segments.}
  \item{parameters}{Matrix of size Kmax*Kmax which elements are the estimated parameters for each segment of the optimal segmentation. If model is Poisson or Normal, the parameter corresponds to the mean of the signal in each segment. If model is Negative binomial, the parameter corresponds to the success-probability of the signal in each segment. If model is normal heteroscedastic, the parameter is the variance assuming known mean.} 
  \item{likelihood}{Vector of size Kmax of resulting negative log-likelihood for each optimal segmentation.}
  \item{overdispersion}{only if model = Negative Binomial, the value of the inverse of overdispersion used for the segmentation}
  \item{mean}{only if model = Normal Heteroscedastic, the value of the mean used for the segmentation}
}



\references{
  Pruned dynamic programming for optimal multiple change-point detection: \url{http://arxiv.org/abs/1004.0887}

  A Generic Implementation of the Pruned Dynamic Programing Algorithm: \url{http://arxiv.org/abs/1204.5564}

  Johnson, Kotz & Kemp: Univariate Discrete Distributions
}

\details{
  \tabular{ll}{
    Package: \tab Segmentor3IsBack\cr
    Type: \tab Package\cr
    Version: \tab 1.1\cr
    Date: \tab 2012-05-10\cr
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
x=c(rpois(N,2.0),rpois(2*N,2.2),rpois(N,1.9));
res=Segmentor(data=x,Kmax=3);  
# Finds the optimal segmentation in 1, 2 and 3 segments with respect to the Poisson model.

y=c(rnbinom(N,prob=0.3,size=0.15),rnbinom(2*N,prob=0.1,size=0.15),rnbinom(N,prob=0.6,size=0.15))
res2=Segmentor(y, model=3,Kmax=10); 
#Finds the optimal segmentation in 1 to 10 segments with respect to the Negative Binomial model.
}


