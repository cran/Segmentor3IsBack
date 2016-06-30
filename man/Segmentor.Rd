\name{Segmentor}
\alias{Segmentor}
\alias{print.Segmentor}
\alias{Segmentor.default}
\title{Segmentor}
\description{
  The functions are used for change-point problems. Given a loss function (Poisson, Normal homoscedastic, Negative Binomial, Normal Heteroscedastic (with given constant mean) or Exponential), the function Segmentor estimates the optimal segmentation with respect to the log-likelihood. The Segmentor gives estimates of the breakpoint locations as well as the loss function parameter of each segment.
}

\usage{
Segmentor(data=numeric(), model=1, Kmax = 15, phi = numeric(), 
m = numeric(), keep=FALSE, bounds=c(0,0), compress = TRUE)
}

\arguments{
  \item{data}{A vector of observations to be segmented. Must have no missing values.}
  \item{model}{Integer between 1 and 4 giving the modelisation of the observed data, 1: poisson (default), 2: normal-homoscedastic, 3: negative binomial, 4: normal-heteroscedastic or 5: exponential}
  \item{Kmax}{The maximum number oe segments wanted for the data. The Segmentor will find all optimal segmentations in 1 to Kmax segments.}
  \item{phi}{Needed only for the Negative Binomial distribution: the value of the inverse of the overdispersion parameter. If the user does not enter a value, the package uses a modified version of Johnson and Kotz's estimator where the mean is replaced by the median. }
  \item{m}{Needed only for Normal Heteroscedastic distribution: the value of the constant mean. If not entered, the function uses the empirical mean of the data.}
  \item{keep}{a boolean stating whether or not to keep Cost and Position matrices}
  \item{bounds}{optional lower and upper bounds on the parameter to segment : will fasten the pruning and hence the algorithm if some values are not allowed (e.g. non negative numbers for Poisson distribution)}
  \item{compress}{A boolean stating whether data should be compressed prior to segmentation}
}

\value{
	\item{data}{The vector of observations to be segmented.} 
  \item{model}{Emission distribution (Poisson, Normal Homoscedastic, Negative Binomial or Normal Heteroscedastic or exponential)} 
  \item{breaks}{Matrix of size Kmax*Kmax of estimated change-point locations for each optimal segmentation in 1 to Kmax segments.}
  \item{parameters}{Matrix of size Kmax*Kmax which elements are the estimated parameters for each segment of the optimal segmentation. If model is Poisson or Normal, the parameter corresponds to the mean of the signal in each segment. If model is Negative binomial, the parameter corresponds to the success-probability of the signal in each segment. If model is normal heteroscedastic, the parameter is the variance assuming known mean.} 
  \item{likelihood}{Vector of size Kmax of resulting negative log-likelihood for each optimal segmentation.}
  \item{Cost}{Matrix of size Kmax x n containing the cost of the segmentation of signal up to point j in i segments} 
  \item{Pos}{Matrix of size Kmax x n containing the last change-point location of the segmentation of signal up to point j in i segments} 
  \item{overdispersion}{only if model = Negative Binomial, the value of the inverse of overdispersion used for the segmentation}
  \item{mean}{only if model = Normal Heteroscedastic, the value of the mean used for the segmentation}
  \item{compression}{The value of the compression factor used (compression>=1)}
}



\references{
  PDPA: Rigaill, G. Pruned dynamic programming for optimal multiple change-point detection: \emph{Submitted} \url{http://arxiv.org/abs/1004.0887}

  PDPA: Cleynen, A. and Koskas, M. and Lebarbier, E. and Rigaill, G. and Robin, S. Segmentor3IsBack (2014): an R package for the fast and exact
segmentation of Seq-data \emph{Algorithms for Molecular Biology}

  overdispersion parameter: Johnson, N. and Kemps, A. and Kotz, S. (2005) Univariate Discrete Distributions: \emph{John Wiley & Sons} 
  
  variance parameter: Hall, P. and Kay, J. and Titterington, D. (1990): Asymptotically optimal difference-based estimation of variance in non-parametric regression \emph{Biometrika}
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
x=c(rpois(N,2.0),rpois(2*N,2.2),rpois(N,1.9));
res=Segmentor(data=x,Kmax=3);  
# Finds the optimal segmentation in 1, 2 and 3 segments with respect to 
#the Poisson model.

y=c(rnbinom(N,prob=0.3,size=0.15),rnbinom(2*N,prob=0.1,size=0.15),
rnbinom(N,prob=0.6,size=0.15),compress=FALSE)
res2=Segmentor(y, model=3,Kmax=10); 
#Finds the optimal segmentation in 1 to 10 segments with respect to 
#the Negative Binomial model, without compression of data.
}


