\name{BestSegmentation}
\alias{BestSegmentation}
\title{BestSegmentation}
\description{
  This function is used to compute the cost of the best segmentation in K segments with change-point t. 
}

\usage{
BestSegmentation(x,K,t=numeric())
}

\arguments{
  \item{x}{An object of class Segmentor returned by function Segmentor}
  \item{K}{The number of segments of the segmentation for which the cost or best segmentation is wanted}
  \item{t}{The position for which the best segmentation with t as change-point is wanted}
}

\value{
  \item{BestSegmentation}{A matrix of size n*K: the cost of the optimal segmentation with jth change-point i} 
}



\references{
  PDPA: Rigaill, G. Pruned dynamic programming for optimal multiple change-point detection: \emph{Submitted} \url{http://arxiv.org/abs/1004.0887}

  PDPA: Cleynen, A. and Koskas, M. and Rigaill, G. A Generic Implementation of the Pruned Dynamic Programing Algorithm: \emph{Submitted} \url{http://arxiv.org/abs/1204.5564}

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
rnbinom(N,size=1.3,prob=0.01));
res=Segmentor(data=x,model=3,Kmax=10, keep=TRUE);  
# Finds the optimal segmentation in up to 10 segments with respect to 
#the negative binomial model.
K<-3
Best<-BestSegmentation(res,K=3,t=3000)
matplot(Best$bestCost, type='l', lty=2)
points(apply(Best$bestCost, 2, which.min), apply(Best$bestCost, 2, min), pch=20, col=1:(K-1))
apply(Best$bestCost, 2, which.min)
getBreaks(res)[K,1:(K-1)]
#computes and plots cost of best segmentation in 3 segments with 
#change-point t, and compares result with change-point estimates.
Best$bestSeg
#returns the optimal segmentation in 3 segments with 3000 as a
#change-point
}


