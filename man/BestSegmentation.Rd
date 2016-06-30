\name{BestSegmentation}
\alias{BestSegmentation}
\title{BestSegmentation}
\description{
  This function is used to compute the cost of the best segmentation in K segments given the position of a change-point, and to return the optimal segmenation having a change-point at location t. 
}

\usage{
BestSegmentation(x,K,t=numeric(),compress=TRUE)
}

\arguments{
  \item{x}{An object of class Segmentor returned by function Segmentor}
  \item{K}{The number of segments of the segmentation for which the cost or best segmentation is wanted}
  \item{t}{The position for which the best segmentation with t as change-point is wanted}
  \item{compress}{A boolean stating whether data should be compressed prior to segmentation}
}

\value{
  \item{bestCost}{A matrix of size n*K: the cost of the optimal segmentation with jth change-point i} 
  \item{bestSeg}{If a t has been specified, a vector of size K+1 containing values of indicating the optimal segmentation with t as a change-point} 
}



\references{
  PDPA: Rigaill, G. Pruned dynamic programming for optimal multiple change-point detection: \emph{Submitted} \url{http://arxiv.org/abs/1004.0887}

  PDPA: Cleynen, A. and Koskas, M. and Lebarbier, E. and Rigaill, G. and Robin, S. Segmentor3IsBack (2014): an R package for the fast and exact
segmentation of Seq-data \emph{Algorithms for Molecular Biology}

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
x=rnbinom(3*N,size=1.3,prob=rep(c(0.7,0.2,0.01),each=N));
res=Segmentor(data=x,model=3,Kmax=10, keep=TRUE);  
# Finds the optimal segmentation in up to 10 segments with respect to 
#the negative binomial model.
K<-3
Best<-BestSegmentation(res,K=3,t=3000,compress=FALSE)
matplot(Best$bestCost, type='l', lty=2)
points(apply(Best$bestCost,2,which.min),apply(Best$bestCost,2,min),pch=20,col=1:(K-1))
apply(Best$bestCost, 2,which.min)
getBreaks(res)[K,1:(K-1)]
#computes and plots cost of best segmentation in 3 segments with 
#change-point t, and compares result with change-point estimates.
Best$bestSeg
#returns the optimal segmentation in 3 segments with 3000 as a
#change-point
Best<-BestSegmentation(res,K=3,t=3000,compress=TRUE)
Best$bestSeg
# this segmentation usually does not make sense because of the
# compress option. t has to be adapted consequently
}


