\name{getKmax}
\alias{getKmax}

\title{ Generic function}
\description{
 Generic function
}
\usage{
getKmax(object)
}
\arguments{
  \item{object}{An object of class Segmentor}
}

\details{
	Returns the slot Kmax of an object of class Segmentor
}
\value{
	An integer of the maximum number of segments considered for the analysis of the dataset.
}
\references{
  PDPA: Rigaill, G. Pruned dynamic programming for optimal multiple change-point detection: \emph{Submitted} \url{http://arxiv.org/abs/1004.0887}

  PDPA: Cleynen, A. and Koskas, M. and Rigaill, G. A Generic Implementation of the Pruned Dynamic Programing Algorithm: \emph{Submitted} \url{http://arxiv.org/abs/1204.5564}
}
\author{Alice Cleynen}


\seealso{
\code{\link{getModel}} }


\examples{
x=new("Segmentor") # new EBS object
getKmax(x) # retrieves the Kmax slot from x
}

\keyword{methods}
\keyword{Kmax}
