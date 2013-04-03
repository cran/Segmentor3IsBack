\name{getOverdispersion}
\alias{getOverdispersion}

\title{ Generic function}
\description{
 Generic function
}
\usage{
getOverdispersion(object)
}
\arguments{
  \item{object}{An object of class Segmentor}
}

\details{
	Returns the slot Overdispersion of an object of class Segmentor
}
\value{
	If model is Negative Binomial, the value of the overdispersion used for the analysis.
}
\references{
  PDPA: Rigaill, G. Pruned dynamic programming for optimal multiple change-point detection: \emph{Submitted} \url{http://arxiv.org/abs/1004.0887}

  PDPA: Cleynen, A. and Koskas, M. and Rigaill, G. A Generic Implementation of the Pruned Dynamic Programing Algorithm: \emph{Submitted} \url{http://arxiv.org/abs/1204.5564}

  overdispersion parameter: Johnson, N. and Kemps, A. and Kotz, S. (2005) Univariate Discrete Distributions: \emph{John Wiley & Sons} 
}
\author{Alice Cleynen}


\seealso{
\code{\link{getMean}} }


\examples{
x=new("Segmentor") # new Segmentor object
getOverdispersion(x) # retrieves the overdispersion slot from x
}

\keyword{methods}
\keyword{Overdispersion}
