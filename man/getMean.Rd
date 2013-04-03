\name{getMean}
\alias{getMean}

\title{ Generic function}
\description{
 Generic function
}
\usage{
getMean(object)
}
\arguments{
  \item{object}{An object of class Segmentor}
}

\details{
	Returns the slot Mean of an object of class Segmentor
}
\value{
	If model is Variance, the value of the mean used for the analysis.
}
\references{
  PDPA: Rigaill, G. Pruned dynamic programming for optimal multiple change-point detection: \emph{Submitted} \url{http://arxiv.org/abs/1004.0887}

  PDPA: Cleynen, A. and Koskas, M. and Rigaill, G. A Generic Implementation of the Pruned Dynamic Programing Algorithm: \emph{Submitted} \url{http://arxiv.org/abs/1204.5564}

  variance parameter: Hall, P. and Kay, J. and Titterington, D. (1990): Asymptotically optimal difference-based estimation of variance in non-parametric regression \emph{Biometrika}
}
\author{Alice Cleynen}


\seealso{
\code{\link{getOverdispersion}} }


\examples{
x=new("Segmentor",model="Variance",mean=2.2) # new Segmentor object
getMean(x) # retrieves the overdispersion slot from x
}

\keyword{methods}
\keyword{Mean}
