\name{getLikelihood}
\alias{getLikelihood}

\title{ Generic function}
\description{
 Generic function
}
\usage{
getLikelihood(object)
}
\arguments{
  \item{object}{An object of class Segmentor}
}

\details{
	Returns the slot likelihood of an object of class Segmentor
}
\value{
	The opposite of the log-likelihood of the optimal segmentation for each value of k
}
\references{
  PDPA: Rigaill, G. Pruned dynamic programming for optimal multiple change-point detection: \emph{Submitted} \url{http://arxiv.org/abs/1004.0887}

  PDPA: Cleynen, A. and Koskas, M. and Lebarbier, E. and Rigaill, G. and Robin, S. Segmentor3IsBack (2014): an R package for the fast and exact
segmentation of Seq-data \emph{Algorithms for Molecular Biology}


}
\author{Alice Cleynen}


\seealso{
\code{\link{getBreaks}} }


\examples{
x=new("Segmentor") # new Segmentor object
getLikelihood(x) # retrieves the likelihood slot from x
}

\keyword{methods}
\keyword{Likelihood}
\keyword{changepoints}
