\name{getCost}
\alias{getCost}

\title{ Generic function}
\description{
 Generic function
}
\usage{
getCost(object)
}
\arguments{
  \item{object}{An object of class Segmentor}
}

\details{
	Returns the slot Cost of an object of class Segmentor
}
\value{
	The matrix of Cost of the segmentation in i segments up to point j
}
\references{
  PDPA: Rigaill, G. Pruned dynamic programming for optimal multiple change-point detection: \emph{Submitted} \url{http://arxiv.org/abs/1004.0887}

  PDPA: Cleynen, A. and Koskas, M. and Lebarbier, E. and Rigaill, G. and Robin, S. Segmentor3IsBack (2014): an R package for the fast and exact
segmentation of Seq-data \emph{Algorithms for Molecular Biology}

}
\author{Alice Cleynen}


\seealso{
\code{\link{getLikelihood}}; \code{\link{getBreaks}} }


\examples{
x=new("Segmentor") # new Segmentor object
getCost(x) # retrieves the Cost slot from x
}

\keyword{methods}
\keyword{Cost}
\keyword{changepoints}
