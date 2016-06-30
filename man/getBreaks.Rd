\name{getBreaks}
\alias{getBreaks}

\title{ Generic function}
\description{
 Generic function
}
\usage{
getBreaks(object)
}
\arguments{
  \item{object}{An object of class Segmentor}
}

\details{
	Returns the slot breaks of an object of class Segmentor
}
\value{
	The ML estimates of the breakpoints of the optimal segmentation for each value of k
}
\references{
  PDPA: Rigaill, G. Pruned dynamic programming for optimal multiple change-point detection: \emph{Submitted} \url{http://arxiv.org/abs/1004.0887}

  PDPA: Cleynen, A. and Koskas, M. and Lebarbier, E. and Rigaill, G. and Robin, S. Segmentor3IsBack (2014): an R package for the fast and exact
segmentation of Seq-data \emph{Algorithms for Molecular Biology}

}
\author{Alice Cleynen}


\seealso{
\code{\link{getLikelihood}} }


\examples{
x=new("Segmentor") # new Segmentor object
getBreaks(x) # retrieves the breaks slot from x
}

\keyword{methods}
\keyword{Breakpoints}
\keyword{changepoints}
