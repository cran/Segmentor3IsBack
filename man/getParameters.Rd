\name{getParameters}
\alias{getParameters}

\title{ Generic function}
\description{
 Generic function
}
\usage{
getParameters(object)
}
\arguments{
  \item{object}{An object of class Segmentor}
}

\details{
	Returns the slot parameters of an object of class Segmentor
}
\value{
	The ML estimates of each segment's parameter of the optimal segmentation for each value of k
}
\references{
  PDPA: Rigaill, G. Pruned dynamic programming for optimal multiple change-point detection: \emph{Submitted} \url{http://arxiv.org/abs/1004.0887}

  PDPA: Cleynen, A. and Koskas, M. and Lebarbier, E. and Rigaill, G. and Robin, S. Segmentor3IsBack: an R package for the fast and exact
segmentation of Seq-data: \emph{Submitted} \url{http://arxiv.org/abs/1204.5564}

}
\author{Alice Cleynen}


\seealso{
\code{\link{getBreaks}} }


\examples{
x=new("Segmentor") # new Segmentor object
getParameters(x) # retrieves the parameters slot from x
}

\keyword{methods}
\keyword{parameters}
\keyword{changepoints}
