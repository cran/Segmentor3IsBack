\name{getData}
\alias{getData}

\title{ Generic function}
\description{
 Generic function
}
\usage{
getData(object)
}
\arguments{
  \item{object}{An object of class Segmentor}
}

\details{
	Returns the slot data of an object of class Segmentor
}
\value{
	The vector of data to be segmented
}
\references{
  PDPA: Rigaill, G. Pruned dynamic programming for optimal multiple change-point detection: \emph{Submitted} \url{http://arxiv.org/abs/1004.0887}

  PDPA: Cleynen, A. and Koskas, M. and Lebarbier, E. and Rigaill, G. and Robin, S. Segmentor3IsBack (2014): an R package for the fast and exact
segmentation of Seq-data \emph{Algorithms for Molecular Biology}

}
\author{Alice Cleynen}


\seealso{
\code{\link{getModel}} }


\examples{
x=new("Segmentor") # new Segmentor object
getData(x) # retrieves the data slot from x
}

\keyword{methods}
\keyword{Data}
\keyword{changepoints}
