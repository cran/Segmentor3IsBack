\name{getCompression}
\alias{getCompression}

\title{ Generic function}
\description{
 Generic function
}
\usage{
getCompression(object)
}
\arguments{
  \item{object}{An object of class Segmentor}
}

\details{
	Returns the slot compression of an object of class Segmentor
}
\value{
	the value of compression that was used when analysing the data.
}
\references{
  PDPA: Rigaill, G. Pruned dynamic programming for optimal multiple change-point detection: \emph{Submitted} \url{http://arxiv.org/abs/1004.0887}

  PDPA: Cleynen, A. and Koskas, M. and Lebarbier, E. and Rigaill, G. and Robin, S. Segmentor3IsBack (2014): an R package for the fast and exact
segmentation of Seq-data \emph{Algorithms for Molecular Biology}

  overdispersion parameter: Johnson, N. and Kemps, A. and Kotz, S. (2005) Univariate Discrete Distributions: \emph{John Wiley & Sons} 
}
\author{Alice Cleynen}


\seealso{
\code{\link{getMean}} }


\examples{
x=new("Segmentor") # new Segmentor object
getCompression(x) # retrieves the compression slot from x
}

\keyword{methods}
\keyword{Compression}
