\name{getModel}
\alias{getModel}

\title{ Generic function}
\description{
 Generic function
}
\usage{
getModel(object)
}
\arguments{
  \item{object}{An object of class Segmentor}
}

\details{
	Returns the slot Model of an object of class Segmentor
}
\value{
	An object of class string returning the model used for the analysis of the dataset.
}
\references{
  PDPA: Rigaill, G. Pruned dynamic programming for optimal multiple change-point detection: \emph{Submitted} \url{http://arxiv.org/abs/1004.0887}

  PDPA: Cleynen, A. and Koskas, M. and Rigaill, G. A Generic Implementation of the Pruned Dynamic Programing Algorithm: \emph{Submitted} \url{http://arxiv.org/abs/1204.5564}
}
\author{Alice Cleynen}


\seealso{
\code{\link{getKmax}} }


\examples{
x=new("Segmentor") # new Segmentor object
getModel(x) # retrieves the Model slot from x
}

\keyword{methods}
\keyword{Model}
