\name{Segmentor-class}
\Rdversion{1.1}
\docType{class}
\alias{Segmentor-class}
\alias{Segmentor-method, show}
\alias{Segmentor-method,getData}
\alias{Segmentor-method,getModel}
\alias{Segmentor-method,getBreaks}
\alias{Segmentor-method,getLikelihood}
\alias{Segmentor-method,getOverdispersion}
\alias{Segmentor-method,getParamters}
\alias{Segmentor-method,getKmax}
\alias{Segmentor-method,getMean}
\alias{Segmentor-method,getCost}
\alias{Segmentor-method,getPos}
\alias{Segmentor-method,getCompression}

\title{Class "Segmentor"}
\description{
A class for Segmentation objects.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Segmentor", ...)}.
  \describe{
    \item{\code{new("Segmentor", ...)}:}{creates a new object with class Segmentor }
  }
}

         
\section{Slots}{
  \describe{
  	\item{\code{data}:}{Object of class \code{"numeric"}, the vector of data }
    \item{\code{model}:}{Object of class \code{"character"}, the assumed distribution of the data }
    \item{\code{Kmax}:}{Object of class \code{"numeric"}, the maximum number of segments considered for the segmentation}
    \item{\code{breaks}:}{Object of class \code{"matrix"}, the breakpoints of the optimal segmentation for each value of k}
    \item{\code{mean}:}{Object of class \code{"numeric"}, if model = Variance, the mean used in the analysis}
    \item{\code{overdispersion}:}{Object of class \code{"numeric"}, if model = Negative Binomial, the overdispersion used in the analysis}
    \item{\code{parameters}:}{Object of class \code{"matrix"}, the ML estimates of the parameter of each segment of the optimal segmentation for each value of k}
    \item{\code{likelihood}:}{Object of class \code{"matrix"}, the likelihood of the segmentation for each value of k}
    \item{\code{Cost}:}{Object of class \code{"matrix"}, the cost of the segmentation up to point j in i segments}
    \item{\code{Pos}:}{Object of class \code{"matrix"}, the last change-point location of the segmentation up to point j in i segments}
    \item{\code{compression}:}{Object of class \code{"numeric"}, the value of the compression obtained prior to segmentation of the data}
  }
}
\section{Methods}{
  \describe{
  	\item{getData}{\code{signature(object = "Segmentor")}: retrieves data slot }
    \item{getModel}{\code{signature(object = "Segmentor")}: retrieves model slot }
    \item{getKmax}{\code{signature(object = "Segmentor")}: retrieves Kmax slot }
    \item{getBreaks}{\code{signature(object = "Segmentor")}: retrieves breaks slot }
    \item{getMean}{\code{signature(object = "Segmentor")}: retrieves mean slot }
    \item{getOverdispersion}{\code{signature(object = "Segmentor")}: retrieves overdispersion slot }
    \item{getLikelihood}{\code{signature(object = "Segmentor")}: retrieves likelihood slot }
    \item{getCost}{\code{signature(object = "Segmentor")}: retrieves Cost slot }
    \item{getParameters}{\code{signature(object = "Segmentor")}: retrieves parameters slot }
    \item{getPos}{\code{signature(object = "Segmentor")}: retrieves Pos slot }
    \item{getCompression}{\code{signature(object = "Segmentor")}: retrieves compression slot }
  }
}

\author{
Alice Cleynen
}

\seealso{
\code{\link{SelectModel}}
}
\keyword{classes}
\keyword{Segmentor}
