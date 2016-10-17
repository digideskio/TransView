
\name{histogram-methods}
\docType{methods}
\alias{histogram-methods}
\alias{histogram}
\alias{histogram,DensityContainer-method}

\title{ Histogram of the read distribution }

\description{
Retrieves the histogram computed by the \code{parseReads} function
}

\usage{
\S4method{histogram}{DensityContainer}(dc)
}


\arguments{

	\item{dc}{
	An object of class \linkS4class{DensityContainer}.
}
}

\details{
The histogram is computed by taking the running average within a window of window size as specified by the argument \code{hwindow} to the function \code{parseReads()}. 
The histogram is only counting local reads within the read density maps and outside of gaps or outside of possible range filters that might be in place.
}

\value{
Returns a numeric vector with the histogram in 1Bp resolution starting from 0.
}
\author{Julius Muller \email{ju-mu@alumni.ethz.ch}}

\keyword{methods}


