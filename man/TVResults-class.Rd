\name{TVResults-class}
\Rdversion{1.1}
\docType{class}

\alias{TVResults-class}
\alias{class:TVResults}
\alias{show,TVResults-method}

\alias{parameters,TVResults-method}
\alias{parameters}
\alias{clusters,TVResults-method}
\alias{clusters}
\alias{cluster_order,TVResults-method}
\alias{cluster_order}
\alias{scores_peaks,TVResults-method}
\alias{scores_peaks}
\alias{scores_rna,TVResults-method}
\alias{scores_rna}
\alias{summaryTV}
\alias{summaryTV,TVResults-method}

\title{Class \code{"TVResults"}}
\description{
Container holding the results of a call to \code{plotTV()}.
}
\section{Objects from the Class}{
Objects are created by the function \code{plotTV()} using an internal constructor.
}


\section{Accessors}{
tvr represents a \code{"TVResults"} instance in the following
  \describe{
	\item{\code{parameters(tvr)}:}{Holds all parameters used to call plotTV}  
	\item{\code{clusters(tvr)}:}{Returns numeric vector with the clsuter of each cluster.}
	\item{\code{cluster_order(tvr)}:}{Ordering of the rows within the original regions passed to plotTV with regard to the clusters.}
	\item{\code{scores_peaks(tvr)}:}{Scores of the peaks. Corresponds to the values within the plot after interpolation and normalization.}
	\item{\code{scores_rna(tvr)}:}{Scores of the transcripts. Corresponds to the values within the plot after interpolation and normalization.}
	\item{\code{summaryTV(tvr)}:}{Returns a data frame with the clustering results of the internal data.}
  }
}


\section{Convenience Methods}{
  \describe{
    \item{plotTVData}{\code{signature(tvr = "TVResults")}: Returns a data frame with summarized clustering results.}
	}
}


\author{Julius Muller \email{ju-mu@alumni.ethz.ch}}
\note{
Not all slots are currently being exported.
}

\seealso{
\code{\link[TransView]{plotTVData-methods}}
}


\examples{
showClass("TVResults")
}
\keyword{classes}
