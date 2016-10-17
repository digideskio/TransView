\name{rmTV}
\alias{rmTV}
\alias{rmTV-methods}
\alias{rmTV,DensityContainer-method}

\title{
Free space occupied by DensityContainer
}
\description{
Free space occupied by DensityContainer
}
\usage{
\S4method{rmTV}{DensityContainer}(dc)
}

\arguments{
	\item{dc}{
		An object of class \linkS4class{DensityContainer}.
	}
}

\value{
None
}
\author{Julius Muller \email{ju-mu@alumni.ethz.ch}}
\examples{
exbam<-dir(system.file("extdata", package="TransView"),full=TRUE,patt="bam$")

#store density maps of the whole sam/bam file in test_data
exden.chip<-parseReads(exbam[2])

rmTV(exden.chip)
}


