\name{tvStats-methods}
\docType{methods}
\alias{tvStats-methods}
\alias{tvStats}
\alias{tvStats,DensityContainer-method}


\title{ DensityContainer accessor function }

\description{
Retrieve important metrics from the outcome of \code{parseReads()} stored in class DensityContainer and its super classes.
}

\usage{
\S4method{tvStats}{DensityContainer}(dc)
}

\arguments{

	\item{dc}{
	An object of class \linkS4class{DensityContainer}.
}
}


\value{
Returns a \code{list} with the slots of the DensityContainer and its super classes.
In detail:
\itemize{

	\item "ex_name": A user provided string to define a name of this dataset
	\item "origin": Filename of the original file
	\item "spliced": Should the class be treated like an RNA-Seq experiment for e.g. plotTV?
	\item "paired": Does the source file contain reads with proper pairs? 
	\item "readthrough_pairs": If TRUE, paired reads will be connected from left to right as one long read.
	\item "filtered": Is there a range filter in place? If yes, slicing should be \bold{only} conducted using the same filter!!
	\item "strands": Which strands were parsed at all. Can be "+", "-" or "both"
	\item "nreads": Total number of reads
	\item "coverage": Total coverage computed by total map mass/(chromosome end - chromosome start). Chromosome length derived from the SAM/BAM header
	\item "maxScore": Maximum read pileup found in file
	\item "lowqual": Amount of reads that did not pass the quality score set by min_quality or were not mapped
	\item "paired_reads": Amount of reads having multiple segments in sequencing
	\item "proper_pairs": Amount of pairs with each segment properly aligned according to the aligner
	\item "collapsed": If maxDups is in place, the reads at the same position and strand exceeding this value will be counted here.
	\item "compression": Size of a gap triggering an index event
	\item "chromosomes": Character string with the chromosomes with reads used for map construction
	\item "filtered":_reads Amount of reads
	\item "pos": Reads used from the forward strand
	\item "neg": Reads used from the reverse strand
	\item "lcoverage": Local coverage which is computed by filtered map mass/covered region
	\item "lmaxScore": Maximum score of the density maps
	\item "size": Size in bytes occupied by the object
}
}
\author{Julius Muller \email{ju-mu@alumni.ethz.ch}}

\keyword{methods}

