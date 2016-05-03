\name{DensityContainer-class}
\Rdversion{1.1}
\docType{class}

\alias{DensityContainer-class}
\alias{class:DensityContainer}
\alias{show,DensityContainer-method}

\alias{ex_name,DensityContainer-method}
\alias{ex_name}
\alias{ex_name<-,DensityContainer-method}
\alias{ex_name<-}
\alias{origin,DensityContainer-method}
\alias{origin}
\alias{spliced,DensityContainer-method}
\alias{spliced<-,DensityContainer-method}
\alias{spliced<-}
\alias{spliced}
\alias{paired,DensityContainer-method}
\alias{paired}
\alias{readthrough_pairs,DensityContainer-method}
\alias{readthrough_pairs}
\alias{filtered,DensityContainer-method}
\alias{filtered}
\alias{strands,DensityContainer-method}
\alias{strands}
\alias{nreads,DensityContainer-method}
\alias{nreads}
\alias{gcoverage,DensityContainer-method}
\alias{gcoverage}
\alias{maxScore,DensityContainer-method}
\alias{maxScore}
\alias{lowqual,DensityContainer-method}
\alias{lowqual}
\alias{paired_reads,DensityContainer-method}
\alias{paired_reads}
\alias{proper_pairs,DensityContainer-method}
\alias{proper_pairs}
\alias{collapsed,DensityContainer-method}
\alias{collapsed}
\alias{compression,DensityContainer-method}
\alias{compression}
\alias{chromosomes,DensityContainer-method}
\alias{chromosomes}
\alias{filtered_reads,DensityContainer-method}
\alias{filtered_reads}
\alias{pos,DensityContainer-method}
\alias{pos}
\alias{neg,DensityContainer-method}
\alias{neg}
\alias{lcoverage,DensityContainer-method}
\alias{lcoverage}
\alias{lmaxScore,DensityContainer-method}
\alias{lmaxScore}
\alias{fmapmass,DensityContainer-method}
\alias{fmapmass}
\alias{data_pointer,DensityContainer-method}
\alias{data_pointer}
\alias{env,DensityContainer-method}
\alias{env}
\alias{size,DensityContainer-method}
\alias{size}
\alias{lsize,DensityContainer-method}
\alias{lsize}
\alias{gsize,DensityContainer-method}
\alias{gsize}


\title{Class \code{"DensityContainer"}}
\description{
Container with the pointer of the actual density maps and a histogram. Inherits from internal classes storing informations about the origin and the details of the results.
}
\section{Objects from the Class}{
Objects are created by the function \code{parseReads()} using an internal constructor.
}
\section{Accessors}{
dc represents a \code{"DensityContainer"} instance in the following
  \describe{
    \item{\code{data_pointer(dc)}:}{A character string pointing to the read density map. 
    It points to a variable in .GlobalEnv which is essentially a list resulting from a call to \code{parseReads}.
    The storage space can be freed with the \code{rmTV} function. }
    \item{\code{ex_name(dc)},\code{ex_name(dc)<-value}:}{Get or set a string to define a name of this data set}
    \item{\code{origin(dc)}:}{Filename of the original file}
    \item{\code{histogram(dc)}:}{A histogram of read pile-ups generated across all read density maps after filtering excluding gaps. }
    \item{\code{env(dc)}:}{The environment which holds the data_pointer target. }
    \item{\code{spliced(dc)},\code{spliced(dc)<-bool}:}{This option will mark the object to be treated like a data set with spliced reads.}
    \item{\code{readthrough_pairs(dc)}:}{If TRUE, paired reads will be connected from left to right and used as one long read.}
    \item{\code{paired(dc)}:}{Does the source file contain reads with proper pairs?  }
    \item{\code{filtered(dc)}:}{Is there a range filter in place? If \code{TRUE}, slicing should be only conducted using the same filter!! }
    \item{\code{strands(dc)}:}{Which strands were parsed at all. Can be "+", "-" or "both" }
    \item{\code{filtered_reads(dc)}:}{FilteredReads class storing information about reads used for read density construction}
	\item{\code{chromosomes(dc)}:}{Character string with the chromosomes used for map construction}
	\item{\code{pos(dc)}:}{Reads used from the forward strand}
	\item{\code{neg(dc)}:}{Reads used from the reverse strand}
	\item{\code{lsize(dc)}:}{Total region covered by reads within the densities returned}
	\item{\code{gsize(dc)}:}{Equals to the sum of the length of all ranges from 0 to the last read per chromosome within the chromosome.}
	\item{\code{lcoverage(dc)}:}{Local coverage within the densities returned which is computed by local mapmass/lsize}
	\item{\code{lmaxScore(dc)}:}{Maximum read pileup within the density maps after filtering}
	\item{\code{fmapmass(dc)}:}{Total map mass after quality filtering present in the file. Equals to filtered_reads*read length}
	\item{\code{nreads(dc)}:}{Total number of reads in the file.}
	\item{\code{coverage(dc)}:}{Total coverage computed by total map mass/(chromosome end - chromosome start). Chromosome length derived from the SAM/BAM header}
	\item{\code{maxScore(dc)}:}{Maximum read pileup found in file after quality filtering}
	\item{\code{lowqual(dc)}:}{Amount of reads that did not pass the quality score set by min_quality or were not mapped}
	\item{\code{paired_reads(dc)}:}{Amount of reads having multiple segments in sequencing}
	\item{\code{proper_pairs(dc)}:}{Amount of pairs with each segment properly aligned according to the aligner}
	\item{\code{collapsed(dc)}:}{If maxDups is in place, the reads at the same position and strand exceeding this value will be counted here.   } 
	\item{\code{size(dc)}:}{Size in bytes occupied by the object.}    

  }
}

\section{Slice Methods}{
  \describe{
    \item{slice1}{\code{signature(dc = "DensityContainer")}: Fetch a slice of read densities. }
    \item{slice1T}{\code{signature(dc = "DensityContainer")}: Recover the structure of a gene from a provided pre-processed GTF and read densities. }
    \item{sliceN}{\code{signature(dc = "DensityContainer", ranges = "data.frame")}: Like slice1 but optimized for repeated slicing. }
    \item{sliceNT}{\code{signature(dc = "DensityContainer", tnames = "character", gtf = "data.frame")}: Like slice1T but optimized for repeated slicing. }
}
}

\section{Convenience Methods}{
  \describe{
    \item{tvStats}{\code{signature(dc = "DensityContainer")}: Returns a list of important metrics about the source file.}
	}
}

\section{Extends}{
Class \code{TransView}, directly.
}
\author{Julius Muller \email{ju-mu@alumni.ethz.ch}}
\note{
Class TotalReads and FilteredReads are not exported but their slots can be fully accessed by several accessors and the \code{tvStats()} method.
}

\seealso{
\code{\link[TransView]{tvStats-methods}},
\code{\link[TransView]{slice1-methods}},
\code{\link[TransView]{sliceN-methods}},
\code{\link[TransView]{histogram-methods}},
\code{\link[TransView]{rmTV-methods}}
}


\examples{
showClass("DensityContainer")
}
\keyword{classes}
