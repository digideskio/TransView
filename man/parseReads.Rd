\name{parseReads}
\alias{parseReads}

\title{
User configurable efficient assembly of read density maps
}
\description{

Generates density maps for further downstream processing. Constructs a \code{DensityContainer}.

}
\usage{
parseReads( filename, spliced=F, read_stranded=0, paired_only=F, readthrough_pairs=F, set_filter=NA, min_quality=0,
		description="NA", extendreads=0, unique_only=F,	max_dups=0, hwindow=1, compression=1, verbose=1 )
}

\arguments{
  \item{filename}{
Character string with the filename of the bam file. The bam file must be sorted according to genomic position.
}
  \item{spliced}{
This option will mark the object to be treated like a data set with spliced reads. Can be switched off also for spliced experiments for special purposes. If \code{TRUE}, switches off \code{extendreads} and \code{readthrough_pairs}. 
}  
  \item{read_stranded}{
0 will read tags from both strands. 1 will skip all tags from the \sQuote{-} strand and -1 will only utilize tags from the \sQuote{-} strand 
}
  \item{paired_only}{
If \code{TRUE}, any reads which are not members of a proper pair according to the 0x0002 FLAG will be discarded. If \code{FALSE} all reads will be used individually.
}
  \item{set_filter}{
Optional \linkS4class{GRanges} object or data.frame with similar structure: data.frame(chromosomes,start,end). Providing this filter will limit density maps to these regions.
}
  \item{min_quality}{
Phred-scaled mapping quality threshold. If 0, all reads will pass this filter.
}
  \item{extendreads}{
If greater 0, this amount of base pairs will be added into the strand direction of each read during density map generation.
} 
  \item{unique_only}{
If TRUE, only unique reads with no multiple alignments will be used. This filter relies on the aligner to use the corresponding flag (0x100).
}
  \item{max_dups}{
If greater 0, maximally this amount of reads are allowed per start position and read direction.
}  
  \item{description}{
An optional character string describing the experiment for labeling purposes.
}
  \item{hwindow}{
A numeric defining the window size used to compute the histogram. This value cannot be bigger than \code{compression}
}
  \item{compression}{
Should be left at the default value. Defines the minimal threshold in base pairs which triggers indexing and collapsing of read free regions. A smaller value leads to faster slicing at the cost of a higher memory footprint. 
}
  \item{readthrough_pairs}{
Currently *experimental*. If \code{TRUE}, \code{parseReads} will attempt to use the region from the left to the right read of the pair for density map assembly. Requires ISIZE to be set within the BAM/SAM file. 
}
  \item{verbose}{
Verbosity level
}

}
\details{

parseReads uses read information of one bam file and scans the entire file read wise. Every read contributes 
to the density track in a user configurable manner. The resulting track will be stored in 
indexed integer vectors within a list. Since each score is stored as a unsigned 16bit integer, the scores can only 
be accessed with one of the slice methods \code{slice1} or \code{sliceN} and not directly. As a consequence of the storage 
format read pile ups greater than \bold{2^16} will be capped and a warning will be issued. 

If memory space is limiting, a filter can be supplied which will limit the density track to these regions. 
Filtered \code{DensityContainer} should only be sliced with the \bold{same} regions used for parsing, since 
all other positions are set to 0 and can produce artificially low read counts.

}

\value{
S4 \code{DensityContainer}
}
\author{Julius Muller \email{ju-mu@alumni.ethz.ch}}
\examples{

exbam<-dir(system.file("extdata", package="TransView"),full=TRUE,patt="bam$")

#store density maps of the whole sam/bam file in test_data
exden.chip<-parseReads(exbam[2],verbose=0)

#display basic information about the content of test.sam 
exden.chip

#all data are easily accessible
test_stat<-tvStats(exden.chip)
test_stat$origin

# histogram of hwindow sized windows
\dontrun{histogram(exden.chip)}

}


