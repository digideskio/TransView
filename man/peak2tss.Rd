\name{peak2tss}
\alias{peak2tss}

\title{
Changes the peak center to the next TSS according to previous annotation 
}

\description{
Sets the peak boundaries of an annotated \linkS4class{GRanges} object with peak locations to TSS centered ranges based on the transcript_id column.
}

\usage{
peak2tss(peaks, gtf, peak_len=500)
}

\arguments{
  \item{peaks}{
An annotated \linkS4class{GRanges} object with a meta data column \sQuote{transcript_id} and \sQuote{exon_id} like e.g. from \code{gtf2gr}.
}
  \item{gtf}{
A \linkS4class{GRanges} object with a meta data column \sQuote{transcript_id} like e.g. from \code{annotatePeaks}.
}
  \item{peak_len}{
The desired total size of the region with the TSS located in the middle.
}
}

\details{
Convenience function to change the peak centers to TSS for e.g. plotting with \code{plotTV}.
}

\value{
A \linkS4class{GRanges} object
}
\author{Julius Muller \email{ju-mu@alumni.ethz.ch}}
\examples{

exgtf<-dir(system.file("extdata", package="TransView"),full=TRUE,patt="gtf.gz$")[2]
fn.macs<-dir(system.file("extdata", package="TransView"),full=TRUE,patt="xls$")

GTF<-gtf2gr(exgtf)
peaks<-macs2gr(fn.macs,psize=500)

peaks.anno<-annotatePeaks(peaks=peaks,gtf=GTF)

peak2tss(peaks.anno, GTF, peak_len=500)

}


