\name{slice1}
\alias{slice1}
\alias{sliceN}
\alias{slice1-methods}
\alias{sliceN-methods}
\alias{slice1,DensityContainer,character,numeric,numeric-method}
\alias{sliceN,DensityContainer-method}

\title{
Slice read densities from a TransView dataset
}
\description{

slice1 returns read densities of a genomic interval. sliceN takes a GRanges object or a data.frame with 
genomic coordinates and returns a list of read densities. 

}
\usage{
\S4method{slice1}{DensityContainer,character,numeric,numeric}(dc, chrom, start, end, control=FALSE, input_method="-",treads_norm=TRUE, nbins=0, bin_method="mean")
\S4method{sliceN}{DensityContainer}(dc, ranges, toRle=FALSE, control=FALSE, input_method="-",treads_norm=TRUE, nbins=0, bin_method="mean")
}

\arguments{
  \item{dc}{
Source \linkS4class{DensityContainer} object
}
  \item{chrom}{
A case sensitive string of the chromosome
}
  \item{start,end}{
Genomic start and end of the slice
}
  \item{ranges}{
A \linkS4class{GRanges} object or a data.frame. 
}  
  \item{toRle}{
The return values will be converted to a \code{RleList}. 
}
  \item{control}{
An optional \linkS4class{DensityContainer} which will used as control and by default subtracted from \code{dc}.
}
  \item{input_method}{
Defines the handling of the optional control \linkS4class{DensityContainer}. \sQuote{-} will subtract the control from the actual data and \sQuote{/} will return log2 fold change ratios with an added pseudo count of 1 read.
} 
  \item{treads_norm}{
If \code{TRUE}, the input densities are normalized to the read counts of the data set. Should not be used if one of the \code{DensityContainer} objects does not contain the whole amount of reads by e.g. placing a filter in \code{parseReads}. 
}
  \item{nbins}{
If all input regions have equal length and nbins greater than 0, all densities will be summarized using the method specified by bin_method into nbins windows of approximately equal size. 
}
  \item{bin_method}{
Character string that specifies the function used to summarize or expand the bins specified by nbins. Valid methods are \sQuote{max}, \sQuote{mean} or \sQuote{median}.
}
}
\details{

slice1 is a fast method to slice a vector of read densities from a \linkS4class{DensityContainer} object.
The vector can be optionally background subtracted. If the query region exceeds chromosome boundaries 
or if an non matching chromosome name will be passed, a warning will be issued and a NULL vector will be returned.

\code{sliceN} returns a list with N regions corresponding to N rows in the \linkS4class{GRanges} object or the data.frame. A list with the 
corresponding read densities will be returned and row names will be conserved. Optionally 
the return values can be converted to a \code{RleList} for seamless integration into the 
\pkg{IRanges} package.

}

\value{
slice1 returns a numeric vector of read densities
sliceN returns a list of read densities and optionally an \code{RleList}
}
\author{Julius Muller \email{ju-mu@alumni.ethz.ch}}

\seealso{
  \itemize{
    \item \code{\link[TransView]{slice1T}}.
    \item \link[TransView]{DensityContainer-class}.
  }
}


\examples{

exbam<-dir(system.file("extdata", package="TransView"),full=TRUE,patt="bam$")
exls<-dir(system.file("extdata", package="TransView"),full=TRUE,patt="xls$")

#store density maps of the whole sam/bam file in test_data
exden.ctrl<-parseReads(exbam[1],verbose=0)
exden.chip<-parseReads(exbam[2],verbose=0)

peaks<-macs2gr(exls,psize=500)

#returns vector of read counts per base pair
slice1(exden.chip,"chr2",30663080,30663580)[300:310]
slice1(exden.ctrl,"chr2",30663080,30663580)[300:310]
slice1(exden.chip,"chr2",30663080,30663580,control=exden.ctrl,treads_norm=FALSE)[300:310]

xout<-sliceN(exden.chip,ranges=peaks)
lapply(xout,function(x)sum(x)/length(x))
xout<-sliceN(exden.ctrl,ranges=peaks)
lapply(xout,function(x)sum(x)/length(x))
xout<-sliceN(exden.chip,ranges=peaks,control=exden.ctrl,treads_norm=FALSE)
lapply(xout,function(x)sum(x)/length(x))


}
