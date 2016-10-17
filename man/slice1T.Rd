\name{slice1T}
\alias{slice1T}
\alias{sliceNT}
\alias{slice1T-methods}
\alias{sliceNT-methods}
\alias{slice1T,DensityContainer,character-method}
\alias{sliceNT,DensityContainer,character-method}

\title{
Slice read densities of whole transcripts from a TransView DensityContainer
}
\description{

slice1T returns read densities of a transcript. sliceNT takes the output of with 
genomic coordinates and returns a list of read densities. 

}
\usage{
\S4method{slice1T}{DensityContainer,character}(dc, tname,  gtf, control=FALSE, input_method="-", concatenate=T, stranded=T, treads_norm=T, nbins=0, bin_method="mean")
\S4method{sliceNT}{DensityContainer,character}(dc, tnames,  gtf, toRle=FALSE, control=FALSE, input_method="-", concatenate=T, stranded=T, treads_norm=T, nbins=0, bin_method="mean")
}

\arguments{
  \item{dc}{
Source \linkS4class{DensityContainer} object
}
  \item{tname,tnames}{
A character string or a character vector with matching identifiers of the provided gtf
}
  \item{gtf}{
A \linkS4class{GRanges} object with a meta data column \sQuote{transcript_id} and \sQuote{exon_id} like e.g. from \code{gtf2gr}. 
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
  \item{concatenate}{
Logical that determines whether exons will be concatenated to one numeric vector (default) or returned as a list of vectors per exon.
}
  \item{stranded}{
If TRUE, the resulting vector will be reversed for reads on the reverse strand.
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
\code{slice1T} and \code{sliceNT} provide a convenient method to access the read densities from a \linkS4class{DensityContainer}
of spliced reads. The transcript structure will be constructed based on the provided gtf information.

slice1T is a fast alternative to sliceNT to slice one vector of read densities corresponding 
to the structure of one transcript and reads can be optionally background subtracted.
If the query region exceeds chromosome boundaries or if an non matching chromosome name will be
passed, a warning will be issued and a NULL vector will be returned.

sliceN slices N regions corresponding to N rows in the range GRanges object. A list with the 
corresponding read densities will be returned and row names will be conserved. Optionally 
the return values can be converted to a \code{RleList} for seamless integration into the 
\pkg{IRanges} package.

}

\value{
slice1T returns a numeric vector of read densities
sliceNT returns a list of read densities and optionally an \code{RleList}
}
\author{Julius Muller \email{ju-mu@alumni.ethz.ch}}
\examples{

library("pasillaBamSubset")

exgtf<-dir(system.file("extdata", package="TransView"),full=TRUE,patt="gtf.gz$")[1]
fn.pas_paired<-untreated1_chr4()

exden.exprs<-parseReads(fn.pas_paired,spliced=TRUE,verbose=0)

GTF.dm3<-gtf2gr(exgtf)

slice1T(exden.exprs,tname="NM_001014688",gtf=GTF.dm3,concatenate=FALSE)

my_genes<-sliceNT(exden.exprs,unique(mcols(GTF.dm3)$transcript_id[101:150]),gtf=GTF.dm3)
lapply(my_genes,function(x)sum(x)/length(x))


}
