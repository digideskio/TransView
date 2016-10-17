\name{annotatePeaks}
\alias{annotatePeaks}

\title{
Associates peaks to TSS
}

\description{
A convenience function to associate the peak center to a TSS or gene body provided by a gtf file. 
}

\usage{
annotatePeaks(peaks, gtf, limit=c(-10e3,10e3), remove_unmatched=T, unifyBy=F, unify_fun="mean", min_genelength=0, reference="tss")
}

\arguments{
  \item{peaks}{
A \linkS4class{GRanges} object. 
}
  \item{gtf}{
A \linkS4class{GRanges} object with a meta data column \sQuote{transcript_id} and \sQuote{exon_id} like e.g. from \code{gtf2gr}. 
}
  \item{limit}{
Maximal distance range for a peak - TSS association in base pairs.  
}
  \item{remove_unmatched}{
If TRUE, only TSS associated peaks will be returned.
}
  \item{unifyBy}{
If a transcript has multiple isoforms, the peak will be associated arbitrarily to the first ID found. In order associate a peak to an isoform with specific characteristics, a \code{DensityContainer} can be provided. The choice of the returned isoform will be made based on unify_fun.
}
  \item{unify_fun}{
A function which will choose the isoform in case of non unique peak - TSS associations. Defaults to the isoform with the highest mean score \code{function(x){mean(x)}}.
}
  \item{min_genelength}{
Genes with a total sum of all exons smaller than this value will not be associated to a peak.
}
  \item{reference}{
If set to \sQuote{tss}, the transcript with the smallest distance from the TSS to the peak center will be returned. If set to \sQuote{gene_body} the transcript with the smallest distance from the gene body (TSS or TES) to the peak center will be returned and the distance will be zero if the peak center is located within the gene body.  
}
}

\details{
Convenience function to annotate a \linkS4class{GRanges} object having one row per peak from e.g. \code{macs2gr}. The resulting peak - TSS associations can be customized by the restricting the distance and resolving multiple matches using unify_fun.
}

\value{
\linkS4class{GRanges} object with row names according to the peak names provided and an added or updated meta data column \sQuote{transcript_id} with the associated transcript IDs and distances.
}
\author{Julius Muller \email{ju-mu@alumni.ethz.ch}}
\examples{

exgtf<-dir(system.file("extdata", package="TransView"),full=TRUE,patt="gtf.gz$")[2]
exls<-dir(system.file("extdata", package="TransView"),full=TRUE,patt="xls$")

GTF<-gtf2gr(exgtf)
peaks<-macs2gr(exls,psize=500)
apeaks<-annotatePeaks(peaks=peaks,gtf=GTF)
apeaks.gb<-annotatePeaks(peaks=peaks,gtf=GTF,reference="gene_body")

}


