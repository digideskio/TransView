\name{meltPeak}
\alias{meltPeak}

\title{
Convenience function which returns a data frame with normalized peak densities suitable for plotting with ggplot2
}
\description{

Returns a data frame with labels and normalized densities of the provided \linkS4class{DensityContainer}

}
\usage{
meltPeak(..., region, control=FALSE, peak_windows = 0, bin_method="mean", rpm=TRUE, smooth=0)
}

\arguments{
  \item{...}{
\linkS4class{DensityContainer} objects
}
  \item{region}{
Can be one entry of the annotated output of annotatePeaks or a \linkS4class{GRanges} object with one entry and with a transcript_id and distance metadata column.
}
  \item{control}{
An optional vector of \linkS4class{DensityContainer} objects, that have to match the order of experiments passed as a first argument. E.g. \code{plotTV(ex1.ChIP,ex2.ChIP,control=c(ex1.Input,ex2.Input)}. The content will be treated as background densities and subtracted from the matching experiment.
}
  \item{peak_windows}{
If set to an integer greater than 0, all binding profiles will be interpolated into this amount of windows by the method specified by \code{bin_method}.
}
  \item{bin_method}{
Specifies the function used to summarize the bins specified by nbins. Possible methods are \sQuote{max}, \sQuote{mean}, \sQuote{median} or \sQuote{approx} for linear interpolation. 
}
  \item{rpm}{
If set to \code{TRUE}, all sample groups will be normalized to Reads Per Million mapped reads after quality filtering according to the filtered_reads slot of the \linkS4class{DensityContainer}. Should not be set in truncated density maps!
}
  \item{smooth}{
If greater than 0, smooth defines the smoother span as described in the function \code{lowess}. This function will be applied to reads or RPM values, depending on \code{rpm} and the results will be stored in the column \sQuote{Smooth}.
}
}
\details{
Convenience function which returns a data frame with one row per BP or, if \code{peak_window} greater than zero, per \code{peak_window}. The label will be taken from the \code{ex_name} slot of the \linkS4class{DensityContainer}. The slot should be set to meaningful names before using this function. All read densities will be normalized to the total map mass and if a control is provided also background subtracted.
}

\value{
\code{data.frame} with 3 columns: \sQuote{NormalizedReads}, \sQuote{Label} and \sQuote{Position}. Optionally a column \sQuote{Smooth} will be appended.
}
\author{Julius Muller \email{ju-mu@alumni.ethz.ch}}
\examples{

exbam<-dir(system.file("extdata", package="TransView"),full=TRUE,patt="bam$")
exls<-dir(system.file("extdata", package="TransView"),full=TRUE,patt="xls$")
exgtf<-dir(system.file("extdata", package="TransView"),full=TRUE,patt="gtf.gz$")[2]
fn.macs<-dir(system.file("extdata", package="TransView"),full=TRUE,patt="xls$")

exden.ctrl<-parseReads(exbam[1],verbose=0)
exden.chip<-parseReads(exbam[2],verbose=0)

peaks<-macs2gr(exls,psize=500)

GTF<-gtf2gr(exgtf)
peaks<-macs2gr(fn.macs,psize=500)
peaks.anno<-annotatePeaks(peaks=peaks,gtf=GTF)

peak1.df<-meltPeak(exden.chip,region=peaks.anno["Peak.1"],bin_method="mean",peak_windows=100,rpm=TRUE)
head(peak1.df)

}


