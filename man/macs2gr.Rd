\name{macs2gr}
\alias{macs2gr}

\title{
Convenience function for MACS output conversion
}
\description{

Parses the output of MACS Peak finding algorithm and returns a \linkS4class{GRanges} object compatible to the down stream functions of TransView

}
\usage{
macs2gr(macs_peaks_xls, psize, amount="all", min_pileup=0, log10qval=0, log10pval=0, fenrichment=0, peak_mid="summit")
}

\arguments{
  \item{macs_peaks_xls}{
Full path to the file ending with \sQuote{_peaks.xls} located in the output folder of a MACS run.
}
  \item{psize}{
An integer setting the total length of the peaks. Setting psize to \sQuote{preserve} will keep the original peak lengths from the output file and override \code{peak_mid}. Note that this is not compatible with \code{plotTV}
}
  \item{amount}{
Amount of peaks returned. If an integer is provided, the returned peaks will be limited to this amount after sorting by pile up score.
}
  \item{min_pileup}{
Minimum pile up.
}
  \item{log10qval}{
Minimal log10 q-value
}
  \item{log10pval}{
Minimal log10 p-value
}
  \item{fenrichment}{
Minimal enrichment.
}
  \item{peak_mid}{
  If set to \sQuote{summit}, the peaks with length \code{psize} will centered on the peak summit. If set to \sQuote{center}, the mid point of start and end will be used.
}
}
\details{
Convenience function parsing the output of a MACS file. Tested with MACS v1.4 and v.2.09
}

\value{
\code{GRanges} object with one row per peak and meta data score, enrichment and log10 pvalue.
}
\author{Julius Muller \email{ju-mu@alumni.ethz.ch}}
\examples{

exls<-dir(system.file("extdata", package="TransView"),full=TRUE,patt="xls$")

peaks<-macs2gr(exls,psize=500)
head(peaks)

}


