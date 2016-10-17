\name{plotTVData}
\alias{plotTVData}
\alias{plotTVData-methods}
\alias{plotTVData,TVResults}
\alias{plotTVData,TVResults-method}


\title{
Summarize plotTV results
}
\description{

plotTVData returns the ordering and clustering results as internally calculated by plotTV. 

}
\usage{
\S4method{plotTVData}{TVResults}(tvr)
}

\arguments{
  \item{tvr}{
A \linkS4class{TVResults} object as returned by plotTV
}
}
\details{

If k-means or manual clustering was performed, row means per cluster will be returned in a data.frame. Otherwise row means over the whole data will be returned.

}

\value{
Returns a \code{data.frame} of the clustering results with five columns: Position, Cluster, Sample, Average_scores and Plot 
}
\author{Julius Muller \email{ju-mu@alumni.ethz.ch}}

\seealso{
  \itemize{
    \item \link[TransView]{TVResults-class}.
  }
}


\examples{

exbam<-dir(system.file("extdata", package="TransView"),full=TRUE,patt="bam$")
exls<-dir(system.file("extdata", package="TransView"),full=TRUE,patt="xls$")

exden.ctrl<-parseReads(exbam[1],verbose=0)
exden.chip<-parseReads(exbam[2],verbose=0)

peaks<-macs2gr(exls,psize=500)

cluster_res<-plotTV(exden.chip,exden.ctrl,regions=peaks,cluster=5,norm_readc=FALSE,showPlot=FALSE)
summaryTV(cluster_res)
tvdata<-plotTVData(cluster_res)


}
