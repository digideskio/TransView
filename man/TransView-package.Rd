
\name{TransView-package}
\alias{TransView-package}
\alias{TransView}
\docType{package}
\title{
Read density map construction and accession. Visualization of
ChIPSeq and RNASeq data sets.
}
\description{
This package provides efficient tools to generate, access and display 
             read densities of sequencing based data sets such as from RNA-Seq and ChIP-Seq.  
}
\details{
\tabular{ll}{
Package: \tab TransView\cr
Type: \tab Package\cr
Version: \tab 1.7.4\cr
URL: \tab http://bioconductor.org/packages/release/bioc/html/TransView.html\cr
License: \tab GPL-3\cr
LazyLoad: \tab yes\cr
Depends: \tab methods,GenomicRanges\cr
Imports: \tab Rsamtools,zlibbioc,gplots,IRanges\cr
Suggests: \tab RUnit,pasillaBamSubset\cr
biocViews: \tab Bioinformatics,DNAMethylation,GeneExpression,Transcription,
Microarray,Sequencing,HighThroughputSequencing,ChIPseq,RNAseq,
Methylseq,DataImport,Visualization,Clustering,MultipleComparisons\cr
LinkingTo: \tab Rsamtools\cr
}

Index:
\preformatted{
DensityContainer-class
                        Class '"DensityContainer"'
TVResults-class         Class '"TVResults"'
TransView-package       The TransView package: Construction and
                        visualisation of read density maps.
annotatePeaks           Associates peaks to TSS
gtf2gr                  GTF file parsing
histogram-methods       Histogram of the read distribution
macs2gr                 Convenience function for MACS output conversion
parseReads              User configurable efficient assembly of read
                        density maps
peak2tss                Changes the peak center to the TSS
plotTV                  Plot and cluster global read densities
plotTVData              Summarize plotTV results
rmTV                    Free space occupied by DensityContainer
slice1                  Slice read densities from a TransView dataset
slice1T                 Slice read densities of whole transcripts from
                        a TransView DensityContainer
tvStats-methods         DensityContainer accessor function
}

Further information is available in the following vignettes:
\tabular{ll}{
\code{TransView} \tab An introduction to TransView (source, pdf)\cr
}
}
\author{
Julius Muller

Maintainer: Julius Muller <ju-mu@alumni.ethz.ch>

}

\examples{
#see vignette
}
\keyword{ package }






