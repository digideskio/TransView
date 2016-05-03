\name{plotTV}
\alias{plotTV}

\title{
Plot and cluster global read densities
}
\description{
Plotting facility for \code{DensityContainer}. 
}
\usage{
plotTV( ..., regions, gtf=NA, scale="global", cluster="none", control = F, peak_windows = 0, ex_windows=100,
		bin_method="mean", show_names=T, label_size=1, zero_alpha=0.5, colr=c("white","blue", "red"), 
		colr_df="redgreen",	colour_spread=c(0.05,0.05), key_limit="auto", key_limit_rna="auto", 
		set_zero="center", rowv=NA,	gclust="peaks", norm_readc=T, no_key=F, stranded_peak=T, 
		ck_size=c(2,1), remove_lowex=0, verbose=1, showPlot=T, name_width=2, pre_mRNA=F)
}

\arguments{
  \item{...}{
Depending on the combination of arguments and limited by the layout up to 20 \linkS4class{DensityContainer} and maximally one \code{matrix} can be supplied. The elements will be plotted in the order they were passed with the expression profiles and the peak profiles on the right hand and the left hand side respectively. The spliced slot determines about the kind of plot. If a \code{matrix} is provided, it will be plotted as a heatmap.
}
  \item{regions}{
\linkS4class{GRanges} object with uniformly sized regions used for plotting or character vector with IDs matching column \sQuote{transcript_id} in the GTF. 
}
  \item{gtf}{
A \linkS4class{GRanges} object with a meta data column \sQuote{transcript_id} and \sQuote{exon_id} like e.g. from \code{gtf2gr}. 
}
  \item{scale}{
A character string that determines the row scaling of the colors. Defaults to \sQuote{global} which results in a global maximum and minimum read value to be plotted across experiments. Alternative is \sQuote{individual} for individual scaling.
}
  \item{cluster}{
Sets the clustering method of the read densities. Defaults to \sQuote{none}. If an integer is passed, kmeans clustering will be performed with \code{cluster} defining the amount of clusters. A colour coded bar will be plotted to the left. For hierarchical clustering the options \sQuote{hc_sp} and \sQuote{hc_pe} for spearman or pearson correlation coefficient based distances respectively, or \sQuote{hc_rm} for distances based on row means are accepted and the results will be displayed as a dendrogram.
}
  \item{control}{
A vector of \linkS4class{DensityContainer} objects, matching the order of experiments passed as a first argument. E.g. \code{plotTV(ex1.ChIP,ex2.ChIP,ex3.RNA_KO,control=c(ex1.Input,ex2.Input,ex3.RNA_WT)}. The content will be treated as background densities and subtracted from the matching experiment.
}
  \item{show_names}{
If \code{TRUE}, peak labels and transcript IDs will be displayed on the left and the right of the plot respectively. 
}
  \item{label_size}{
Font size of the row and axis labels.
}
  \item{zero_alpha}{
Determines the alpha level of the line indicating the zero point within the peaks.
}
  \item{colr}{
A vector containing the 3 colors used for the lowest, middle and highest values respectively.  
}
  \item{colr_df}{
Determines the color in case a \code{matrix} is provided and uses \code{greenred(100)} from \pkg{gplots} by default. If changed, the arguments should be formatted analogous to \code{colr}.   
}
  \item{colour_spread}{
sets the distance of the maximum and minimum value to the saturation levels of the plot. The first value for the left side (Peak profiles) and the right for the expression plots. Can be used to adjust the contrast.
}
  \item{key_limit}{
If left at the default, the upper and lower saturation levels the peak profile colour keys will be automatically determined based on colour_spread. Can be manually overridden by a numeric vector with upper and lower levels.
}
  \item{key_limit_rna}{
If left at the default, the upper and lower saturation levels the transcript profile colour keys will be automatically determined based on colour_spread. Can be manually overridden by a numeric vector with upper and lower levels.
}
  \item{set_zero}{
if set to an integer, it determines the zero point of the x axis below the plot. E.g. a value of 250 will scale the x-axis of a 500bp peak from -250 to +250.
}
  \item{rowv}{
If a numeric vector is provided, no clustering will be performed and all rows will be ordered based on the values of this vector. Alternatively a TVResults object can be provided to reproduce previous k-means clustering. 
}  
  \item{peak_windows}{
If set to an integer greater than 0, all binding profiles will be interpolated into this amount of windows by the method specified by \code{bin_method}.
}
  \item{ex_windows}{
An integer that determines the amount of points at which the read densities of an expression experiment will get interpolated by the method specified by \code{bin_method}.
}  
  \item{bin_method}{
Specifies the function used to summarize the bins specified by nbins. Possible methods are \sQuote{max}, \sQuote{mean}, \sQuote{median} or \sQuote{approx} for linear interpolation. 
}
  \item{gclust}{
If \code{cluster} is not set to \sQuote{none}, this character string determines the cluster group. If set to \sQuote{expression} or \sQuote{peaks}, only the expression profile or peak profile data sets will be used to perform the clustering respectively. All data sets passed will be reordered based on the results of the clustering. If set to \sQuote{both}, all data sets will be treated as one matrix and clustered altogether.  
}
  \item{norm_readc}{
If set to \code{TRUE}, all sample groups will be normalized based on the map mass which is defined here as all mapped reads after quality filtering multiplied by their individual read length.
}
  \item{no_key}{
If \code{TRUE}, no color keys will be displayed.  
}
  \item{stranded_peak}{
If \code{TRUE} and strand informations are provided in \code{regions}, peak profiles will flipped if located on the negative strand.
}  
  \item{ck_size}{
  Determines the size of the colour key in the form \code{c(height,width)}
}  
  \item{remove_lowex}{
  Numeric that sets the threshold for the average read density per base pair for expression data sets. Transcripts not passing will be filtered out and a message will be displayed.
}  
  \item{verbose}{
Verbosity level
}  
  \item{showPlot}{
If \code{FALSE}, plotting will be suppressed and only the \linkS4class{TVResults} will be returned. 
}  
  \item{name_width}{
Determines the width of the space for the peak and gene names.
}  
  \item{pre_mRNA}{
All expression data will be plotted from the start of the first exon to the end of the last exon including all introns.
}  
}
\details{

Plots a false color image using the \code{image} function similar to \code{heatmap.2} of \pkg{gplots} but based on read densities. 
There are 2 different kind of plots, that can be combined or plotted individually: expression profiles and peak profiles. 
\itemize{
  \item "Peak profile plots": Peak profiles are plotted if a \linkS4class{DensityContainer} instance is supplied with the spliced slot set to \code{FALSE}. The image consists of color coded, optionally total read normalized read pileups as a stacked false color image with one peak per row. The size of the peaks is soleley relying on the genomic range passed with \code{peaks}. If strand information is available through \code{peaks}, all peaks on the reverse strand will be reversed.
  \item "Transcript profile plots": If the spliced slot of the respective \linkS4class{DensityContainer} is set to \code{TRUE}, an expression profile will be plotted. First, each expression profile will be normalized to the total amount of reads of the source BAM/SAM file and reduced to \code{ex_windows} as calculated by the \code{approx} function. The optional clustering will then be performed and subsequently all expression profiles will be scaled across rows so that each row has a mean of zero and standard deviation of one.
  \item "Heatmap": Instead of a \linkS4class{DensityContainer} with spliced set to \code{TRUE}, one \code{matrix} can be provided. The data will be scaled analogous to \sQuote{Expression profile plots} and plotted as a heatmap using the \code{image} command.
  \item "Mixed plots": If \linkS4class{DensityContainer} instances with spliced slot set to \code{TRUE} or a \code{matrix} are combined with \linkS4class{DensityContainer} with the spliced slot set to \code{FALSE}, the peak profiles will be plotted on the left and the expression plots will be plotted on the right. The \code{gclust} argument determines the clustered groups. 
}
}

\value{
Returns a \linkS4class{TVResults} class object with the results of the clustering.
}
\author{Julius Muller \email{ju-mu@alumni.ethz.ch}}
\examples{

exbam<-dir(system.file("extdata", package="TransView"),full=TRUE,patt="bam$")
exls<-dir(system.file("extdata", package="TransView"),full=TRUE,patt="xls$")

exden.ctrl<-parseReads(exbam[1],verbose=0)
exden.chip<-parseReads(exbam[2],verbose=0)

peaks<-macs2gr(exls,psize=500)

cluster_res<-plotTV(exden.chip,exden.ctrl,regions=peaks,cluster=5,norm_readc=FALSE,showPlot=FALSE)
summary(cluster_res)

}


