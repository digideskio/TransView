
#' Stores summary of all reads that were used for constructing the read density map before eventual further filtering by a range filter
#' @slot compression Size of a gap triggering an index event
#' @slot chromosomes Character string with the chromosomes with reads used for map construction
#' @slot filtered_reads Amount of reads
#' @slot pos Reads used from the forward strand
#' @slot neg Reads used from the reverse strand
#' @slot lcoverage Local coverage which is computed by filtered map mass/covered region
#' @slot lmaxScore Maximum score of the density maps
#' @slot fmapmass total map mass after filtering
#' @slot lsize Equals the toal amount of base pairs covered excluding empty regions

#' @author Julius Muller
#' @export
setClass("FilteredReads",
		representation(
				compression="numeric",
				chromosomes="character",
				filtered_reads="numeric",
				pos="numeric",
				neg="numeric",
				lcoverage="numeric",
				lmaxScore="numeric",
				fmapmass="numeric",
				lsize="numeric"
		)
)

#' Storage of information from file parsing of all reads not necessarily used for construction of the read density map
#' @slot nreads Total number of reads
#' @slot gcoverage Total coverage computed by total map mass/(chromosome end - chromosome start). Chromosome length derived from the SAM/BAM header
#' @slot maxScore Maximum read pileup found in file
#' @slot lowqual Amount of reads that did not pass the quality score set by min_quality or were not mapped
#' @slot paired_reads Amount of reads having multiple segments in sequencing
#' @slot proper_pairs Amount of pairs with each segment properly aligned according to the aligner
#' @slot collapsed If maxDups is in place, the reads at the same position and strand exceeding this value will be counted here.
#' @slot gsize Equals to the sum of the length of all ranges from 0 to the last read per chromosome.
#' @author Julius Muller
#' @export
setClass("TotalReads",
		representation(
				nreads="numeric",
				gsize="numeric",
				gcoverage="numeric",
				maxScore="numeric",
				lowqual="numeric",
				paired_reads="numeric",
				proper_pairs="numeric",
				collapsed="numeric"

		)
)


#' Virtual S4 class storing all information about the experiment and parser settings and inheriting read information classes 
#' @slot ex_name A user provided string to define a name of this dataset
#' @slot origin Filename of the original file
#' @slot paired Does the source file contain reads with proper pairs? 
#' @slot spliced Should the class be treated like an RNA-Seq experiment for e.g. plotTV?
#' @slot readthrough_pairs If TRUE, paired reads will be connected from left to right as one long read.
#' @slot filtered Is there a range filter in place? If yes, slicing should be only conducted using the same filter!!
#' @slot strands Which strands were parsed at all. Can be "+", "-" or "both"
#' @slot total_reads TotalReads class with information about the all reads in the source file
#' @slot filtered_reads FilteredReads class storing information about reads used for read density construction
#' @author Julius Muller
#' @export
setClass("TransView",
		representation(
				ex_name="character",
				origin="character",
				spliced="logical",
				paired="logical",
				filtered="logical",
				readthrough_pairs="logical",
				strands="character",
				total_reads="TotalReads",
				filtered_reads="FilteredReads",
				"VIRTUAL"
		)
)


#' Container with the pointer of the actual density maps and a histogram
#' @slot data_pointer A private character that should not be accessed or altered directly. It points to a variable in DensityContainer:::env which is essentially a list with all results from file parsing and can be deleted with the rmTV() function.
#' @slot histogram A histogram of read pileups generated across all read density maps after filtering excluding gaps. 
#' @slot size size of the object in bytes 
#' @slot env The environment which will keep the data_pointer target 
#' @author Julius Muller
#' @export
setClass("DensityContainer",
		contains="TransView",
		representation(
				data_pointer="character",
				histogram="numeric",
				size="numeric",
				env="environment"
	)
)

#' Container with the results of a plotTV call. Includes all visual results including clusters and scores
#' @slot parameters Holds all parameters used to call plotTV
#' @slot clusters Clustering results
#' @slot cluster_order Ordering of the rows with regard to the clusters
#' @slot scores_peaks Scores of the peaks. Corresponds to the values within the plot after interpolation and normalization.
#' @slot scores_rna Scores of the transcripts. Corresponds to the values within the plot after interpolation and normalization.
#' @slot scores_matrix Scores of the matrix. Corresponds to the values within the plot after interpolation and normalization.
#' @author Julius Muller
#' @export
setClass("TVResults",
		representation(
				parameters="list",
				ptv_order="data.frame",
				scores_peaks="list",
				scores_rna="list"
		)
)

#Initialize / not necessary!
setMethod("initialize","DensityContainer",
		function(.Object,
				ex_name, compression,paired,RNA_Seq,filtered,
				origin,sequences,gcoverage,
				lcoverage,maxScore,lmaxScore,
				nreads,filtered_reads, lowqual,spliced,strands,data_pointer,
				...){
			.Object<-callNextMethod()  			
		}
)


setValidity("DensityContainer",
		function(object){
			TRUE	
		}
)
