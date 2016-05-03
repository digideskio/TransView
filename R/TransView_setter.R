

setTV<-function(dc,filename,seqs,stats,pileup_hist,call_args,spliced,ex_name,data_pointer,filt_status,env){
			dc@env<-env
			dc@origin<-filename 
			dc@spliced<-spliced
			dc@readthrough_pairs<-ifelse(call_args[1],TRUE,FALSE)
			if(call_args[3]==0)dc@strands<-"both"
			else if(call_args[3]==-1)dc@strands<-"-"
			else if(call_args[3]==1)dc@strands<-"+"
			dc@filtered<-filt_status
			dc@ex_name<-ex_name
			
			dc@total_reads@nreads<-round(stats[1])  
			dc@total_reads@gcoverage<-stats[2]
			dc@total_reads@maxScore<-round(stats[4])
			dc@total_reads@lowqual<-stats[6]
			dc@total_reads@collapsed<-stats[8]
			dc@total_reads@paired_reads<-stats[9]
			dc@total_reads@proper_pairs<-stats[10]
			dc@total_reads@gsize<-stats[15]
			dc@paired<-ifelse(dc@total_reads@proper_pairs,TRUE,FALSE)
			
			dc@filtered_reads@lcoverage<-stats[3]
			dc@filtered_reads@lmaxScore<-stats[5]
			dc@filtered_reads@filtered_reads<-stats[7]
			dc@filtered_reads@compression<-call_args[8]
			dc@filtered_reads@pos<-stats[11]
			dc@filtered_reads@neg<-stats[12]      
			dc@filtered_reads@chromosomes<-seqs[!sapply(seqs,function(x) x=="")]
			dc@filtered_reads@chromosomes<-dc@filtered_reads@chromosomes[-grep("_*ind|Statistics|Histogram",dc@filtered_reads@chromosomes)]
			dc@filtered_reads@fmapmass<-stats[13] 
			dc@filtered_reads@lsize<-stats[14] 
			
			dc@histogram<-pileup_hist[2:length(pileup_hist)]
			names(dc@histogram)<-1:(length(pileup_hist)-1)
			
			dc@data_pointer<-data_pointer
			dc@size<-as.numeric(object.size(get(dc@data_pointer,dc@env ))+object.size(dc))/1048576
			
			return(dc)
}

#' Fill the DensityContainer class
#' @param dc DensityContainer instance
#' @param filename The filename of the source file
#' @param seqs Chromosomes
#' @param stats List returned from construct_dc
#' @param hist A numeric vector containing a histogram of read pileups generated across all read density maps after filtering excluding gaps.
#' @param call_args The arguments partially set by the user and passed to construct_dc
#' @param spliced Should the class be treated like an RNA-Seq experiment for e.g. plotTV?
#' @param ex_name A user provided string to define a name of this dataset
#' @param data_pointer A private character that should not be accessed or altered directly. It points to a variable in .GlobalEnv which is essentially a list with all results from file parsing and can be deleted with the rm.data() function.
#' @param filt_status Is there a range filter in place? If yes, slicing should be only conducted using the same filter!!
#' @returnType DensityContainer
#' @return dc
#' @author Julius Muller
#' @export
setMethod(".setTV", signature(dc="DensityContainer",filename="character",seqs="character",stats="numeric",pileup_hist="numeric",call_args="numeric",
				spliced="logical",ex_name="character",data_pointer="character",filt_status="logical"),setTV)


setReplaceMethod("ex_name", "DensityContainer",
		function(dc, value) {
			dc@ex_name <- value
			dc
		}
) 

setReplaceMethod("spliced", "DensityContainer",
		function(dc, value) {
			dc@spliced <- value
			ifelse(is.logical(dc),dc,stop("spliced must be TRUE or FALSE"))
		}
) 




setTVResults<-function(tvr,parameters,ptv_order,scores_peaks,scores_rna){
	tvr@parameters<-parameters
	tvr@ptv_order<-ptv_order
	tvr@scores_peaks<-scores_peaks
	tvr@scores_rna<-scores_rna
	
	return(tvr)
}

#' Fill the DensityContainer class
#' @param parameters Holds all parameters used to call plotTV
#' @param ptv_order data.frame with the clustering results and the ordering
#' @param scores_peaks Scores of the peaks. Corresponds to the values within the plot after interpolation and normalization.
#' @param scores_rna Scores of the transcripts. Corresponds to the values within the plot after interpolation and normalization.
#' @returnType DensityContainer
#' @return dc
#' @author Julius Muller
#' @export
setMethod(".setTVResults", signature(tvr="TVResults",parameters="list",ptv_order="data.frame",scores_peaks="list",scores_rna="list"),setTVResults)



