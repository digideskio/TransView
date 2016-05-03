

#' Constructor function
#' @param filename 
#' @param spliced 
#' @param read_stranded 
#' @param paired_only 
#' @param readthrough_pairs 
#' @param set_filter 
#' @param min_quality 
#' @param description 
#' @param extendreads 
#' @param unique_only 
#' @param max_dups 
#' @param hwindow 
#' @param compression 
#' @param verbose 
#' @returnType 
#' @return 
#' @author Julius Muller
#' @export
parseReads<-function( filename, spliced=F, read_stranded=0, paired_only=F, readthrough_pairs=F, set_filter=NA, min_quality=0, description="NA", extendreads=0, unique_only=F,	max_dups=0, hwindow=1, compression=1, verbose=1 ){
	if(hwindow>compression){hwindow=compression;warning(sprintf("hwindow has been reset to %d",compression))}
	if(compression<1 || compression>1e+6)stop("Argument 'compression' must be an integer between 1 and 1e+6")
	if(hwindow<0)stop(sprintf("Argument 'hwindow' must be a positive integer"))
	if(!is.logical(spliced))stop("Argument 'spliced' must be TRUE or FALSE")
	if(!is.logical(readthrough_pairs))stop("Argument 'readthrough_pairs' must be TRUE or FALSE")
	if(!is.logical(unique_only))stop("Argument 'unique_only' must be TRUE or FALSE")
	if(!is.logical(paired_only))stop("Argument 'paired_only' must be TRUE or FALSE")
	if(read_stranded!=0 && read_stranded!=1 && read_stranded!=-1){
		stop("Argument 'read_stranded' must be either 1 for '+', -1 for '-' or 0 for all")
	}
	if(readthrough_pairs && extendreads)stop("readthrough_pairs and extendreads cannot be simultaneously switched on!")
	if(max_dups<0 || max_dups>1000)stop("Argument 'max_dup' must be a positive integer between 0 (=off) and 1000")
	if(extendreads<0 || extendreads>10000)stop("Argument 'extendreads' must be a positive integer between 1 and 10000")
	if(min_quality<0 || min_quality>255)stop("Argument 'min_mapq' must be a positive integer between 0 and 255")
	
	
	if(verbose>1){
		message(sprintf("Filename: %s",filename))
		message(sprintf("CallArgs: %d %d %d %d %d %d %d %d %d %d %d",readthrough_pairs,paired_only,read_stranded,min_quality,max_dups,extendreads,hwindow,compression,verbose,unique_only))
	}
	
	if(class(set_filter)[1] %in% c("GRanges","data.frame")){
		if(class(set_filter)[1] == "GRanges")set_filter<-as.data.frame(set_filter)
		set_filter<-.prepare_flist(set_filter)
		if(length(set_filter[[1]])<2 && !is.numeric(set_filter))stop("The filter contains not enough columns. Structure: Chromosome-Start-End-[Strand]")
		if(length(unique(names(set_filter)))>1000)stop("The filter contains too many sequences / unique chromosome IDs (>1000).")
	} else set_filter<-0
	#1 readthrough_pairs, 2 paired, 3 strandedR, 4 min_mapq, 5 collapseR, 6 fragmentExtension, 7 hwindow, 8 compRession, 9 verbose, 10 unique_only
	#set RNA Seq specific options
	callargs<-c(readthrough_pairs,paired_only,read_stranded,min_quality,max_dups,extendreads,hwindow,compression,verbose,unique_only)
	.TransViewEnv <- new.env(hash=TRUE, parent=emptyenv())
	
	score_pointer<-paste(".TV",make.names(filename),sample(1:10e6,size=1),sep="_")#internal storage of densities
	
	assign(sprintf("%s",score_pointer),.Call("construct_dc",filename,callargs,set_filter,PACKAGE = "TransView"),envir=.TransViewEnv)
	
	if(is.null(.TransViewEnv[[sprintf("%s",score_pointer)]]))stop(paste(filename,"could not be parsed!"))
	
	seqs<-as.character(names(get(sprintf("%s",score_pointer),envir=.TransViewEnv)))
	tvhist<-as.numeric(get(sprintf("%s",score_pointer),envir=.TransViewEnv)$Histogram)
	tvhist<-tvhist[1:which.max(cumsum(tvhist!=0))]
	stats<-as.numeric(get(sprintf("%s",score_pointer),envir=.TransViewEnv)$Statistics)
	filt_status<-ifelse(length(set_filter)>1,T,F)
	
	.setTV(new("DensityContainer"),filename,seqs,stats,tvhist,callargs,spliced,description,score_pointer,filt_status,.TransViewEnv)
	
}


