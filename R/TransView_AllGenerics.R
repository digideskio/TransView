
#Accessors
setGeneric("tvStats",function(dc) standardGeneric("tvStats"))
setGeneric("histogram",function(dc) standardGeneric("histogram"))
setGeneric("ex_name", function(dc) standardGeneric("ex_name"))
setGeneric("origin", function(dc) standardGeneric("origin"))
setGeneric("spliced", function(dc) standardGeneric("spliced"))
setGeneric("paired", function(dc) standardGeneric("paired"))
setGeneric("readthrough_pairs", function(dc) standardGeneric("readthrough_pairs"))
setGeneric("filtered", function(dc) standardGeneric("filtered"))
setGeneric("strands", function(dc) standardGeneric("strands"))
setGeneric("nreads", function(dc) standardGeneric("nreads"))
setGeneric("gcoverage", function(dc) standardGeneric("gcoverage"))
setGeneric("maxScore", function(dc) standardGeneric("maxScore"))
setGeneric("lowqual", function(dc) standardGeneric("lowqual"))
setGeneric("paired_reads", function(dc) standardGeneric("paired_reads"))
setGeneric("proper_pairs", function(dc) standardGeneric("proper_pairs"))
setGeneric("collapsed", function(dc) standardGeneric("collapsed"))
setGeneric("compression", function(dc) standardGeneric("compression"))
setGeneric("chromosomes", function(dc) standardGeneric("chromosomes"))
setGeneric("filtered_reads", function(dc) standardGeneric("filtered_reads"))
setGeneric("pos", function(dc) standardGeneric("pos"))
setGeneric("neg", function(dc) standardGeneric("neg"))
setGeneric("lcoverage", function(dc) standardGeneric("lcoverage"))
setGeneric("lmaxScore", function(dc) standardGeneric("lmaxScore"))
setGeneric("fmapmass", function(dc) standardGeneric("fmapmass"))
setGeneric("data_pointer", function(dc) standardGeneric("data_pointer"))
setGeneric("env", function(dc) standardGeneric("env"))
setGeneric("size", function(dc) standardGeneric("size"))
setGeneric("lsize", function(dc) standardGeneric("lsize"))
setGeneric("gsize", function(dc) standardGeneric("gsize"))

setGeneric("parameters", function(tvr) standardGeneric("parameters"))
setGeneric("clusters", function(tvr) standardGeneric("clusters"))
setGeneric("cluster_order", function(tvr) standardGeneric("cluster_order"))
setGeneric("scores_peaks", function(tvr) standardGeneric("scores_peaks"))
setGeneric("scores_rna", function(tvr) standardGeneric("scores_rna"))
setGeneric("summaryTV", function(tvr) standardGeneric("summaryTV"))

#Setters
setGeneric("spliced<-", function(dc, value) standardGeneric("spliced<-"))
setGeneric("ex_name<-", function(dc, value) standardGeneric("ex_name<-"))

#Slice DensityContainer
setGeneric("slice1",function(dc,chrom,start,end,control=FALSE,input_method="-", treads_norm=TRUE, nbins=0, bin_method="mean") standardGeneric("slice1"))
setGeneric("sliceN",function(dc,ranges,toRle=FALSE,control=FALSE,input_method="-", treads_norm=TRUE, nbins=0, bin_method="mean") standardGeneric("sliceN"))


#Slice Transcripts from DensityContainer
setGeneric("slice1T",function(dc,tname,gtf,control=FALSE,input_method="-",concatenate=T,stranded=T, treads_norm=T, nbins=0, bin_method="mean") standardGeneric("slice1T"))
setGeneric("sliceNT",function(dc,tnames,gtf,toRle=FALSE,control=FALSE,input_method="-",concatenate=T,stranded=T, treads_norm=T, nbins=0, bin_method="mean") standardGeneric("sliceNT"))


#Populate DensityContainer after parseReads 
setGeneric(".setTV",function(dc,filename,seqs,stats,pileup_hist,call_args,spliced,ex_name,data_pointer,filt_status,env) standardGeneric(".setTV"))

#Populate TVResults after plotTV
setGeneric(".setTVResults",function(tvr,parameters,ptv_order,scores_peaks,scores_rna) standardGeneric(".setTVResults"))

#Return dTVResults with results of each cluster
setGeneric("plotTVData",function(tvr) standardGeneric("plotTVData"))


