
#' Getter Methods
#' @param dc 
#' @returnType 
#' @return 
#' @author Julius Muller
#' @export
setMethod("ex_name", signature("DensityContainer"), function(dc) dc@ex_name)
setMethod("origin", signature("DensityContainer"), function(dc) dc@origin)
setMethod("spliced", signature("DensityContainer"), function(dc) dc@spliced)
setMethod("paired", signature("DensityContainer"), function(dc) dc@paired)
setMethod("readthrough_pairs", signature("DensityContainer"), function(dc) dc@readthrough_pairs)
setMethod("filtered", signature("DensityContainer"), function(dc) dc@filtered)
setMethod("strands", signature("DensityContainer"), function(dc) dc@strands)
setMethod("nreads", signature("DensityContainer"), function(dc) dc@total_reads@nreads)
setMethod("gcoverage", signature("DensityContainer"), function(dc) dc@total_reads@gcoverage)
setMethod("maxScore", signature("DensityContainer"), function(dc) dc@total_reads@maxScore)
setMethod("lowqual", signature("DensityContainer"), function(dc) dc@total_reads@lowqual)
setMethod("paired_reads", signature("DensityContainer"), function(dc) dc@total_reads@paired_reads)
setMethod("proper_pairs", signature("DensityContainer"), function(dc) dc@total_reads@proper_pairs)
setMethod("collapsed", signature("DensityContainer"), function(dc) dc@total_reads@collapsed)
setMethod("compression", signature("DensityContainer"), function(dc) dc@filtered_reads@compression)
setMethod("chromosomes", signature("DensityContainer"), function(dc) dc@filtered_reads@chromosomes)
setMethod("filtered_reads", signature("DensityContainer"), function(dc) dc@filtered_reads@filtered_reads)
setMethod("pos", signature("DensityContainer"), function(dc) dc@filtered_reads@pos)
setMethod("neg", signature("DensityContainer"), function(dc) dc@filtered_reads@neg)
setMethod("lcoverage", signature("DensityContainer"), function(dc) dc@filtered_reads@lcoverage)
setMethod("lmaxScore", signature("DensityContainer"), function(dc) dc@filtered_reads@lmaxScore)
setMethod("fmapmass", signature("DensityContainer"), function(dc) dc@filtered_reads@fmapmass)
setMethod("data_pointer", signature("DensityContainer"), function(dc) dc@data_pointer)
setMethod("env", signature("DensityContainer"), function(dc) dc@env)
setMethod("size", signature("DensityContainer"), function(dc) dc@size)
setMethod("histogram", signature("DensityContainer"), function(dc) dc@histogram)
setMethod("lsize", signature("DensityContainer"), function(dc) dc@lsize)
setMethod("gsize", signature("DensityContainer"), function(dc) dc@gsize)

setMethod("parameters", signature("TVResults"), function(tvr) tvr@parameters)
setMethod("clusters", signature("TVResults"), function(tvr) tvr@ptv_order$Cluster)
setMethod("cluster_order", signature("TVResults"), function(tvr) tvr@ptv_order$NewPosition)
setMethod("scores_peaks", signature("TVResults"), function(tvr) tvr@scores_peaks)
setMethod("scores_rna", signature("TVResults"), function(tvr) tvr@scores_rna)
setMethod("summaryTV", signature("TVResults"), function(tvr) tvr@ptv_order[order(tvr@ptv_order$NewPosition),])

#' 
#' @param dc 
#' @returnType 
#' @return 
#' @author Julius Muller
#' @export
.tvStats=function(dc) {
	return(
		list(
			ex_name=ex_name(dc),
			origin=origin(dc),
			spliced=spliced(dc),
			paired=paired(dc),
			readthrough_pairs=readthrough_pairs(dc),
			filtered=filtered(dc),
			strands=strands(dc),
			nreads=nreads(dc),
			gcoverage=gcoverage(dc),
			maxScore=maxScore(dc),
			lowqual=lowqual(dc),
			paired_reads=paired_reads(dc),
			proper_pairs=proper_pairs(dc),
			collapsed=collapsed(dc),
			compression=compression(dc),
			chromosomes=chromosomes(dc),
			filtered_reads=filtered_reads(dc),
			pos=pos(dc),
			neg=neg(dc),
			lcoverage=lcoverage(dc),
			lmaxScore=lmaxScore(dc),
			fmapmass=fmapmass(dc)
		)
	)
}

#' 
#' @param dc 
#' @returnType 
#' @return 
#' @author Julius Muller
#' @export
setMethod("tvStats",
		signature(dc="DensityContainer"),
		definition=.tvStats
)



#' 
#' @param dc 
#' @returnType 
#' @return 
#' @author Julius Muller
#' @export
rmTV<-function(dc){
	assign(data_pointer(dc),1,envir=env(dc))
	rm(list=ls(all.names = TRUE,envir=env(dc)),envir=env(dc))
	invisible(gc(verbose=FALSE))
}



