
#' 
#' @param dc 
#' @param tnames 
#' @param gtf 
#' @param toRle 
#' @param control 
#' @param input_method 
#' @param concatenate 
#' @param stranded 
#' @param treads_norm 
#' @returnType 
#' @return 
#' @author Julius Muller
#' @export

.sliceNT<-function(dc,tnames,gtf,toRle=FALSE,control=FALSE,input_method="-",concatenate=T,stranded=T,treads_norm=T, nbins=0, bin_method="mean") {
	if(!spliced(dc))warning("dc does not contain spliced reads.\nFor ChIP-Seq slicing use sliceN.")
	if(!is.logical(control))env2<-env(control)
	env1<-env(dc)
	oritnames<-tnames
	tnames<-unique(tnames)
	
	if(length(nbins)==1 && is.numeric(nbins) && nbins>0){
		if(!(bin_method %in% c("mean","median","max")))stop("bin_method must be either 'mean','median' or 'max'");
		if(!concatenate)stop("if nbins is > 0, concatenate must be TRUE!");
	} else if(nbins!=0) {
		stop("nbins must be a positive integer");
	}	
	
	if(class(gtf)[1] == "GRanges"){
		if(!("transcript_id" %in% colnames(values(gtf))))stop("Column 'transcript_id' missing in gtf metadata")
		transc<-gtf[which(values(gtf)$transcript_id %in% tnames)]
		transc<-as.data.frame(transc)[,c("seqnames","start","end","strand","transcript_id")]
		transc$seqnames<-as.character(transc$seqnames)
		transc$transcript_id<-as.character(transc$transcript_id)
		transc$strand<-as.character(transc$strand)
	}else if(class(gtf)[1] == "data.frame"){
		if(!("transcript_id" %in% colnames(gtf)))stop("Column 'transcript_id' missing in gtf")
		transc<-gtf[which(gtf[,'transcript_id'] %in% tnames),1:5]
		
		if(!(class(transc[,2])[1] %in% c("numeric","integer")) || !(class(transc[,3])[1] %in% c("numeric","integer"))){
			mode(transc[,2])<-"integer"
			mode(transc[,3])<-"integer"
			if(any(is.na(transc[,2])) || any(is.na(transc[,3])))stop("column 2 [start] and column 3 [end] in gtf must be numeric")
		}
		
	}else{stop("gtf must be of class GRanges or data.frame")}
	
	if(length(tnames)!=length(unique(transc$transcript_id)))warning(sprintf("%d geneIDs could not be located in GTF\n",length(unique(tnames))-length(unique(transc$transcript_id))))
	if(dim(transc)[1]==0)stop("No matching identifiers found in ",tnames)
	tnames<-tnames[which(tnames %in% unique(transc$transcript_id))]
	cnames<-transc[,5]
	negnames<-rownames(transc[which(transc[,4]=="-"),])
	orinames<-rownames(transc)
	
	seqes<-vector("list", length(transc[,1]))
	names(seqes)<-orinames
	chroms<-unique(transc[,1])
	dc_chroms<-chromosomes(dc)
	nom_chroms<-setdiff(chroms,dc_chroms)
	is_chroms<-intersect(chroms,dc_chroms)
	if(min(transc$end-transc$start)<1)stop("Exons with length < 1 detected")
	for(chrom in is_chroms){
		chrl<-paste(chrom,"_lind",sep="");chrg<-paste(chrom,"_gind",sep="");
		uslices<-transc[which(transc[,1] == chrom),c(2,3)]
		uslices<-uslices[order(uslices[,1]),]
		slicenames<-rownames(uslices)
		dp<-data_pointer(dc)
		tlist<-.Call("slice_dc",env1[[dp]][[chrg]],env1[[dp]][[chrl]],env1[[dp]][[chrom]],uslices[,1],uslices[,2],0,bin_method,PACKAGE = "TransView")
		names(tlist)<-slicenames
		if(!is.logical(control)){#check input
			if(class(control)[1]!="DensityContainer")stop("Input must be of class 'DensityContainer'")
			dp2<-data_pointer(control)
			input_dense<-.Call("slice_dc",env2[[dp2]][[chrg]],env2[[dp2]][[chrl]],env2[[dp2]][[chrom]],uslices[,1],uslices[,2],0,bin_method,PACKAGE = "TransView")
			if(treads_norm){
				norm_fact<-nreads(dc)/nreads(control)#read normalization factor
				input_dense<-lapply(input_dense,"*",norm_fact)
			}
			
			if(input_method=="-"){
				tlist<-Map("-",tlist,input_dense)
				tlist<-lapply(tlist,function(x)ifelse(x<0,0,x))
			}else if(input_method=="/"){
				tlist<-Map("+",tlist,1)#add one pseudoread to avoid div by zero
				input_dense<-Map("+",input_dense,1)#add one pseudoread to avoid div by zero
				tlist<-Map(log2,Map("/",tlist,input_dense))
			} else stop("input_method must be either '+' or '/'")
		}
		seqes[slicenames]<-tlist
	}
	if(!is.logical(control))message(sprintf("Normalization factor: %.2f",norm_fact))
	seqes<-seqes[orinames]
	if(stranded){
		seqes[negnames]<-lapply(seqes[negnames],rev)
		#unegnames<-unique(sapply(strsplit(negnames,"\\."),"[",1))
		unegnames<-unique(transc[which(rownames(transc) %in% negnames),"transcript_id"])
		sortind<-1:length(orinames)
		invisible(lapply(unegnames,function(x){y<-which(transc$transcript_id==x);sortind[y]<<-sortind[rev(y)]}))
		seqes<-seqes[sortind]
	}
	gc()
	if(length(nom_chroms)>0)warning(sprintf("The following %d chromosome(s) were not found within the DensityContainer:\n %s",length(nom_chroms),paste(nom_chroms,collapse="|")))
	if(concatenate){
		b<-vector("list", length(unique(transc[,5])))
		names(b)<-unique(transc[,5])
		non_un_names<-transc[which(rownames(transc) %in% orinames),5]
		invisible(lapply(names(b),function(rn){b[[rn]]<<-unlist(seqes[which(non_un_names == rn)],use.names=F);NULL}))
		
		if(nbins)b<-.Call("approx_window",nbins,b,bin_method,PACKAGE = "TransView")
		
		if(toRle)RleList(b[oritnames])
		else b[oritnames]
	} else {
		if(toRle)RleList(seqes)
		else seqes
	}
}


#' Slices read densities from a DensityContainer dc
#' @param dc 
#' @param tnames 
#' @param gtf 
#' @returnType 
#' @return 
#' @author Julius Muller
#' @export
setMethod("sliceNT", signature(dc="DensityContainer",tnames="character"), .sliceNT)


#' Slice read densities from a DensityContainer dc
#' @param dc 
#' @param tname 
#' @param gtf 
#' @param control 
#' @param input_method 
#' @param concatenate 
#' @param stranded 
#' @param treads_norm 
#' @returnType 
#' @return 
#' @author Julius Muller
#' @export
.slice1T<-function(dc,tname,gtf,control=FALSE,input_method="-",concatenate=T,stranded=T,treads_norm=T, nbins=0, bin_method="mean") {
	if(!spliced(dc))warning("Input dc does not contain spliced reads.\nFor ChIP-Seq slicing use slice1.")
	if(!is.logical(control))env2<-env(control)
	env<-env(dc)
	
	if(length(nbins)==1 && is.numeric(nbins) && nbins>0){
		if(!(bin_method %in% c("mean","median","max")))stop("bin_method must be either 'mean','median' or 'max'");
		if(!concatenate)stop("if nbins is > 0, concatenate must be TRUE!");
	} else if(nbins!=0) {
		stop("nbins must be a positive integer");
	}	
	
	if(class(gtf)[1] == "GRanges"){
		if(!("transcript_id" %in% colnames(values(gtf))))stop("Column 'transcript_id' missing in gtf metadata")
		transc<-gtf[which(values(gtf)$transcript_id == tname)]
		transc<-as.data.frame(transc)[,c("seqnames","start","end","strand")]
	}else if(class(gtf)[1] == "data.frame"){
		if(!("transcript_id" %in% colnames(gtf)))stop("Column 'transcript_id' missing in gtf")
		transc<-gtf[which(gtf[,'transcript_id']  == tname),1:4]
		if(!(class(transc[,2])[1] %in% c("numeric","integer")) || !(class(transc[,3])[1] %in% c("numeric","integer"))){
			mode(transc[,2])<-"integer"
			mode(transc[,3])<-"integer"
			if(any(is.na(transc[,2])) || any(is.na(transc[,3])))stop("column 2 [start] and column 3 [end] in gtf must be numeric")
		}
	}else{stop("gtf must be of class GRanges or data.frame")}

	if(dim(transc)[1]==0){
		transc<-gtf[which(tolower(gtf$transcript_id) == tolower(tname)),1:4]
		if(length(dim(transc)[1])==0)stop(paste(tname,"not found in 'transcript_id' column of the GTF"))
	}
	chrom<-transc[1,1]
	strand<-transc[1,4]
	chrl<-paste(chrom,"_lind",sep="");chrg<-paste(chrom,"_gind",sep="");
	slicenames<-rownames(transc)
	starts<-transc[,"start"]
	ends<-transc[,"end"]
	
	tlist<-.Call("slice_dc",env[[data_pointer(dc)]][[chrg]],env[[data_pointer(dc)]][[chrl]],env[[data_pointer(dc)]][[chrom]],starts,ends,0,bin_method,PACKAGE = "TransView")
	
	if(all(is.na(tlist))){
		warning(sprintf("%s was not found in the DensityContainer:\n%s",chrom))
	}else if(!is.logical(control)){#check input
		subname<-data_pointer(control)
		
		input_dense<-.Call("slice_dc",env2[[subname]][[chrg]],env2[[subname]][[chrl]],env2[[subname]][[chrom]],starts,ends,0,bin_method,PACKAGE = "TransView")
		
		if(treads_norm){
			norm_fact<-nreads(dc)/nreads(control)#read normalization factor
			input_dense<-lapply(input_dense,"*",norm_fact)
		}
		
		if(input_method=="-"){
			tlist<-Map("-",tlist,input_dense)
			tlist<-lapply(tlist,function(x)ifelse(x<0,0,x))
		}else if(input_method=="/"){
			tlist<-Map("+",tlist,1)#add one pseudoread to avoid div by zero
			input_dense<-Map("+",input_dense,1)#add one pseudoread to avoid div by zero
			tlist<-Map(log2,Map("/",tlist,input_dense))
		} else stop("input_method must be either '+' or '/'")
	}
	names(tlist)<-slicenames
	if(stranded && strand=="-"){
		tlist<-lapply(tlist,rev)
		tlist<-tlist[order(names(tlist))]
	}
	if(concatenate){
		tlist<-unlist(tlist,use.names=F)
		if(nbins)tlist<-.Call("approx_window",nbins,list(tlist),bin_method,PACKAGE = "TransView")[[1]]
	}
	return(tlist)
}

#' Slice read densities from a DensityContainer dc
#' @param dc 
#' @returnType 
#' @return 
#' @author Julius Muller
#' @export
setMethod("slice1T", signature(dc="DensityContainer",tname="character"), .slice1T)


