

#' 
#' @param dc 
#' @param ranges 
#' @param toRle 
#' @param control 
#' @param input_method 
#' @param treads_norm 
#' @returnType 
#' @return 
#' @author Julius Muller
#' @export
.sliceN<-function(dc,ranges,toRle=FALSE,control=FALSE,input_method="-",treads_norm=TRUE, nbins=0, bin_method="mean") {
	env1<-env(dc)
	if(!is.logical(control))env2<-env(control)
	if(class(ranges)[1] == "GRanges"){
		ranges<-as.data.frame(ranges)
		ranges$seqnames<-as.character(ranges$seqnames)
		ranges$strand<-as.character(ranges$strand)
	}
	if(class(ranges)[1] == "data.frame" & length(names(ranges))<3)stop("The ranges must have at least 3 columns with the structure: Chromosome-Start-End-...")
	orinames<-rownames(ranges)
	nrns<-paste("P",1:length(orinames),sep="")
	rownames(ranges)<-nrns
	if(!(class(ranges[,2])[1] %in% c("numeric","integer")) || !(class(ranges[,3])[1] %in% c("numeric","integer"))){
		mode(ranges[,2])<-"integer"
		mode(ranges[,3])<-"integer"
		if(any(is.na(ranges[,2])) || any(is.na(ranges[,3])))stop("column 2 [start] and column 3 [end] in ranges must be numeric")
	}
	ranges<-ranges[order(ranges[,1],ranges[,2]),,drop=F]
	
	if(length(nbins)==1 && is.numeric(nbins) && nbins>0){
		if(!(bin_method %in% c("mean","median","max")))stop("bin_method must be either 'mean','median' or 'max'");
		#if(length(unique(ranges[,3]-ranges[,2]))>1)stop("All ranges must have equal length if nbins is greater than 0")
	} else if(nbins!=0) {
		stop("nbins must be a positive integer");
	}	
	
	seqes<-vector("list", length(ranges[,1]))
	names(seqes)<-rownames(ranges)
	chroms<-unique(ranges[,1])
	dc_chroms<-chromosomes(dc)
	nom_chroms<-setdiff(chroms,dc_chroms)
	is_chroms<-intersect(chroms,dc_chroms)
	if(!is.logical(control)){
		if(class(control)[1]!="DensityContainer")stop("Input must be of class 'DensityContainer'")
	}
	
	for(chrom in is_chroms){
		chrl<-paste(chrom,"_lind",sep="");chrg<-paste(chrom,"_gind",sep="");
		uslices<-ranges[which(ranges[,1] == chrom),c(2,3)]
		slicenames<-rownames(uslices)
		dp<-data_pointer(dc)
		tlist<-.Call("slice_dc",env1[[dp]][[chrg]],env1[[dp]][[chrl]],env1[[dp]][[chrom]],uslices[,1],uslices[,2],nbins,bin_method,PACKAGE = "TransView")
		
		if(!is.logical(control)){#check input
			#starts<-ranges[which(ranges[,1] == chrom),2]#re set as starts might have changed in slice_dc
			dp2<-data_pointer(control)
			input_dense<-.Call("slice_dc",env2[[dp2]][[chrg]],env2[[dp2]][[chrl]],env2[[dp2]][[chrom]],uslices[,1],uslices[,2],nbins,bin_method,PACKAGE = "TransView")
			if(treads_norm){
				norm_fact<-fmapmass(dc)/fmapmass(control)#read normalization factor
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
	if(length(nom_chroms)>0)warning(sprintf("The following %d chromosome(s) were not found within the DensityContainer:\n %s",length(nom_chroms),paste(nom_chroms,collapse="|")))
	if(!is.logical(control) && treads_norm)message(sprintf("Normalization factor: %.2f",norm_fact))
	seqes<-seqes[nrns]
	names(seqes)<-orinames
	if(toRle)RleList(seqes)
	else seqes
}


#' 
#' @param dc 
#' @param ranges 
#' @returnType 
#' @return 
#' @author Julius Muller
#' @export
setMethod("sliceN", signature(dc="DensityContainer"), .sliceN)





#' 
#' @param dc 
#' @param chrom 
#' @param start 
#' @param end 
#' @param control 
#' @param input_method 
#' @param treads_norm 
#' @returnType 
#' @return 
#' @author Julius Muller
#' @export
.slice1<-function(dc,chrom,start,end,control=FALSE,input_method="-",treads_norm=TRUE, nbins=0, bin_method="mean") {
	env<-env(dc)
	if(start<1)start=1
	if(start>end)stop("end smaller than start")

	if(!(class(start)[1] %in% c("numeric","integer")) || !(class(end)[1] %in% c("numeric","integer"))){
		mode(start)<-"integer"
		mode(end)<-"integer"
		if(is.na(start) || is.na(end))stop("start and end must be numeric")
	}

	if(length(nbins)==1 && is.numeric(nbins) && nbins>0){
		if(!(bin_method %in% c("mean","median","max")))stop("bin_method must be either 'mean','median' or 'max'");
	} else if(nbins!=0) {
		stop("nbins must be a positive integer");
	}	
	
	chrl<-paste(chrom,"_lind",sep="");chrg<-paste(chrom,"_gind",sep="");
	numvec<-.Call("slice_dc",env[[data_pointer(dc)]][[chrg]],env[[data_pointer(dc)]][[chrl]],env[[data_pointer(dc)]][[chrom]],start,end,nbins,bin_method,PACKAGE = "TransView")[[1]];
	
	if(!is.logical(control)){#check input
		env<-env(control)
		if(class(control)[1]!="DensityContainer")stop("Input must be of class 'DensityContainer'")
		subname<-data_pointer(control)
		input_dense<-.Call("slice_dc",env[[subname]][[chrg]],env[[subname]][[chrl]],env[[subname]][[chrom]],start,end,nbins,bin_method,PACKAGE = "TransView")[[1]]
		
		if(all(is.na(input_dense))){
			warning(sprintf("%s was not found in the DensityContainer:\n%s",chrom))
			return(input_dense)
		}
		
		if(treads_norm){
			norm_fact<-fmapmass(dc)/fmapmass(dc)#read normalization factor
			input_dense<-input_dense*norm_fact
		}
		
		if(input_method=="-"){numvec<-round(numvec-input_dense)
		}else if(input_method=="/"){
			numvec<-log2((numvec+1)/(input_dense+1))
		}else stop("input_method must be either '+' or '/'")
	}
	return(numvec)
}


#' 
#' @param dc 
#' @returnType 
#' @return 
#' @author Julius Muller
#' @export
setMethod("slice1", signature(dc="DensityContainer",chrom="character",start="numeric",end="numeric"), .slice1)



