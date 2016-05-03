
#' 
#' @param cl 
#' @returnType 
#' @return 
#' @author Julius Muller
#' @export
.is.dc<-function(cl){
	ifelse(class(cl)[1] == "DensityContainer",TRUE,FALSE)
}

#' Annotates a list of peaks with a gtf data frame
#'
#' @name sliceNT
#' @docType methods
#' @author Julius Muller
.prepare_flist<-function(df){
	nl<-list()
	seqs<-unique(df[,1])
	for(seq in seqs){
		pv<-df[which(df[,1] %in% seq),]
		ir<-IRanges(start=pv[,2]-1,end=pv[,3])#closed interval
		ir<-reduce(ir)
		nl[[seq]]<-sort(c(start(ir),end(ir)))
	}
	return(nl)
}


#' 
#' @param x 
#' @returnType 
#' @return 
#' @author Julius Muller
#' @export
.row_z_score<-function(x){
	rm <- rowMeans(x)
	x <- sweep(x, 1, rm)
	sx <- apply(x, 1, sd)+0.00001#pseudocount to avoid sd=0 -> NA
	x <- sweep(x, 1, sx, "/")
	return(x)
}


#' Determine lib directory within package source directory
#'
#' @name .lib_fname
#' @docType methods
#' @author Julius Muller
.lib_fname<-function(qlib){
	if(.Platform$OS.type=="windows"){
		if(R.Version()$arch=="x86_64")return(file.path(system.file(package=qlib),"libs","x64",paste0(qlib, .Platform$dynlib.ext)))
		else if(R.Version()$arch=="i386")return(file.path(system.file(package=qlib),"libs","i386",paste0(qlib, .Platform$dynlib.ext)))
	} else return(file.path(system.file(package=qlib),"libs",paste0(qlib, .Platform$dynlib.ext)))
}

#' Splits a list of densities into wsize windows with a unique value calculated by window_fun
#'
#' @name .gene2window
#' @docType methods
#' @author Julius Muller
.gene2window<-function(dlist,wsize,window_fun="median"){
	if(window_fun=="median"){lapply(dlist,function(d){y<-split(d,ceiling(seq_along(d)/(length(d)/wsize)));sapply(y,function(x){median(x)})}) 
	}else if(window_fun=="mean"){lapply(dlist,function(d){y<-split(d,ceiling(seq_along(d)/(length(d)/wsize)));sapply(y,function(x){mean(x)})}) 
	}else if(window_fun=="approx")sapply(lapply(dlist,approx,n=wsize), "[", "y")
}

#' 
#' @param macs_peaks_xls 
#' @param psize 
#' @param amount 
#' @param min_pileup 
#' @param log10qval 
#' @param log10pval 
#' @param fenrichment 
#' @param peak_mid 
#' @returnType 
#' @return 
#' @author Julius Muller
#' @export
macs2gr<-function(macs_peaks_xls,psize=500,amount="all",min_pileup=0,log10qval=0,log10pval=0,fenrichment=0,peak_mid="summit"){
	
	head_required_v1_41p<-c("chr","start","end","abs_summit","pileup","log10_pval","enrichment")
	head_required_v2_09p<-c("chr","start","end","abs_summit","pileup","log10_pval","enrichment","log10_qval")
	
	fh<-file(macs_peaks_xls,"r")
	lh<-0
	while (tolower(substr(readLines(fh, n=1),1,3))!="chr"){lh=lh+1}
	close(fh)
	peaks<-read.delim(macs_peaks_xls,header=T,as.is=T,skip=lh,stringsAsFactors =F)
	colnames(peaks)[grep("pvalue",colnames(peaks))]<-"log10_pval"
	colnames(peaks)[grep("qvalue",colnames(peaks))]<-"log10_qval"
	colnames(peaks)[grep("tags",colnames(peaks))]<-"pileup"
	colnames(peaks)[grep("fold_enrichment",colnames(peaks))]<-"enrichment"
	if(all(length(grep("summit",colnames(peaks)))>0,length(grep("abs_summit",colnames(peaks)))==0)){
		colnames(peaks)[grep("summit",colnames(peaks))]<-"abs_summit"
		peaks$abs_summit<-peaks$abs_summit+peaks$start
	}
	version<-0
	if(all(head_required_v2_09p %in% colnames(peaks))){
		version<-2.09;cat("Version 2.09+ detected\n")
		if(log10qval>0)peaks<-peaks[which(peaks[,"log10_qval"]>log10qval),]
	}else if(all(head_required_v1_41p %in% colnames(peaks))){
		version<-1.41;cat("Version 1.4x detected\n")
		if(all(log10qval!=0,length(grep("qvalue",colnames(peaks)))==0))stop("No qvalue present in this MACS file. Please set log10qval to 0")
	}else{stop("MACS _peaks.xls file [tested on version 1.4x or 2.09] not detected. \nFor support of untested MACS versions please contact the maintainer")}
	
	if(fenrichment>0)peaks<-peaks[which(peaks$enrichment>=fenrichment),]
	if(log10pval>0){peaks<-peaks[which(peaks$log10_pval>=log10pval),]}
	if(min_pileup>0)peaks<-peaks[which(peaks$pileup>=min_pileup),]
	
	if(psize!="preserve"&!is.numeric(psize))stop("Please provide psize or set psize to 'preserve'")
	
	cat(paste(nrow(peaks),"peaks matching\n"))
	peaks<-peaks[order(peaks[,"pileup"],decreasing=T),]
	
	if(psize=="preserve"){
		pstarts<-peaks$start
		pends<-peaks$end
		psize<-1
	}else{
		if(peak_mid=="summit"){
			pstarts<-peaks$abs_summit
			pends<-peaks$abs_summit
		} else if(peak_mid=="center"){
			pstarts<-pends<-rowMeans(peaks[,c("start","end")])
		}else stop(paste("peak_mid='",peak_mid,"' not implemented",sep=""))
	}
	
	gr<-GRanges(ranges=IRanges(start=pstarts,end=pends),strand="*",seqnames=peaks$chr)+((psize-1)/2)
	mcols(gr)<-cbind(pileup=peaks$pileup,enrichment=peaks$enrichment,log10_pval=peaks$log10_pval)
	names(gr)<-paste("Peak",1:nrow(peaks),sep=".")
	
	if(amount!="all" & amount<length(gr))gr<-gr[1:amount]
	
	gr
}

#' Annotates a list of peaks with a gtf data frame
#' @param peaks 
#' @param gtf 
#' @param limit 
#' @param remove_unmatched 
#' @param unifyBy 
#' @param unify_fun 
#' @param min_genelength 
#' @returnType 
#' @return 
#' @author Julius Muller
#' @export
annotatePeaks<-function(peaks, gtf, limit=c(-10e3,10e3), remove_unmatched=T, unifyBy=F, unify_fun="mean", min_genelength=0,reference="tss"){
	
	if(!reference%in%c("gene_body","tss"))stop("reference must be one of: gene_body or tss")
	if(class(peaks)[1]!="GRanges")stop("peaks must be of class 'GRanges'")
	if(class(gtf)[1]!="GRanges")stop("gtf must be of class 'GRanges'")
	if(!all(c("exon_id","transcript_id") %in% names(values(gtf))))stop("gtf must have columns 'exon_id' and 'transcript_id'")
	if(length(limit)==1)limit<-c(-limit,limit)
	
	if(min_genelength>0){
		tlens<-sum(width(split(gtf, gtf$transcript_id)))>=min_genelength
		gtf<-gtf[!gtf$transcript_id %in% names(tlens[which(tlens==F)])]
	}
	
	if(reference=="gene_body"){
		gn2tr<-mcols(gtf[gtf$exon_id==1,])[,c("transcript_id","gene_id")]
		gn2tr2<-gn2tr$gene_id;names(gn2tr2)<-gn2tr$transcript_id
		strans<-split(gtf,gtf$transcript_id,drop=F)
		transr<-unlist(reduce(strans,min.gapwidth=99999999))
		transr$gene_id<-gn2tr2[names(transr)]
		stopifnot(!any(duplicated(names(transr))))
		transr$tss<-NA;transr$tes<-NA;
		transr$tss[as.character(strand(transr)) == "+"]<-start(transr[strand(transr) == "+"])
		transr$tss[as.character(strand(transr)) == "-"]<-end(transr[strand(transr) == "-"])
		transr$tes[as.character(strand(transr)) == "+"]<-end(transr[strand(transr) == "+"])
		transr$tes[as.character(strand(transr)) == "-"]<-start(transr[strand(transr) == "-"])    
	}else if(reference=="tss"){
		transr<-gtf[gtf$exon_id==1,]
		end(transr[strand(transr) == "+"])<-start(transr[strand(transr) == "+"])
		start(transr[strand(transr) == "-"])<-end(transr[strand(transr) == "-"])
		transr$tss<-start(transr)  
		names(transr)<-transr$transcript_id
	}
	
	start(transr[strand(transr) == "+"])<-start(transr[strand(transr) == "+"])+limit[1]
	end(transr[strand(transr) == "+"])<-end(transr[strand(transr) == "+"])+limit[2]
	
	start(transr[strand(transr) == "-"])<-start(transr[strand(transr) == "-"])-limit[2]
	end(transr[strand(transr) == "-"])<-end(transr[strand(transr) == "-"])-limit[1]
	
	closest_ref_all<-suppressWarnings(findOverlaps(peaks,transr,type="any",select="all"))
	closest_ids<-names(transr[subjectHits(closest_ref_all)])
	gtf<-gtf[gtf$transcript_id %in% closest_ids]
	
	peaks$transcript_id<-as.character(NA);peaks$distance<-as.integer(NA)
	if("gene_id" %in% names(values(gtf)))peaks$gene_id<-as.character(NA)
	
	transr<-transr[subjectHits(closest_ref_all)]
	peakind<-unique(queryHits(closest_ref_all))
	
	if(reference=="gene_body"){
		start(transr)[as.character(strand(transr)) == "+"]<-transr[strand(transr) == "+"]$tss
		end(transr)[as.character(strand(transr)) == "+"]<-transr[strand(transr) == "+"]$tes
		
		start(transr)[as.character(strand(transr)) == "-"]<-transr[strand(transr) == "-"]$tes
		end(transr)[as.character(strand(transr)) == "-"]<-transr[strand(transr) == "-"]$tss
	}else if(reference=="tss"){
		start(transr)<-transr$tss
		end(transr)<-transr$tss
	}
	
	peaksm<-peaks[peakind]
	peakmids<-mid(ranges(peaksm))
	start(peaksm)<-peakmids
	end(peaksm)<-peakmids
	
	closest_ref<-suppressWarnings(nearest(peaksm, transr, select = "arbitrary"))
	
	if("gene_id" %in% names(values(gtf)))peaks[peakind]$gene_id<-transr[closest_ref]$gene_id
	peaks[peakind]$transcript_id<-names(transr[closest_ref])
	
	mstrands<-as.character(strand(transr[closest_ref]))
	
	if(reference=="gene_body"){
		peaks$distance_tss<-as.integer(NA);peaks$distance_tes<-as.integer(NA);
		if("+" %in% mstrands){
			pdist.tss<-peakmids[mstrands=="+"]-start(transr[closest_ref])[mstrands=="+"]
			pdist.tes<-end(transr[closest_ref])[mstrands=="+"]-peakmids[ mstrands=="+"]
			peaks[peakind]$distance[mstrands=="+"]<-sapply(1:length(pdist.tss),function(dii){mdis<-c(pdist.tss[dii],pdist.tes[dii]);mdis[mdis>0]<-0;ifelse(min(mdis)<0,min(mdis),0)})
			peaks[peakind]$distance_tss[mstrands=="+"]<-pdist.tss;peaks[peakind]$distance_tes[mstrands=="+"]<-pdist.tes;
		}
		if("-" %in% mstrands){
			ndist.tss<-end(transr[closest_ref])[mstrands=="-"]-peakmids[ mstrands=="-"]
			ndist.tes<-peakmids[ mstrands=="-"]-start(transr[closest_ref])[mstrands=="-"]
			peaks[peakind]$distance[mstrands=="-"]<-sapply(1:length(ndist.tss),function(dii){mdis<-c(ndist.tss[dii],ndist.tes[dii]);mdis[mdis>0]<-0;ifelse(min(mdis)<0,min(mdis),0)})
			peaks[peakind]$distance_tss[mstrands=="-"]<-ndist.tss;peaks[peakind]$distance_tes[mstrands=="-"]<-ndist.tes;
		}
	}else if(reference=="tss"){
		if("+" %in% mstrands)peaks[peakind]$distance[mstrands=="+"]<-peakmids[ mstrands=="+"]-start(transr[closest_ref])[mstrands=="+"]
		if("-" %in% mstrands)peaks[peakind]$distance[mstrands=="-"]<-start(transr[closest_ref])[mstrands=="-"]-peakmids[ mstrands=="-"]
	}
	
	resm<-table(!is.na(peaks$distance))
	message(sprintf("Successful annotation within %gkBps to %gkBps: %d, %d peaks without hit",limit[1]/1000,limit[2]/1000,resm["TRUE"],resm["FALSE"]))
	
	if(class(unifyBy)[1]=="DensityContainer"){
		if(!spliced(unifyBy))stop("unifyBy can only be used with RNA-Seq DensityContainer")
		
		transr$tss_group<-paste(as.character(seqnames(transr)),start(transr),sep="|")
		peak_df<-data.frame(transcript_id=peaks$transcript_id,gene_id=peaks$gene_id,stringsAsFactors=F)
		peak_df$tss_group<-as.character(NA)
		peak_df[peakind,]$tss_group<-transr[closest_ref]$tss_group
		
		amb_list<-split(transr,transr$tss_group)
		amb_list<-amb_list[elementNROWS(amb_list)>1]
		transr<-transr[which(transr$tss_group %in% names(amb_list) & transr$tss_group %in% peak_df$tss_group)]
		
		sob<-sliceNT(unifyBy,transr$transcript_id,gtf,concatenate=T,stranded=T)  
		if(is.character(unify_fun) && unify_fun=="mean"){transr$score<-sapply(sob,function(x){mean(x)})
		}else if(is.character(unify_fun) && unify_fun=="median"){transr$score<-sapply(sob,function(x){median(x)})
		}else{transr$score<-sapply(sob,function(x){unify_fun(x)})}
		
		xdf<-as.data.frame(transr,row.names=1:length(transr))[,c("tss_group","transcript_id","score")]
		xdf$transcript_id<-as.character(xdf$transcript_id)
		invisible(lapply(unique(transr$tss_group),function(gid){ndf<-xdf[xdf$tss_group==gid,];ndf<-ndf[which.max(ndf$score),,drop=T];peak_df[which(peak_df$tss_group == gid),]$transcript_id<<-ndf$transcript_id}))
		xrec<-table(peak_df$transcript_id==peaks$transcript_id)
		message(sprintf("%d transcript_id values changed based on unifyBy",xrec["FALSE"]))
		peaks$transcript_id<-peak_df$transcript_id
	}
	
	if(remove_unmatched)peaks<-peaks[!is.na(peaks$distance)]
	peaks
}


#' Converts a gtf file to a data frame with one row per exon
#' @param gtf_file 
#' @param chromosomes 
#' @param refseq_nm 
#' @returnType 
#' @return 
#' @author Julius Muller
#' @export
gtf2gr<-function(gtf_file,chromosomes=NA,refseq_nm=F, gtf_feature=c("exon"),transcript_id="transcript_id",gene_id="gene_id"){
	GTF <- read.table(gtf_file, sep="\t", header=F, stringsAsFactors=F, colClasses=c(rep("character", 3), rep("numeric", 2), rep("character", 4)))[,c(9,1,4,5,7,3)]
	colnames(GTF)<-c("transcript_id","chr","start","end","strand","type")
	GTF<-GTF[which(GTF$type %in% gtf_feature),]
	GTF<-GTF[which(GTF$end-GTF$start > 1),]
	GTF<-GTF[,-which(colnames(GTF)=="type")]
	if(nrow(GTF)==0)stop("No matching exon entries found in GTF")
	if(substr(GTF[1,2],1,3)!="chr")GTF[,2]<-paste("chr",GTF[,2],sep="")# for chromosome IDs without chr
	if(is.character(chromosomes))GTF<-GTF[which(tolower(GTF[,2]) %in% tolower(chromosomes)),]
	if(nrow(GTF)==0){cat("Chromosome names not matching GTF\n");stop()}
	GTF<-GTF[which(GTF$strand %in% c("+","-")),]
	transpos<-grep(transcript_id,strsplit(GTF[1,1],";")[[1]])
	gpos<-grep(gene_id,strsplit(GTF[1,1],";")[[1]])
	if(length(gpos)>0)GTF$gene_id<-gsub(paste("^ *",gene_id," ",sep=""),"",sapply(strsplit(GTF[,1],";"),"[",gpos))
	GTF[,1]<-gsub(paste("^ *",transcript_id," ",sep=""),"",sapply(strsplit(GTF[,1],";"),"[",transpos))
	if(refseq_nm)GTF<-GTF[which(substring(GTF[,1],1,2)=="NM"),]
	
	GTF_se<-GTF[which(GTF$strand == "+"),]
	GTF_se<-GTF_se[order(GTF_se$chr,GTF_se$transcript_id,GTF_se$start),]
	GTF_se<-cbind(GTF_se,ExonNr=ave(1:length(GTF_se$transcript_id),GTF_se$transcript_id, FUN=rank ))#rank 30% slower compared to make.unique
	
	GTF_as<-GTF[which(GTF$strand == "-"),]
	GTF_as<-GTF_as[order(GTF_as$chr,GTF_as$transcript_id,GTF_as$start,decreasing=T),]
	GTF_as<-cbind(GTF_as,ExonNr=ave(1:length(GTF_as$transcript_id),GTF_as$transcript_id, FUN=rank ))
	
	GTF<-rbind(GTF_se,GTF_as)
	GTF<-GTF[order(GTF$chr,GTF$transcript_id,GTF$start),]
	cat(paste(dim(GTF)[1],"rows matching\n"))  
	
	gr<-GRanges(seqnames=GTF$chr,ranges=IRanges(start=GTF$start,end=GTF$end),strand=GTF$strand,transcript_id=GTF$transcript_id,exon_id=GTF$ExonNr )
	gr$gene_id<-GTF$gene_id
	names(gr)<-paste(GTF$transcript_id,GTF$ExonNr,sep=".")
	
	gr  
}

#' Convert a vector of IDs matching the transcript_id column of the GTF into a TSS centred data.frame
#' @param peaks
#' @param gtf GRanges output from gtf2gr()
#' @param peak_len The desired total size of the region with the TSS located in the middle.
#' @returnType GRanges
#' @return GRanges
#' @author Julius Muller
#' @export
peak2tss<-function(peaks, gtf, peak_len=500){
	if(class(peaks)[1]!="GRanges")stop("peaks must be of class 'GRanges'")
	if(class(gtf)[1]!="GRanges")stop("gtf must be of class 'GRanges'")
	if(!all(c("exon_id","transcript_id") %in% names(values(gtf))))stop("gtf must have columns 'exon_id' and 'transcript_id'")
	if(!"transcript_id" %in% names(values(peaks)))stop("peaks must be annotated and have colum 'transcript_id'. Please run annotatePeaks first")
	ids<-values(peaks)$transcript_id
	ids_fnd<-which(ids %in% values(gtf)$transcript_id)
	if(length(ids_fnd)<length(ids)){
		warning(sprintf("%d ids not found in GTF and removed",length(ids)-length(ids_fnd)))
		peaks<-peaks[ids_fnd]
		ids<-values(peaks)$transcript_id
	}
	mingtf<-gtf[which(values(gtf)$transcript_id %in% ids),]
	names(mingtf)<-paste(values(mingtf)$transcript_id,values(mingtf)$exon_id,sep=".")
	mids<-paste(ids,"1",sep=".")
	mingtf<-as.data.frame(mingtf[which(values(mingtf)$exon_id ==1 ),])
	tss<-mingtf[which(mingtf$strand == "+"),c("seqnames","start","strand"),drop=F]
	names(tss)<-c("chr","TSS","strand")
	ntss<-mingtf[which(mingtf$strand == "-"),c("seqnames","end","strand"),drop=F]
	names(ntss)<-c("chr","TSS","strand")
	tss<-rbind(tss,ntss)
	tss<-tss[mids,]
	ranges(peaks)<-IRanges(start=tss$TSS-(peak_len/2),end=tss$TSS+(peak_len/2),names=names(peaks))
	peaks
}

#' Convenience function meltPeaks(), which returns a data frame with normalized peak densities suitable for plotting with ggplot2 
#' @param ... 
#' @param region
#' @param control 
#' @param peak_windows 
#' @param bin_method
#' @param rpm
#' @param smooth
#' @returnType data.frame
#' @return 
#' @author Julius Muller
#' @export
meltPeak<-function (..., region ,control=FALSE, peak_windows = 0, bin_method="mean", rpm=TRUE, smooth=0)
{
	argList<-list(...)
	
	if(class(region)[1]!="GRanges" || !all(c("transcript_id","distance") %in% names(mcols(region)))){stop("region must be the annotated output of annotatePeaks or a GRanges object with transcript_id and distance metadata columns")
	}
	if(!is.logical(control) && length(control)!=length(argList))stop("If control is provided, it must match the amount of experiments.")
	if(length(region)!=1)stop("meltPeak cann only process one peak at a time")
	ttl<-sapply(argList,ex_name)
	
	if(is.numeric(peak_windows) && peak_windows>=0){usize<-ifelse(peak_windows==0,width(region),peak_windows)
	}else(stop("peak_windows must be a positive integer"))
	
	argc <- 0;mpops<-approx(start(region):end(region),n=usize)$y
	if(rpm){plotdf<-data.frame(RPM=c(rep(rep(0,usize),length(argList))),Position=c(rep(mpops,length(argList))),Label=unlist(lapply(ttl,rep,usize)),stringsAsFactors=F)
	}else{plotdf<-data.frame(Reads=c(rep(rep(0,usize),length(argList))),Position=c(rep(mpops,length(argList))),Label=unlist(lapply(ttl,rep,usize)),stringsAsFactors=F)}
	rlab<-colnames(plotdf)[1]
	if(smooth)plotdf$Smooth<-0
	for (arg in argList) {
		if (!.is.dc(arg)) stop("Data sets must be of any number of class 'DensityContainer'")
		argc <- argc + 1
		if (!is.logical(control)){
			if(!.is.dc(control[[argc]]))stop("Input must be of class 'DensityContainer'")
			ctrl<-control[[argc]]
		}else{ctrl<-F}
		
		if(peak_windows>0){
			if(bin_method=="approx"){
				dsts <- slice1(arg, chrom=as.character(seqnames(region)), start=start(region), end=end(region), control = ctrl, treads_norm = T)
				dsts<-approx(dsts,n=peak_windows)$y
			}else{dsts<-slice1(arg, chrom=as.character(seqnames(region)), start=start(region), end=end(region), control=ctrl, treads_norm=T, nbins=peak_windows, bin_method=bin_method)}
		}else{dsts <- slice1(arg, chrom=as.character(seqnames(region)), start=start(region), end=end(region), control = ctrl, treads_norm = T)}      
		plotdf[(1+(argc-1)*usize):(argc*usize),rlab] <- if(rpm) dsts/(filtered_reads(arg)/10^6) else dsts
		plotdf[(1+(argc-1)*usize):(argc*usize),"Label"] <- ex_name(arg)
		if(smooth>0){
			mloe <- lowess(plotdf[(1+(argc-1)*usize):(argc*usize),"Position"],plotdf[(1+(argc-1)*usize):(argc*usize),rlab], f =smooth)
			plotdf[(1+(argc-1)*usize):(argc*usize),"Smooth"] <- mloe$y
			plotdf[(1+(argc-1)*usize):(argc*usize),"Position"] <- mloe$x
		}
	}
	return(plotdf)
}

