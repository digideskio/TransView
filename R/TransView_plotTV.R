

#' Plots overview of peak locations
#' @param ... 
#' @param peaks 
#' @param gtf 
#' @param scale 
#' @param cluster 
#' @param control 
#' @param peak_windows 
#' @param show_names 
#' @param label_size 
#' @param zero_alpha 
#' @param colr 
#' @param colr_df 
#' @param colour_spread 
#' @param key_limit 
#' @param set_zero 
#' @param rowv 
#' @param ex_windows 
#' @param gclust 
#' @param norm_readc 
#' @param no_key 
#' @param stranded_peak 
#' @param ck_size 
#' @param remove_lowex 
#' @param verbose 
#' @returnType 
#' @return 
#' @author Julius Muller
#' @export
plotTV<-function ( ..., regions, gtf=NA, scale="global", cluster="none", control = F, peak_windows = 0, ex_windows=100,
			bin_method="mean", show_names=T, label_size=1, zero_alpha=0.5, colr=c("white","blue", "red"), 
			colr_df="redgreen",	colour_spread=c(0.05,0.05), key_limit="auto", key_limit_rna="auto", 
			set_zero="center", rowv=NA,	gclust="peaks", norm_readc=T, no_key=F, stranded_peak=T, 
			ck_size=c(2,1), remove_lowex=0, verbose=1, showPlot=T, name_width=2, pre_mRNA=F)
{
	
	argList<-list(...)
	ttlRNA<-c();ttl<-c()
	tcvg<-c();rcvg<-c();hmapc<-0;
	
	if(class(gtf)[1]!="GRanges"){
		if(class(gtf)[1]!="logical" || !is.na(gtf))stop("gtf must be of class 'GRanges'")
	}
	if(class(regions)[1]!="GRanges"){
		if(class(regions)[1]=="character"){
			trefs<-length(unique(regions))
			regions<-as.data.frame(gtf[which(mcols(gtf)$transcript_id %in% regions),])
			tpeaks<-length(unique(regions$transcript_id))
			if(tpeaks!=trefs){
				if(tpeaks==0){stop("No identifier in column 'transcript_id' of the gtf is matching the regions!")
				}else{warning(paste(trefs-tpeaks,"transcript_id's not found in GTF"))}
			}
			rg<-split(regions,f=regions$transcript_id)
			rg<-lapply(rg,function(x){c(as.character(x$seqnames)[1],min(x$start),max(x$end),sum(x$width),as.character(x$strand)[1],as.character(x$transcript_id)[1])})
			regions<-as.data.frame(do.call("rbind",rg),stringsAsFactors=F)
			colnames(regions)<-c("seqnames","start","end","width","strand","transcript_id")
		}else if(class(regions)[1]!="logical" || !is.na(regions))stop("regions must be of class 'GRanges' or 'character'")
	}else{
		regions<-as.data.frame(regions)
		regions$seqnames<-as.character(regions$seqnames)
		regions$strand<-as.character(regions$strand)
		tpeaks<-nrow(regions)
		if(!"transcript_id" %in% colnames(regions))regions$transcript_id<-NA
	}
	
	ptv_order<-data.frame("Original"=1:nrow(regions),"Peak"=rownames(regions),"Transcript"=regions$transcript_id,"Cluster"=rep(1,nrow(regions)),"NewPosition"=1:nrow(regions),stringsAsFactors=F)
	
	if(tpeaks<2)stop("At least 2 rows have to be present in regions")
	if(length(ex_windows)!=1 || ex_windows<1)stop("ex_windows needs to be greater than 0")
	if((cluster!="none" & cluster<2) | length(cluster)>1)stop("cluster must be a numeric > 1") 
	if(is.numeric(rowv) || class(rowv)[[1]]=="TVResults"){
		if(cluster!="none"){stop("rowv cannot be used if cluster is set")
		}else cluster<-1
	}else if(class(rowv)[[1]]!="TVResults"){if(!all(is.na(rowv)))stop("rowv has to be a numeric vector or of class TVResults")}
	if(!is.character(key_limit) & (!is.numeric(key_limit) | length(key_limit)!=2))stop("key_limit must be a numeric vector of length 2")
	if(!is.character(key_limit_rna) & (!is.numeric(key_limit_rna) | length(key_limit_rna)!=2))stop("key_limit_rna must be a numeric vector of length 2")
	if(!(gclust %in% c("expression","peaks","both")))stop("Argument gclust must be either 'expression','peaks' or 'both'")
	if(!(scale %in% c("global","individual")))stop("Argument scale must be either 'global' or 'individual'")
	if(!is.numeric(set_zero) && set_zero!="center")stop("set_zero must be numeric")
	
	for (arg in argList) {
		
		if (!.is.dc(arg)){
			if (length( dim(arg)) == 2 && is.matrix(arg)){
				if(hmapc)stop("Only one matrix allowed!")
				hmap<-arg
				hmapc<-1
				ttlRNA<-c(ttlRNA,"Matrix")
				if(nrow(hmap)!=tpeaks)stop("A matrix has to have the same amount of rows like regions")
				next
			}else stop("Data sets must be of any number of class 'DensityContainer' and maximally one matrix")
		} 
		if(spliced(arg)){
			if(class(gtf)[1]!="GRanges")stop("Expression data detected but no GTF found! Please re-run with gtf2gr")
			rcvg <- c(rcvg, ifelse(norm_readc,fmapmass(arg),1))
			if(!any(colnames(regions) %in% "transcript_id"))stop("'transcript_id' column is missing in regions. RNA-Seq can not be associated to regions.")
			regions$transcript_id<-as.character(regions$transcript_id)
			ttlRNA<-c(ttlRNA,ex_name(arg))
		}else{
			tcvg <- c(tcvg, ifelse(norm_readc,fmapmass(arg),1))
			ttl<-c(ttl,ex_name(arg))
		}
	}
	
	if (hmapc && length(rcvg)>0)stop("Expression data can be passed as one matrix or one or many DensityContainer but not both!")
	
	if(!is.logical(control) && length(control)!=(length(argList)-hmapc))stop("If control is provided, it must match the amount of experiments.")
	
	if(!is.null(tcvg))nvec <- tcvg/min(tcvg)
	if(!is.null(rcvg))rvec <- rcvg/min(rcvg)
	argc <- 0; argcRNA <- 0
	plotmat<-list();scalevec<-list();key_limits<-list();usize<-NULL
	plotmatRNA<-list()
	if(hmapc){
		plotmatRNA<-hmap
	}
	scalevecRNA<-list();key_limitsRNA<-list()
	
	if(verbose>0)message("Fetching densities...")
	
	for (arg in argList) {
		if(.is.dc(arg)){
			if (!is.logical(control)){
				if(!.is.dc(control[[argc+argcRNA+1]]))stop("Input must be of class 'DensityContainer'")
				ctrl<-control[[argc+argcRNA+1]]
			}else{ctrl=F}
			if(spliced(arg)){
				argcRNA <- argcRNA+1

				if(pre_mRNA){
					transc<-gtf[which(mcols(gtf)$transcript_id %in% regions$transcript_id)]
					transc<-as.data.frame(transc)[,c("seqnames","start","end","strand","transcript_id","exon_id")]
					mtrans<-do.call(rbind,lapply(as.character(unique(transc$transcript_id)),function(x){z<-transc[which(transc$transcript_id==x),];maxc<-z[paste0(x,".",max(z$exon_id)),]$end;minc<-z[paste0(x,".",1),]$start ;c(as.character(z[1,]$seqnames),min(minc,maxc),max(minc,maxc),as.character(z[1,]$strand))}))
					rownames(mtrans)<-as.character(unique(transc$transcript_id))
					tregions<-data.frame(chr=mtrans[,1],start=as.numeric(mtrans[,2]),end=as.numeric(mtrans[,3]),strand=mtrans[,4])
					tregions<-tregions[regions$transcript_id,]
					if(bin_method=="approx"){
						dsts <- sliceN(arg, tregions, control = ctrl,treads_norm = norm_readc)
						dsts<-.gene2window(dsts,ex_windows,window_fun="approx")
					}else{dsts <- sliceN(arg, tregions, control = ctrl,treads_norm = norm_readc,nbins=ex_windows,bin_method=bin_method)}
					if("strand" %in% colnames(regions)){#flip negative strand at the tss
						oord<-which(regions$strand == "-")
						dsts[oord]<-lapply(dsts[oord],rev)
					}
				}else{
					if(bin_method=="approx"){
						dsts <- sliceNT(arg, gtf=gtf, tnames=regions$transcript_id, control = ctrl,treads_norm = norm_readc,concatenate=T,stranded=T)
						dsts<-.gene2window(dsts,ex_windows,window_fun="approx")
					}else{dsts <- sliceNT(arg, gtf=gtf, tnames=regions$transcript_id, control = ctrl,treads_norm = norm_readc,concatenate=T,stranded=T,nbins=ex_windows,bin_method=bin_method)}
				}

				dsts<- do.call(rbind, dsts)/rvec[argcRNA]#matrix for fast plotting normalized by total reads
				plotmatRNA[[argcRNA]] <- dsts
				
			}else{
				argc <- argc + 1
				
				if(argc==1){
					usize<-unique(regions[,3]-regions[,2])+1
					if(length(usize)>1)stop("If non spliced DensityContainer are passed, all regions must have equal length")
				}
				
				if(peak_windows>0){
					if(bin_method=="approx"){
						dsts <- sliceN(arg, regions, control = ctrl,treads_norm = norm_readc)
						dsts<-.gene2window(dsts,peak_windows,window_fun="approx")
					}else{dsts<-sliceN(arg, regions, control = ctrl,treads_norm = norm_readc,nbins=peak_windows,bin_method=bin_method)}#.Call("approx_window",peak_windows,dsts,bin_method,PACKAGE = "TransView")}
				}else{
					dsts <- sliceN(arg, regions, control = ctrl,treads_norm = norm_readc)
					peak_windows<-usize
				}      

				if(stranded_peak && "strand" %in% colnames(regions)){#flip negative strand at the tss
					oord<-which(regions$strand == "-")
					dsts[oord]<-lapply(dsts[oord],rev)
				}
				
				dsts<- do.call(rbind, dsts)/nvec[argc]
				plotmat[[argc]] <- dsts
				
				scalevec[[argc]]<-as.vector(plotmat[[argc]])#for keys and scale finding
				kh<-quantile( scalevec[[argc]][scalevec[[argc]]>mean(scalevec[[argc]])[1]],prob=c(colour_spread[1],1-colour_spread[1]))
				key_limits[[argc]]<-c(floor(min(kh)),ceiling(max(kh)))
			}
		}
	}
	### ###
	
	if(argcRNA)ptv_order$Mean_dens_ex<-rowMeans(do.call("cbind",plotmatRNA))#apply(do.call("cbind",plotmatRNA),1,median)
	
	#### Remove not expressed ####
	if(remove_lowex>0 && argcRNA>0){
		
		lowex<-which(ptv_order$Mean_dens_ex<remove_lowex)
		if(length(lowex)>0){
			regions<-regions[-lowex,]
			ptv_order<-ptv_order[-lowex,]
			ptv_order$Original<-1:nrow(ptv_order)
			ptv_order$NewPosition<-1:nrow(ptv_order)
			message(sprintf("%d genes did not pass the expression threshold",length(lowex)))
			for(x in 1:argcRNA)plotmatRNA[[x]] <- plotmatRNA[[x]][-lowex,]
			if(argc>0){
				for(x in 1:argc){
					plotmat[[x]] <- plotmat[[x]][-lowex,]
					kh<-quantile( scalevec[[x]][scalevec[[x]]>mean(scalevec[[x]])[1]],prob=c(colour_spread[1],1-colour_spread[1]))
					key_limits[[argc]]<-c(floor(min(kh)),ceiling(max(kh)))
				}
			}
		}
		tpeaks<-nrow(regions)
	}
	
	### ###
	
	
	
	#### CLUSTER ####
	
	if(cluster!="none"){
		if(verbose>0)message("Clustering...")
		if(argc && !hmapc && !argcRNA){cob<-do.call("cbind",plotmat)
		}else if(argcRNA && !argc && !hmapc){cob<-do.call("cbind",plotmatRNA)                   
		}else if(hmapc && !argcRNA && !argc){cob<-plotmatRNA
		}else if(gclust=="peaks"){cob<-do.call("cbind",plotmat)
		}else if(gclust=="expression" && !hmapc){cob<-do.call("cbind",plotmatRNA)
		}else if(gclust=="expression" && hmapc){ cob<-plotmatRNA                                         
		}else if(gclust=="both" && !hmapc){ cob<-cbind(do.call("cbind",lapply(plotmat,.row_z_score)),do.call("cbind",lapply(plotmatRNA,.row_z_score)))   
		}else if(gclust=="both"){ cob<-cbind(do.call("cbind",lapply(plotmat,.row_z_score)),.row_z_score(plotmatRNA))
		}
		
		cob<-.row_z_score(cob)#do kmeans clustering only on z scores
		if(is.numeric(cluster)){
			ptv_order$Cluster<-kmeans(cob,cluster)$cluster
			ptv_order$NewPosition<-order(ptv_order$Cluster)
		}else if(substr(cluster,1,3)=="hc_"){
			if(substr(cluster,4,5)=="sp"){dend<-as.dendrogram(hclust(as.dist(1-cor(t(cob), method="spearman"))))
			}else if(substr(cluster,4,5)=="pe"){dend<-as.dendrogram(hclust(as.dist(1-cor(t(cob), method="pearson"))))
			}else if(substr(cluster,4,5)=="rm"){dend<-as.dendrogram(hclust(dist(rowMeans(cob))))
			}else stop("Clustering not implemented:",cluster) 
			ptv_order$NewPosition<-order.dendrogram(dend)
		}else stop("Clustering not implemented:",cluster) 
	}
	
	### ###
	
	if(verbose>0 && showPlot)message("Plotting...")
	
	#### Re scale expression ####
	
	if(argcRNA>0){
		cob<-do.call("cbind",plotmatRNA)
		cob<-.row_z_score(cob)
		vcob<-as.vector(cob)#vector for keys and scale finding
		xa<-quantile( vcob[abs(vcob)>mean(vcob)],prob=c(colour_spread[2],1-colour_spread[2]))
		maa<-abs(xa[which(abs(xa)==max(abs(xa)))[1]])
		for(x in 1:argcRNA){
			key_limitsRNA[[x]]<-c(floor(-maa),ceiling(maa))
			plotmatRNA[[x]] <- cob[,(1+(x-1)*ex_windows):(ex_windows+((x-1)*ex_windows))]
			scalevecRNA[[x]]<-as.vector(plotmatRNA[[x]])#for keys and scale finding
		}
	} else if (hmapc){
		plotmatRNA<-.row_z_score(plotmatRNA)
		scalevecRNA<-as.vector(plotmatRNA)#vector for keys and scale finding
		xa<-quantile( scalevecRNA[abs(scalevecRNA)>mean(scalevecRNA)],prob=c(colour_spread[2],1-colour_spread[2]))
		maa<-abs(xa[which(abs(xa)==max(abs(xa)))[1]])
		key_limitsRNA<-c(floor(-maa),ceiling(maa))
	}
	### ###
	
	
	#### LAYOUT ####
	
	if(showPlot){
		op <- par(no.readonly = TRUE)
		on.exit(par(op))
	}
	lhei <- c(ck_size[1], 8)
	
	lwid<-rep(10,argc+argcRNA+hmapc)
	lmax<-(2*(argc+argcRNA+hmapc))
	uorder<-1:(argc+argcRNA+hmapc);lorder<-(argc+argcRNA+hmapc+1):lmax
	if(show_names){
		if(argc>0){
			lwid<-c(name_width,lwid)#reserve some space for peak names in first column
			uorder<-c(0,uorder)#dont plot upper panel
			lmax<-lmax+1
			lorder<-(argc+argcRNA+1):lmax
		}
		if(argcRNA>0){
			lwid<-c(lwid,name_width)#reserve some space for gene names in last column
			uorder<-c(uorder,0)#dont plot upper panel
			lmax<-lmax+1
			lorder<-(argc+argcRNA+1):lmax
		} else if(hmapc){
			lwid<-c(lwid,name_width)#reserve some space for gene names in last column
			uorder<-c(uorder,0)#dont plot upper panel
			lmax<-lmax+1
			lorder<-(argc+hmapc+1):lmax
		}
	}
	
	if(cluster!="none"){
		uorder<-c(0,uorder)
		lmax<-lmax+1
		lorder<-(argc+argcRNA+hmapc+1):lmax
		if(substr(cluster,1,3)=="hc_"){
			lwid<-c(3,lwid)
		} else {#names are plotted last so positions have to be flipped for kmeans
			lwid<-c(1,lwid)
			if(argc>1 && show_names){
				lorder[1:2]<-rev(lorder[1:2])
				lwid[1:2]<-rev(lwid[1:2])
			}
		}
	}  
	lmat <- rbind(uorder, lorder)
	below<-rep(0,length(uorder))
	below[which(uorder>0)]<-(max(lorder)+1):(max(lorder)+length(which(uorder>0)))
	lmat<-rbind(lmat,below)
	lhei<-c(lhei,0.1)
	
	if(no_key){
		lhei<-lhei[2:3]
		mu<-max(uorder)
		lorder[lorder>0]<-lorder[lorder>0]-mu
		below[below>0]<-below[below>0]-mu
		lmat<-rbind(lorder,below)
	}
	
	if(showPlot)layout(lmat, widths = lwid, heights = lhei)
	
	
#layout.show(nf) 
	### ###
	
	
	
	#### PLOT KEYS AND SET SCALES ####
	
	if(scale=="global"){
		if(argc>0){
			gmin<-min(unlist(key_limits))
			gmax<-max(unlist(key_limits))
			for (x in 1:argc)key_limits[[x]]<-c(floor(gmin),ceiling(gmax))
		}
		if(argcRNA>0){
			gmin<-min(unlist(key_limitsRNA))
			gmax<-max(unlist(key_limitsRNA))
			for (x in 1:argcRNA)key_limitsRNA[[x]]<-c(floor(gmin),ceiling(gmax))
		}
	}
	
	shrink<-function(x,ext){x[x < ext[1]] <- ext[1];x[x > ext[2]] <- ext[2];x}
	rowMax<-function(x){apply(x,1,max)}
	if(!no_key & showPlot)par(mar = c(1.5, 3/ck_size[2], 1.5,3/ck_size[2]), cex = 0.45)#c(bottom, left, top, right)
	col<-list();breaks<-list()
	
	for (argn in 1:argc) {
		if(argc==0)break
		rmax<-min(rowMax(plotmat[[argn]]))
		
		if(!is.character(key_limit)){key_limits[[argn]]<-key_limit
		}else if(key_limits[[argn]][1]>rmax && rmax>=1 && scale!="global")key_limits[[argn]][1]<-rmax-1
		
		plotmat[[argn]]<-shrink(plotmat[[argn]],key_limits[[argn]])
		if(length(colr)<3){col[[argn]] <- colorpanel(key_limits[[argn]][2]-key_limits[[argn]][1], low=colr[1],high=colr[2])
		}else{col[[argn]] <- colorpanel(key_limits[[argn]][2]-key_limits[[argn]][1],colr[1],colr[2],colr[3])}
		if(!no_key & showPlot){
			scalevec[[argn]]<-shrink(scalevec[[argn]],key_limits[[argn]])
			breaks <- seq(key_limits[[argn]][1], key_limits[[argn]][2], length = length(col[[argn]])+1)
			image( matrix(breaks, ncol = 1), col = col[[argn]],breaks=breaks, xaxt = "n", yaxt = "n",ylim=c(0,1))
			
			hc <- hist(scalevec[[argn]], plot = F, breaks = breaks)$counts
			ktitle<-sprintf("Reads %s > %d","",key_limits[[argn]][1])
			kpeak<-if(length(hc)>1)2*max(hc[2:(length(hc)-1)]) else hc[1]
			if(hc[1]>kpeak){#correct for very skewed distributions
				ktitle<-sprintf("Reads %s > %d","",key_limits[[argn]][1]+1)
				hc[1]<-0
			}
			if(hc[length(hc)]>kpeak){#correct for very skewed distributions
				ktitle<-sprintf("Reads %s > %d",paste(key_limits[[argn]][2]-1,"> x"),key_limits[[argn]][1]+1)
				hc[length(hc)]<-0
			}
			if(ttl[argn]!="NA")title(ttl[argn])
			mtext(side = 1, ktitle, line = 1.2,cex=label_size*0.5)
			hy <- c(hc, hc[length(hc)])
			lines(seq(0,1,length.out=length(breaks)), hy/max(hy) * 0.95, lwd = 2, type = "s", col = "cyan")
			lv <- round(seq(key_limits[[argn]][1], key_limits[[argn]][2], length = 5))
			axis(1, at = seq(0,1,length.out=length(lv)), labels = lv, line = -1, tick = 0,font=2,cex.axis=label_size*0.7)
			axis(2, at = pretty(hy)/max(hy) * 0.8, pretty(hy),cex.axis=label_size*0.5)
			mtext(side = 2, "Count", line = 2,cex=label_size*0.5)#font=1
		}
	}
	colres<-100
	colRNA <- colorpanel(colres, "blue","white", "red")
	breaksRNA<-list()
	for (argn in 1:argcRNA) {
		if(argcRNA==0)break
		
		rmax<-min(rowMax(plotmatRNA[[argn]]))
		
		if(!is.character(key_limit_rna)){key_limitsRNA[[argn]]<-key_limit_rna
		}else if(key_limitsRNA[[argn]][1]>rmax && rmax>=1 && scale!="global")key_limitsRNA[[argn]][1]<-rmax-1
		plotmatRNA[[argn]]<-shrink(plotmatRNA[[argn]],key_limitsRNA[[argn]])
		if(!no_key & showPlot){
			scalevecRNA[[argn]]<-shrink(scalevecRNA[[argn]],key_limitsRNA[[argn]])
			breaks <- seq(key_limitsRNA[[argn]][1], key_limitsRNA[[argn]][2], length = colres+1)
			image( matrix(breaks, ncol = 1), col = colRNA,breaks=breaks, xaxt = "n", yaxt = "n",ylim=c(0,1))
			
			hc <- hist(scalevecRNA[[argn]], plot = F, breaks = breaks)$counts
			kpeak<-3*max(hc[2:(length(hc)-1)])
			if(ttlRNA[argn]!="NA")title(ttlRNA[argn])
			mtext(side = 1, "Z-Score", line = 1.2,cex=label_size*0.5)
			hy <- c(hc, hc[length(hc)])
			lines(seq(0,1,length.out=length(breaks)), hy/max(hy) * 0.95, lwd = 2, type = "s", col = "cyan")
			lv <- round(seq(key_limitsRNA[[argn]][1], key_limitsRNA[[argn]][2], length = 5))
			axis(1, at = seq(0,1,length.out=length(lv)), labels = lv, line = -1, tick = 0,font=2,cex.axis=label_size*0.7)
			axis(2, at = pretty(hy)/max(hy) * 0.8, pretty(hy),cex.axis=label_size*0.7)
			mtext(side = 2, "Count", line = 2,cex=label_size*0.5)#font=1
		}
	}
	
	if(hmapc){# heatmap
		if(!is.character(key_limit_rna))key_limitsRNA<-key_limit_rna
		plotmatRNA<-shrink(plotmatRNA,key_limitsRNA)
		colres<-100
		if(colr_df=="redgreen")colRNA=greenred(100)#looks better for microarrays?
		if(!no_key & showPlot){
			scalevecRNA<-shrink(scalevecRNA,key_limitsRNA)
			breaks <- seq(key_limitsRNA[1], key_limitsRNA[2], length = colres+1)
			image( matrix(breaks, ncol = 1), col = colRNA,breaks=breaks, xaxt = "n", yaxt = "n",ylim=c(0,1))
			
			hc <- hist(scalevecRNA, plot = F, breaks = breaks)$counts
			kpeak<-3*max(hc[2:(length(hc)-1)])
			
			mtext(side = 1, "Z-Score", line = 1.2,cex=label_size*0.5)
			hy <- c(hc, hc[length(hc)])
			lines(seq(0,1,length.out=length(breaks)), hy/max(hy) * 0.95, lwd = 2, type = "s", col = "cyan")
			lv <- round(seq(key_limitsRNA[1], key_limitsRNA[2], length = 5))
			axis(1, at = seq(0,1,length.out=length(lv)), labels = lv, line = -1, tick = 0,font=2,cex.axis=label_size*0.7)
			axis(2, at = pretty(hy)/max(hy) * 0.8, pretty(hy),cex.axis=label_size*0.7)
			mtext(side = 2, "Count", line = 2,cex=label_size*0.5)#font=1
		}
	}
	
	### ###
	
	
	#### PLOT CLUSTER ####
	if(showPlot)par(mar=c(0,0,1,0))#c(bottom, left, top, right)
	if(is.numeric(rowv)){
		if(length(rowv)!=tpeaks)stop("rowv has to have the same row number as the data")
		ptv_order$NewPosition<-order(rowv)
		ptv_order$Cluster<-rowv
		cluster<-length(unique(rowv))
	}else if(class(rowv)[[1]]=="TVResults"){
		ptv_order$NewPosition<-order(clusters(rowv))
		ptv_order$Cluster<-clusters(rowv)
		ptv_order$KClust_color<-ptv_order$KClust_color[order(ptv_order$Cluster)]
		cluster<-length(unique(ptv_order$Cluster))
	}
	if(cluster!="none"){
		if(is.numeric(cluster)){
			if(!"KClust_color" %in% colnames(ptv_order)){
				kcol<-rainbow(cluster)
				ptv_order$KClust_color<-NA
				for(k in 1:cluster)ptv_order$KClust_color[which(ptv_order$Cluster==k)]<-kcol[k]
			}
			if(showPlot)image(rbind(seq(1,0,-1/(tpeaks-1))),col = ptv_order$KClust_color[ptv_order$NewPosition], axes = F)
		}else if(substr(cluster,1,3)=="hc_"){
			if(showPlot)plot(dend, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none",edgePar = list(lwd=1))
		}
	}
	### ###
	
	
	#### PEAK NAMES ####
	if(argc>0 && show_names && showPlot){
		par(mar=c(0,0,1,0))#c(bottom, left, top, right)
		plot(x=c(0.5,0.5),y=c(0,tpeaks),type="n",axes=F,ylim=c(-0.5,tpeaks-0.5),yaxs="i")
		text(.5,seq((tpeaks-1),0), labels =  ptv_order$Peak[ptv_order$NewPosition], xpd = F,cex = label_size)
	}
	### ###
	
	#### PLOT MAIN FIGURE ####
	if(showPlot)par(mar = c(0, 1, 1, 1))#c(bottom, left, top, right)
	rotate = function(mat) t(mat[nrow(mat):1,,drop=FALSE])
	for (argn in 1:argc) {
		if(argc==0)break
		if(showPlot)image(rotate(plotmat[[argn]][ptv_order$NewPosition,]), col = col[[argn]], useRaster = T, axes = F)
		if(set_zero=="center")set_zero<-usize/2
		if(showPlot)lines(c(set_zero/usize,set_zero/usize),c(0,1),col=rgb(0,0,0,alpha=zero_alpha),lwd=3,lty=1)
	}
	for (argn in 1:argcRNA) {
		if(argcRNA==0)break
		if(showPlot)image(rotate(plotmatRNA[[argn]][ptv_order$NewPosition,]), col = colRNA, useRaster = T, axes = F)
	}
	if(hmapc && showPlot)image(rotate(plotmatRNA[ptv_order$NewPosition,]), col = colRNA, useRaster = T, axes = F)
	### ###
	
	#### GENE NAMES ####
	if(show_names && showPlot){
		if(argcRNA>0){
			par(mar=c(0,0,1,1))#c(bottom, left, top, right)
			plot(x=c(0.5,0.5),y=c(0,tpeaks),type="n",axes=F,ylim=c(-0.5,tpeaks-0.5),yaxs="i")
			text(y=seq((tpeaks-1),0), .5, labels = ptv_order$Transcript[ptv_order$NewPosition], xpd = F,cex = label_size)
		}else if(hmapc){
			par(mar=c(0,0,1,1))#c(bottom, left, top, right)
			plot(x=c(0.5,0.5),y=c(0,tpeaks),type="n",axes=F,ylim=c(-0.5,tpeaks-0.5),yaxs="i")
			text(y=seq((tpeaks-1),0), .5, labels = ptv_order$Transcript[ptv_order$NewPosition], xpd = F,cex = label_size)
		}
	}
	### ###
	
	#### PLOT AXIS LABELS ####
	if(showPlot){
		par(mar = c(0, 0, 0, 0))#c(bottom, left, top, right)
		
		for (argn in 1:argc) {
			if(argc==0)break
			plot.new()
			text(seq(0,1,length.out=5),0.5, labels = round(seq(-set_zero,usize-set_zero,usize/4)), cex=label_size, font=2)
		}
		for (argn in 1:argcRNA) {
			if(argcRNA==0)break
			plot.new()
			text(c(0,1),0.5, labels = c("5'","3'"), cex=label_size, font=2)
		}
		if(hmapc){
			plot.new()
			corec<-(1/ncol(plotmatRNA))/2.2
			text(seq(0+corec,1-corec,length.out=ncol(plotmatRNA)),.5, labels = colnames(plotmatRNA), cex=label_size, font=2)
		}
	}
	### ###
	
	
	
	#### Finally return the results ####
	
	params<-list("colnames_peaks"=ttl,"colnames_expression"=ttlRNA,scale=scale, cluster=cluster, control = control, peak_windows = peak_windows, show_names=show_names, label_size=label_size, 
			zero_alpha=zero_alpha, colr=colr, colr_df=colr_df,colour_spread=colour_spread, key_limit=key_limit, 
			key_limit_rna=key_limit_rna, rowv=rowv,ex_windows=ex_windows, gclust=gclust, norm_readc=norm_readc, no_key=no_key, 
			stranded_peak=stranded_peak, ck_size=ck_size, remove_lowex=remove_lowex, verbose=verbose, pre_mRNA=pre_mRNA ,"Matrices"=hmapc,"Expression_data"=argcRNA,"Peak_data"=argc,showPlot=showPlot)
	
	ptv_order$NewPosition<-order(ptv_order$NewPosition)
	if(!is.list(plotmatRNA))plotmatRNA<-list(plotmatRNA)
	invisible(.setTVResults(new("TVResults"),params,ptv_order,plotmat,plotmatRNA))
	
	### ###	
	
	
}






