

#' 
#' @param tvr
#' @returnType 
#' @return 
#' @author Julius Muller
#' @export
.plotTVData<-function(tvr) {
	params<-parameters(tvr)
	ptv_order<-summaryTV(tvr)
	argcRNA<-params[["Expression_data"]]
	argc<-params[["Peak_data"]]
	hmapc<-params[["Matrices"]]
	ttl<-params[["colnames_peaks"]]
	ttlRNA<-params[["colnames_expression"]]
	cluster<-max(unique(ptv_order$Cluster))
	kml<-data.frame();kml2<-data.frame()
	if(argcRNA && hmapc==0){
		plotmatRNA<-lapply(scores_rna(tvr),function(x){x[ptv_order$Original,]})
		usize<-ncol(plotmatRNA[[1]])
		knames<-unlist(lapply(paste("Cluster",1:cluster,sep=""),rep,usize))
		rnames<-unlist(lapply(paste("Column",1:argcRNA,sep="_"),rep,usize*cluster))
		
		kml<-data.frame(Position=rep(rep(1:usize,argcRNA),cluster),Cluster=rep(knames,argcRNA),Sample=rnames,Average_scores=rep(NA,usize*cluster*argcRNA),Plot=rep("Transcript",usize*cluster*argcRNA),stringsAsFactors=F)
		clustcoords<-lapply(1:cluster,function(x){which(ptv_order$Cluster==x)})
		for(x in 1:argcRNA){
			kml[which(kml$Sample==paste("Column",x,sep="_")),"Average_scores"]<-unlist(lapply(1:cluster,function(z){if(length(clustcoords[[z]])==1) plotmatRNA[[x]][clustcoords[[z]],] else colMeans(plotmatRNA[[x]][clustcoords[[z]],])}))
		}
		if(!any(is.na(ttlRNA)) & length(unique(ttlRNA))==length(ttlRNA)){
			for(x in 1:argcRNA)kml[which(kml$Sample==paste("Column",x,sep="_")),"Sample"]<-ttlRNA[x];
		}
	}
	if(hmapc){
		plotmatRNA<-scores_rna(tvr)[[1]][ptv_order$Original,]
		usize<-ncol(plotmatRNA)
		knames<-unlist(lapply(paste("Cluster",1:cluster,sep=""),rep,usize))
		rnames<-rep("Matrix",usize*cluster)
		clustcoords<-lapply(1:cluster,function(x){which(ptv_order$Cluster==x)})
		kscores<-unlist(lapply(1:cluster,function(z){if(length(clustcoords[[z]])==1) plotmatRNA[clustcoords[[z]],] else colMeans(plotmatRNA[clustcoords[[z]],])}))
		kml<-data.frame(Position=rep(1:usize,cluster),Cluster=knames,Sample=rnames,Average_scores=kscores,Plot=rep("Matrix",usize*cluster),stringsAsFactors=F)
	}
	if(argc){
		plotmat<-lapply(scores_peaks(tvr),function(x){x[ptv_order$Original,]}) 
		usize<-ncol(plotmat[[1]])
		mmid<-round(seq( from=-(usize/2),to=usize/2,length.out=usize))
		
		knames<-unlist(lapply(paste("Cluster",1:cluster,sep=""),rep,usize))
		rnames<-unlist(lapply(paste("Column",1:argc,sep="_"),rep,usize*cluster))
		
		kml2<-data.frame(Position=rep(rep(mmid,argc),cluster),Cluster=rep(knames,argc),Sample=rnames,Average_scores=rep(NA,usize*cluster*argc),Plot=rep("Peak",usize*cluster*argc),stringsAsFactors=F)
		clustcoords<-lapply(1:cluster,function(x){which(ptv_order$Cluster==x)})
		for(x in 1:argc){
			kml2[which(kml2$Sample==paste("Column",x,sep="_")),"Average_scores"]<-unlist(lapply(1:cluster,function(z){if(length(clustcoords[[z]])==1) plotmat[[x]][clustcoords[[z]],] else colMeans(plotmat[[x]][clustcoords[[z]],])}))
		}
		if(!any(is.na(ttl)) & length(unique(ttl))==length(ttl)){
			for(x in 1:argc)kml2[which(kml2$Sample==paste("Column",x,sep="_")),"Sample"]<-ttl[x];
		}
	}
	return(rbind(kml,kml2))
	
}


#' 
#' @param tvr
#' @returnType 
#' @return 
#' @author Julius Muller
#' @export
setMethod("plotTVData", signature(tvr="TVResults"), .plotTVData)







