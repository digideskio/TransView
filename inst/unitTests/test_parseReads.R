
library("pasillaBamSubset")

exbam<-dir(system.file("extdata", package="TransView"),full=T,patt="bam$")
filename<-system.file("extdata", "test.sam", package="TransView")

test_parseReads <- function(){
	test_data<-parseReads(filename, spliced=TRUE,verbose=0,hwindow=20,compression=20)
	test_stat<-tvStats(test_data)
	checkEqualsNumeric(c(1, NA, NA), histogram(test_data)[1:3])
	checkEqualsNumeric(1.34783, round(test_stat$gcoverage,digits=5))
}


exls<-dir(system.file("extdata", package="TransView"),full=T,patt="xls$")
peaks<-macs2gr(exls,psize=500)
peak<-as.data.frame(peaks[18,])
peak$seqnames<-as.character(peak$seqnames)

test_parseVariations <- function(){
	
	exden.ctrl<-parseReads(exbam[1],verbose=0,hwindow=20,compression=20)
	exden.chip<-parseReads(exbam[2],verbose=0,hwindow=20,compression=20)
	test_stat<-tvStats(exden.chip)
	checkIdentical(c("1"=49, "2"=40, "3"=21, "4"=28, "5"=23, "6"=21, "7"=20, "8"=14, "9"= 7,"10"=4), histogram(exden.chip)[1:10])
	x1<-slice1(exden.chip,peak$seqnames,peak$start,peak$end)[300:310]
	x2<-slice1(exden.ctrl,peak$seqnames,peak$start,peak$end)[300:310]
	x3<-slice1(exden.chip,peak$seqnames,peak$start,peak$end,control=exden.ctrl,treads_norm=F)[300:310]
	checkIdentical(0.00023, round(test_stat$gcoverage,digits=5))
	checkEqualsNumeric(c(77 ,74 ,74 ,73 ,71 ,71 ,69 ,72 ,76 ,74 ,75),x1)
	checkEqualsNumeric(c( 7, 7, 7 ,7, 7, 7, 5, 5, 5, 5, 5),x2)
	checkEqualsNumeric(c(70, 67, 67, 66, 64, 64, 64, 67, 71, 69, 70),x3)
}


test_parseVariations2 <- function(){
	
	exden.ctrl<-parseReads(exbam[1],verbose=0,hwindow=20,compression=20)
	exden.chip<-parseReads(exbam[2],min_quality=10,extendreads=10,read_stranded=-1,max_dups=4,verbose=0,hwindow=20,compression=20)
	test_stat<-tvStats(exden.chip)
	checkEqualsNumeric(c(64,47,32,30,19,18,12,10,8,9), histogram(exden.chip)[1:10])
	x1<-slice1(exden.chip,peak$seqnames,peak$start,peak$end)[300:310]
	x2<-slice1(exden.ctrl,peak$seqnames,peak$start,peak$end)[300:310]
	x3<-slice1(exden.chip,peak$seqnames,peak$start,peak$end,control=exden.ctrl,treads_norm=T,input_method="/")[300:310]
	checkIdentical(0.0001112, round(test_stat$gcoverage,digits=7))
	checkEqualsNumeric(c(84,81,80,80,76,76,73,74,74,73,76),x1)
	checkEqualsNumeric(c( 7, 7, 7 ,7, 7, 7, 5, 5, 5, 5, 5),x2)
	checkIdentical(c(3.40939,3.35755,3.33985,3.33985,3.26679,3.26679,3.62449,3.64386,3.64386,3.62449,3.68182),round(x3,digits=5))
}

exgtf<-dir(system.file("extdata", package="TransView"),full=T,patt="gtf.gz$")
GTF.mm9<-gtf2gr(exgtf[2])
GTF.dm3<-gtf2gr(exgtf[1])
fn.pas_paired<-untreated1_chr4()

test_parse_T <- function(){
	exden.exprs<-parseReads(fn.pas_paired,spliced=T,min_quality=30,read_stranded=1,max_dups=5,verbose=0)
	x1<-slice1T(exden.exprs,mcols(GTF.dm3)$transcript_id[50],gtf=GTF.dm3,concatenate=T,stranded=T)[2735:2745]
	checkEqualsNumeric(c(34,35,35,36,36,36,38,38,38,38,37),x1)
	#GTF.dm3[which(names(GTF.dm3) == "NM_001014696.2")]
	x2<-sliceNT(exden.exprs,unique(mcols(GTF.dm3)$transcript_id[101:150]),gtf=GTF.dm3,concatenate=F,stranded=F)$"NM_001014696.2"  [31:40]
	checkEqualsNumeric(c( 8,  8, 11, 12, 12, 13, 13, 13, 12, 13),x2)
	
	x3<-sliceNT(exden.exprs,unique(mcols(GTF.dm3)$transcript_id[101:150]),gtf=GTF.dm3,concatenate=F,stranded=T)$"NM_001014696.2"  [31:40]
	checkEqualsNumeric(c( 26, 26, 26, 25, 25, 25, 23, 23, 23, 23),x3)
	
}

test_annotate<- function(){
	apeaks<-annotatePeaks(peaks=peaks,gtf=GTF.mm9)
	checkIdentical(mcols(apeaks)$transcript_id[5],"NM_011031")
	
	peaks2<-shift(peaks,10000)
	apeaks2<-annotatePeaks(peaks=peaks2,limit=c(-5000,2000),gtf=GTF.mm9,reference="tss",remove_unmatched=F)
	apeaks2.gb<-annotatePeaks(peaks=peaks2,limit=c(-5000,2000),gtf=GTF.mm9,reference="gene_body",remove_unmatched=F)
	checkEqualsNumeric(apeaks2.gb[which(apeaks2.gb$transcript_id=="NM_013686")]$distance,-1364)
	checkTrue(all(is.na(apeaks2$transcript_id)))
}

test_melt_peak_plotTVData<- function(){
	exden.ctrl<-parseReads(exbam[1],verbose=0,hwindow=20,compression=20)
	exden.chip<-parseReads(exbam[2],min_quality=10,extendreads=10,read_stranded=-1,max_dups=4,verbose=0,hwindow=20,compression=20)
	ex_name(exden.ctrl)<-"Test1"
	ex_name(exden.chip)<-"Test2"
	cluster_res<-plotTV(exden.chip,exden.ctrl,regions=peaks,norm_readc=FALSE,showPlot=FALSE,verbose=0,cluster="hc_pe")
	checkEquals(as.character(na.omit(unlist(parameters(cluster_res)))),c("Test2","Test1","global","hc_pe", "FALSE", "500","TRUE", "1", "0.5","white","blue" ,"red", "redgreen","0.05", "0.05","auto", "auto" , "100", "peaks","FALSE" , "FALSE" ,"TRUE","2", "1", "0","0","FALSE", "0", "0","2","FALSE" )  )
	xdf<-summaryTV(cluster_res)
	checkEqualsNumeric(c(15,20,5,18,1,8,3,12),xdf$Original[3:10])
	tvdata<-plotTVData(cluster_res)
	checkEqualsNumeric(c(17,17,18,18,19,19),round(tvdata$Average_scores[200:205]))
	apeaks<-annotatePeaks(peaks=peaks,gtf=GTF.mm9)
	peak5.df<-meltPeak(exden.chip,exden.chip,region=apeaks["Peak.5"],bin_method="mean",peak_windows=800,rpm=F)
	checkEqualsNumeric(peak5.df[which(peak5.df$Label=="Test2" & peak5.df$Position>53914536 & peak5.df$Position<53914538),]$Reads,c(29,28,27,29,28,27))
}


