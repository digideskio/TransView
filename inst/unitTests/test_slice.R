

filename<-system.file("extdata", "test.sam", package="TransView")
test_data2<-parseReads(filename, spliced=TRUE,verbose=0)

test_slice1 <- function(){
	#returns vector of read counts per base pair
	resv<-as.numeric(slice1(test_data2,"chr1",1,50))
	checkIdentical(c(1,1,3,3,3,3,3,3,2,3,3,3,1,2,2,1,0,0,0,0,0,0,1,1,1,1,1,0,0,1),resv[7:36])
}

#generate named pseudo intervals
positions<-data.frame(c(rep("chr1",3),"chr2"),c(11,31,61,21),c(21,41,71,31))
rownames(positions)<-paste(rep("interval_"),seq(1,4),sep="")

test_sliceN <- function(){
	
	#returns named list with read counts per base pair
	suppressWarnings(resl<-sliceN(test_data2,positions))
	
	checkIdentical(as.numeric(resl[["interval_1"]]), c(3, 3, 3, 3, 2, 3, 3, 3, 1, 2, 2))
	checkIdentical(as.numeric(resl[["interval_2"]]), c(1, 1, 1, 0, 0, 1, 2, 2, 2, 2, 1))
	checkIdentical(as.numeric(resl[["interval_3"]]), rep(0,11))
	checkTrue(is.null(resl[["interval_4"]]))
}


test_sliceN_filter <- function(){
	#set a filter and check the same regions	
	test_data3<-parseReads(filename,verbose=0,set_filter=data.frame(c("chr1","chr2"),c(5,21),c(200,31)))
	
	#returns named list with read counts per base pair
	suppressWarnings(resl<-sliceN(test_data3,positions))
	
	checkIdentical(as.numeric(resl[["interval_1"]]), c(3, 3, 3, 3, 2, 3, 3, 3, 1, 2, 2))
	checkIdentical(as.numeric(resl[["interval_2"]]), c(1, 1, 1, 0, 0, 1, 2, 2, 2, 2, 1))
	checkIdentical(as.numeric(resl[["interval_3"]]), rep(0,11))
	checkTrue(is.null(resl[["interval_4"]]))
}
