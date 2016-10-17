

#### DensityContainer ####

#' prints out dc information
#'
#' @name show
#' @docType methods
#'
#'
setMethod("show", 
		signature(object="DensityContainer"), 
		definition=function(object) {
			cat("class: DensityContainer")
			cat("\n  Experiment:",ex_name(object))
			cat("\n  Source:",origin(object))
			cat("\n  Spliced:",spliced(object))
			cat("\n  Paired:",paired(object))
			cat("\n  Filtered:",filtered(object))
			cat("\n  Reads in file:",nreads(object))
			cat("\n  Reads used:",filtered_reads(object))
			cat("\n  Coverage:",gcoverage(object))    
			cat("\n  Local Coverage:",lcoverage(object))
			cat("\n  Max Score:",maxScore(object))
			cat("\n  Local Max Score:",lmaxScore(object))
			cat("\n  Low Quality / Unmapped:",lowqual(object))
			cat("\n  Strands:",strands(object))
			cat("\n  Memory usage [MB]:",round(size(object),digits=2))
			cat("\n Available Slots:\n")
			cat(slotNames(object),"\n")
			cat(slotNames(object@filtered_reads),"\n")
			cat(slotNames(object@total_reads),"\n")
			cat("Chromosomes:",paste(chromosomes(object),collapse="|"))
			cat("\n")
		}
)


#' prints out tvr information
#'
#' @name show
#' @docType methods
#'
#'
setMethod("show", 
		signature(object="TVResults"), 
		definition=function(object) {
			cat("class: TVResults\n")
			print(format(parameters(object)))
			cat("\n")
		}
)
