
#pragma once


#include <stdint.h>
#include <R.h>
#include <Rdefines.h>

#define verbose 0 //4 prints out each value for each position within scores!!
#define MAGIC "BAM\001"
#define lMAGIC 4
#define pedantic 0 //extra check for boundary violations
#define max(a,b)  ((a < b) ?  (b) : (a))
#define printf Rprintf
#define umaxof(t) (((0x1ULL << ((sizeof(t) * 8ULL) - 1ULL)) - 1ULL) | \
                    (0xFULL << ((sizeof(t) * 8ULL) - 4ULL))) // determine maximum size of any *unsigned* type

typedef uint16_t usersize;//uint16_t with 65k max reads per position should be sufficient for sequencing

SEXP construct_dc(SEXP bamfilenameR, SEXP aRgvals, SEXP filterList);


/** @struct filter
 *  @brief Stores the intervals used to shrink the total genome size
 * 	@var filter::sequence_name
 *  Chromosome name which has to fit the chromosome name found in the source file
 * 	@var filter::seqn
 *  Amount of chromosomes that are stored in sequence_name
 *  @var filter::values
 *  Interval of type start1-end1-start2-end2-...-startn-endn
 *  @var filter::value_length
 *  Amount of intervals found in values
 */
typedef struct filter{
	char **sequence_name;
	uint32_t seqn;
	int32_t **values;
	int32_t *value_length;
}filter_t;

/** @struct chromosome_size
 *  @brief Stores basic information about the chromosome dimensions
 * 	@var chromosome_size::min_scorespace
 *  Pre allocated space needed to store all scores
 * 	@var chromosome_size::min_indexspace
 *  Predicted amount of indexes
 *  @var chromosome_size::min_pos
 *  Minimal genomic position to be screened
 *  @var chromosome_size::max_pos
 *  Maximal genomic position to be screened
 */
typedef struct chromosome_size{
	int min_scorespace;
	int min_indexspace;
	int min_pos;
	int max_pos;
}chromosome_size_t;

/*! \var struct user_arguments_t
\brief Stores the arguments passed from R through the main function to C
Details.
*/
/** @struct user_arguments
 *  @brief Stores the arguments passed from R through the main function to C
 * 	@var user_arguments::bamfilename
 *  Filename. At the moment only bam files are supported
 * 	@var user_arguments::COMPRESSION
 *  A positive integer to determine the minimum length of the sequences that should be skipped and indexed
 * 	@var user_arguments::READTHROUGH
 *  bool. read from start to end and 0 take whole read whithout CIGAR splice info
 *  @var user_arguments::PAIRED
 *  Use only proper pairs for read density map assembly
 *  @var user_arguments::EXTEND
 *  A positive integer defining the number of basepairs that should be extended in the direction of the read. only works with 'SPLICED'=0
 *  @var user_arguments::TMAPQ
 *  Threshold for MAPQ qualities
 *  @var user_arguments::FILTER
 *  if a list is provided from R in the form list("chr1"=c(100,1000,2000,2200.... start,end)) only these regions will be passed back to R after scanning.
 *  @var user_arguments::STRANDED
 *  0 for taking both strands into account. -1 takes only reverse starnds and +1 only forward strands. Relies on the bam FLAG
 *  @var user_arguments::COLLAPSE
 *  A positive integer defining the maximum amount of reads per position with the same direction.
 *  @var user_arguments::HWINDOW
 *  Window size of the histogram
 *  @var user_arguments::VERBOSE
 *  verbose levels. Basically 1 for normal and 0 for off
 *  @var user_arguments::UNIQUE
 *  Will only be used if 1 and 0x100 flag is set in the file.
 */
typedef struct user_arguments{
	const char *bamfilename;
	int COMPRESSION;
	int READTHROUGH;
	int PAIRED;
	int EXTEND;
	int TMAPQ;
	int FILTER;
	int STRANDED;
	int COLLAPSE;
	int HWINDOW;
	int VERBOSE;
	int UNIQUE;
}user_arguments_t;

/** @struct seq_block
 *  @brief Stores information returned after reading a whole chromosome block from the source file
 * 	@var seq_block::total_reads
 *  Total amount of reads found in block
 * 	@var seq_block::filtered_reads
 *  Amount of reads found in block after applying quality filters
 * 	@var seq_block::mapmass
 *  Mapmass found in block
 * 	@var seq_block::maxScore
 *  Stores maximum score found so far in file
 *  @var seq_block::lowqual
 *  Amount of reads not passing the quality threshold in block.
 *  @var seq_block::collapsed
 *  Read count of reads collapsed according to the user argument
 *  @var seq_block::ppairs
 *  Amount of proper pairs
 *  @var seq_block::paired
 *  Amount of reads in a pair
 *  @var seq_block::pos_strand
 *  Amount of reads on the forward strand
 *  @var seq_block::neg_strand
 *  Amount of reads on the reverse strand
 *  @var seq_block::file_status
 *  File status of the reader after parsing
 *  @var seq_block::chrom_index_next
 *  The index of the next chromosome as found during linear scanning
 */
typedef struct seq_block{
	long long unsigned int mapmass;
	uint32_t maxScore,total_reads,filtered_reads;
	uint32_t lowqual,collapsed,ppairs,paired,pos_strand,neg_strand;
	int32_t file_status;
	int32_t chrom_index_next;
}seq_block_t;


/** @struct global_densities
 *  @brief Stores the accumulated informations about all reads parsed
 * 	@var global_densities::list
 *  R list containing all density vectors
 * 	@var global_densities::list_names
 *  The list names of each vector
 * 	@var global_densities::mapmass
 *  Total amount of reads times read length found in all blocks
 * 	@var global_densities::lmapmass
 *  Total amount of reads times read length found in all blocks outside of gaps
 * 	@var global_densities::total_reads
 *  Total amount of reads found in all blocks
 * 	@var global_densities::gsize
 *  Total amount of base pairs covered in all blocks
 * 	@var global_densities::lsize
 *  Total amount of base pairs covered in all blocks outside of gaps
 *  @var global_densities::lowqual
 *  Amount of reads not passing the quality threshold.
 * 	@var global_densities::lmaxScore
 *  Stores maximum score in all blocks outside of gaps
 * 	@var global_densities::maxscore
 *  Stores maximum score in all blocks
 * 	@var global_densities::upcounter
 *  Amount of R elements on the stack to be freed
 * 	@var global_densities::total_elements
 *  Total elements in file header
 * 	@var global_densities::filtered_reads
 *  Amount of reads found in block after applying quality filters
 *  @var global_densities::collapsed
 *  Read count of reads collapsed according to the user argument
 *  @var global_densities::ppairs
 *  Amount of proper pairs
 *  @var global_densities::paired
 *  Amount of reads in a pair
 *  @var global_densities::pos_strand
 *  Amount of reads on the forward strand
 *  @var global_densities::neg_strand
 *  Amount of reads on the reverse strand
 *  @var global_densities::histogramp
 *  Pointer to the histogram vector
 */
typedef struct global_densities{
	SEXP list,list_names;
	long long unsigned int mapmass,lmapmass,total_reads,gsize,lsize;
	uint32_t lowqual,lmaxScore,maxscore,upcounter,total_elements;
	uint32_t filtered_reads,collapsed,ppairs,paired,pos_strand,neg_strand;
	uint32_t *histogramp;
}global_densities_t;



