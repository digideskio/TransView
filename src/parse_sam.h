
#pragma once

#include "samtools/sam.h"
#include "construct_dc.h"
#include <stddef.h>

#define bam1_ppair(b) (((b)->core.flag&BAM_FPROPER_PAIR) != 0)
#define bam1_pair(b) (((b)->core.flag&BAM_FPAIRED) != 0)
#define bam1_revpair(b) (((b)->core.flag&BAM_FREVERSE) != 0)
#define bam1_mrevpair(b) (((b)->core.flag&BAM_FMREVERSE) != 0)
#define bam1_firstr(b) (((b)->core.flag&BAM_FREAD1) != 0)
#define bam1_secondr(b) (((b)->core.flag&BAM_FREAD2) != 0)
#define bam1_unmapped(b) (((b)->core.flag&BAM_FUNMAP) != 0)
#define bam1_notprimary(b) (((b)->core.flag&BAM_FSECONDARY) != 0)
#define bam1_multimap(b) (((b)->core.flag&100) != 0)

seq_block_t seq_density(usersize *cptr,uint32_t *lstart_ind,uint32_t *lend_ind,user_arguments_t *user_args,chromosome_size_t *cs,samfile_t *bam_file,int *first_call);
samfile_t * open_samtools(const char * name);

void close_bamfile(samfile_t *fp);


/*! \var struct buffered_read_t
\brief Stores the values of a read if the next chromosome is reached during linear file parsing
Details.
*/
/** @struct buffered_read
 *  @brief Stores the values of a read if the next chromosome is reached during linear file parsing
 * 	@var buffered_read::chrom_index
 *  holds the chromosome ID from the bam header
 * 	@var buffered_read::cigar
 *  Pointer to the current CIGAR string
 * 	@var buffered_read::n_cigar
 *  Number of CIGAR operations
 *  @var buffered_read::pos
 *  holds the position in BPs 0 indexed
 *  @var buffered_read::l_seq
 *  Length of the sequence in BPs without gaps
 *  @var buffered_read::tlen
 *  Length of paired reads as stated by core.isize
 *  @var buffered_read::genomic_end
 *  End of current read in the genome
 *  @var buffered_read::revcomp
 *  0 for + strand -1 for reverse strand
 *  @var buffered_read::proper_pair
 *  whether current read is in a proper pair
 *  @var buffered_read::mapq
 *  MAPQ quality
 *  @var buffered_read::written
 *  0 for information has to be written to cptr. The struct will be ignored if written is 1
 */
typedef struct buffered_read{
	int32_t chrom_index;
	uint32_t *cigar,n_cigar;
	int32_t pos,l_seq,tlen,genomic_end;
	int8_t revcomp,proper_pair;
	int32_t mapq;
	uint8_t written;
} buffered_read_t;

/*! \var struct read_metrics_t
\brief Stores the results from the quality assessment by quality_check()
*/
/** @struct read_metrics
 *  @brief Stores the results from the quality assessment by quality_check()
 * 	@var read_metrics::revcomp
 *  Is the read on the reverse strand
 * 	@var read_metrics::skip
 *  If TRUE, the read did not pass quality assessment.
 * 	@var read_metrics::read_length
 *  read length
 *  @var read_metrics::genomic_end
 *  End of current read in the genome
 */
typedef struct read_metrics{
	uint8_t revcomp,skip;
	int read_length,genomic_end;
} read_metrics_t;

/**
* @brief Updates the cptr until l_seq and adds one count per read per position
*
* @param cptr Major pointer holding the scores for one chromosome
* @param l_seq The length of the read
* @return maximal score
* @details Updates cptr and updates also the global maximal score count.
* @todo nothing
*/
static inline void write_density_ungapped(usersize *cptr,int32_t l_seq,uint32_t *maxscore){
	int k=1;//set k to 1 to achieve half opened intervals -> 110 - 100 will be 10 Bps!
	for(;k<=l_seq;++k){//go through every bp and write
		if(*cptr<UINT16_MAX){
			++*cptr;
			*maxscore=max(*cptr,*maxscore);
		}
		++cptr;
	}
}

/**
* @brief Updates the cptr based on a CIGAR string and adds one count per read per position
*
* @param cptr Major pointer holding the scores for one chromosome
* @param cigar Pointer to CIGAR array digested by samtools
* @param ncigarops Amount of cigar operations
* @return maxscore
* @details Updates cptr and updates also the global maximal score count.
* @todo nothing
*/
static inline void write_density_gapped(usersize *cptr,uint32_t *cigar,uint16_t ncigarops,uint32_t *maxscore){
	int op,oplen,k,i=0;
	for(;i<ncigarops;++i){
		op = cigar[i] & BAM_CIGAR_MASK;
		oplen=cigar[i] >> BAM_CIGAR_SHIFT;
		switch(op){
			case BAM_CMATCH://0
			case 7://7
				for(k=1;k<=oplen;++k){//set k to 1 to achieve half opened intervals -> 110 - 100 will be 10 Bps!
					if(*cptr<UINT16_MAX){
						++*cptr;
						*maxscore=max(*cptr,*maxscore);
					}
					++cptr;
				}
				break;
			case BAM_CINS://1
			case BAM_CPAD://6
			case BAM_CSOFT_CLIP://4
			case BAM_CHARD_CLIP://5
				break;
			case BAM_CDEL://2
			case BAM_CREF_SKIP://3
			case 8:
				cptr+=oplen;
				break;
			default:warning("Illegal CIGAR operation: %d\n",op);break;
		}
	}
}


