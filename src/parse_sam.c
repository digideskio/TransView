
/*
 * parse_sam.c
 *
 *  Created on: May 29, 2012
 *      Author: Julius Muller
 */
/**
 *  @file	parse_sam.c
 *  @brief	Contains and imports all samtools related IO functions. Relies on libbam.a
 */

#include "parse_sam.h"
#include <string.h>


#define MAX_NCIGAR 1000 // guess of the maximal cigar size for preallocation

/**
* @brief Open the file handler
*
* @param name Filename
* @return void
* @details sets the global bam_file or file variable to the start of the file
* @todo nothing
*/
samfile_t * open_samtools(const char * filename) {
    const char *dot = strrchr(filename, '.'),*mode;
    if(dot && !strcmp(dot,".bam")){
    	mode="rb";
    }else if(dot && !strcmp(dot,".sam")){
    	mode="r";
    }else{
    	warning("File ending not '.sam' or '.bam'");
    	return(0);
    }
	return(samopen(filename,mode,0));
}


/**
* @brief Close the file handler
*
* @return void
* @details Closes samfile_t
* @note
*/
void close_bamfile(samfile_t *fp){
	samclose(fp);
}

/**
* @brief Updates a buffered_read with information from bamread
*
* @param bufread Buffered reads
* @param bamread Read info stored in the samtools struct bam1_t
* @param rm Results from quality_check()
* @return void
* @details Copies information from bamread and destroys bamread after passing
* @todo nothing
*/
void store_read(buffered_read_t *bufread, bam1_t *bamread,read_metrics_t *rm){

	bufread->chrom_index=bamread->core.tid;

	if(sizeof(uint32_t)*bamread->core.n_cigar>MAX_NCIGAR)bufread->cigar=(uint32_t *)Realloc(bufread->cigar,bamread->core.n_cigar,uint32_t);

	memcpy(bufread->cigar,bam1_cigar(bamread),sizeof(uint32_t)*bamread->core.n_cigar);

	bufread->l_seq=rm->read_length;
	bufread->tlen=bamread->core.isize;
	bufread->mapq=bamread->core.qual;
	bufread->n_cigar=bamread->core.n_cigar;
	bufread->pos=bamread->core.pos;// this program uses 1 based positions
	bufread->revcomp=rm->revcomp;
	bufread->proper_pair=bam1_ppair(bamread);
	bufread->written=0;
	bufread->genomic_end=rm->genomic_end;
	bam_destroy1(bamread);
}

/**
* @brief Free space occupied by buffered_read_t and bam1_t
*
* @param fr Pointer to buffered_read_t
* @param current_read pointer to a bam1_t struct
* @return void
* @details Free space occupied by buffered_read_t and bam1_t
* @todo nothing
*/
void free_samio(buffered_read_t *fr,bam1_t *current_read){
	bam_destroy1(current_read);
	if(fr->cigar)Free(fr->cigar);
	fr->written=0;
}


/**
* @brief Major quality check point for each read
*
* @param rm Empty read_metrics_t to be updated with the results
* @param temp_read The read to be assessed
* @param user_args User arguments to be considered during assessment
* @param bresults Block wide results of the parsing
* @param lpos Current position
* @return void
* @details Major quality check point for each read
* @todo nothing
*/
void quality_check(read_metrics_t *rm,bam1_t *temp_read,user_arguments_t *user_args,seq_block_t *bresults,int lpos){
	static int pos_dupcounter=0,neg_dupcounter=0;
	rm->skip=0;
	rm->read_length=0;

	rm->genomic_end= bam_calend(&temp_read->core,bam1_cigar(temp_read));

	/* Determine read length */
	if(bam1_pair(temp_read)){
		++bresults->paired;
		if (bam1_ppair(temp_read))++bresults->ppairs;
	}

	++bresults->total_reads;
	if(temp_read->core.qual < user_args->TMAPQ || bam1_unmapped(temp_read)){
		++bresults->lowqual;
		rm->skip=1;
		return;
	}

	if(user_args->UNIQUE && bam1_multimap(temp_read)){
		rm->skip=1;
		return;
	}

	if(!user_args->PAIRED){
		rm->revcomp=bam1_strand(temp_read);
		rm->read_length=bam_cigar2qlen(&temp_read->core,bam1_cigar(temp_read));
	} else if (bam1_ppair(temp_read) && !bam1_notprimary(temp_read)){
		rm->revcomp=bam1_revpair(temp_read);
		if(!user_args->READTHROUGH){rm->read_length=bam_cigar2qlen(&temp_read->core,bam1_cigar(temp_read));//sets the read length only!!
		}else if(temp_read->core.isize!=0 ){
				if((bam1_firstr(temp_read)&&!bam1_revpair(temp_read))||(bam1_secondr(temp_read)&&bam1_mrevpair(temp_read))){
					rm->read_length=temp_read->core.isize;
				} else {
					rm->skip=1;
					return;
				}
		} else{
			warning("ISIZE not set in SAM/BAM file. Re-run without using the readthrough_pairs option\n");
			rm->skip=-4;
			return;
		}
	} else{
		rm->skip=1;
		return;
	}

	if(!rm->read_length){
		rm->read_length=temp_read->core.l_qseq;
		if(!rm->read_length){
			warning("Read length neither found in core.isize=%d, core.l_qseq=%d or cigar=%d!\n",temp_read->core.isize,temp_read->core.l_qseq,bam1_cigar(temp_read));
			rm->skip=-4;
			return;
		}
	}
	/* END */

	if(user_args->STRANDED!=0){
		if((user_args->STRANDED==-1 && !rm->revcomp) || (user_args->STRANDED==1 && rm->revcomp)){
			rm->skip=1;return;
		}
	}

	if(user_args->COLLAPSE>0){
		if(lpos==temp_read->core.pos){
			if(!rm->revcomp)++pos_dupcounter;
			else ++neg_dupcounter;
			if(pos_dupcounter>=user_args->COLLAPSE || neg_dupcounter>=user_args->COLLAPSE){
				++bresults->collapsed;
				rm->skip=1;
				return;
			}
		}else{
			pos_dupcounter=0;
			neg_dupcounter=0;
		}
	}

	if(!rm->skip){
		rm->revcomp ? ++bresults->neg_strand : ++bresults->pos_strand;
		++bresults->filtered_reads;
		bresults->mapmass+=rm->read_length;
	}
}

/**
* @brief Debugging
*/
void print_readinfo(seq_block_t *bresults,bam1_t *current_read,read_metrics_t *rm,samfile_t *bam_file){
	printf("\nREADING %d\n",bresults->total_reads);
	printf("Chrom %s\n",bam_file->header->target_name[current_read->core.tid]);
	printf("Pos %d\n",current_read->core.pos);
	printf("Len %d -> END: %d\n",rm->read_length,rm->read_length+current_read->core.pos);
	printf("REVCOMP: %d\n",rm->revcomp);
	printf("SKIP: %d\n",rm->skip);
	printf("mapq: %d\n",*bam1_qual(current_read));
	printf("flag %d\n",current_read->core.flag);
	printf("n_cigar_op %d\n",current_read->core.n_cigar);
	printf("NAME %s\n",bam1_qname(current_read));
}



/**
* @brief Major file/read parsing function
*
* @param cptr Pointer to calloc initialized scores. Uses samtools API as opposed to seq_density -> faster
* @param databl_start Pointer to genomic start positions of data blocks
* @param databl_end Pointer to genomic end positions of data blocks
* @param user_args Structure holding all user input
* @param cs Pointer to chromosome_size structure holding chromosome dimensions
* @param bam_file Pointer to samfile_t structure holding the file handler
* @param first_call Is this the first call to this function?
* @return Structure seq_block with results of the read scanning such as mapmass and file status
* @details Major function looping through every read in file until a new chromosome starts. Updates cptr scores and returns statistics on reads
* @note Very long. Should be refactored?
* @todo
*/
seq_block_t seq_density(usersize *cptr, uint32_t *databl_start, uint32_t *databl_end, user_arguments_t *user_args,
		chromosome_size_t *cs, samfile_t *bam_file, int *first_call)
/*reads all reads from bam until a new chromosome is reached and sets first read infos stored in the struct fr.
 * returns a flag about the parse status and writes the scores to cptr*/
{
	int abs_gen_end;
	uint32_t nindex=0,*ind_start,*ind_end,max_index=cs->min_indexspace;
	usersize * cptr_beg=cptr;
	seq_block_t bresults={0};
	read_metrics_t rm={0};
	static buffered_read_t fr={0};
	if(*first_call)fr.cigar=(uint32_t *)Calloc(MAX_NCIGAR,uint32_t);

	int lpos=fr.pos,last_end=fr.pos+fr.l_seq;
	uint32_t BUFFERLIMIT=cs->max_pos+1;//there are no more than this many items calloc'ed
	uint32_t INDEXLIMIT=cs->min_indexspace;
	ind_start=databl_start;ind_end=databl_end;//note the start of pointers

	while(1){//heavy main loop will go through this n=[amount of reads] times!
		bam1_t * current_read=bam_init1();
		/* ######### SCAN BLOCK ################## */
		bresults.file_status=samread(bam_file,current_read);
		if (bresults.file_status == -1){// EOF
			if (*first_call){// EOF
				warning("No compatible read found for this settings!\n");
				if(!bresults.paired && user_args->PAIRED)warning("No proper pair [flag 2] found in this file. \nSet 'paired_only' to FALSE\n");
				bresults.file_status = -10;
				return bresults;
			}
			#if verbose==1
			printf("EOF detected -> %d read(s) screened!\n",bresults.total_reads);
			#endif
			cs->max_pos=last_end;//set maximal absolut position
			*(databl_end++)=cs->max_pos;//completes the missing end entry to be on par with start
			cs->min_scorespace+=*(databl_end-1)-*(databl_start-1)+1;

			cs->min_indexspace=nindex;
			*databl_end=0;*databl_start=0;//end flag
			databl_start=ind_start;databl_end=ind_end;//reset pointer to the beginning
			free_samio(&fr,current_read);
			bresults.file_status = 0;
			return bresults;
		}else if(bresults.file_status<-1){
			warning("File truncated!\n");
			if (!*first_call)free_samio(&fr,current_read);
			return bresults; // truncated
		}

		/* MAIN QUALITY CHECKPOINT */
		quality_check(&rm,current_read,user_args,&bresults,lpos);
		#if verbose==4
		print_readinfo(&bresults,current_read,&rm,bam_file);
		#endif
		if(rm.skip<0){
			bresults.file_status = rm.skip;
			return bresults;
		}else if(rm.skip)goto SKIP_READ;
		abs_gen_end=rm.genomic_end+user_args->EXTEND;
		/* ######### END SCAN ######################*/

		/*################  BOUNDARY CHECK ######################*/

		if(fr.chrom_index!=current_read->core.tid || *first_call){//Did we reach the end of the chromosome? If yes save data and return!
			if (!*first_call){
				cs->max_pos=last_end;//set maximal absolut position
				*databl_end++=cs->max_pos;//completes the missing end entry to be on par with start
				cs->min_scorespace+=*(databl_end-1)-*(databl_start-1)+1;
				cs->min_indexspace=nindex;
				*databl_end=0;*databl_start=0;//end flag
				databl_start=ind_start;databl_end=ind_end;//reset pointer to the beginning
			}
			*first_call=0;
			bresults.file_status=1;
			bresults.chrom_index_next=current_read->core.tid;
			store_read(&fr,current_read,&rm);
			return bresults;

		} else if(lpos>current_read->core.pos){
			bresults.file_status=-5;//something wrong with the positions
			warning("Last position>current position. File doesn't seem to be sorted!\n");
			free_samio(&fr,current_read);
			return bresults;
		} else if(BUFFERLIMIT<abs_gen_end || abs_gen_end < 0){
			//skip read if sequence out of bounce
			//possibly bad header with wrong chromosome margins or EXTEND too large!
			warning("BUFFER only %d\n But POS: %d cur_seq_len: %d EXTEND: %d -> %d \n GLOBAL %d\n",
					BUFFERLIMIT,current_read->core.pos,rm.read_length,user_args->EXTEND,
					current_read->core.pos+rm.read_length+user_args->EXTEND,abs_gen_end);
			#if pedantic==1
			bresults.file_status=-4;
			return bresults;
			#endif
			bam_destroy1(current_read);
			continue;
		}
		/*#####################################  WRITE  ###############*/

		cptr+=current_read->core.pos;//align pointer to current position
		if(user_args->READTHROUGH){
			write_density_ungapped(cptr,rm.read_length,&bresults.maxScore);
		}else{
			write_density_gapped(cptr,bam1_cigar(current_read),current_read->core.n_cigar,&bresults.maxScore);//minor speed panelty to ungapped
		}

		if(user_args->EXTEND>0){

			if(rm.revcomp){
				abs_gen_end-=user_args->EXTEND;
				if(cptr-user_args->EXTEND>cptr_beg)cptr-=user_args->EXTEND;
				else goto NOEXTEND;
			}else{
				cptr=cptr_beg;
				cptr+=rm.genomic_end;
			}
			int k=0;
			for(;k<user_args->EXTEND;++k){++*cptr;++cptr;}
		}
		NOEXTEND:

		cptr=cptr_beg;//jump back to the beginning of the chromosome using the helper

		if(!fr.written){//flush the first read of the chromosome stored in the struct fr
			cptr+=fr.pos;

			if(user_args->READTHROUGH){write_density_ungapped(cptr,fr.tlen,&bresults.maxScore);
			}else{write_density_gapped(cptr,fr.cigar,fr.n_cigar,&bresults.maxScore);}

			if(user_args->EXTEND>0){
				if(fr.revcomp){
					if(cptr-user_args->EXTEND>cptr_beg)cptr-=user_args->EXTEND;
					else goto FRNOEXTEND;
				}else{
					cptr=cptr_beg;
					cptr+=fr.genomic_end;
				}
				int k=0;
				for(;k<user_args->EXTEND;++k){++*cptr;++cptr;}
			}
			FRNOEXTEND:

			cptr=cptr_beg;//jump back to the beginning of the chromosome using the helper
			*databl_start++=fr.pos;++nindex;//genomic coordinate where the data block starts
			cs->min_pos=fr.pos;
			cs->min_scorespace=0;//will only be set based on the index
			cs->min_indexspace=nindex;
			cs->max_pos=abs_gen_end;
			lpos=fr.pos;
			last_end=fr.genomic_end;//load last read info of first read

			fr.written=1;//indicate that information has been used
		}

		/*################  INDEXING ######################*/

		cptr=cptr_beg;//jump back to the beginning of the chromosome using the helper
		if((current_read->core.pos-last_end)>=user_args->COMPRESSION){//check whether there was a large block without any data -> Triggers a jump on the sequence!
			nindex++;//add one index
			if(max_index<=nindex){//check whether we are already over the allocated index space
				printf("Index space found: %d > Index space allocated: %d !!",nindex,max_index);
				error("Error in indexing allocation detected!");
			}
			*databl_end++=last_end;//genomic coordinate where the data block ends | lags one index position behind start in the main loop!
			cs->min_scorespace+=*(databl_end-1)-*(databl_start-1)+1;//+1 because start=end is still one Bp!
			*databl_start++=current_read->core.pos;//genomic coordinate where the data block starts beginning with the current read
		}
		lpos=current_read->core.pos;
		last_end= user_args->READTHROUGH ? max(current_read->core.pos+rm.read_length,last_end) : max(abs_gen_end,last_end);
		SKIP_READ:
		if(nindex>=INDEXLIMIT && !first_call)error("Index overflow!\n");
		bam_destroy1(current_read);
	}//end of bam index parsing skip tag section

	bresults.file_status=-4;
	return bresults;//can never happen
}
