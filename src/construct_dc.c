/*
 * construct_dc.c
 *
 *  Created on: May 29, 2012
 *      Author: Julius Muller
 */

/**
 *  @file	construct_dc.c
 *  @brief	Main entry point for R by calling dyn.load("construct_dc.so");.Call("construct_dc",...);
 */


#include <R.h>
#include <Rdefines.h>
#include <signal.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "construct_dc.h"
#include "parse_sam.h"
#include "visuals.h"

/**
* @brief Returns end-start of a member of ft
*
* @param filter_index Index in ft
* @param ft filter struct
* @return end - start of the filter
* @todo nothing
*/
uint32_t filter_coverage(int filter_index, filter_t * ft){//calculates the bps covered by an overlapping sorted filter
	int fcount=0,fcvg=0,totfs;
	totfs=ft->value_length[filter_index];
	for(;fcount<totfs;fcount+=2){
		fcvg+=ft->values[filter_index][fcount+1]-ft->values[filter_index][fcount]+1;
	}
	return(fcvg);
}

/**
* @brief Sets optionally provided filter which is passed from R as a list
*
* @param filterList passed from R in the form list("chr1"=c(100,200,3000,3010...start,end...))
* @param ft pointer to struct filter
* @return upc -> unprotect counter to remember SEXPs to unprotect
* @details This function processes the filter passed from R and makes it compatible to work on the densities
* @note Filters are used after the density has been copied
* @todo nothing.
*/
int set_filter(SEXP filterList, filter_t * ft){

	SEXP fnames = getAttrib(filterList, R_NamesSymbol);
	SEXP fvals;
	int flen=LENGTH(fnames),upc=0;
	ft->sequence_name = (char**)Calloc(flen,char*);
	ft->values = (int32_t **)Calloc(flen,int32_t*);
	ft->value_length = (int32_t *)Calloc(flen,int32_t);
	PROTECT(fnames = AS_CHARACTER(fnames));upc++;
	int j=0;
	for(;j<flen;j++){
		PROTECT(fvals=AS_INTEGER(VECTOR_ELT(filterList, j)));upc++;
		ft->sequence_name[j] = Calloc(strlen(CHAR(STRING_ELT(fnames, j))), char);
		strcpy(ft->sequence_name[j], CHAR(STRING_ELT(fnames, j)));
		ft->values[j] = INTEGER_POINTER(fvals);
		ft->value_length[j] = LENGTH(VECTOR_ELT(filterList, j));
		if(ft->value_length[j]<2)error("Filter must have the form: list('chr1'=c(100,200,3000,3010,start,end,...),...");
	}
	ft->seqn=flen;
	return(upc);
}

/**
* @brief Free resources occupied by filter
*
* @param ft Pointer to filter
* @return void
* @todo nothing
*/
void destroy_filter(filter_t * ft){
	int j=0;
	for(;j<ft->seqn;j++){
		if(ft->sequence_name[j])Free(ft->sequence_name[j]);
	}
	if(ft->sequence_name)Free(ft->sequence_name);
	if(ft->value_length)Free(ft->value_length);
	if(ft->values)Free(ft->values);
}


/**
* @brief Match Filter and chromosome
*
* @param seqi chromosome string pointer
* @param ft stored struct holding the filter chromosomes.
* @return Returns index of which chromosome in ft matches seqi. -1 for no match
* @details Compares passed chromosome string and looks up ft struct
* @todo nothing
*/
int seq_match(char * seqi, filter_t * ft){//returns index of filter values ft->values
	int i=0;
	for(;i<ft->seqn;i++){
		if(!strcmp(seqi,ft->sequence_name[i]))return(i);
	}
	return(-1);
}

/**
* @brief Free resources of one chromosome vector
*
* @param cptr  Pointer to cpt vector holding the scores to be destroyed
* @param lstart_ind The lstart_ind to be destroyed
* @param lend_ind The lend_ind to be destroyed
* @return void
* @details Free resources of one chromosome vector
* @note
* @todo nothing
*/
void destroy_scores(usersize *cptr,uint32_t *lstart_ind, uint32_t *lend_ind){
	Free(cptr);
	Free(lstart_ind);
	Free(lend_ind);
}

/**
* @brief Fill global_densities_t struct
*
* @param gd global_densities_t to be filled
* @param bresults Pointer to seq_block_t struct holding the results from parsing
* @return void
* @details Fill global_densities_t struct
* @note
* @todo nothing
*/
void copy2globals(global_densities_t *gd,seq_block_t *bresults){
	gd->total_reads+=bresults->total_reads;
	gd->mapmass+=bresults->mapmass;
	gd->lowqual+=bresults->lowqual;
	gd->maxscore=max(gd->maxscore,bresults->maxScore);
	gd->collapsed+=bresults->collapsed;
	gd->filtered_reads+=bresults->filtered_reads;
	gd->neg_strand+=bresults->neg_strand;
	gd->pos_strand+=bresults->pos_strand;
	gd->paired+=bresults->paired;
	gd->ppairs+=bresults->ppairs;
}

/**
* @brief Generate score vector and dispatch to parse_sam.c
*
* @param gd global_densities_t to be filled
* @param user_args Pointer to seq_block_t struct holding the results from parsing
* @param bam_file samfile_t structure holding the header
* @param ft Filter previously set by set_filter()
* @return void
* @details Loops through chromosomes in the header and generates a score vector for each. Dispatch vector to parse_sam.c and digests and filters the return values.
* @note
* @todo nothing
*/
void write_density(global_densities_t * gd,user_arguments_t * user_args,samfile_t *bam_file,filter_t * ft){
	int *l_indp,*gen_indp,scount_ind,filter_index=-1,first_call=1;
    uint32_t chrom_tid=0,tid_vector=0,maxscores=0,cc=0,*lstart_ind,*lend_ind,*proc_chrom;
	usersize *cptr,*scoresp;//should be user defined

	time_t start,stop;
	seq_block_t bresults;
	SEXP scores,l_ind,gen_ind;
	chromosome_size_t cs={0};//here are all the limits for a chromosome

	proc_chrom=Calloc(bam_file->header->n_targets,uint32_t);
	bresults=seq_density(NULL,NULL,NULL,user_args,&cs,bam_file,&first_call);//fetch first read

	copy2globals(gd,&bresults);
	chrom_tid=bresults.chrom_index_next;
	if(user_args->VERBOSE>0)printStatus(bam_file->header->target_name[chrom_tid], &cc, bam_file->header->n_targets);

	while(bresults.file_status>0 && chrom_tid<bam_file->header->n_targets){//loop through chromosomes and fetch densities
		if(user_args->VERBOSE>0)printStatus(bam_file->header->target_name[chrom_tid], &cc, bam_file->header->n_targets);

    	time(&start);

		cs.min_pos=0;//start position first read
		cs.max_pos=bam_file->header->target_len[chrom_tid]+1;//+1 for 0 based sources

		cs.min_scorespace=cs.max_pos+user_args->EXTEND;//expect maximally this amount of values in resulting vector - will be refined after scanning the bam
		gd->gsize+=cs.min_scorespace;//calculate total amount of genomic region covered [spanning first read to last read]
		//the average read size is set to only 30 to be on the safe side
		cs.min_indexspace = cs.min_scorespace/(user_args->COMPRESSION+30);
		if(cs.min_indexspace<50)cs.min_indexspace=1000;//leave some space for testing

		#if verbose==1
			printf("\nSTATUS %d\n",bresults.file_status);
			printf("TID+1 %d <= TE %d\n",chrom_tid+1,bam_file->header->n_targets);
			printf("Index Preallocated for %s: %d Scores Preallocated: %d | TID: %d CC: %d\n",bam_file->header->target_name[chrom_tid],cs.min_indexspace,cs.max_pos,chrom_tid,cc);
		#endif
		cptr = Calloc(cs.max_pos+1,usersize);
		lstart_ind=Calloc(cs.min_indexspace,uint32_t);
		lend_ind=Calloc(cs.min_indexspace,uint32_t);

		bresults=seq_density(cptr,lstart_ind,lend_ind,user_args,&cs,bam_file,&first_call);//check if EOF -> status==0)

		if(user_args->FILTER){//skip chromosome in case its not found in filter sequences
			filter_index=seq_match(bam_file->header->target_name[chrom_tid],ft);//is the current chromosome in the filter list? If yes return filter index
			if(filter_index<0){
				time(&stop);
				chrom_tid=bresults.chrom_index_next;
				copy2globals(gd,&bresults);
				destroy_scores(cptr,lstart_ind,lend_ind);
				continue;
			}

			maxscores=bam_file->header->target_len[chrom_tid];//maximal chromosomal position
			cs.min_scorespace=filter_coverage(filter_index,ft);//override scorespace by bps covered by the filter
			cs.min_indexspace=ft->value_length[filter_index]/2;//override index space allocation
		}
		proc_chrom[chrom_tid]=1;//memorize which chromosomes where processed

		if(bresults.file_status<-5)warning("Error. Return value: %d\n",bresults.file_status);
		if(bresults.file_status==-5)warning("bam file doesn\'t appear to be sorted!");
		if(bresults.file_status<0){
			destroy_scores(cptr,lstart_ind,lend_ind);
			break;
		}

		#if verbose==1
		printf("SCAN MEM ALLOCATED: %ldMB\n",((sizeof(usersize)*cs.max_pos)+(2*cs.min_indexspace*sizeof(uint32_t)))/1000000);
		printf("RETURN VALUE %d\n",bresults.file_status);
		printf("R VECTOR ALLOCATION SIZE/LENGTH: %dMB / %d\n",(sizeof(usersize)*cs.min_scorespace)/1000000,sizeof(usersize)*cs.min_scorespace);
		printf("R VECTOR INDEX ALLOCATION SIZE/LENGTH: %dMB / %d\n",(4*(cs.min_indexspace+1))/1000000,cs.min_indexspace+1);
		printf("%d LOW QUALITY READS SKIPPED\n",bresults.lowqual);
		#endif

		cs.min_scorespace++;
		cs.min_scorespace=cs.min_scorespace/(4/sizeof(usersize));//Only half the space needed with uint16_t
		PROTECT(scores = NEW_INTEGER(cs.min_scorespace++));gd->upcounter++;//initialize compressed scores
		PROTECT(l_ind = NEW_INTEGER(cs.min_indexspace+1));gd->upcounter++;//initialize genomic index with starting positions of 0 blocks
		PROTECT(gen_ind = NEW_INTEGER(cs.min_indexspace+1));gd->upcounter++;//initialize linear index with list indexes of block starts

		scoresp =(usersize*) INTEGER_POINTER(scores);

		l_indp = INTEGER_POINTER(l_ind);
		gen_indp = INTEGER_POINTER(gen_ind);

		/*Send results to R*/
		uint32_t bstart=0,bend=0,bpos,scount=0,fcount=0,windowt=0,histc=0;
		for(scount_ind=0;scount_ind<cs.min_indexspace;scount_ind++){//go through the index
			if(!user_args->FILTER){
				bstart=lstart_ind[scount_ind];//start of data block in scores
				bend=lend_ind[scount_ind];//end of data block in scores
				if(bstart>bend){
					printf("--START--> %d --END--> %d --DIV--> %d\n",bstart,bend,bend-bstart);
					error("--POSSIBLE INDEX ERROR--> BEND-BSTART<0!\n");
				}
			}else{//intersect filter and index
				bstart=ft->values[filter_index][fcount++];//C index 0 based
				bend=ft->values[filter_index][fcount++];//C index 0 based
				if(bstart>maxscores){
					gd->lsize+=bend-bstart;
					if(fcount>3)bend=ft->values[filter_index][fcount-3];//just in case the end of the chromosome is reached -> use last end as index end
					else bend=bstart;
					break;
				}else if(bend>maxscores){
					gd->lsize+=bend-bstart;
					bend=maxscores;
				}
			}
			gen_indp[scount_ind]=bstart+1;//C index 0 based -> genomic position +1
			l_indp[scount_ind]=scount;//R index 1 based -> but accession will be in C
			for(bpos=bstart;bpos<=bend;bpos++){//go through every base pair
				#if pedantic==1
					if(scount>cs.min_scorespace*(4/sizeof(usersize)) || scount_ind>cs.min_indexspace){error("INDEXING FAULT\n");
					}
				#endif
				scoresp[scount++]=*(cptr+bpos);
				gd->lmaxScore=max(*(cptr+bpos),gd->lmaxScore);
				gd->lmapmass+=*(cptr+bpos);
				if(user_args->HWINDOW>1){
					histc++;
					windowt+=*(cptr+bpos);
					if(histc%user_args->HWINDOW==0){
						++gd->histogramp[windowt/user_args->HWINDOW];
						windowt=0;
					}
				} else if(user_args->HWINDOW==1)++gd->histogramp[*(cptr+bpos)];

			#if verbose==2
				printf("FINALLY VALUE2@ %d -> %d || POS: %d\n",scount,scoresp[scount-1],bpos);
			#endif
			}
			if(user_args->HWINDOW>1)++gd->histogramp[windowt/user_args->HWINDOW];//empty last window | prerequisite: HWINDOW<=COMPRESSION!
			histc=0;windowt=0;
			gd->lsize+=bend-bstart;//calculate amount of bps covered outside of blocks

		}

		gen_indp[scount_ind]=user_args->FILTER ? bend+1 : cs.max_pos+1;//override read based cs.maxpos to filter end in case of FILTER. +1 to indicate start of next block
		l_indp[scount_ind]=scount;
		#if pedantic==1
			if(scount>cs.min_scorespace*(4/sizeof(usersize)))error("Uninitialized variables in score vector detected: COUNT %d >  SPACE %d\n",scount,cs.min_scorespace*(4/sizeof(usersize)));
			if(scount_ind!=cs.min_indexspace)error("Uninitialized variables in index vector detected: %d <> %d\n",scount_ind,cs.min_indexspace);
		#endif


    	char lindname[strlen(bam_file->header->target_name[chrom_tid])+5];
    	sprintf(lindname,"%s",bam_file->header->target_name[chrom_tid]);
    	strcat(lindname,"_lind");

    	char gindname[strlen(bam_file->header->target_name[chrom_tid])+5];
    	sprintf(gindname,"%s",bam_file->header->target_name[chrom_tid]);
    	strcat(gindname,"_gind");

    	SET_STRING_ELT(gd->list_names,tid_vector,mkChar(bam_file->header->target_name[chrom_tid]));
    	SET_STRING_ELT(gd->list_names,tid_vector+gd->total_elements,mkChar(gindname));
    	SET_STRING_ELT(gd->list_names,tid_vector+(gd->total_elements*2),mkChar(lindname));
		SET_VECTOR_ELT(gd->list, tid_vector, scores);
		SET_VECTOR_ELT(gd->list, tid_vector+gd->total_elements, gen_ind);
		SET_VECTOR_ELT(gd->list, tid_vector+(gd->total_elements*2), l_ind);

		tid_vector++;

		time(&stop);

		#if verbose>0
		printf("\nAbout %.0f seconds. %d reads processed on %s\n", difftime(stop, start),bresults.total_reads,bam_file->header->target_name[chrom_tid]);
		printf("l_ind set at EL# %d g_ind at %d and scores at %d \n\n",tid_vector,tid_vector+gd->total_elements,tid_vector+(gd->total_elements*2));
		#endif
		chrom_tid=bresults.chrom_index_next;
		copy2globals(gd,&bresults);
		destroy_scores(cptr,lstart_ind,lend_ind);
    }
	if(user_args->VERBOSE>0 && bresults.file_status>=0){
		if(cc<=bam_file->header->n_targets){
			printf("\nWarning: the following chromosomes have no reads\n");
			int x=0,y=0;
			for(;x<bam_file->header->n_targets;x++){
				if(!proc_chrom[x]){
					printf("%s ",bam_file->header->target_name[x]);
					++y;
					if(y%5==0)printf("\n");
				}
			}
			printf("\n");
		}
	}
	Free(proc_chrom);
	if(bresults.file_status<0)gd->total_reads=0;
	return;
}



/**
* @brief Main entry point for R
*
* @param bamfilenameR Filename of read container
* @param aRgvals Vector containing the user arguments
* @param filterList passed from R in the form list("chr1"=c(100,200,3000,3010...start,end...))
* @return R list in the form list("chr1"=c(1,2,1,2,1,1,...),"chr1_gind"=c(1100,1200...),"chr1_lind"=c(0,112,...),"chrX"=NA,...) and a Statistics vector
* @details All chromosome of the filter or all chromosomes in the file header will be scanned and passed to an R list
* @note
* @todo high_cov not yet implemented.
*/
SEXP construct_dc(SEXP bamfilenameR, SEXP aRgvals, SEXP filterList) {

	double *statsp;//resulting statistics in the order "total reads" "coverage" "local coverage" "max score"
	uint32_t upcounter=0,i=0;
	time_t tstart,tstop;
	global_densities_t gd={0};
    user_arguments_t user_args;
    filter_t ft;
    SEXP histogram,stats;
    int *argvalsp;

    signal(SIGINT,SIG_DFL);//make this thing stop on CTRL+C
    time(&tstart);

    /* Set user defined values */

    PROTECT(aRgvals=AS_INTEGER(aRgvals));upcounter++;
    if(LENGTH(aRgvals)!=10)error("Invalid amount of arguments - arguments[%d] / should be %d!\n",LENGTH(aRgvals),9);
    argvalsp=INTEGER_POINTER(aRgvals);
    user_args.bamfilename = STRING_VALUE(bamfilenameR);
    user_args.READTHROUGH = argvalsp[0];//bool. read from start to end and 0 take whole read whithout CIGAR splice info
    user_args.PAIRED = argvalsp[1];
    user_args.STRANDED = argvalsp[2];//Set to 1 / -1 it will use only forward / reverse reads respectively. 0 means all reads are processed
    user_args.TMAPQ = argvalsp[3];//Minimum MPAQ score. Lower scored reads will be skipped
    user_args.COLLAPSE = argvalsp[4];
    user_args.EXTEND = argvalsp[5];//extend each read in its direction by this amount of BPs
    user_args.HWINDOW = argvalsp[6];
    user_args.COMPRESSION = argvalsp[7];//minimum BPs needed between data blocks to collapse the gap and index it
    user_args.VERBOSE = argvalsp[8];
    user_args.UNIQUE = argvalsp[9];


    /* Try to open the file */
    samfile_t *bam_file;
    bam_file=open_samtools(user_args.bamfilename);
	if(!bam_file){
		warning("sam/bam file not found!\n");
		UNPROTECT(upcounter);
	    return(R_NilValue);
	}

    if(user_args.HWINDOW>user_args.COMPRESSION){
    	warning("HWINDOW has to be smaller than COMPRESSION! HWINDOW updated to %d\n",user_args.COMPRESSION);
    	user_args.HWINDOW=user_args.COMPRESSION;
    }

	PROTECT(histogram = NEW_INTEGER(UINT16_MAX));upcounter++;//initialize compressed scores
	gd.histogramp = (uint32_t*) INTEGER_POINTER(histogram);
	for(i = 0; i < UINT16_MAX; i++) gd.histogramp[i] = 0;

	gd.total_elements=bam_file->header->n_targets;//one vector per chromosome needed
	/* ####  CHECK IF THERE IS AN ACTIVE FILTER IN PLACE */
    user_args.FILTER=isNewList(filterList) ? 1 : 0;
    if(user_args.FILTER){
    	upcounter+=set_filter(filterList,&ft);
    	gd.total_elements=ft.seqn;//overwrite total elements if filter is passed, since one density is returned per slice
    }

	// Creating a list with vector elements as many as sequences plus a character string vector:
    PROTECT(gd.list = allocVector(VECSXP, (gd.total_elements*3)+2));upcounter++;//3x for the two indexes and scores per chromosome
    PROTECT(gd.list_names = allocVector(STRSXP,(gd.total_elements*3)+2));upcounter++;//+1 for statistics vector +1 for the histogram

	/* PASS EVERYTHING */
	write_density(&gd,&user_args,bam_file,&ft);
	if(!gd.total_reads)goto NO_READS_FOUND;
	// 1 total_reads  2 gcoverage  3 lcoverage  4 maxscore  5 lmaxscore  6 lowqual  7 filtered  8 collapsed  9 paired  10 proper_pairs 11 pos  12 neg 13 fmapmass 14 lsize 15 gsize
	SET_STRING_ELT(gd.list_names,gd.total_elements*3,mkChar("Statistics"));
	PROTECT(stats = NEW_NUMERIC(15));upcounter++;
	statsp = NUMERIC_POINTER(stats);
	*statsp++=(double)gd.total_reads;
	*statsp++=(double)gd.mapmass/(double)gd.gsize;
	*statsp++=(double)gd.lmapmass/(double)gd.lsize;
	*statsp++=(double)gd.maxscore;
	*statsp++=(double)gd.lmaxScore;
	*statsp++=(double)gd.lowqual;
	*statsp++=(double)gd.filtered_reads;
	*statsp++=(double)gd.collapsed;
	*statsp++=(double)gd.paired;
	*statsp++=(double)gd.ppairs/2;
	*statsp++=(double)gd.pos_strand;
	*statsp++=(double)gd.neg_strand;
	*statsp=(double)gd.mapmass;
	*statsp++=(double)gd.lsize;
	*statsp++=(double)gd.gsize;


	if(gd.lmaxScore>=umaxof(usersize)-1){
		warning("\nThe maximum pile up is exceeding the maximal value of UINT16_MAX=%d. Reads have been capped to that value.\nConsider to rerun using the maxDups option!\n",UINT16_MAX);
	}

	SET_VECTOR_ELT(gd.list,gd.total_elements*3, stats);

	SET_STRING_ELT(gd.list_names,(gd.total_elements*3)+1,mkChar("Histogram"));
	SET_VECTOR_ELT(gd.list,(gd.total_elements*3)+1,histogram);

    setAttrib(gd.list, R_NamesSymbol, gd.list_names);

    NO_READS_FOUND:
    time(&tstop);
	if(user_args.VERBOSE>0)printf("About %.0f seconds passed. %llu reads processed \n", difftime(tstop, tstart),gd.total_reads);
	close_bamfile(bam_file);
	if(user_args.FILTER)destroy_filter(&ft);
    UNPROTECT(upcounter+gd.upcounter);
    if(!gd.total_reads)return(R_NilValue);
    else return(gd.list);
}


