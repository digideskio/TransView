/*
 * visuals.c
 *
 *  Created on: May 29, 2012
 *      Author: Julius Muller
 */

/**
 *  @file	visuals.c
 *  @brief	Provides some simple visual functions for TransView
 */

#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include "visuals.h"

/**
* @brief Simple progress bar
*
* @param x completed rounds
* @param n total rounds
* @param r resolution of bar
* @param w width of bar
* @return void
* @details Writes progress bar to the screen
* @note
* @todo Nothing
*/
void progress_bar(int x, int n, int r, int w)
{
	int i;
	printf("\r"); // Move to the first column
    if ( x % (n/r) != 0 ) return;// Calculuate the ratio of complete-to-incomplete.
    float ratio = x/(float)n;
    int c= ratio * w;
    printf("%3d%% [", (int)(ratio*100) );// Show the percentage completed.
    for (i=0; i<c; i++)printf("=");// Show the load bar.
    for (i=c; i<w; i++)printf(" ");
    printf("]");
    R_FlushConsole();
    R_CheckUserInterrupt();
#ifdef _Win32
  R_ProcessEvents();
#endif
}



/**
* @brief Simple line displaying processed chromosomes
*
* @param cur_chrom Currently processed Chromosome
* @param total_chrom Total amount of chromosomes
* @param chrom_numb Current chromosome index
* @return void
* @details Writes processed chromosome to the screen
* @note
* @todo Nothing
*/
void printStatus(char *cur_chrom, uint32_t *chrom_numb, uint32_t total_chrom)
{
	int i=0;
	printf("\r"); // Move to the first column
    printf("Reading %s (%d of %d)",cur_chrom,(*chrom_numb)++,total_chrom);
    for (; i<20; i++)printf(" ");//overwrite previous longer chromosome names
    if(total_chrom+1==*chrom_numb)printf("\n");//if the end is reached, go to the next line
    R_FlushConsole();//not necessary?
    R_CheckUserInterrupt();/* Minor Memory leak -> At this point only struct filter_t, the bam_header and proc_chrom are allocated. */
#ifdef _Win32
  R_ProcessEvents();
#endif
}
