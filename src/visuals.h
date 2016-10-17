
#pragma once

#include <stdint.h>

#define printf Rprintf

void progress_bar(int x, int n, int r, int w);
void printStatus(char *cur_chrom, uint32_t *chrom_numb, uint32_t total_chrom);
