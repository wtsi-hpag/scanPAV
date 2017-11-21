/**************************************************************************
 *
 *-------------------------------------------------------------------------
 * Project Notes : various fasta and fastq manipulation routines
 *
 *
 *-------------------------------------------------------------------------
 #######################################################################
 # This software has been created by Genome Research Limited (GRL).    # 
 # GRL hereby grants permission to use, copy, modify and distribute    # 
 # this software and its documentation for non-commercial purposes     # 
 # without fee at the user's own risk on the basis set out below.      #
 # GRL neither undertakes nor accepts any duty whether contractual or  # 
 # otherwise in connection with the software, its use or the use of    # 
 # any derivative, and makes no representations or warranties, express #
 # or implied, concerning the software, its suitability, fitness for   #
 # a particular purpose or non-infringement.                           #
 # In no event shall the authors of the software or GRL be responsible # 
 # or liable for any loss or damage whatsoever arising in any way      # 
 # directly or indirectly out of the use of this software or its       # 
 # derivatives, even if advised of the possibility of such damage.     #
 # Our software can be freely distributed under the conditions set out # 
 # above, and must contain this copyright notice.                      #
 #######################################################################
 *
 *  Author : Zemin Ning and Adam Spargo 
 *
 *	Copyright (C) 2005-2007 by Genome Research Limited, All rights reserved.
 *
 **************************************************************************/


#include <math.h>
#include <values.h>
#include <stdio.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>
#include "fasta.h"

int countQthresh (fasta *seg, int qthresh)
{
	int j;
	fasta *segp;
	char *b;
	int n = 0;
	segp = seg;
	if(segp->finished) return(segp->length);
	b = (char*)segp->qual;
	for(j=segp->length;--j>=0;b++)
		if(*b >= qthresh)
			n++;
	return(n);
}

int extractFastq(char *fname, char *pdata, B64_long Size_pdata) {
  FILE *fil; 
  char *dp;
  B64_long dataSize;
  B64_long i;
  int nSeg;
  int state;

  if((fil = fopen(fname,"r")) == NULL) printf("error: file not found\n");
  fseek( fil, 0, SEEK_END);
  dataSize = ftell( fil);
  rewind( fil);
  if(dataSize == 0) {
    fclose(fil);
    printf("error: no data in file\n");
  }
  memset(pdata,'\0',Size_pdata);
  for(i=0;i<dataSize;i+=1024*1000000) {
    B64_long size = 1024*1000000;
    if(dataSize-i < size) size = dataSize-i;
    if(fread(pdata+i,sizeof(char), size, fil) != size) {
      fclose(fil);
      printf("error: file is too small\n");
    }
  }
  fclose(fil);
  pdata[dataSize] = 0;
  nSeg = 0;
  dp = pdata;
  if(*dp != '>' && *dp != '@') fprintf(stderr,"Corrupt input file %s *dp=>%c< [1].\n", fname, *dp);

  /****************************************************************************
   * Use simple 4-state automaton.
   * -----------------------------
   * state 0 - character is part of sequence header.
   * state 1 - character is part of sequence bases.
   * state 2 - character is part of quality header.
   * state 3 - character is part of quality values.
   ***************************************************************************/
  state=0;
  for(i=dataSize;--i>=0L;dp++) {
    if(state == 0) {
      if(*dp == '@' || *dp == '>') {
        nSeg++;
        continue;
      }
      else if(*dp == '\n') {
        if(i == 0L) printf("Corrupt input file [2].\n");
        state = 1;
        continue;
      }
      else continue;
    }
    else if(state == 1) {
      if((*dp == '\n')) {
        if(i == 0L) break;
        else if((*(dp+1) == '@') || (*(dp+1) == '>')) {
          state=0;
          continue;
        }
        else if(*(dp+1) == '+') {
          state=2;
          continue;
        }
        else continue;
      }
      else continue;
    }
    else if(state == 2) {
      if((*dp == '\n')) {
        if(i == 0L) printf("Corrupt input file [3].\n");
        state = 3;
        continue;
      }
      else continue;
    }
    else if(state == 3) {
      if((*dp == '\n')) {
        if(i == 0L) break;
        else if(*(dp+1) == '@' || *(dp+1) == '>') {
          state=0;
          continue;
        }
        else continue;
      }
      else continue;

    }
    else printf("Corrupt input file [4].\n");
  }
  if(nSeg == 0) printf("error: no segments found\n");
  return(nSeg);
}




fasta *decodeFastq(char *fname,int *nContigs,B64_long *tB,char* pdata,B64_long Size_pdata,fasta *segg)
{
  FILE *fil;
  char *dp, *cp;
  B64_long dataSize, totalbases=0;
  B64_long i;
  int nSeg;
  fasta *seg,*segp; 
  int state, isfastq = 0;
        
  if((fil = fopen(fname,"r")) == NULL) printf("error: file not found\n");
  fseek( fil, 0, SEEK_END);
  dataSize = ftell( fil);
  rewind( fil); 
  if(dataSize == 0) {
    fclose(fil);
    printf("error decodeFastq: no data in file\n");
  } 
  memset(pdata,'\0',Size_pdata);
  for(i=0;i<dataSize;i+=1024*1000000) {
    B64_long size = 1024*1000000;
    if(dataSize-i < size) size = dataSize-i; 
    if(fread(pdata+i,sizeof(char), size, fil) != size) {
      fclose(fil);
      printf("%s\n",fname);
      printf("error: file is too small\n");
    }     
  }       
  fclose(fil);
  pdata[dataSize] = 0;
  nSeg = 0;
  dp = pdata;
    
  if     (*dp == '@') isfastq = 0;
  else if(*dp == '>') isfastq = 1;
  else                printf("Corrupt input file [5].\n");
        
  /****************************************************************************
   * Use simple 4-state automaton.
   * -----------------------------
   * state 0 - character is part of sequence header.
   * state 1 - character is part of sequence bases.
   * state 2 - character is part of quality header.
   * state 3 - character is part of quality values.
   *
   * This should be used for reading data too - current code is unstable.
   *
   ***************************************************************************/
  state=0;
  for(i=dataSize;--i>=0L;dp++) {
    if(state == 0) {
      if(*dp == '@' || *dp == '>') {
        nSeg++;
        continue;
      }
      else if(*dp == '\n') {
        if(i == 0L) printf("Corrupt input file [6].\n");
        state = 1;
        continue;
      }
      else continue;
    }
    else if(state == 1) {
      if((*dp == '\n')) {
        if(i == 0L) break;
        else if((*(dp+1) == '@') || (*(dp+1) == '>')) {
          state=0;
          continue;
        }
        else if(*(dp+1) == '+') {
          state=2;
          continue;
        }
        else continue;
      }
      else continue;
    }
    else if(state == 2) {
      if((*dp == '\n')) {
        if(i == 0L) printf("Corrupt input file [7].\n");
        state = 3;
        continue;
      }
      else continue;
    }
    else if(state == 3) {
      if((*dp == '\n')) {
        if(i == 0L) break;
        else if(*(dp+1) == '@' || *(dp+1) == '>') {
          state=0;
          continue;
        }
        else continue;
      }
      else continue;

    }
    else printf("Corrupt input file [8].\n");
  }
  if(nSeg == 0) printf("no segments found\n");
  cp = dp = pdata;
  seg = segg;
  segp = seg;
  while(dp < pdata+dataSize) {
    if(*dp != '@' && *dp != '>') {
      dp++;
      continue;
    }
    else if ( dp > pdata && dp[-1] != '\n' ) {
      dp++;
      continue;
    }
    else {
      int tmp;
      int k;
      dp++;
      segp->name = (char *)(long long int)(cp-pdata);
      while((*cp++ = tmp = *dp++) != ' ' && tmp != '\t' && tmp != '\n');
      cp[-1] = 0;
      if(tmp != '\n') {
        segp->name2 = (char *)(B64_long)(cp-pdata);
        while((*cp++ = *dp++) != '\n');
        cp[-1] = 0;
      }
      else segp->name2 = NULL;
      segp->data = (char *)(B64_long)(cp-pdata);
      for(k=0;(tmp = *dp) &&
              !((tmp == '+' || tmp == '>') &&
              dp[-1] == '\n');dp++) {
        if(tmp == '\n') continue;
        *cp++ = tmp;
        k++;
      }
      segp->length = k;
      segp->finished = 1;
      segp->qual = NULL;
      if(isfastq) {
        segp->finished = 1;
        segp->qual = NULL;
      }
      else {
        segp->finished = 0;
        k = 0;
        if(tmp == '+') {
          int cp_set = 0;
          int comp = 1;
          while(*dp++ != '\n');
          for(i=0;i<20;i++) {
            if(dp[i] == 0) break;
            if(dp[i] == ' ') {comp = 0; break;}
          }
          segp->qual = (char *)(B64_long)(cp-pdata);
          *cp = 0;
          if(comp) for(k=0;
                       (tmp = *dp) && (k < segp->length);
                       dp++) {
//             if(tmp == '\n') continue;
             if(tmp == '\n') break;
             *cp++ = tmp-041;
             k++;
          }
          else for(k=0;
                   (tmp = *dp) && (k < segp->length);
                   dp++) {
            if(tmp == '\n') break;
            if(tmp >= '0' && tmp <= '9') {
              if(*cp > 9) continue;
              *cp *= 10;
              *cp += tmp-'0';
              cp_set = 1;
              continue;
            }
            if(cp_set) {
              cp++;
              *cp = 0;
              cp_set = 0;
              k++;
            }
          }
        }
        if((segp-seg + 1 < nSeg) && !(*(dp+1) == '@' && *dp == '\n')) {
          fprintf(stderr, "Warning: %-20.20s has %d bases and %d quals, skipping entry\n",pdata+((B64_long)segp->name),segp->length,k);
          cp = pdata+((B64_long)segp->name);
          continue;
        }
      }
      segp++;
      if(segp-seg > nSeg) break;
     }
   }
   nSeg = segp-seg;
   for(i=nSeg,segp=seg;--i>=0;segp++) {
     segp->name = pdata+((B64_long)segp->name);
     if(segp->name2) segp->name2 = pdata+((B64_long)segp->name2);
     segp->data = pdata+((B64_long)segp->data);
     segp->qual = pdata+((B64_long)segp->qual);
     if(segp->finished) segp->qual = NULL;
     totalbases+=segp->length;
   }
   *nContigs = nSeg;
   *tB = totalbases;
   return(seg);
}


void fastaLC (fasta *seg, int nSeg)
{
	int i,j;
	fasta *segp;
	segp = seg;
	for(i=0;i<nSeg;i++,segp++)
	{
		char *b;
		b = (char*)segp->data;
		for(j=segp->length;--j>=0;b++)
			*b = tolower(*b);
	}
}

void fastaUC (fasta *seg, int nSeg)
{
	int i,j;
	fasta *segp;
	segp = seg;
	for(i=0;i<nSeg;i++,segp++)
	{
		char *b;
		b = (char*)segp->data;
		for(j=segp->length;--j>=0;b++)
			*b = toupper(*b);
	}
}


int reverseCompliment(fasta *seq, fasta *rseq)
{
	int i;
	int size=((seq->data != 0)+(seq->qual != 0))*seq->length;
	char *tp,*dp;

	rseq->length = 0;
	if(size <= 0 || seq->data == NULL)
		return(1);

	//printf("reverseCompliment: malloc %d",(size+10*strlen(seq->name)+100) * sizeof(char));
	//fflush(stdout);
	if((rseq->data = (char *)malloc((size+10*strlen(seq->name)+40) * sizeof(char))) == NULL)
	{
		printf("ERROR reverseCompliment: malloc\n");
		return(1);
	}
	//printf("  DONE\n");
	//fflush(stdout);

	rseq->name = rseq->data+size;
	strcpy(rseq->name,seq->name);
	rseq->path = seq->path;
	rseq->length = seq->length;
	rseq->finished = seq->finished;
	rseq->qual = NULL;
	if(seq->qual != NULL)
		rseq->qual = rseq->data+rseq->length;
	
	dp = rseq->data;
	tp = seq->data+seq->length;
	for(i=seq->length;--i>=0;)
	{
		int tmp = *--tp;
		if     (tmp == 'T') *dp++ = 'A';
		else if(tmp == 'G') *dp++ = 'C';
		else if(tmp == 'C') *dp++ = 'G';
		else if(tmp == 'A') *dp++ = 'T';
		else                *dp++ = tmp;
	}
	if(seq->qual != 0)
	{
		dp = rseq->qual;
		tp = seq->qual+seq->length;
		for(i=seq->length;--i>=0;)
			*dp++ = *--tp;
	}
	return(0);
}

int duplicateRead(fasta *seq, fasta *dseq)
{
	int i;
	int size=((seq->data != 0)+(!seq->finished))*seq->length;
	char *tp,*dp;

	dseq->length = 0;
	dseq->data = NULL;
	if(size <= 0 || seq->data == NULL)
		return(1);

	//printf("duplicateRead: malloc %d",(size+10*strlen(seq->name)+100) * sizeof(char));
	//fflush(stdout);
	if((dseq->data = (char *)malloc((size + 10*strlen(seq->name) + 40 )* sizeof(char))) == NULL)
	{
		printf("ERROR reverseCompliment: malloc\n");
		return(1);
	}
	//printf("  DONE\n");
	//fflush(stdout);

	dseq->name = dseq->data+size;
	strcpy(dseq->name,seq->name);
	dseq->path = seq->path;
	dseq->SCFname = NULL;
	if(seq->SCFname != NULL)
	{
		dseq->SCFname = dseq->data+size+strlen(seq->name) + 1;
		strcpy(dseq->SCFname,seq->SCFname);
	}
	dseq->length = seq->length;
	dseq->finished = seq->finished;
	dseq->qual = NULL;
	if(!seq->finished)
		dseq->qual = dseq->data+dseq->length;

	dp = dseq->data;
	tp = seq->data;
	for(i=seq->length;--i>=0;)
		*dp++ = *tp++;

	if(!seq->finished && seq->qual != NULL)
	{
		dp = dseq->qual;
		tp = seq->qual;
		for(i=seq->length;--i>=0;)
			*dp++ = *tp++;
	}
	return(0);
}

