/***********************************************************************\
 ****************************************************************************
 *                                                                          *
 *  Copyright (C) 2017  Genome Research Ltd.                                *
 *                                                                          *
 *  Author: Zemin Ning (zn1@sanger.ac.uk)                                   *
 *                                                                          *
 *  This file is part of scanPAV  pipeline.                                 *
 *                                                                          *
 *  Scaff10x is a free software: you can redistribute it and/or modify it   *
 *  under the terms of the GNU General Public License as published by the   *
 *  Free Software Foundation, either version 3 of the License, or (at your  *
 *  option) any later version.                                              *
 *                                                                          *
 *  This program is distributed in the hope that it will be useful, but     *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of              *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU        *
 *  General Public License for more details.                                *
 *                                                                          *
 *  You should have received a copy of the GNU General Public License along *
 *  with this program.  If not, see <http://www.gnu.org/licenses/>.         *
 *                                                                          *
 ****************************************************************************

****************************************************************************/

 
#include <math.h>
#include <values.h>
#include <stdio.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>
#include "fasta.h"


static int *list,*readlength,*out_list,*ctg_list,*rd_group;
static long *hist;
static long n_Entry;
static int K2_LEN=17;
static int mhist=10;
static int mhist2=2;
static int breakmod=0;
static int NEDGE=250;
static int nmatch=30;
static int nmatch2=30;
static int N_SET=2000;
static int N2_SET=12000;
static char *SCGname=NULL;//SingleCopyGenome
static long baseSize;
static int n_grouped;
static int fastq_flag=1;
static int n_maxblock;
static int n_blockreads;
static int n_block;
static int Max_N_NameBase=60;

/*   Subroutine to sort the DNA sequences into a matrix table   */
/* =============================================  */
void main(int argc, char **argv)
/* =============================================  */
{
     long i,j,k;
     long ps=0,ns=0,IntSeg=0,IntBase=0,lii;
     fasta *seq;
     int nSeq;
     long totalBases;
     long mhistc = 0;
     long mhistcc = 0;
     unsigned long *sarr, *sarrp;
     int nsorts = 1024;
//     int tb = 0;
     int args,ac;
     int m_st,m_ed,step_len,pmod;
     int nseq = 0,seqc,n_patch;
     int sshift,rshift=6;
     unsigned long nmask;
     int id_read[64],num_sect,n_reads,max_edge; //id_read will never exceed mhist
     int rd_st,rd_ed,n_input;
     unsigned int **imatrix(long nrl,long nrh,long ncl,long nch);
     char line[500] = {0},outName[60]={0},outFast[60]={0},syscmd[200]={0};
     char *ptr,base[20],zero[20]={0},line2[500]={0};
     FILE *fpMate,*fpOutname,*fpOutfast;
     int rd,read_pair[200];
     long gtBases=0,nclip=0;
     void Phusion_Stage(char **argv, int argc, int args);

     if(argc < 2)
     {
         printf("Usage: %s <assembly.fasta> <assembly_gapsized.fasta> <assembly_gapsized.agp>\n",argv[0]);
         exit(1);
     }
/*   sort all the names of genes or name entries   */
     args=1;
     for(i=1;i<argc;i++)
     {
       if(!strcmp(argv[i],"-kmer"))
       {
	 i++;
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-set2"))
       {
         sscanf(argv[++i],"%d",&N2_SET);
         args=args+2;
       } 
     }
     nseq = 0;
     Phusion_Stage(argv,argc,args);

     printf("All jobs finished\n");
     return;
}


/*   Subroutine to sort the DNA sequences into a matrix table   */
/* =============================================  */
void Phusion_Stage(char **argv, int argc, int args)
/* =============================================  */
{
     long i,j,k,iseq;
     fasta *seqp;
     fasta *seq;
     long totalBases,total_len;
     int nSeq,num_Ns,i_contig;
     int qthresh=23;
     int ac,stopflag,n_contig,offset,n_Ns,offset_gap;
     char outName[Max_N_NameBase],name_tag[10],namep_tag[10];
     void ArraySort_String(int n,char Pair_Name[][Max_N_NameBase],int *brr);
     char **cmatrix(long nrl,long nrh,long ncl,long nch);
     FILE *fp,*fpOutfast,*fpOutfast2,*namef;
     fasta *segg;
     long Size_pdata,Size_q_pdata;
     int num_seqque;
     char *pdata,*st;

     printf("www: \n");
     if((fp=fopen(argv[args],"rb"))==NULL) error("Cannot open file\n");
     fseek(fp, 0, SEEK_END);
     Size_q_pdata = ftell(fp) + 1;
     fclose(fp);
     if((pdata=(char*)calloc(Size_q_pdata,sizeof(char)))==NULL)
       error("calloc pdata\n");
     num_seqque = extractFastq(argv[args],pdata,Size_q_pdata);
     if((segg=(fasta*)calloc((num_seqque+1),sizeof(fasta)))==NULL)
       error("calloc segg\n");
     if((seq=decodeFastq(argv[args],&num_seqque,&totalBases,pdata,Size_q_pdata,segg))==NULL)
       error("no query data found.\n");
     nSeq = num_seqque;
     printf("Number of shotgun reads  %d \n",nSeq);

     fastaUC(seq,nSeq);

     if((fpOutfast = fopen(argv[args+1],"w")) == NULL)
     {
       printf("Unable to open file for fastq out\n");
       exit(1);
     }
     if((fpOutfast2 = fopen(argv[args+2],"w")) == NULL)
     {
       printf("Unable to open file for fastq out\n");
       exit(1);
     }

     n_contig = 0;
     seqp=seq;
     offset=0;
     total_len = 0;
     memset(namep_tag,'\0',10);
     strncpy(namep_tag,"NNNNNN",5);
     memset(name_tag,'\0',10);
     strncpy(name_tag,seqp->name,5);
     for(iseq=0;iseq<nSeq;iseq++)
     {
        int kk,rc,slength,start=0,nline=0,outlen=0,n_base,olen=60,offset2=0;
        char *st;

        offset_gap = 0;
        seqp=seq+iseq;
        slength=seqp->length;
	memset(name_tag,'\0',10);
	strncpy(name_tag,seqp->name,5);
        if(strcmp(name_tag,namep_tag) != 0)
	{
	  n_contig = 0;
          printf("%s %s %d\n",namep_tag,name_tag,n_contig);
          memset(namep_tag,'\0',10);
	  strncpy(namep_tag,name_tag,5);
	}

        i_contig = 1;
        n_base=0;
        n_Ns=0;
        fprintf(fpOutfast,">%s\n",seqp->name);
//        fprintf(fpOutfast2,"%ld %s\n",iseq,seqp->name);
        for(j=0;j<slength;j++) 
        {
           kk=j+1;
           n_Ns=1;
           if((seqp->data[j]!='N')&&(seqp->data[j+1]=='N')&&(seqp->data[j+2]!='N')&&(j<(slength-2)))
              n_Ns=10;
           while((seqp->data[kk]=='N')&&(seqp->data[j]=='N')&&(kk<slength))
           {
             kk++;
             n_Ns++;
           }
           if((j==0)&&(n_Ns>=2))
           {
             offset2=n_Ns;
           }
           if(n_Ns>=2)
           {
             offset=kk-1;
             outlen=offset-n_Ns+1-start-offset2;
             i_contig++;
             nline=outlen/olen;
             total_len = total_len+outlen;
             st=seqp->data+start+offset2;
             i_contig++;
             for(i=0;i<outlen;i++)
                fprintf(fpOutfast,"%c",seqp->data[start+offset2+i]);
             if(n_Ns >= 200)
             {
               n_Ns = 200;
             }
             else if(n_Ns == 2)
             {
               n_Ns = 10;
             }
             if(outlen>0)
             {
               fprintf(fpOutfast2,"gapsm: %s %d %ld %d %d\n",seqp->name,outlen,iseq,offset_gap,start+offset2);
               offset_gap = offset_gap + outlen+n_Ns;
             }
             if((outlen > 0)&&((slength - kk) > 1))
             {
               for(i=0;i<n_Ns;i++)
                  fprintf(fpOutfast,"%c",'N');
             }
             if((n_Ns> 0)&&((slength - kk)< 1))
          printf("offset: %d %d %d %s\n",n_Ns,kk,slength,seqp->name);
             start=offset+1;
             n_Ns=0;
             n_base=0;
             offset2=0;
             n_contig++;
           }
           j=kk-1; 
        }

        num_Ns=n_Ns-1;
        outlen=slength-start-n_Ns+1;
        if(outlen >=10)
        {
          i_contig++;
          nline=outlen/olen;
          total_len = total_len+outlen;
          fprintf(fpOutfast2,"gapsm: %s %d %ld %d %d\n",seqp->name,outlen,iseq,offset_gap,start+offset2);
//          printf("offset: %d %d %d %ld\n",start,outlen,n_Ns,total_len);
          {
            for(i=0;i<outlen;i++)
               fprintf(fpOutfast,"%c",seqp->data[start+i]);
          }
        }
        fprintf(fpOutfast,"\n");
     }
     fclose(fpOutfast);
     if(seq){
         free(seq->name);
         free(seq);
         seq = NULL;
     }

}


/* creat an int matrix with subscript ange m[nrl...nrh][ncl...nch]  */
unsigned int     **imatrix(long nrl,long nrh,long ncl,long nch)
{
        long i, nrow=((nrh-nrl+1)/100 + 1)*100,ncol=nch-ncl+1;
        unsigned int  **m;
	long nri,nrn=nrow/100;

        /* allocate pointers to rows        */
        if((m=(unsigned int **)calloc(nrow,sizeof(int*)))==NULL)
        {
           printf("error imatrix: calloc error No. 1 \n");
           return(NULL);
        }

        /* allocate rows and set pointers to them        */
	/* allocate in 100 batches to use freed memory */
	nrl = 0;
	for(nri=0;nri<100;nri++,nrl+=nrn) {
           if((m[nrl]=(unsigned int *)calloc(nrn*ncol,sizeof(int)))==NULL)
           {
              printf("error imatrix: calloc error No. 2 \n");
              return(NULL);
           }

           for(i=1;i<nrn;i++)
              m[i+nrl]=m[i+nrl-1]+ncol;
	}
       /* return pointer to array of pointers to rows   */
        return m;
}

#define SWAP(a,b) temp=(a);(a)=b;(b)=temp;

/*   Subroutine to sort an array arr[0,...,n-1] into ascending order while
     making the corresponding reaarangement of the array brr[0,...,n-1]
     by the use of Quicksort (Sedgwick, R. 1978, Communications o fthe ACM,
     vol. 21, pp. 847-857) also see Numerical Recipes in C                  */  



/*   to swap the string arrays           */
/* ============================================= */
void s_swap(char Pair_Name[][Max_N_NameBase], int i, int j)
/* ============================================= */
{
     char temp[Max_N_NameBase];

     strcpy(temp,Pair_Name[i]);
     strcpy(Pair_Name[i],Pair_Name[j]);
     strcpy(Pair_Name[j],temp);
}


/*   to sort the string array in order          */
/* ============================================= */
void ArraySort_String(int n, char Pair_Name[][Max_N_NameBase], int *brr)
/* ============================================= */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int temp,MIN=7;
     char p[Max_N_NameBase];

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             strcpy(p,Pair_Name[j]);
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(strcmp(Pair_Name[i],p)<=0) break;
                strcpy(Pair_Name[i+1],Pair_Name[i]);
                brr[i+1]=brr[i];
             }
             strcpy(Pair_Name[i+1],p);
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          s_swap(Pair_Name,k,m+1);
          SWAP(brr[k],brr[m+1]);

          if(strcmp(Pair_Name[m],Pair_Name[ir])>0)
          {
            s_swap(Pair_Name,m,ir);
            SWAP(brr[m],brr[ir]);
          }

          if(strcmp(Pair_Name[m+1],Pair_Name[ir])>0)
          {
            s_swap(Pair_Name,m+1,ir);
            SWAP(brr[m+1],brr[ir]);
          }

          if(strcmp(Pair_Name[m],Pair_Name[m+1])>0)
          {
            s_swap(Pair_Name,m,m+1);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          strcpy(p,Pair_Name[m+1]);
          b=brr[m+1];
          for(;;)
          {
             do i++; while (strcmp(Pair_Name[i],p)<0);
             do j--; while (strcmp(Pair_Name[j],p)>0);
             if(j<i) break;
             s_swap(Pair_Name,i,j);
             SWAP(brr[i],brr[j]);
          }
          strcpy(Pair_Name[m+1],Pair_Name[j]);
          strcpy(Pair_Name[j],p);
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/* creat a char matrix with subscript ange c[nrl...nrh][ncl...nch]  */
char     **cmatrix(long nrl,long nrh,long ncl,long nch)
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        char  **m;

        /* allocate pointers to rows        */
        if((m=(char **)calloc(nrow,sizeof(char*)))==NULL)
        {
           printf("error cmatrix: calloc error No. 1 \n");
           return(NULL);
        }
        m+=0;
        m-=nrl;

        /* allocate rows and set pointers to them        */
        if((m[nrl]=(char *)calloc(nrow*ncol,sizeof(char)))==NULL)
        {
           printf("error cmatrix: calloc error No. 2 \n");
           return(NULL);
        }
        m[nrl]+=0;
        m[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           m[i]=m[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return m;
}





