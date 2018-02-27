/***********************************************************************\
 ****************************************************************************
 *                                                                          *
 *  Copyright (C) 2017  Genome Research Ltd.                                *
 *                                                                          *
 *  Author: Zemin Ning (zn1@sanger.ac.uk)                                   *
 *                                                                          *
 *  This file is part of scanPAV pipeline.                                 *
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
#include <unistd.h>
#include <sys/types.h> 
#include <sys/socket.h>
#include <sys/wait.h>
#include <sys/signal.h>
#include <errno.h>
#include "fasta.h"

#define MAXLINE 4096
#define ENDS_EXTRA 0
#define PADCHAR '-'
#define Max_N_NameBase 200
#define Max_N_Pair 100
static char **cell_name,**S_Name;
static char **rdname;
static int *map_score,*hit_refst,*hit_refed,*hit_rcdex,*hit_quest,*hit_queed;
static int *readIndex,*ctg_index,*cell_index,*read_mask;
static int *cigar_list,*contig_list,*contig_head;
static B64_long *cigar_head,sBase;
static char *dataline,*cigar_line;
static int *hit_score,*hit_iddex,*indel_st,*indel_size,*hit_offset;

/* SSAS default parameters   */
static int IMOD=0;
static B64_long line_len=0;
static int num_insreads=0;
static int num_delreads=0;
static int num_hitreads=0;
static int num_reads=0;
static int num_cline=0;
static int set_shift = 10;
static int insert_flag = 0;
static int file_flag = 0;
static int delete_flag = 0;
static int n_contigs = 0;
static char strain_name[100];

typedef struct
{
       int foffset;
       int fsindex;
} SIO;

fasta *sub;


int main(int argc, char **argv)
{
    FILE *fp,*namef,*fpOutfast,*fpOutfast2;
    int i,j,nSeq,args,i_contig,idt,n_ctg;
    char line[2000]={0},ctgnameout[Max_N_NameBase],tmptext[Max_N_NameBase],tmptext2[Max_N_NameBase];
    char **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
    fasta *seq,*seq2;
    char file_name[2000] = {"presence.fasta"}; 
    void ArraySort_Mix(int n, B64_long *arr, int *brr);
    void Read_Pairs(char **argv,int args,int nLib,int nSeq);
    void Align_Process(char **argv,int args,int nRead);
    void Cigar_Filter(char **argv,int args,int nRead);
    fasta *segg,*seqp,*seqp2;
    B64_long Size_q_pdata;
    int num_seqque,size_range[5],*ctg_list,*ctg_locist,*ctg_locied,*ctg_locist2,*ctg_locied2;
    char *pdata;


    memset(strain_name,'\0',100);
    if(argc < 2)
    {
      printf("Usage: %s -file file_name <input_contig_fastq> <input_match_file> <input_tag_file> <output_vcf/bed_file> <output_fasta_file>\n",argv[0]);
      exit(1);
    }

    nSeq=0;
    args=1;
    for(i=1;i<argc;i++)
    {
       if(!strcmp(argv[i],"-mod"))
       {
         sscanf(argv[++i],"%d",&IMOD); 
         args=args+2;
       }
       else if(!strcmp(argv[i],"-insertion"))
       {
         sscanf(argv[++i],"%d",&insert_flag);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-file"))
       {
         sscanf(argv[++i],"%s",file_name);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-shift"))
       {
         sscanf(argv[++i],"%d",&set_shift);
         args=args+2;
       }
    }

    if((fp=fopen(argv[args],"rb"))==NULL) printf("Cannot open file\n");
      fseek(fp, 0, SEEK_END);
    Size_q_pdata = ftell(fp) + 1;
    fclose(fp);
    if((pdata=(char*)calloc(Size_q_pdata,sizeof(char)))==NULL)
      printf("calloc pdata\n");
    num_seqque = extractFastq(argv[args],pdata,Size_q_pdata);
    if((segg=(fasta*)calloc((num_seqque),sizeof(fasta)))==NULL)
      printf("calloc segg\n");
    if((seq=decodeFastq(argv[args],&num_seqque,&sBase,pdata,Size_q_pdata,segg))==NULL)
      printf("no query data found.\n");
    nSeq = num_seqque;
    fastaUC(seq,nSeq);

    if((fp = fopen(argv[args+1],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }
    i_contig = 0;
    while(!feof(fp))
    {
      fgets(line,2000,fp);
      if(feof(fp)) break;
      i_contig++;
    }
    fclose(fp);


    S_Name=cmatrix(0,nSeq+1,0,Max_N_NameBase);

    i_contig = i_contig + nSeq;
    if((ctg_list= (int *)calloc(i_contig,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: Align_Process - ctg_list\n");
      exit(1);
    }
    if((ctg_locist= (int *)calloc(i_contig,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: Align_Process - ctg_head\n");
      exit(1);
    }
    if((ctg_locied= (int *)calloc(i_contig,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: Align_Process - ctg_loci\n");
      exit(1);
    }
    if((ctg_locist2= (int *)calloc(i_contig,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: Align_Process - ctg_head\n");
      exit(1);
    }
    if((ctg_locied2= (int *)calloc(i_contig,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: Align_Process - ctg_loci\n");
      exit(1);
    }

    if((namef = fopen(argv[args+1],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }
    i_contig = 0;
    while(fscanf(namef,"%s %s %d %d %d %d %d %d",tmptext,tmptext,&idt,&ctg_list[i_contig],&ctg_locist[i_contig],&ctg_locied[i_contig],&ctg_locist2[i_contig],&ctg_locied2[i_contig])!=EOF)
    {
      i_contig++;
    }
    fclose(namef);


/*  read the tag file         */
    if((namef = fopen(argv[args+2],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }
    i=0;
    while(fscanf(namef,"%s %s %s %s",tmptext,tmptext,tmptext,S_Name[i])!=EOF)
    {
        i++;
    }
    fclose(namef);

/* output the VCF file   */
    if((fpOutfast2 = fopen(argv[args+3],"w")) == NULL)
    {
      printf("ERROR main:: reads group file: %s \n",argv[args]);
      exit(1);
    }
/* output the fasta file   */
    if((fpOutfast = fopen(argv[args+4],"w")) == NULL)
    {
      printf("ERROR main:: reads group file: %s \n",argv[args]);
      exit(1);
    }

    fprintf(fpOutfast2,"##fileformat=VCFv4.2\n");
    fprintf(fpOutfast2,"##FILTER=<ID=PASS,Description=%c%s%c>\n",'"',"All filters passed",'"');
    fprintf(fpOutfast2,"##reference=file:%s\n",file_name);
    fprintf(fpOutfast2,"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=%c%s>%c\n",'"',"Type of structural variant",'"');
    fprintf(fpOutfast2,"##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=%c%s%c>\n",'"',"Difference in length between REF and ALT alleles",'"');
    fprintf(fpOutfast2,"##INFO=<ID=END,Number=1,Type=Integer,Description=%c%s%c>\n",'"',"End position of the variant described in this record",'"');
    fprintf(fpOutfast2,"##INFO=<ID=STRANDS,Number=.,Type=String,Description=%c%s%c>\n",'"',"Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)",'"');
    fprintf(fpOutfast2,"##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=%c%s%c>\n",'"',"Imprecise structural variation",'"');
    fprintf(fpOutfast2,"##ALT=<ID=DEL,Description=%c%s%c>\n",'"',"Deletion",'"');
    fprintf(fpOutfast2,"##ALT=<ID=DUP,Description=%c%s%c>\n",'"',"Duplication",'"');
    fprintf(fpOutfast2,"##ALT=<ID=INV,Description=%c%s%c>\n",'"',"Inversion",'"');
    fprintf(fpOutfast2,"##ALT=<ID=DUP:TANDEM,Description=%c%s%c>\n",'"',"Tandem duplication",'"');
    fprintf(fpOutfast2,"##ALT=<ID=INS,Description=%c%s%c>\n",'"',"Insertion of novel sequence",'"');
    fprintf(fpOutfast2,"##ALT=<ID=CNV,Description=%c%s%c>\n",'"',"Copy number variable region",'"');
    fprintf(fpOutfast2,"##FORMAT=<ID=GT,Number=1,Type=String,Description=%c%s%c>\n",'"',"Genotype",'"');
    fprintf(fpOutfast2,"##source=scanPAV\n");


    for(i=0;i<nSeq;i++)
    {
       seqp = seq+i; 
       fprintf(fpOutfast2,"##contig=<ID=%s,length=%d>\n",S_Name[i],seqp->length);
    }
    fprintf(fpOutfast2,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

    n_ctg = 0;
    for(i=0;i<i_contig;i++)
    {
       int k,seq_st,seq_ed,rc,seq_len,st2,ed2,cut_st,cut_ed;
       int trash_flag = 0;

       
       seqp= seq + ctg_list[i];
       seq_st = ctg_locist[i];
       seq_ed = ctg_locied[i];
       if(seq_ed > seqp->length)
         seq_ed = seqp->length;
       if(seq_st > seqp->length)
       {
         seq_st = 0;
         seq_ed = 0;
       }
       if((seq_ed > 0)&&(seq_st >= 0)&&((seq_ed - seq_st)>= 500))
       {
         int nline = (seq_ed-seq_st)/60;
         seq_len = seq_ed-seq_st;

         cut_st = seq_st;
         cut_ed = seq_ed;
         if(seqp->data[seq_ed-1] != 'N')
           cut_ed = seq_ed;
         else
         {
           for(rc=(seq_ed-1);rc>seq_st;rc--)
           {
              if(seqp->data[rc]!='N')
              {
                idt = rc+1;
                rc = 0;
              }
           }
           cut_ed = idt; 
         }
         if(seqp->data[seq_st] != 'N')
           cut_st = seq_st;
         else
         {
           for(rc=seq_st;rc<seq_ed;rc++)
           {
              if(seqp->data[rc]!='N')
              {
                idt = rc;
                rc = seq_ed;
              }
           }
           cut_st = idt; 
         }

         seq_len = cut_ed - cut_st;
         nline = (cut_ed-cut_st)/60;

         fprintf(fpOutfast2,"%s\t%d\t.\t%c\t<DEL>",S_Name[ctg_list[i]],ctg_locist2[i]+1,'N');
         fprintf(fpOutfast2,"\t50\tPASS\tSVTYPE=DEL;SVLEN=-%d;END=%d\n",seq_len,ctg_locied2[i]+1);  

         fprintf(fpOutfast,">PAV_%s_%09d_%08d\n",S_Name[ctg_list[i]],ctg_locist2[i]+1,seq_len);
         for(k=0;k<nline;k++)
         {
            for(j=0;j<60;j++)
               fprintf(fpOutfast,"%c",seqp->data[k*60+j+cut_st]);
            fprintf(fpOutfast,"\n");
         }
         for(j=0;j<(seq_len-(nline*60));j++)
            fprintf(fpOutfast,"%c",seqp->data[nline*60+j+cut_st]);
         if((seq_len%60)!=0)
           fprintf(fpOutfast,"\n");
         n_ctg++;
       }
    }
    fclose(fpOutfast2);
    fclose(fpOutfast);

    return EXIT_SUCCESS;

}
/* end of the main */


#define SWAP(a,b) temp=(a);(a)=b;(b)=temp;

/*   Subroutine to sort an array arr[0,...,n-1] into ascending order while
     making the corresponding reaarangement of the array brr[0,...,n-1]
     by the use of Quicksort (Sedgwick, R. 1978, Communications o fthe ACM,
     vol. 21, pp. 847-857) also see Numerical Recipes in C                  */  

/* =============================== */
void ArraySort_Long(int n, B64_long *arr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,NSTACK=50,istack[NSTACK];
     B64_long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
             }
             arr[i+1]=a;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
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


/* =============================== */
void ArraySort_Int(int n, int *arr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
             }
             arr[i+1]=a;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
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


/* =============================== */
void ArraySort_Mix(int n, B64_long *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     B64_long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
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

/* =============================== */
void ArraySort_Int2(int n, int *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
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

/*   function to sort an array into a decreasing order:  a>b>c>....    */  
/* =============================== */
void ArraySort2_Int2(int n, int *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]>=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]<arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]<arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]<arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]>a);
             do j--; while (arr[j]<a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
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

/* =============================== */
void ArraySort_Mix3(int n, B64_long *arr, int *brr, int *crr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,c,NSTACK=50,istack[NSTACK];
     B64_long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             c=crr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
                crr[i+1]=crr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
             crr[i+1]=c;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);
          SWAP(crr[k],crr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
            SWAP(crr[m],crr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
            SWAP(crr[m+1],crr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
            SWAP(crr[m],crr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          c=crr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
             SWAP(crr[i],crr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          crr[m+1]=crr[j];
          crr[j]=c;
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


/*   to swap the string arrays           */
/* ============================================= */
void s_swap(char **Pair_Name, int i, int j)
/* ============================================= */
{
     char temp[Max_N_NameBase];

     strcpy(temp,Pair_Name[i]);
     strcpy(Pair_Name[i],Pair_Name[j]);
     strcpy(Pair_Name[j],temp);
}


/*   to sort the string array in order          */
/* ============================================= */
void ArraySort_String(int n, char **Pair_Name, int *brr)
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


/* creat an int matrix with subscript ange m[nrl...nrh][ncl...nch]  */
int     **mmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch)
{
        B64_long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        int  **m;

        
        /* allocate pointers to rows        */
        if((m=(int **)calloc(nrow,sizeof(int*)))==NULL)
        {
           printf("error imatrix: calloc error No. 1 \n");
           return(NULL);
        }
        m+=0;
        m-=nrl;

        /* allocate rows and set pointers to them        */
        if((m[nrl]=(int *)calloc(nrow*ncol,sizeof(int)))==NULL)
        {
           printf("error imatrix: calloc error No. 2 \n");
           return(NULL);
        }
        m[nrl]+=0;
        m[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           m[i]=m[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return m;
}

/* creat an int matrix with subscript ange m[nrl...nrh][ncl...nch]  */
int     **imatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch)
{
        B64_long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        int  **m;

        /* allocate pointers to rows        */
        if((m=(int **)calloc(nrow,sizeof(int*)))==NULL)
        {
           printf("error imatrix: calloc error No. 1 \n");
           return(NULL);
        }
        m+=0;
        m-=nrl;

        /* allocate rows and set pointers to them        */
        if((m[nrl]=(int *)calloc(nrow*ncol,sizeof(int)))==NULL)
        {
           printf("error imatrix: calloc error No. 2 \n");
           return(NULL);
        }
        m[nrl]+=0;
        m[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           m[i]=m[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return m;
}

/* creat char matrix with subscript ange cm[nrl...nrh][ncl...nch]  */
char    **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch)
{
        B64_long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        char **cm;

        /* allocate pointers to rows        */
        if((cm=(char **)calloc(nrow,sizeof(char*)))==NULL)
        {
           printf("error cmatrix: calloc error No. 1 \n");
           return(NULL);
        }
        cm+=0;
        cm-=nrl;

        /* allocate rows and set pointers to them        */
        if((cm[nrl]=(char *)calloc(nrow*ncol,sizeof(char)))==NULL)
        {
           printf("error cmatrix: calloc error No. 2 \n");
           return(NULL);
        }
        cm[nrl]+=0;
        cm[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           cm[i]=cm[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return cm;
}

