/****************************************************************************
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
/****************************************************************************/

 
#include <math.h>
#include <values.h>
#include <stdio.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include "fasta.h"

#define GT '>'
#define GT4 (((((GT<<8)+GT)<<8)+GT)<<8)+GT

#define ENDS_EXTRA 0
#define PADCHAR '-'
#define MAX_N_BRG 50000 
#define MAX_N_ROW 40000 
#define Max_N_NameBase 400 
#define Max_N_NameBase2 400 
#define Max_N_Pair 100
static char bindir[2000];
static char tmpdir[2000];
static char **S_Name;
static int *insert_siz,*insert_dev,*core_list,*ctg_list,*ctg_head,*read2contig;
static int *readIndex;

/* SSAS default parameters   */
static int n_group=0;
static int m_score=550;
static int num_reads=10;
static int num_links=8;
static int gap_size=20000;
static int len_block=2500;
static int len_edges=50000;
static int len_matrx=2000;
static int min_edge = 5;
static int uplinks = 50;
static int run_align = 1;
static int min_len = 3000;

typedef struct
{
       int foffset;
       int fsindex;
} SIO;


int main(int argc, char **argv)
{
    int i,nSeq,args;
    char *st,*ed;
    char **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
    void ArraySort_String2(int n,char **Pair_Name,int *brr);
    fasta *seq;
    FILE *fp,*namef,*namef2;
    int size_file;
    int n_nodes,n_gap1,n_gap2,n_debug,num_sigma;
    void decodeReadpair(int nSeq);
    void HashFasta_Head(int i, int nSeq);
    void HashFasta_Table(int i, int nSeq);
    void Search_SM(fasta *seq,int nSeq);
    void Assemble_SM(int arr,int brr);
    void Read_Index(int nSeq,char *namefile);
    void Read_Group(fasta *seq,int nSeq,int nRead,int cindex);
    void File_Output(int aaa);
    void Memory_Allocate(int arr);
    char tempa[2000],tempc[2000],syscmd[2000],file_presence[2000],file_absence[2000],file_scanPAV[2000],workdir[2000];
    char file_read1[900],file_read2[900],samname[500],toolname[500];
    int systemRet = system (syscmd);
    int systemChd = chdir(tmpdir);
    pid_t pid;

    seq=NULL;
    system("ps aux | grep scanPAV; date");

    if(argc < 2)
    {
         printf("Usage: %s -nodes 30 -align bwa/smalt -score 550 <presence_assembly.fasta> <absence_assembly.fasta> <presence_absence.fasta>\n",argv[0]);
         printf("       nodes  (30)    - number of CPUs requested\n");
         printf("       score  (550)   - Smith-Waterman alignment score\n");
         printf("       align  (bwa)   - use BWA as the aligner\n");
         exit(1);
    }

    n_nodes = 20;
    n_debug = 0;

    memset(toolname,'\0',500);
    strcpy(toolname,"bwa");
    nSeq=0;
    args=1;
    for(i=1;i<argc;i++)
    {
       if(!strcmp(argv[i],"-gap"))
       {
         sscanf(argv[++i],"%d",&gap_size); 
         args=args+2;
       }
       else if(!strcmp(argv[i],"-debug"))
       {
         sscanf(argv[++i],"%d",&n_debug);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-align"))
       {
         sscanf(argv[++i],"%s",toolname);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-nodes"))
       {
         sscanf(argv[++i],"%d",&n_nodes);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-score"))
       {
         sscanf(argv[++i],"%d",&m_score);
         args=args+2;
       }
    }

    pid = getpid();
    memset(tempa,'\0',2000);
    if (!getcwd(tempa, sizeof(tempa)))
    {
      exit(1);
    } 
    memset(tmpdir,'\0',2000);
    memset(workdir,'\0',2000);
    sprintf(tmpdir,"%s/",tempa);
    sprintf(workdir,"%s/tmp_rununik_%d/",tempa,pid);

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"mkdir %s",workdir);
    if(system(syscmd) == -1)
    {
      printf("System command error:\n");
    }
//    system (syscmd);

    if(chdir(workdir) == -1)
    {
      printf("System command error: chdir\n");
    }
     
    st = argv[0];
    ed = strrchr(argv[0],'/');
    memset(tempc,'\0',2000);
    strncpy(tempc,argv[0],ed-st);
    memset(bindir,'\0',2000);
    sprintf(bindir,"%s/scanPAV-bin",tempc);


    memset(file_presence,'\0',2000);
    memset(file_absence,'\0',2000);
    memset(file_scanPAV,'\0',2000);
    sprintf(file_presence,"%s/%s",tempa,argv[args]);
    sprintf(file_absence,"%s/%s",tempa,argv[args+1]);
    sprintf(file_scanPAV,"%s/%s",tempa,argv[args+2]);

    if((namef = fopen(file_presence,"r")) == NULL)
    {
      printf("File not in the working directory!\n");
      if((namef = fopen(argv[args],"r")) == NULL)
      {
        printf("File %s not found and please copy it to your working directory!\n",argv[args]);
        exit(1);
      }
      else
      {
        memset(file_presence,'\0',2000);
        strcpy(file_presence,argv[args]);
        printf("Input presence file: %s\n",file_presence);
      }
    }
    else
    {
      printf("Input presence file: %s\n",file_presence);
    } 

    if((namef = fopen(file_absence,"r")) == NULL)
    {
      printf("File not in the working directory!\n");
      if((namef = fopen(argv[args+1],"r")) == NULL)
      {
        printf("File %s not found and please copy it to your working directory!\n",argv[args+1]);
        exit(1);
      }
      else
      {
        memset(file_absence,'\0',2000);
        strcpy(file_absence,argv[args+1]);
        printf("Input absence file: %s\n",file_absence);
      }
    }
    else
    {
      printf("Input absence file: %s\n",file_absence);
    } 

    printf("Calling %s %s %s %s\n",argv[args],argv[args+1],argv[args+2],toolname);
   
    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/scanPAV_gapsize %s genome-gapsize.fasta > try.out",bindir,file_presence);
    if(system(syscmd) == -1)
    {
//      printf("System command error:\n);
    }

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/scanPAV_fasta -name presence -len 10 genome-gapsize.fasta presence.fasta presence.tag > try.out",bindir);
    if(system(syscmd) == -1)
    {
//      printf("System command error:\n);
    }

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/scanPAV_fasta -name absence -len 10 %s absence.fasta absence.tag > try.out",bindir,file_absence);
    if(system(syscmd) == -1)
    {
//      printf("System command error:\n);
    }

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/scanPAV_shred -name presence -length 1000 presence.fasta presence.shred > try.out",bindir);
    if(system(syscmd) == -1)
    {
//      printf("System command error:\n);
    }

    if(strcmp(toolname,"bwa") == 0)
    {
      printf("Calling bwa:\n");
      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/bwa index absence.fasta > try.out",bindir);
      if(system(syscmd) == -1)
      {
//      printf("System command error:\n);
      }

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/bwa mem -t %d absence.fasta presence.shred > align.sam",bindir,n_nodes);
      if(system(syscmd) == -1)
      {
//      printf("System command error:\n);
      }

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"cat align.sam | awk '%s' | egrep -v '^@' > align-1.dat","($2<100)&&($5>=0){print $1,$2,$3,$4,$5,$14}");
      printf("%s\n",syscmd);
      if(system(syscmd) == -1)
      {
//      printf("System command error:\n);
      }

    }
    else if(strcmp(toolname,"smalt") == 0)
    {
      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/smalt index -k 13 -s 13 hash_absence absence.fasta  > try.out",bindir);
      if(system(syscmd) == -1)
      {
//    printf("System command error:\n);
      }

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"%s/smalt map -f samsoft -m %d -n %d -o align.sam -O hash_absence presence.shred > try.out",bindir,200,n_nodes);
      if(system(syscmd) == -1)
      {
//    printf("System command error:\n);
      }

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"cat align.sam | awk '%s' | egrep -v '^@' > align-1.dat","($2<100)&&($5>=0){print $1,$2,$3,$4,$5,$13}");
      printf("%s\n",syscmd);
      if(system(syscmd) == -1)
      {
//      printf("System command error:\n);
      }

    }
    else
    {
      printf("No aligner is assigned and input -align bwa or -align smalt\n");
      exit(1);
    }

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/scanPAV_bwasam -score %d align-1.dat align-2.dat > try.out",bindir,m_score);
    if(system(syscmd) == -1)
    {
//  printf("System command error:\n);
    }

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/scanPAV_bwatag absence.tag align-2.dat align-3.dat > try.out",bindir);
    if(system(syscmd) == -1)
    {
//    printf("System command error:\n);
    }

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/scanPAV_clean align-3.dat align.clean1 > try.out",bindir);
    if(system(syscmd) == -1)
    {
//    printf("System command error:\n);
    }

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/scanPAV_bwagap align.clean1 align.dat > try.out",bindir);
    if(system(syscmd) == -1)
    {
//    printf("System command error:\n);
    }

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/scanPAV_clean align.dat align.clean > try.out",bindir);
    if(system(syscmd) == -1)
    {
//      printf("System command error:\n);
    }

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/scanPAV_process align.clean | egrep match:  > align.match",bindir);
    if(system(syscmd) == -1)
    {
//      printf("System command error:\n);
    }

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"%s/scanPAV_seqout presence.fasta align.match presence-absence.fasta > try.out",bindir);
    if(system(syscmd) == -1)
    {
//      printf("System command error:\n);
    }

    memset(syscmd,'\0',2000);
    sprintf(syscmd,"mv presence-absence.fasta %s",file_scanPAV);
    if(system(syscmd) == -1)
    {
//      printf("System command error:\n);
    }

    if(n_debug == 0)
    {
      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf * > /dev/null");
      if(system(syscmd) == -1)
      {
//      printf("System command error:\n);
      }

      chdir(tmpdir);

      memset(syscmd,'\0',2000);
      sprintf(syscmd,"rm -rf %s > /dev/null",workdir);
      if(system(syscmd) == -1)
      {
//      printf("System command error:\n);
      }
    }
    return EXIT_SUCCESS;

}
/* end of the main */



/*   subroutine to sort out read pairs    */
/* =============================== */
void Read_Index(int nSeq, char *namefile)
/* =============================== */
{
     int i,j,nseq;
     int i_reads,n_reads,c_reads,insertsize=0;
     FILE *namef;
     char *ptr;
     char line[500];
     char **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);


     if((namef = fopen(namefile,"r")) == NULL)
     {
       printf("ERROR main:: reads group file \n");
       exit(1);
     }
     nseq = 0;
     i_reads = 0;
     n_reads = 0;
     c_reads = 0;
     while(!feof(namef))
     {
       if(fgets(line,500,namef) == NULL)
            printf("Data input file problem!\n");
       if(feof(namef)) break;
       nseq++;
     }
     fclose(namef);

     nseq = 2*nseq;
     if((ctg_head = (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_head\n");
       exit(1);
     }
     if((ctg_list = (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ctg_list\n");
       exit(1);
     }
     if((core_list = (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - core_list\n");
       exit(1);
     }
     if((read2contig = (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - read2contig\n");
       exit(1);
     }

     nseq = nseq*3;
     S_Name=cmatrix(0,nseq,0,Max_N_NameBase);
     if((insert_siz = (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - insert\n");
       exit(1);
     }
     if((insert_dev = (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - insert\n");
       exit(1);
     }
 
     if((namef = fopen(namefile,"r")) == NULL)
     {
       printf("ERROR main:: reads group file \n");
       exit(1);
     }

     j = 0;
     insertsize = 0;
     while(!feof(namef))
     {
       int nPair=0,len=0;
       char line2[500],line3[500],base[500];

       if(fgets(line,500,namef) == NULL)
            printf("Data input file problem!\n");
       if(feof(namef)) break;
       strcpy(line2,line);
       strcpy(line3,line);
       insertsize = 0;
       if((strncmp(line,"readnames",9))==0)
       {
         i=0;
         c_reads = 0;
         for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
         {
            if(i==5)
            {
                memset(base,'\0',500);
//                len=strlen(ptr);
//                strncpy(base,ptr,len-1);
                strcat(base,ptr);
                c_reads = atoi(base);
            }
         }
//       printf("creads: %d %d\n",c_reads,n_reads);
         if(n_group>0)
           ctg_list[n_group-1]=n_reads;
         n_group++;
         n_reads = 0;
       }
       else
       {      
         for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),nPair++)
         {
         }
         i=0;
         for(ptr=strtok(line3," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
         {
            if(nPair>1)
            {
              if(i==(nPair-2))
              {
                memset(base,'\0',500);
                strcat(base,ptr);
                insertsize = atoi(base);
              }
            }
         }
         i=0;
         j=0;
         for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
         {
            if(i==0)
            {

              len=strlen(ptr);
              if(nPair==1)
                strncpy(S_Name[i_reads],ptr,len-1);
              else
                strncpy(S_Name[i_reads],ptr,len);
//       printf("reads: %d %d %d %d %s\n",j,i_reads,insertsize,c_reads,S_Name[i_reads]);
              i_reads++;
              j++;
            }
//            else if(insertsize<50000)
            else if((insertsize<50000)&&(c_reads>4))
            {
              if(i==1)
              {
                len=strlen(ptr);
                strncpy(S_Name[i_reads],ptr,len);
                i_reads++;
                j++;
              }
              else if((i==2)&&(i<(nPair-2)))
              {
                len=strlen(ptr);
                strncpy(S_Name[i_reads],ptr,len);
                i_reads++;
                j++;
              }
/*              else if(i==(nPair-2))
              {
                memset(base,'\0',500);
                strcat(base,ptr);
                for(k=0;k<j;k++)
                   insert_siz[i_reads-k-1] = atoi(base);
              }
              else if(i==(nPair-1))
              {
                memset(base,'\0',500);
                strcat(base,ptr);
                for(k=0;k<j;k++)
                   insert_dev[i_reads-k-1] = atoi(base);
              }   */
            }
         }
         n_reads = n_reads+j;
       }
       c_reads++;
     }
     fclose(namef);
     printf("contig: %d %d\n",n_reads,i_reads);
     ctg_list[n_group-1]=n_reads;
     ctg_head[0]=0;
     for(i=1;i<n_group;i++)
        ctg_head[i]=ctg_head[i-1]+ctg_list[i-1];
      
}

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

/*   to swap the string arrays           */
/* ============================================= */
void s_swap2(char **Pair_Name, int i, int j)
/* ============================================= */
{
     char temp[Max_N_NameBase2];

     strcpy(temp,Pair_Name[i]);
     strcpy(Pair_Name[i],Pair_Name[j]);
     strcpy(Pair_Name[j],temp);
}


/*   to sort the string array in order          */
/* ============================================= */
void ArraySort_String2(int n, char **Pair_Name, int *brr)
/* ============================================= */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int temp,MIN=7;
     char p[Max_N_NameBase2];

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
          s_swap2(Pair_Name,k,m+1);
          SWAP(brr[k],brr[m+1]);

          if(strcmp(Pair_Name[m],Pair_Name[ir])>0)
          {
            s_swap2(Pair_Name,m,ir);
            SWAP(brr[m],brr[ir]);
          }

          if(strcmp(Pair_Name[m+1],Pair_Name[ir])>0)
          {
            s_swap2(Pair_Name,m+1,ir);
            SWAP(brr[m+1],brr[ir]);
          }

          if(strcmp(Pair_Name[m],Pair_Name[m+1])>0)
          {
            s_swap2(Pair_Name,m,m+1);
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
             s_swap2(Pair_Name,i,j);
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

