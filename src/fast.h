
typedef struct
{
        char *name;
        char *name2;
        char *path;
        char *SCFname;
        int  length;
        char *data;
        char *qual;
        int  finished;
} fasta;

typedef struct
{
        char *name;
        char *name2;
        char *path;
        char *SCFname;
        int  length;
        char *data;
        char *qual;
        int  finished;
} fastq;


fasta *decodeFasta (char *fname, int *nContigs);
fasta *decodeFastq (char *fname, int *nContigs, long *tB, int qT);
void fastaLC (fasta *seg, int nSeg);
void fastaUC (fasta *seg, int nSeg);
