
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

#define B64_long long int
//#define B64_long long long int
fasta *decodeFastq (char *fname, int *nContigs, B64_long *tB, char* pdata, B64_long Size_pdata,fasta *segg);
int extractFastq(char *fname, char *pdata, B64_long Size_pdata);
fasta *splitFastq ( fasta *iseg, int inSeg, int tB, int *nReads, int length, int step);
void fastaLC (fasta *seg, int nSeg);
void fastaUC (fasta *seg, int nSeg);

