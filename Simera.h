/* This software is copyright © by Ben Nichols 2016. */
/* Permission is granted for anyone to copy, use or modify this software for the purposes of research or education provided that this copyright notice is retained and note is given of any modifications. */

//Simera Header File

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

#define TRUE  1
#define FALSE 0

#define MATCH 0
#define DEL	  1
#define INS   2

typedef struct s_Data
{	
	char** aszSeq;
	
	char** aszFrag;
	
	char** aszID;
	
	int    nMaxLen;
	
	int    nSeq;
	
	int	   nNextID;
	
	int*   anID;
	
	int*   anParentA;
	
	int*   anParentB;
	
	int*   anLen;
	
	int*   anBreak;
	
	int64_t*   anAbund;
	
	int64_t    nTotAbund;
	
	int64_t	 nMaxAbund;
	
	int	   nSize;
	
	int*   anNmera;
	
	double* adWeight;
	
	double* adProb;
	
} t_Data;

#define MAX_LINE_LENGTH 65536
#define INIT_SIZE 1024

void PCR(gsl_rng *ptRNG, t_Data *ptData, t_Data *ptPrimer, t_Data *ptChimeras, double dLambda, int64_t nPrimAbund, int nRounds, int nSeed, int nTotal);

void Chimera(t_Data *ptData, t_Data *ptFrag, int64_t pool, char* Seq);

void initData(t_Data *ptData);

void freeData(t_Data *ptData);

void incSeq(int id, int64_t inc, t_Data *ptData);

int existSeq(char* seq, t_Data *ptData);

void addSeq(char* seq, char* fragment, char* name, int64_t abund, int chimera, int parentA, int parentB, int frag, int bp, double weight, double prob, t_Data *ptData);

void addSeq2(char* seq, char* fragment, char* name, int64_t abund, int chimera, int parentA, int parentB, int frag, int bp, double weight, double prob, t_Data *ptData);

void emptyData(t_Data *ptData);

void removeSeq(size_t nID, int64_t abund, t_Data *ptData);

void refreshData(t_Data *ptData);

void printData(t_Data *ptData);

void printChim(t_Data *ptData);

void printGood(t_Data *ptData);

void printParents(t_Data *ptData);

void printBP(t_Data *ptData);

void printAbund(t_Data *ptData);

char* substring(char* str, size_t begin, size_t len);

void readData(char* szInputFile, t_Data *ptData);

void SimsFWD(char *SeqA, char *SeqB, int *nPos, int *nMax, int nBase);

void SimsREV(char *SeqA, char *SeqB, int *nPos, int *nMax, int nBase);

double SimsPrimerFWD(char *SeqA, char *SeqB);

double SimsPrimerREV(char *SeqA, char *SeqB);

void Sample(gsl_rng *ptRNG, t_Data *ptIn, t_Data *ptOut, int64_t nChoose);

void writeData(t_Data *ptData, char* fileName);

void writeChim(t_Data *ptData, char* fileName);

void writeGood(t_Data *ptData, char* fileName);

void writeParents(t_Data *ptData, char* fileName);

void writeBP(t_Data *ptData, char* fileName);

void writeAbund(t_Data *ptData, char* fileName);

void writeInputInfo(int nSeed, char* szInput, char* szVersion, int nRounds, double dLambda, int nSamp, char* PrimerFWD, char* PrimerREV, char* fileName);

void writeSummary(t_Data *ptData, char* fileName);

char* Rev_Comp(char* SeqA);