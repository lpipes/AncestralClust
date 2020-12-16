/*
 * global.h
 */
#include "needleman_wunsch.h"
#include "hashmap.h"
#ifndef _GLOBAL_
#define _GLOBAL_

#define FASTA_MAXLINE 6000
#define MAXNAME 30
#define DISTMAX 30.0
#define MINBL 0.00001
#define MAXBL 2.0
#define MAXNUMBEROFCLUSTERS 500
#define MAXNUMBEROFKSEQS 10000
#define MAXNUMINCLUSTER 10000
#define PADDING 500
#define MIN_REQ_SSIZE 81920
#define STATESPACE 20 /*number of categories in approximation of gamma distribution for Ne. Must be at least 4 because some of the memory is used for the nucleotide model*/
#define NUMCAT 1/*number of categories in the discretization of the gamma for the nucleotide substituion model*/
#define MAXNUMBEROFINDINSPECIES 500 /*maximum number of individuals belonging to a species*/
#define type_of_PP double
//#define MIN_REQ_SSIZE 83886080
typedef struct node{
	int down;
	int up[2];
	double bl;
	int nd;
	char* name;
	int nodeToCut;
	int depth;
	double distance;
	double** likenc;
	double** posteriornc;
}node;

typedef struct Options{
	char output_directory[2000];
	char output_file[2000];
	char fastaToAssign[2000];
	char fasta[1000];
	char taxonomy[1000];
	int number_of_clusters;
	int number_of_kseqs;
	int slash;
	int default_directory;
	int numthreads;
	int numberOfLinesToRead;
	int largest_cluster;
	int hasTaxFile;
	int clstr_format;
	int output_fasta;
}Options;

typedef struct nw_alignment{
	nw_aligner_t *nw;
	alignment_t *aln;
	scoring_t *scoring;
}nw_alignment;

typedef struct resultsStruct{
	char** accession;
	//char** assigned;
	int numassigned;
	int* clusterNumber;
	double average;
	char** savedForNewClusters;
	int number_of_clusters;
	int* clusterSizes;
	char buffer[99999];
}resultsStruct;

typedef struct mystruct{
	int start;
	int end;
	int num_threads;
	double average;
	int number_of_clusters;
	int number_of_kseqs;
	int largest_cluster;
	int* clusterSize;
	//char*** clusterNames;
	char*** clusterSeqs;
	int* chooseK;
	nw_alignment* nw_struct;
	int* fasta_specs;
	char** seqNames;
	char** sequences;
	char** taxonomy;
	int threadnumber;
	//struct hashmap seqsToCompare;
	//struct hashmap assignedSeqs;
	//char** assignedSeqs;
	resultsStruct *str;
	int numAssigned;
	char buffer[9999];
}mystruct;

typedef struct msa{
	struct msa_seq** sequences;
	int** sip;
	int* nsip;
	int* plen;
	int numseq;
	int num_profiles;
	int alloc_numseq;
	int aligned;
	int letter_freq[128];
	int L;
}msa;

typedef struct msa_seq{
	char* name;
	char* seq;
	uint8_t* s;
	int* gaps;
	int len;
	int alloc_len;
}msa_seq;

extern char*** clusters;
extern struct hashmap map;
extern double LRVECnc[4][4], RRVECnc[4][4], RRVALnc[4], PMATnc[2][4][5];
extern double LRVEC[STATESPACE][STATESPACE], RRVEC[STATESPACE][STATESPACE], RRVAL[STATESPACE], PMAT1[STATESPACE][STATESPACE], PMAT2[STATESPACE][STATESPACE];
extern double parameters[10];
extern int COUNT2;//counting the number of tiems the likelihood function is called.
extern int COUNT; //this one counts how many times the likelihood function has been called
extern double *statevector, *UFCnc, *localpi, **templike_nc;
extern double Logfactorial[MAXNUMBEROFINDINSPECIES];
extern double parameters[10];
extern node** treeArr;
#endif
