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
#define MAXNUMBEROFCLUSTERS 500
#define MAXNUMBEROFKSEQS 10000
#define MAXNUMINCLUSTER 10000
#define PADDING 500

typedef struct node{
	int down;
	int up[2];
	double bl;
	int nd;
	char* name;
	int nodeToCut;
	int depth;
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
}Options;

typedef struct nw_alignment{
	nw_aligner_t *nw;
	alignment_t *aln;
	scoring_t *scoring;
}nw_alignment;

typedef struct resultsStruct{
	char** accession;
	char** assigned;
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
	char*** clusterNames;
	char*** clusterSeqs;
	int* chooseK;
	int* fasta_specs;
	char** seqNames;
	char** sequences;
	char** taxonomy;
	int threadnumber;
	//struct hashmap seqsToCompare;
	//struct hashmap assignedSeqs;
	char** assignedSeqs;
	resultsStruct *str;
	int numAssigned;
	char buffer[9999];
}mystruct;

extern char*** clusters;
extern struct hashmap map;
#endif
