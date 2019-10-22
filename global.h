/*
 * global.h
 */
#ifndef _GLOBAL_
#define _GLOBAL_

#define FASTA_MAXLINE 5000
#define MAXNAME 30
#define DISTMAX 30.0
#define MINBL 0.00001
#define MAXNAME 30
#define MAXNUMBEROFCLUSTERS 100
#define MAXNUMBEROFKSEQS 1000

typedef struct node{
	int down;
	int up[2];
	double bl;
	int nd;
	char* name;
	int nodeToCut;
}node;

typedef struct Options{
	char output_directory[2000];
	char fasta[1000];
	char taxonomy[1000];
	int number_of_clusters;
	int number_of_kseqs;
}Options;

#endif
