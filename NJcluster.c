#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include "needleman_wunsch.h"
#include "global.h"
#include "hashmap.h"

void setNumSeq(FILE* fasta, int* fasta_specs){
	char buffer[FASTA_MAXLINE];
	int numSeq = 0;
	int maxLine = 0;
	int maxName = 0;
	int i=0;
	while( fgets(buffer,FASTA_MAXLINE,fasta) != NULL ){
		if ( buffer[0] == '>' ){
			numSeq++;
			int nameLength = 0;
			for(i=0; buffer[i]!='\0'; i++){
				nameLength++;
			}
			if (nameLength > maxName ){
				maxName = nameLength;
			}
		}else{
			int seqLength = 0;
			for(i=0; buffer[i]!='\0'; i++){
				seqLength++;
			}
			if ( seqLength > maxLine ){
				maxLine = seqLength;
			}
		}
	}
	fasta_specs[0] = numSeq; //set number of sequences
	fasta_specs[1] = maxLine; //set maxline
	fasta_specs[2] = maxName; //set max name length
}

/*void allocateMemForArr(int* fasta_specs, char** seqNames, char** sequences){
	int i;
	seqNames = (char **)malloc(fasta_specs[0]*sizeof(char *));
	sequences = (char **)malloc(fasta_specs[0]*sizeof(char *));
	for(i=0; i<fasta_specs[0]; i++){
		seqNames[i]=(char *)malloc(fasta_specs[2]*sizeof(char));
		sequences[i]=(char *)malloc(fasta_specs[1]*sizeof(char));
	}
}*/

void readInFasta(FILE* fasta, char** seqNames, char** sequences, struct hashmap map){
	char buffer[FASTA_MAXLINE];
	int i=0;
	int j=0;
	while( fgets(buffer,FASTA_MAXLINE,fasta) != NULL ){
		if ( buffer[0] == '>' ){
			for(j=1;buffer[j]!='\n';j++){
				seqNames[i][j-1]=buffer[j];
			}
			seqNames[i][j-1]='\0';
		}else{
			for(j=0;buffer[j]!='\n';j++){
				sequences[i][j]=buffer[j];
			}
			sequences[i][j-1]='\0';
			hashmap_put(&map,seqNames[i],sequences[i]);
			i++;
		}
	}
}

int generateRandom(int numSeq){
	int i;
	int num = (rand() % (numSeq+1));
	return num;
}

void print_distance_matrix(double** distMat, char*** clusters, int whichCluster, int clusterSize){
	int i,j;
	for(i=0; i<clusterSize; i++){
		printf("row %i %s\t",i,clusters[whichCluster][i]);
		for(j=0; j<clusterSize; j++){
			printf("%lf\t",distMat[i][j]);
		}
		printf("\n");
	}
}
int populate_DATA(char* query1, char* query2, int** DATA, int alignment_length, int* mult){
	
	int i;
	//printf("query1: %s\n",query1);
	for(i=0; i<alignment_length; i++){
		//printf("%d %c",i,query1[i]);
		if (query1[i] == 'a' || query1[i] == 'A'){
			DATA[0][i]=0;
		}else if ( query1[i] == 'c' || query1[i] == 'C'){
			DATA[0][i]=1;
		}else if ( query1[i] == 'g' || query1[i] == 'G'){
			DATA[0][i]=2;
		}else if ( query1[i] == 't' || query1[i] == 'u' || query1[i] == 'T' || query1[i] == 'U'){
			DATA[0][i]=3;
		}else if ( query1[i] == 'n' || query1[i] == '-' || query1[i] == 'N'){
			DATA[0][i]=-1;
		}else if ( query1[i] == '\n'){
		}else{
			printf("\nBAD BASE in sequence 1 base %i: %c\n",i,query1[i]);  exit(-1);
		}
	}
	//printf("query2: %s\n",query2);
	for(i=0; i<alignment_length; i++){
		if (query2[i] == 'a' || query2[i] == 'A'){
			DATA[1][i]=0;
		}else if ( query2[i] == 'c' || query2[i] == 'C'){
			DATA[1][i]=1;
		}else if ( query2[i] == 'g' || query2[i] == 'G'){
			DATA[1][i]=2;
		}else if ( query2[i] == 't' || query2[i] == 'u' || query2[i] == 'T' || query2[i] == 'U'){
			DATA[1][i]=3;
		}else if ( query2[i] == 'n' || query2[i] == '-' || query2[i] == 'N'){
			DATA[1][i]=-1;
		}else if (query2[i] == '\n'){
		}else{
			printf("\nBAD BASE in sequence 2 base %i: %c\n",i,query2[i]);  exit(-1);
		}
	}
	alignment_length = sortseq(alignment_length, 2, DATA, mult);
	/*for(i=0; i<alignment_length; i++){
		printf("%d",DATA[0][i]);
	}
	printf("\n");
	for(i=0; i<alignment_length; i++){
		printf("%d",DATA[1][i]);
	}
	printf("\n");*/
}
void Get_dist_JC(int alignment_length, double **M, int** DATA, int* mult, int index1, int index2){
	int i;
	double rawd, dmax;
	int numdiff, numrealsites;
	numdiff=0;
	numrealsites=0;
	dmax = 0.75*(1.0 - exp(-4.0*DISTMAX/3.0)); //precalculate instead
	for( i=0; i<alignment_length; i++){
		if (DATA[0][i]> -1 && DATA[1][i] > -1){
			numrealsites=numrealsites+mult[i];
			if (DATA[0][i]!=DATA[1][i]){
				numdiff=numdiff+mult[i];
			}
		}
	}
	if (numrealsites < 1) {printf("Wacky distance in function 'Get_dist_JC' between sequence %i and %i\n",index1,index2); exit(-1);}
	if (numdiff==0) {
		M[index1][index2] = 0.0;
	}else{
		if (numdiff == 0) {
			M[index1][index2] = 0.0;
		}else{
			rawd = (double)numdiff/(double)numrealsites;
			if (rawd >= dmax) M[index1][index2] = DISTMAX; //printf("upper bound hit on distances in Get_dist_JC (%i %i)\n",numdiff,numrealsites);}'
			else M[index1][index2]= - 0.75*log(1.0 -4.0*rawd/3.0);
		}
	}
}
int sortseq(int alignment_length, int number_of_sequences, int** seq, int* mult){
	int i,j;
	for(i=0; i<alignment_length; i++){
		mult[i]=1;
	}
	for(i=0; i<alignment_length; i++){
		for(j=i+1; j<alignment_length; j++){
			if(seqid(i,j,number_of_sequences,seq)==1){
				elimfrom(j,alignment_length,number_of_sequences,seq);
				j--;
				alignment_length--;
				mult[i]++;
			}
		}
	}
	return alignment_length;
}
int seqid(int i, int j, int number_of_sequences, int** seq){
	int n;
	for(n=0; n<number_of_sequences; n++){
		if( seq[n][i] != seq[n][j]){
			return 0;
		}
	}
	return 1;
}
void elimfrom(int site, int alignment_length, int number_of_sequences, int **seq){
	int i,j;
	for(i=0; i<number_of_sequences; i++){
		for(j=site;j<alignment_length-1;j++){
			seq[i][j] = seq[i][j+1];
		}
	}
}
int NJ(node** tree, double** distMat,char*** clusters, int clusterSize,int whichTree){
	double *r, minval, D, u1;
	int child1, child2, i, j, n, min, pair[2], newnode, *nodenum;
	r = malloc(clusterSize*(sizeof(double)));
	nodenum = malloc(clusterSize*(sizeof(int)));
	n=clusterSize;
	int* index = (int *)malloc(clusterSize*sizeof(int));
	for(i=0; i<clusterSize; i++){
		nodenum[i]=i;
		tree[whichTree][i].up[0]=tree[whichTree][i].up[1]=-1;
		index[i]=i;
	}
	//recurse until tree is found
	do{
		for (i=0; i<n; i++){
			r[i]=0;
			for (j=0; j<i; j++){
				r[i] = r[i] + distMat[j][i];
			}
			for (j=i+1; j<n; j++){
				r[i] = r[i] + distMat[i][j];
			}
			r[i] = r[i]/(double)(n-2);
		}
		//find minival value of D
		minval = DISTMAX;
		for (i=0; i<n; i++){
			for (j=i+1; j<n; j++){
				if ((D=distMat[i][j] -(r[i]+r[j])) < minval){
					minval=D;
					pair[0]=i;
					pair[1]=j;
				}
			}
		}
		if (minval==DISTMAX) {printf("error in 'NJ' (%lf, %i)",minval,n); exit(-1);}
		newnode=2*clusterSize-n;
		child1 = nodenum[pair[0]];
		child2 = nodenum[pair[1]];
		tree[whichTree][newnode].up[0]=child1;
		tree[whichTree][newnode].up[1]=child2;
		tree[whichTree][child1].down=newnode;
		tree[whichTree][child2].down=newnode;
		if (tree[whichTree][child1].up[0]==-1){
			strcpy(tree[whichTree][child1].name,clusters[whichTree][index[pair[0]]]);
		}
		if (tree[whichTree][child2].up[1]==-1){
			strcpy(tree[whichTree][child2].name,clusters[whichTree][index[pair[1]]]);
		}
		// HERE WE DISALLOW NEGATIVE BRANCHLENGTHS!  IS THIS THE BEST THING TO DO?
		if ((u1 = (distMat[pair[0]][pair[1]] +r[pair[0]]-r[pair[1]])/2.0) < MINBL){
			tree[whichTree][child1].bl = MINBL;
		}else{
			tree[whichTree][child1].bl = u1;
		}
		if ((tree[whichTree][child2].bl = distMat[pair[0]][pair[1]]-u1) < MINBL)/*(DM[pair[0]][pair[1]] +r[pair[1]]-r[pair[0]])/2.0;*/{
			tree[whichTree][child2].bl = MINBL;
		}
		for (i=0; i<n; i++){
			if (i != pair[0] && i != pair[1]){
				if (i<pair[0]) {
					distMat[i][n] = (distMat[i][pair[0]]+distMat[i][pair[1]] - distMat[pair[0]][pair[1]])/2.0;
				}else if (i<pair[1]){
					distMat[i][n] = (distMat[pair[0]][i]+distMat[i][pair[1]] - distMat[pair[0]][pair[1]])/2.0;
				}else{
					distMat[i][n] = (distMat[pair[0]][i]+distMat[pair[1]][i] - distMat[pair[0]][pair[1]])/2.0;
				}
			}
		}
		for (i=pair[0]; i<pair[1]-1;i++){
			nodenum[i]=nodenum[i+1];
		}
		for (i=pair[1]-1; i<n-2;i++){
			nodenum[i]=nodenum[i+2];
		}
		nodenum[n-2]=newnode;
		updatematrix(distMat, pair[1], n+1);
		updatematrix(distMat, pair[0], n);
		updateChooseK(index, pair[1], clusterSize);
		updateChooseK(index, pair[0], clusterSize);
		n--;
	} while (n>2);
	newnode=2*clusterSize-2;
	child1 = nodenum[0];
	child2 = nodenum[1];
	tree[whichTree][newnode].up[0]=child1;
	tree[whichTree][newnode].up[1]=child2;
	tree[whichTree][newnode].down=-1;
	tree[whichTree][newnode].bl=-1.0;
	tree[whichTree][child1].down=newnode;
	tree[whichTree][child2].down=newnode;
	if (tree[whichTree][child1].up[0]==-1 && strcmp(tree[whichTree][child1].name,"internal")==0){
		strcpy(tree[whichTree][child1].name,clusters[whichTree][index[0]]);
	}
	if (tree[whichTree][child2].up[0]==-1 && strcmp(tree[whichTree][child2].name,"internal")==0){
		strcpy(tree[whichTree][child1].name,clusters[whichTree][index[0]]);
	}	
	if (distMat[0][1] > MINBL*2.0)  // HERE WE DISALLOW NEGATIVE BRANCHLENGTHS!  IS THIS THE BEST THING TO DO?
		tree[whichTree][child1].bl = tree[whichTree][child2].bl = distMat[0][1]/2.0;
	else tree[whichTree][child1].bl = tree[whichTree][child2].bl = MINBL;
	//for (i=0; i<opt.number_of_kseqs+1; i++){
	//	free(distMat[i]);
	//}
	//free(distMat);
	free(nodenum);
	free(r);
	free(index);
	return(newnode);
}
void updateChooseK(int* chooseK, int index, int clusterSize){
	int i;
	int tmp[clusterSize];
	for(i=0; i<clusterSize; i++){
		tmp[i] = chooseK[i];
	}
	for(i=0; i<clusterSize-1; i++){
		if ( i>=index ){
			chooseK[i]=tmp[i+1];
		}
	}
}
//this algorithm takes an nxn upper triangular matrix and converts it to an (n-1)x(n-1)
//upper triangular matrix by deleting row k and column k
void updatematrix(double **matrix, int k, int n){
	int i,j;
	for (i=0; i<n; i++){
		for (j=i+1; j<n; j++){
			if (i>=k && j>=k){
				matrix[i][j]=matrix[i+1][j+1];
			}else if (i>=k) {
				matrix[i][j]=matrix[i+1][j];
			}else if (j>=k) {
				matrix[i][j]=matrix[i][j+1];
			}
		}
	}
}
void printtree(node** tree, int whichTree, int clusterSize){
	int i;
	for(i=0; i<2*clusterSize-1; i++){
		printf("%i: up: %i %i, down: %i, bl: %f, nd: %d, name: %s\n",i,tree[whichTree][i].up[0],tree[whichTree][i].up[1],tree[whichTree][i].down,tree[whichTree][i].bl/*totsites*/,tree[whichTree][i].nd,tree[whichTree][i].name);
	}
}
void sortArray(double* branchLengths, int* indexArray, Options opt){
	int i,j, tmp2;
	double tmp;
	//Bubble Sort
	for(i=0; i<2*opt.number_of_kseqs-1; i++){
		indexArray[i]=i;
	}
	for(i=0; i<2*opt.number_of_kseqs; i++){
		for(j=i+1;j<2*opt.number_of_kseqs-1;j++){
			if(branchLengths[i] < branchLengths[j]){
				tmp=branchLengths[i];
				branchLengths[i] = branchLengths[j];
				branchLengths[j] = tmp;
				tmp2=indexArray[i];
				indexArray[i] = indexArray[j];
				indexArray[j]=tmp2;
			}
		}
	}
	/*for(i=0; i<2*opt.number_of_kseqs-1; i++){
		printf("branchLengths[%d]: %lf\n",i,branchLengths[i]);
	}*/
}
int get_number_descendants(node** tree, int node, int whichTree){
	if (tree[whichTree][node].up[0]==-1)  return (tree[whichTree][node].nd=1);
	else return (tree[whichTree][node].nd=(get_number_descendants(tree,tree[whichTree][node].up[0],whichTree)+get_number_descendants(tree,tree[whichTree][node].up[1],whichTree)));
}
void printdescendants(node** tree, int node, char*** clusters, int clusterNumber, int whichTree, Options opt){
	int child1 = tree[whichTree][node].up[0];
	int child2 = tree[whichTree][node].up[1];
	int i;
	if (tree[whichTree][node].nodeToCut==1){
		return;
	}
	if (tree[whichTree][node].up[0]==-1 && tree[whichTree][node].up[1]==-1){
		int count=0;
		for(i=opt.number_of_kseqs-1; i>=0; i--){
			if (clusters[clusterNumber][i][0]=='\0'){
				count=i;
			}
		}
		printf("count %d node %d: %s\n",count,node,tree[whichTree][node].name);
		strcpy(clusters[clusterNumber][count],tree[whichTree][node].name);
	}else{
		printdescendants(tree,child1,clusters,clusterNumber,whichTree,opt);
		printdescendants(tree,child2,clusters,clusterNumber,whichTree,opt);
	}
}
void clearDescendants(node** tree, int node, int whichTree){
	int child1 = tree[whichTree][node].up[0];
	int child2 = tree[whichTree][node].up[1];
	if (tree[whichTree][node].up[0] == -1 && tree[whichTree][node].up[1] == -1){
		tree[whichTree][node].nd = 0;
		return;
	}else{
		tree[whichTree][node].nd = 0;
		clearDescendants(tree,child1,whichTree);
		clearDescendants(tree,child2,whichTree);
	}
}
int countNumInCluster(char*** clusters, int index, Options opt){
	int i,j;
	for(i=0; i<opt.number_of_kseqs; i++){
		if(clusters[index][i][0]!='\0'){
			j=i;
		}
	}
	return j+1;
}
void allocateMemForTreeArr(int numberOfClusters, int* clusterSize, node** treeArr, Options opt){
	int i,j;
	treeArr=(node **)malloc(numberOfClusters*sizeof(node *));
	for(i=0; i<numberOfClusters; i++){
		treeArr[i]=malloc((2*opt.number_of_kseqs-1)*sizeof(node));
		for(j=0; j<2*opt.number_of_kseqs-1; j++){
			treeArr[i][j].name = (char *)malloc(MAXNAME*sizeof(char));
		}
	}
}
void createDistMat(char** sequences, double** distMat, int clusterSize, int longestSeq){
	int i,j;
	nw_aligner_t *nw;
	alignment_t *aln;
	scoring_t *scoring;
	nw=needleman_wunsch_new();
	aln=alignment_create(longestSeq);
	scoring=malloc(sizeof(scoring_t));
	int match=2;
        int mismatch=-1;
        int gap_open=-3;
        int gap_extend=-1;
        int num_of_mismatches=0;
        int num_of_indels = 0;
        bool no_start_gap_penalty=true;
        bool no_end_gap_penalty=true;
        bool no_gaps_in_a = false;
        bool no_gaps_in_b = false;
        bool no_mismatches = false;
        bool case_sensitive=false;
        scoring_init(scoring, match, mismatch, gap_open, gap_extend, no_start_gap_penalty, no_end_gap_penalty, no_gaps_in_a, no_gaps_in_b, no_mismatches, case_sensitive);
	for (i=0; i<clusterSize; i++){
		for (j=i+1; j<clusterSize; j++){
			needleman_wunsch_align(sequences[i],sequences[j], scoring, nw, aln);
			int alignment_length = strlen(aln->result_a);
			//printf("alignment length: %d\n",alignment_length);
			//printf("RESULT_A: %s\n",aln->result_a);
			//printf("RESULT_B: %s\n",aln->result_b);
			int** DATA = (int **)malloc(2*sizeof(int *));
			int k;
			for(k=0; k<2; k++){
				DATA[k] = (int *)malloc(alignment_length*sizeof(int));
			}
			int* mult = (int *)malloc(alignment_length*sizeof(int));
			alignment_length = populate_DATA(aln->result_a,aln->result_b,DATA,alignment_length,mult);
			Get_dist_JC(alignment_length,distMat,DATA,mult,i,j);
			free(mult);
			free(DATA[0]);
			free(DATA[1]);
			free(DATA);
		}
	}
	free(scoring);
	alignment_free(aln);
	needleman_wunsch_free(nw);
}
void updateNumberOfDescendants(node** tree, int node, int descendants, int whichTree){
	int parent = tree[whichTree][node].down;
	if ( tree[whichTree][node].down == -1 ){
		return;
	}
	if ( tree[whichTree][node].nd == 0 ){
		return;
	}
	tree[whichTree][node].nd=tree[whichTree][node].nd-descendants;
	updateNumberOfDescendants(tree,parent,descendants,whichTree);
}
double calculateAverageDist(char** seqsA, char** seqsB, int sizeA, int sizeB, int longestSeq, double*** distMat, int count, double* calculateAvg){
	int i,j;
	nw_aligner_t *nw;
	alignment_t *aln;
	scoring_t *scoring;
	nw=needleman_wunsch_new();
	aln=alignment_create(longestSeq);
	scoring=malloc(sizeof(scoring_t));
	int match=2;
        int mismatch=-1;
        int gap_open=-3;
        int gap_extend=-1;
        int num_of_mismatches=0;
        int num_of_indels = 0;
        bool no_start_gap_penalty=true;
        bool no_end_gap_penalty=true;
        bool no_gaps_in_a = false;
        bool no_gaps_in_b = false;
        bool no_mismatches = false;
        bool case_sensitive=false;
        scoring_init(scoring, match, mismatch, gap_open, gap_extend, no_start_gap_penalty, no_end_gap_penalty, no_gaps_in_a, no_gaps_in_b, no_mismatches, case_sensitive);
	for(i=0; i<sizeA; i++){
		for(j=0; j<sizeB; j++){
			needleman_wunsch_align(seqsA[i],seqsB[j],scoring,nw,aln);
			int alignment_length = strlen(aln->result_a);
			//printf("alignment length: %d\n",alignment_length);
			//printf("RESULT_A: %s\n",aln->result_a);
			//printf("RESULT_B: %s\n",aln->result_b);
			int** DATA = (int **)malloc(2*sizeof(int *));
			int k;
			for(k=0; k<2; k++){
				DATA[k] = (int *)malloc(alignment_length*sizeof(int));
			}
			int* mult = (int *)malloc(alignment_length*sizeof(int));
			alignment_length = populate_DATA(aln->result_a,aln->result_b,DATA,alignment_length,mult);
			Get_dist_JC(alignment_length,distMat[count],DATA,mult,i,j);
			free(mult);
			free(DATA[0]);
			free(DATA[1]);
			free(DATA);
		}
	}
	double totalDist=calculateAvg[0];
	double numberOfPairs=calculateAvg[1];
	for(i=0; i<sizeA; i++){
		for(j=0; j<sizeB; j++){
			totalDist = totalDist + distMat[count][i][j];
			numberOfPairs++;
		}
	}
	calculateAvg[0]=totalDist;
	calculateAvg[1]=numberOfPairs;
	free(scoring);
	alignment_free(aln);
	needleman_wunsch_free(nw);
}
double findShortestDist(char** clusterSeqs, char* seq, int clusterSize, int longestSeq, nw_alignment* nw_struct, double** distMat, int** DATA, int* mult){
	int i,j;
	/*nw_aligner_t *nw;
	alignment_t *aln;
	scoring_t *scoring;
	nw=needleman_wunsch_new();
	aln=alignment_create(longestSeq);
	scoring=malloc(sizeof(scoring_t));
	int match=2;
        int mismatch=-1;
        int gap_open=-3;
        int gap_extend=-1;
        int num_of_mismatches=0;
        int num_of_indels = 0;
        bool no_start_gap_penalty=true;
        bool no_end_gap_penalty=true;
        bool no_gaps_in_a = false;
        bool no_gaps_in_b = false;
        bool no_mismatches = false;
        bool case_sensitive=false;
        scoring_init(scoring, match, mismatch, gap_open, gap_extend, no_start_gap_penalty, no_end_gap_penalty, no_gaps_in_a, no_gaps_in_b, no_mismatches, case_sensitive);
	double** distMat = (double**)malloc(clusterSize*sizeof(double *));
	for(i=0;i<clusterSize; i++){
		distMat[i] = (double*)malloc(sizeof(double));
	}
	*/
	for(i=0; i<clusterSize; i++){
		needleman_wunsch_align(clusterSeqs[i],seq,nw_struct->scoring,nw_struct->nw,nw_struct->aln);
		int alignment_length = strlen(nw_struct->aln->result_a);
		//int** DATA = (int **)malloc(2*sizeof(int *));
		//int k;
		//for(k=0; k<2; k++){
		//	DATA[k] = (int *)malloc(alignment_length*sizeof(int));
		//}
		//int* mult = (int *)malloc(alignment_length*sizeof(int));
		alignment_length = populate_DATA(nw_struct->aln->result_a,nw_struct->aln->result_b,DATA,alignment_length,mult);
		Get_dist_JC(alignment_length,distMat,DATA,mult,i,0);
	}
	double shortestDist = 1;
	for(i=0; i<clusterSize; i++){
		//printf("distMat[%d][0]: %lf\n",i,distMat[i][0]);
		if (shortestDist > distMat[i][0]){
			shortestDist = distMat[i][0];
		}
		//free(distMat[i]);
	}
	//free(distMat);
	//printf("shortest distance: %lf\n",shortestDist);
	//free(scoring);
	//alignment_free(aln);
	//needleman_wunsch_free(nw);
	return shortestDist;
}
void makeNewCluster(char*** clusters,int number_of_clusters, char* sequence_to_add, int* chooseK, int seqIndex, int number_in_initial_clust, struct hashmap seqsToCompare){
	strcpy(clusters[number_of_clusters][0],sequence_to_add);
	strcpy(clusters[0][number_in_initial_clust],sequence_to_add);
	hashmap_put(&seqsToCompare,sequence_to_add,1);
	chooseK[number_in_initial_clust] = seqIndex;
}
void addToCluster(char*** clusters, char* sequence_to_add, int clusterNumber, int clusterSize){
	strcpy(clusters[clusterNumber][clusterSize],sequence_to_add);
}
void printClusters(char*** clusters, int number_of_clusters, int* clusterSize, char** seqNames, char** sequences, int number_of_total_seqs, char** taxonomy, Options opt){
	FILE *clusterFile, *clusterTaxFile;
	char fileName[FASTA_MAXLINE];
	char *directory = strdup(opt.output_directory);
	int i, j, k;
	for(i=1; i<number_of_clusters; i++){
		if (opt.slash==0 && opt.default_directory==0){
			snprintf(fileName,FASTA_MAXLINE,"%s/%d.fasta",directory,i);
		}else if (opt.slash==1){
			snprintf(fileName,FASTA_MAXLINE,"%s%d.fasta",directory,i);
		}else if (opt.default_directory==1){
			snprintf(fileName,FASTA_MAXLINE,"%d.fasta",i);
		}
		clusterFile = fopen(fileName, "w");
		if (clusterFile == NULL){ printf("Error opening cluster file!"); exit(1); }
		if (opt.slash==0 && opt.default_directory==0){
			snprintf(fileName,FASTA_MAXLINE, "%s/%d_taxonomy.txt",directory,i);
		}else if (opt.slash==1){
			snprintf(fileName,FASTA_MAXLINE, "%s%d_taxonomy.txt",directory,i);
		}else if (opt.default_directory==1){
			snprintf(fileName,FASTA_MAXLINE, "%d_taxonomy.txt",i);
		}
		clusterTaxFile = fopen(fileName, "w");
		if (clusterTaxFile == NULL){ printf("Error opening cluster taxonomy file!"); exit(1); }
		for(j=0; j<clusterSize[i]; j++){
			for(k=0; k<number_of_total_seqs; k++){
				if (strcmp(clusters[i][j],seqNames[k])==0){
					fprintf(clusterFile,">%s\n",seqNames[k]);
					fprintf(clusterFile,"%s\n",sequences[k]);
					fprintf(clusterTaxFile,"%s\t%s\n",seqNames[k],taxonomy[k]);
				}
			}
		}
		fclose(clusterFile);
		fclose(clusterTaxFile);
	}
}
void readInTaxFile(FILE* taxonomyFile,char** taxonomy){
	char buffer[FASTA_MAXLINE];
	char *accession, *lineTaxonomy;
	int i=0;
	while( fgets(buffer,FASTA_MAXLINE,taxonomyFile) != NULL){
		accession = strtok(buffer,"\t");
		lineTaxonomy = strtok(NULL,"\n");
		assert(strlen(accession) <= MAXNAME);
		strcpy(taxonomy[i],lineTaxonomy);
		i++;
	}
}
void print_taxonomy(char*** clusters, char** seqNames, char** taxonomy, int clusterNumber, int clusterSize, int* chooseK, Options opt){
	int i, j;
	printf("cluster %d:\n",clusterNumber);
	for(i=0; i<clusterSize; i++){
		for(j=0; j<opt.number_of_kseqs; j++){
			if (strcmp(clusters[clusterNumber][i],seqNames[chooseK[j]])==0){
				printf("%s\t%s\n",seqNames[chooseK[j]],taxonomy[chooseK[j]]);
			}
		}
	}
}
void freeDistMat(double** distMat, int size){
	int i;
	for(i=0; i<size; i++){
		free(distMat[i]);
	}
	free(distMat);
}
void freeDistMatAvg(double*** distMatAvg, int* clusterSize, int number_of_clusters){
	int i, j;
	for(i=0; i<number_of_clusters; i++){
		for(j=0; j<clusterSize[i+1]; j++){
			free(distMatAvg[i][j]);
		}
		free(distMatAvg[i]);
	}
	free(distMatAvg);
}
void freeSeqsInClusterAvg(char*** seqsInClusterAvg, Options opt){
	int i,j;
	for(i=0; i<opt.number_of_clusters-1; i++){
		for(j=0; j<opt.number_of_kseqs; j++){
			free(seqsInClusterAvg[i][j]);
		}
		free(seqsInClusterAvg[i]);
	}
	free(seqsInClusterAvg);
}
void freeTreeArr(node** tree,int* clusterSizes){
	int i,j;
}
void freeClusters(char*** clusters,Options opt){
	int i,j;
	for(i=0; i<MAXNUMBEROFCLUSTERS; i++){
		for(j=0; j<opt.number_of_kseqs; j++){
			free(clusters[i][j]);
		}
		free(clusters[i]);
	}
	free(clusters);
}
void freeSequences(char** sequences, int number_of_seqs){
	int i;
	for(i=0;i<number_of_seqs; i++){
		free(sequences[i]);
	}
	free(sequences);
}
void freeSeqsInCluster(char** seqsInCluster, Options opt){
	int i;
	for(i=0; i<opt.number_of_kseqs; i++){
		free(seqsInCluster[i]);
	}
	free(seqsInCluster);
}
nw_alignment* initialize_nw(int longestSeq){
	nw_alignment* nw_struct = malloc(sizeof(nw_alignment));
	nw_struct->nw=needleman_wunsch_new();
	nw_struct->aln=alignment_create(2*longestSeq); //MAX alignment length is 2 times the longest sequence
	nw_struct->scoring=malloc(sizeof(scoring_t));
	int match=2;
        int mismatch=-1;
        int gap_open=-3;
        int gap_extend=-1;
        int num_of_mismatches=0;
        int num_of_indels = 0;
        bool no_start_gap_penalty=true;
        bool no_end_gap_penalty=true;
        bool no_gaps_in_a = false;
        bool no_gaps_in_b = false;
        bool no_mismatches = false;
        bool case_sensitive=false;
        scoring_init(nw_struct->scoring, match, mismatch, gap_open, gap_extend, no_start_gap_penalty, no_end_gap_penalty, no_gaps_in_a, no_gaps_in_b, no_mismatches, case_sensitive);
	return nw_struct;	
}
void free_nw(nw_alignment* nw_struct){
	free(nw_struct->scoring);
	alignment_free(nw_struct->aln);
	needleman_wunsch_free(nw_struct->nw);
	free(nw_struct);
}
int findLargestCluster(int* clusterSizes, int numberOfClusters){
	int i;
	int largest = clusterSizes[1];
	int largestCluster = 1;
	for(i=1; i<numberOfClusters; i++){
		if ( clusterSizes[i] > largest ){
			largest = clusterSizes[i];
			largestCluster = i;
		}
	}
	return largestCluster;
}
void allocateMemForDistMat(int* clusterSizes, int largestCluster, double** distMat){
	int i;
	distMat = (double**)malloc(clusterSizes[largestCluster]*sizeof(double *));
	for(i=0; i<clusterSizes[largestCluster]; i++){
		distMat[i] = (double*)malloc(sizeof(double));
	}
}
void allocateMemForAlign(int** DATA, int longest_seq, int* mult){
	DATA = (int **)malloc(2*sizeof(int *));
	int k;
	for(k=0; k<2; k++){
		DATA[k] = (int *)malloc(longest_seq*sizeof(int));
	}
	mult = (int *)malloc(longest_seq*sizeof(int));
}
void freeMemForAlign(int** DATA, int longest_seq, int* mult){
	free(DATA[0]);
	free(DATA[1]);
	free(DATA);
	free(mult);
}
void freeMemForDistMat(int* clusterSizes, int largestCluster, double** distMat){
	int i;
	for(i=0; i<clusterSizes[largestCluster]; i++){
		free(distMat[i]);
	}
	free(distMat);
}
int main(int argc, char **argv){
	Options opt;
	opt.number_of_clusters = 10;
	opt.number_of_kseqs=50;
	opt.slash=0;
	opt.default_directory=1;
	strcpy(opt.output_directory,"");
	parse_options(argc, argv, &opt);
	struct hashmap map;
	FILE* fasta_for_clustering;
	//if (( fasta_for_clustering = fopen("/space/s1/lenore/crux_db/crux_db2/new_blast/GazF1_GazR1/GazF1_GazR1_db_filtered/GazF1_GazR1_fasta_and_taxonomy/GazF1_GazR1_.rmAmbig.fasta","r")) == (FILE *) NULL ) fprintf(stderr,"FASTA file could not be opened.\n");
	if (( fasta_for_clustering = fopen(opt.fasta,"r")) == (FILE *) NULL ) fprintf(stderr,"FASTA file could not be opened.\n");
	int number_of_sequences = 0;
	int* fasta_specs = (int *)malloc(4*sizeof(int));
	setNumSeq(fasta_for_clustering,fasta_specs);
	if (opt.number_of_clusters > fasta_specs[0] || opt.number_of_clusters < 2){
		printf("please enter a number of clusters > 1 and less than the number of sequences provided\n");
		exit(1);
	}
	fasta_specs[3] = opt.number_of_clusters+1; //NUMBER OF CLUSTERS
	printf("Number of sequences: %d\n",fasta_specs[0]);
	//printf("Longest sequence: %d\n",fasta_specs[1]);
	//printf("Longest name: %d\n",fasta_specs[2]);
	printf("Number of clusters: %d\n",fasta_specs[3]-1);
	fclose(fasta_for_clustering);
	char **seqNames;
	char **sequences;
	int i,j;
	seqNames = (char **)malloc(fasta_specs[0]*sizeof(char *));
	sequences = (char **)malloc(fasta_specs[0]*sizeof(char *));
	for(i=0; i<fasta_specs[0]; i++){
		seqNames[i]=(char *)malloc(fasta_specs[2]*sizeof(char));
		sequences[i]=(char *)malloc(fasta_specs[1]*sizeof(char));
	}
	hashmap_init(&map, hashmap_hash_string, hashmap_compare_string, fasta_specs[0]);
	if (( fasta_for_clustering = fopen(opt.fasta,"r")) == (FILE *) NULL ) fprintf(stderr,"FASTA file could not be opened.\n");
	readInFasta(fasta_for_clustering,seqNames,sequences,map);
	fclose(fasta_for_clustering);
	char **taxonomy = (char **)malloc(fasta_specs[0]*sizeof(char *));
	for(i=0; i<fasta_specs[0]; i++){
		taxonomy[i]=(char *)malloc(FASTA_MAXLINE*sizeof(char));
	}
	FILE* taxonomyFile;
	//if (( taxonomyFile = fopen("/space/s1/lenore/crux_db/crux_db2/new_blast/GazF1_GazR1/GazF1_GazR1_db_filtered/GazF1_GazR1_fasta_and_taxonomy/GazF1_GazR1_taxonomy.sorted.txt","r")) == (FILE *) NULL ) fprintf(stderr,"TAXONOMY file could not be opened.\n");
	if (( taxonomyFile = fopen(opt.taxonomy,"r")) == (FILE *) NULL ) fprintf(stderr,"TAXONOMY file could not be opened.\n");
	readInTaxFile(taxonomyFile,taxonomy);
	fclose(taxonomyFile);
	struct hashmap seqsToCompare;
	hashmap_init(&seqsToCompare, hashmap_hash_string, hashmap_compare_string, fasta_specs[0]);
	int* chooseK = (int *)malloc(MAXNUMBEROFKSEQS*sizeof(int));
	for(i=0; i<MAXNUMBEROFKSEQS; i++){
		chooseK[i]=-1;
	}
	struct timespec tstart={0,0}, tend={0,0};
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	printf("Choosing %d sequences at random...\n",opt.number_of_kseqs);
	int choosing=0;
		srand(time(0));
	while(choosing==0){
		int random_number = generateRandom(fasta_specs[0]-1);
		int index=0;
		for(i=opt.number_of_kseqs; i>=0; i--){
			if (chooseK[i]==-1){
				index=i;
			}
		}
		int dontadd =0;
		for(i=0; i<index; i++){
			if (chooseK[i]==random_number){
				dontadd = 1;
			}
		}
		if (dontadd==0){
			chooseK[index] = random_number;
			hashmap_put(&seqsToCompare,seqNames[random_number],1);
		}
		if (chooseK[opt.number_of_kseqs-1]!=-1){
			choosing=1;
		}
	}
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %lf seconds\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	double** distMat = (double **)malloc((opt.number_of_kseqs+1)*sizeof(double *));
	for(i=0; i<opt.number_of_kseqs+1; i++){
		distMat[i] = (double *)malloc((opt.number_of_kseqs+1)*sizeof(double));
	}
	for(i=0; i<opt.number_of_kseqs+1; i++){
		for(j=0; j<opt.number_of_kseqs+1; j++){
			distMat[i][j] = 0;
		}
	}
	char** seqsInCluster = (char **)malloc(MAXNUMBEROFKSEQS*sizeof(char *));
	for(i=0; i<MAXNUMBEROFKSEQS; i++){
		seqsInCluster[i] = (char *)malloc(fasta_specs[1]*sizeof(char));
		//printf("seqsInCluster[%d]: %s\n",i,seqsInCluster[i]);
	}
	for(i=0; i<opt.number_of_kseqs; i++){
		strcpy(seqsInCluster[i],sequences[chooseK[i]]);
	}
	printf("Creating distance matrix...\n");
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	createDistMat(seqsInCluster,distMat,opt.number_of_kseqs,fasta_specs[1]);
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %lf seconds\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	char*** clusters = (char***)malloc(MAXNUMBEROFCLUSTERS*sizeof(char **));
	for(i=0; i<MAXNUMBEROFCLUSTERS; i++){
		clusters[i]=(char **)malloc((fasta_specs[0])*sizeof(char *));
		for(j=0; j<fasta_specs[0]; j++){
			clusters[i][j]=(char *)malloc(MAXNAME*sizeof(char));
			clusters[i][j][0]='\0';
		}
	}
	for(i=0; i<opt.number_of_kseqs; i++){
		strcpy(clusters[0][i],seqNames[chooseK[i]]);
	}
	int* clusterSize = (int *)malloc(fasta_specs[3]*sizeof(int));
	clusterSize[0]=opt.number_of_kseqs;
	//print_distance_matrix(distMat,clusters,0,opt.number_of_kseqs);
	node **tree = malloc((fasta_specs[3]+1)*sizeof(node *));
	for(i=0; i<fasta_specs[3]+1; i++){
		tree[i]=(node *)malloc((2*opt.number_of_kseqs-1)*sizeof(node));
		for(j=0; j<2*opt.number_of_kseqs-1; j++){
			tree[i][j].name=(char *)malloc(MAXNAME*sizeof(char));
			strcpy(tree[i][j].name,"internal");
			tree[i][j].nodeToCut=0;
		}
	}
	int root = NJ(tree,distMat,clusters,opt.number_of_kseqs,0);
	get_number_descendants(tree,root,0);
	//printtree(tree,0,opt.number_of_kseqs);
	double* branchLengths = (double *)malloc((2*opt.number_of_kseqs-1)*sizeof(double));
	for(i=0; i<2*opt.number_of_kseqs-1; i++){
		branchLengths[i] = tree[0][i].bl;
	}
	int* indexArray = (int *)malloc((2*opt.number_of_kseqs-1)*sizeof(int));
	sortArray(branchLengths,indexArray,opt);
	printf("Longest Branch: %lf node: %d\n",branchLengths[0],indexArray[0]);
	int numberOfNodesToCut=0;
	for(i=0; i<2*opt.number_of_kseqs-1; i++){
		if (tree[0][indexArray[i]].nd > 4 && numberOfNodesToCut < opt.number_of_clusters-1){
			printf("cutting at node %d\n",indexArray[i]);
			numberOfNodesToCut++;
			printdescendants(tree,indexArray[i],clusters,numberOfNodesToCut,0,opt);
			tree[0][indexArray[i]].nodeToCut=1;
			clearDescendants(tree,tree[0][indexArray[i]].up[0],0);
			clearDescendants(tree,tree[0][indexArray[i]].up[1],0);
			updateNumberOfDescendants(tree,indexArray[i],tree[0][indexArray[i]].nd,0);
			printf("updating tree\n");
		}
	}
	free(indexArray);
	free(branchLengths);
	numberOfNodesToCut++;
	printf("number of nodes to cut: %d\n",numberOfNodesToCut);
	printdescendants(tree,root,clusters,numberOfNodesToCut,0,opt);
	for(i=0; i<numberOfNodesToCut; i++){
		clusterSize[i]=countNumInCluster(clusters,i,opt);
		printf("number in cluster%d: %d\n",i,clusterSize[i]);
	}
	//node** treeArr;
	//allocateMemForTreeArr(numBranchesToFind,clusterSize,treeArr);
	int k,l,m;
	int count=0;
	double averageDist;
	char*** seqsInClusterAvg = (char ***)malloc((opt.number_of_clusters-1)*sizeof(char **));
	for(i=0; i<opt.number_of_clusters-1; i++){
		seqsInClusterAvg[i] = (char **)malloc(opt.number_of_kseqs*sizeof(char *));
		for(j=0; j<opt.number_of_kseqs; j++){
			seqsInClusterAvg[i][j] = (char *)malloc(fasta_specs[1]*sizeof(char));
		}
	}
	double* calculateAvg = (double *)malloc(2*sizeof(double));
	calculateAvg[0]=0;
	calculateAvg[1]=0;
	double*** distMatAvg = (double ***)malloc(numberOfNodesToCut*sizeof(double **));
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	printf("Calculating pairwise average...\n");
	for(i=1; i<numberOfNodesToCut; i++){
		for(j=0; j<clusterSize[i]; j++){
			for(k=0; k<opt.number_of_kseqs; k++){
				if (strcmp(clusters[i][j],clusters[0][k])==0){
					strcpy(seqsInCluster[j],sequences[chooseK[k]]);
					strcpy(seqsInClusterAvg[i-1][j],sequences[chooseK[k]]);
				}
			}
		}
		for(l=0; l<opt.number_of_kseqs+1; l++){
			for(m=0; m<opt.number_of_kseqs+1; m++){
				distMat[l][m]=0;
			}
		}
		if (i>1){
			distMatAvg[count] = (double **)malloc(clusterSize[i-1]*sizeof(double *));
			for(l=0; l<clusterSize[i-1]; l++){
				distMatAvg[count][l] = (double *)malloc(clusterSize[i]*sizeof(double));
			}
			averageDist = calculateAverageDist(seqsInClusterAvg[i-2],seqsInClusterAvg[i-1],clusterSize[i-1],clusterSize[i],fasta_specs[1],distMatAvg,count,calculateAvg);
			count++;
		}
		createDistMat(seqsInCluster,distMat,clusterSize[i],fasta_specs[1]);
		//print_distance_matrix(distMat,clusters,i,clusterSize[i]);
		root = NJ(tree,distMat,clusters,clusterSize[i],i);
		get_number_descendants(tree,root,i);
		printf("cluster %d tree\n",i);
		//printtree(tree,i,clusterSize[i]);
		print_taxonomy(clusters,seqNames,taxonomy,i,clusterSize[i],chooseK,opt);
	}
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %lf seconds\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	freeSeqsInClusterAvg(seqsInClusterAvg,opt);
	freeDistMat(distMat,opt.number_of_kseqs+1);
	freeDistMatAvg(distMatAvg,clusterSize,count);
	freeTreeArr(tree,clusterSize);
	double average = calculateAvg[0]/calculateAvg[1];
	free(calculateAvg);
	printf("average pairwise: %lf\n",average);
	double shortest_distance=1;
	double distance=0;
	int closestCluster=0;
	int next=0;
	int update_initial_cluster = opt.number_of_kseqs;
	int* newClusterSizes = (int *)malloc(MAXNUMBEROFCLUSTERS*sizeof(int));
	for(i=0; i<numberOfNodesToCut; i++){
		newClusterSizes[i] = clusterSize[i];
	}
	int largest_cluster = findLargestCluster(clusterSize,numberOfNodesToCut);
	nw_alignment* nw_struct;
	nw_struct = initialize_nw(fasta_specs[3]);
	double** distMat2;
	distMat2 = (double**)malloc((clusterSize[largest_cluster]+1)*sizeof(double *));
	for(i=0; i<clusterSize[largest_cluster]+1; i++){
		distMat2[i] = (double*)malloc(sizeof(double));
	}
	int** DATA;
	int* mult;
	DATA = (int **)malloc(2*sizeof(int *));
	for(i=0; i<2; i++){
		DATA[i] = (int *)malloc(2*fasta_specs[1]*sizeof(int)); //MAX alignment length (2 times the longest sequence)
	}
	mult = (int *)malloc(2*fasta_specs[1]*sizeof(int)); //MAX alignment length (2 times the longest sequence)
	//allocateMemForDistMat(clusterSize,largest_cluster,&distMat2);
	//allocateMemForAlign(&DATA,fasta_specs[1],&mult);
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	printf("Starting to assign sequences to clusters\n");
	for(i=0; i<fasta_specs[0]; i++){
		shortest_distance=1;
		//for(j=0; j<opt.number_of_kseqs; j++){
		//	if (i==chooseK[j]){
		//		next=1;
		//	}
		//}
		if ( 1==(int)hashmap_get(&seqsToCompare,seqNames[i]) ){
			next=1;
		}
		if (next != 1){
			for(j=1; j<numberOfNodesToCut; j++){
				for(k=0; k<clusterSize[j]; k++){
					//for(l=0; l<update_initial_cluster; l++){
						//if (strcmp(clusters[j][k],clusters[0][l])==0){
						//	strcpy(seqsInCluster[k],sequences[chooseK[l]]);
						//}
					//}
					if ( 1==(int)hashmap_get(&seqsToCompare,clusters[j][k]) ){
						strcpy(seqsInCluster[k],(char *)hashmap_get(&map,clusters[j][k]));
					}
				}
				distance=findShortestDist(seqsInCluster,sequences[i],clusterSize[j],fasta_specs[3],nw_struct,distMat2,DATA,mult);
				if (distance < shortest_distance){
					shortest_distance = distance;
					closestCluster=j;
				}
			}
			printf("%s: cluster: %d distance: %lf tax: %s\n",seqNames[i],closestCluster,shortest_distance,taxonomy[i]);
			if (shortest_distance > average){
				printf("making new cluster for %s number of clusters now %d\n",seqNames[i],numberOfNodesToCut);
				makeNewCluster(clusters,numberOfNodesToCut,seqNames[i],chooseK,i,update_initial_cluster,seqsToCompare);
				newClusterSizes[numberOfNodesToCut]=1;
				clusterSize[numberOfNodesToCut]=1;
				update_initial_cluster++;
				numberOfNodesToCut++;
			}else{
				addToCluster(clusters,seqNames[i],closestCluster,newClusterSizes[closestCluster]);
				newClusterSizes[closestCluster]++;
			}
		}
		next=0;
	}
	free(DATA[0]);
	free(DATA[1]);
	free(DATA);
	free(mult);
	for(i=0; i<clusterSize[largest_cluster];i++){
		free(distMat2[i]);
	}
	free(distMat2);
	free_nw(nw_struct);
	free(chooseK);
	freeSeqsInCluster(seqsInCluster,opt);
	printClusters(clusters,numberOfNodesToCut,newClusterSizes,seqNames,sequences,fasta_specs[0],taxonomy,opt);
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Finished! Took %lf seconds\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	free(clusterSize);
	//freeMemForAlign(DATA,fasta_specs[1],mult);
	//freeMemForDistMat(clusterSize,largest_cluster,distMat2);
	free(newClusterSizes);
	freeClusters(clusters,opt);
	freeSequences(seqNames,fasta_specs[0]);
	freeSequences(sequences,fasta_specs[0]);
	freeSequences(taxonomy,fasta_specs[0]);
	free(fasta_specs);
	hashmap_destroy(&map);
	hashmap_destroy(&seqsToCompare);
}
