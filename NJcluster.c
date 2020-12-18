#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <pthread.h>
#include <string.h>
#include <dirent.h>
#include <regex.h>
#include "needleman_wunsch.h"
#include "global.h"
#include "hashmap.h"
#include "WFA/affine_wavefront_align.h"

//struct hashmap map;
char*** clusters;
//char** seqNames;
//char** sequences;
pthread_mutex_t lock;
double LRVECnc[4][4], RRVECnc[4][4], RRVALnc[4], PMATnc[2][4][5];
double LRVEC[STATESPACE][STATESPACE], RRVEC[STATESPACE][STATESPACE], RRVAL[STATESPACE], PMAT1[STATESPACE][STATESPACE], PMAT2[STATESPACE][STATESPACE];
double parameters[10] = {0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
int COUNT2;//counting the number of tiems the likelihood function is called.
int COUNT; //this one counts how many times the likelihood function has been called
int localnode;
double *localpi, *statevector, *UFCnc, **templike_nc;
double Logfactorial[MAXNUMBEROFINDINSPECIES];
double currentestimate[10]; //DELETEME
node** treeArr;
char** rootSeqs;
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
void assignDepth(node** tree, int node0, int node1, int depth){
	if( node0 != -1 && node1 != -1){
		tree[0][node0].depth = depth;
		tree[0][node1].depth = depth;
		assignDepth(tree,tree[0][node0].up[0], tree[0][node0].up[1],depth+1);
		assignDepth(tree,tree[0][node1].up[0], tree[0][node1].up[1],depth+1);
	}
}
void readInFasta(FILE* fasta, char** seqNames, char** sequences){
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
			//for(j=0; j<100; j++){
			for(j=0;buffer[j]!='\n';j++){
				sequences[i][j]=buffer[j];
			}
			sequences[i][j]='\0';
			//hashmap_put(&map,seqNames[i],sequences[i]);
			i++;
		}
	}
}
/*int readInFastaToAssign(FILE* fasta, char** seqNames, char** sequences, struct hashmap mapToAssign){
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
			char sequence[FASTA_MAXLINE];
			for(j=0;buffer[j]!='\n';j++){
				sequence[j]=buffer[j];
			}
			sequence[j-1]='\0';
			strcpy(sequences[i],sequence);
			hashmap_put(&mapToAssign,seqNames[i],sequence);
			i++;
		}
	}
	return i;
}*/
void readInClusterSeqs(char** filesForReading, int number_of_clusters){
	char buffer[FASTA_MAXLINE];
	int i=0;
	int j=0;
	FILE *fasta;
	for(i=0; i<number_of_clusters; i++){
		int k=0;
		if (( fasta= fopen(filesForReading[i],"r")) == (FILE *) NULL ) fprintf(stderr,"FASTA file could not be opened.\n");
		while( fgets(buffer,FASTA_MAXLINE,fasta) != NULL ){
			if ( buffer[0] == '>' ){
				for(j=1;buffer[j]!='\n';j++){
					clusters[i][k][j-1]=buffer[j];
				}
				clusters[i][k][j-1]='\0';
			}else{
				char sequences[FASTA_MAXLINE];
				for(j=0;buffer[j]!='\n';j++){
					sequences[j]=buffer[j];
				}
				sequences[j-1]='\0';
				//hashmap_put(&map,clusters[i][k],sequences);
				k++;
			}
		}
	}
}
int generateRandom(int numSeq){
	int i;
	int num = (rand() % (numSeq+1));
	return num;
}

void print_distance_matrix(double** distMat, int whichCluster, int clusterSize){
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
	//if (numrealsites < 1) {printf("Wacky distance in function 'Get_dist_JC' between sequence %i and %i\n",index1,index2); exit(-1);}
	if (numrealsites < 5) {
		printf("Wacky distance in function 'Get_dist_JC' between sequence %i and %i\n",index1,index2); M[index1][index2]=1; return;}
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
void elimfrom(int site, int alignment_length, int number_of_sequences, int **seq){
	int i,j;
	for(i=0; i<number_of_sequences; i++){
		for(j=site;j<alignment_length-1;j++){
			seq[i][j] = seq[i][j+1];
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
int NJ(node** tree, double** distMat,int clusterSize,int whichTree, int whichTree2){
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
			strcpy(tree[whichTree][child1].name,clusters[whichTree2][index[pair[0]]]);
		}
		if (tree[whichTree][child2].up[1]==-1){
			strcpy(tree[whichTree][child2].name,clusters[whichTree2][index[pair[1]]]);
		}
		// HERE WE DISALLOW NEGATIVE BRANCHLENGTHS!  IS THIS THE BEST THING TO DO?
		tree[whichTree][child2].distance = distMat[pair[0]][pair[1]];
		if ((u1 = (distMat[pair[0]][pair[1]] +r[pair[0]]-r[pair[1]])/2.0) < MINBL){
			tree[whichTree][child1].bl = MINBL;
		}else{
			tree[whichTree][child1].bl = u1;
		}
		if ((tree[whichTree][child2].bl = distMat[pair[0]][pair[1]]-u1) < MINBL)/*(DM[pair[0]][pair[1]] +r[pair[1]]-r[pair[0]])/2.0;*/{
			tree[whichTree][child2].bl = MINBL;
			tree[whichTree][child2].distance = distMat[pair[0]][pair[1]];
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
		strcpy(tree[whichTree][child1].name,clusters[whichTree2][index[0]]);
	}
	if (tree[whichTree][child2].up[0]==-1 && strcmp(tree[whichTree][child2].name,"internal")==0){
		strcpy(tree[whichTree][child1].name,clusters[whichTree2][index[0]]);
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
void printtree(node** tree, int whichTree, int clusterSize){
	int i;
	for(i=0; i<2*clusterSize-1; i++){
		printf("%i: up: %i %i, down: %i, bl: %f, nd: %d, name: %s\n",i,tree[whichTree][i].up[0],tree[whichTree][i].up[1],tree[whichTree][i].down,tree[whichTree][i].bl/*totsites*/,tree[whichTree][i].nd,tree[whichTree][i].name);
	}
}
void sortArray(double* branchLengths, int* indexArray, int kseqs){
	int i,j, tmp2;
	double tmp;
	//Bubble Sort
	for(i=0; i<2*kseqs-1; i++){
		indexArray[i]=i;
	}
	for(i=0; i<2*kseqs; i++){
		for(j=i+1;j<2*kseqs-1;j++){
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
void printdescendants(node** tree, int node, int clusterNumber, int whichTree, int kseqs){
	int child1 = tree[whichTree][node].up[0];
	int child2 = tree[whichTree][node].up[1];
	int i;
	if (tree[whichTree][node].nodeToCut==1){
		return;
	}
	if (tree[whichTree][node].up[0]==-1 && tree[whichTree][node].up[1]==-1){
		int count=0;
		for(i=kseqs-1; i>=0; i--){
			if (clusters[clusterNumber][i][0]=='\0'){
				count=i;
			}
		}
		printf("count %d node %d: %s\n",count,node,tree[whichTree][node].name);
		strcpy(clusters[clusterNumber][count],tree[whichTree][node].name);
	}else{
		printdescendants(tree,child1,clusterNumber,whichTree,kseqs);
		printdescendants(tree,child2,clusterNumber,whichTree,kseqs);
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
int countNumInCluster(int index, int kseqs){
	int i,j;
	for(i=0; i<kseqs; i++){
		if(clusters[index][i][0]!='\0'){
			j=i;
		}
	}
	return j+1;
}
void allocateMemForTreeArr(int numberOfClusters, int* clusterSize, node** treeArr, int kseqs){
	int i,j;
	treeArr=(node **)malloc(numberOfClusters*sizeof(node *));
	for(i=0; i<numberOfClusters; i++){
		treeArr[i]=malloc((2*kseqs-1)*sizeof(node));
		for(j=0; j<2*kseqs-1; j++){
			treeArr[i][j].name = (char *)malloc(MAXNAME*sizeof(char));
		}
	}
}
void createDistMat_WFA(char** seqsInCluster, double** distMat, int clusterSize, int* fasta_specs){
	int i,j;
	affine_penalties_t affine_penalties = {
		.match = 0,
		.mismatch =4,
		.gap_opening = 6,
		.gap_extension = 2,
	};
	//char* seq1 = (char *)malloc(2*fasta_specs[1]*sizeof(char));
	//char* seq2 = (char *)malloc(2*fasta_specs[1]*sizeof(char));
	for (i=0; i<clusterSize; i++){
		for (j=i+1; j<clusterSize; j++){
			mm_allocator_t* const mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);
			char* seq1 = strdup(seqsInCluster[i]);
			char* seq2 = strdup(seqsInCluster[j]);
			affine_wavefronts_t* affine_wavefronts = affine_wavefronts_new_complete(strlen(seq1),strlen(seq2),&affine_penalties,NULL,mm_allocator);
			affine_wavefronts_align(affine_wavefronts,seq1,strlen(seq1),seq2,strlen(seq2));
			const int score = edit_cigar_score_gap_affine(&affine_wavefronts->edit_cigar,&affine_penalties);
			//fprintf(stderr,"  PATTERN  %s\n",seq1);
			//fprintf(stderr,"  TEXT     %s\n",seq2);
			//fprintf(stderr,"  SCORE COMPUTED %d\t",score);
			//edit_cigar_print_pretty(stderr,seq1,strlen(seq1),seq2,strlen(seq2),&affine_wavefronts->edit_cigar,mm_allocator);
			char* const pattern_alg = mm_allocator_calloc(mm_allocator,strlen(seq1)+strlen(seq2)+1,char,true);
			char* const text_alg = mm_allocator_calloc(mm_allocator,strlen(seq1)+strlen(seq2)+1,char,true);
			int k, alg_pos =0, pattern_pos= 0, text_pos =0;
			for (k=affine_wavefronts->edit_cigar.begin_offset;k<affine_wavefronts->edit_cigar.end_offset;++k) {
				switch (affine_wavefronts->edit_cigar.operations[k]) {
					case 'M':
						if (seq1[pattern_pos] != seq2[text_pos]) {
							pattern_alg[alg_pos] = seq1[pattern_pos];
							text_alg[alg_pos++] = seq2[text_pos];
						}else{
							pattern_alg[alg_pos] = seq1[pattern_pos];
							text_alg[alg_pos++] = seq2[text_pos];
						}
						pattern_pos++; text_pos++;
						break;
					case 'X':
						if (seq1[pattern_pos] != seq2[text_pos]) {
							pattern_alg[alg_pos] = seq1[pattern_pos++];
							text_alg[alg_pos++] = seq2[text_pos++];
						}else{
							pattern_alg[alg_pos] = seq1[pattern_pos++];
							text_alg[alg_pos++] = seq2[text_pos++];
						}
						break;
					case 'I':
						pattern_alg[alg_pos] = '-';
						text_alg[alg_pos++] = seq2[text_pos++];
						break;
					case 'D':
						pattern_alg[alg_pos] = seq1[pattern_pos++];
						text_alg[alg_pos++] = '-';
						break;
					default:
						break;
				}
			}
			k=0;
			while (pattern_pos < strlen(seq1)) {
				pattern_alg[alg_pos+k] = seq1[pattern_pos++];
				++k;
			}
			while (text_pos < strlen(seq2)) {
				text_alg[alg_pos+k] = seq2[text_pos++];
				++k;
			}
			int alignment_length = strlen(pattern_alg);
			int** DATA = (int **)malloc(2*sizeof(int *));
			for(k=0; k<2; k++){
				DATA[k] = (int *)malloc(alignment_length*sizeof(int));
			}
			int* mult = (int *)malloc(alignment_length*sizeof(int));
			alignment_length = populate_DATA(pattern_alg,text_alg,DATA,alignment_length,mult);
			Get_dist_JC(alignment_length,distMat,DATA,mult,i,j);
			free(mult);
			free(DATA[0]);
			free(DATA[1]);
			free(DATA);
			affine_wavefronts_delete(affine_wavefronts);
			mm_allocator_delete(mm_allocator);
			free(seq1);
			free(seq2);
		}
	}
}
void createDistMat(char** seqsInCluster, double** distMat, int clusterSize, int* fasta_specs){
	int i,j;
	nw_aligner_t *nw;
	alignment_t *aln;
	scoring_t *scoring;
	nw=needleman_wunsch_new();
	aln=alignment_create(2*fasta_specs[1]);
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
			needleman_wunsch_align(seqsInCluster[i],seqsInCluster[j], scoring, nw, aln);
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
void calculateAverageDist_WFA(char** seqsA, char** seqsB, int sizeA, int sizeB, int longestSeq, double** distMat, double* calculateAvg){
	int i,j;
	affine_penalties_t affine_penalties = {
		.match = 0,
		.mismatch =4,
		.gap_opening = 6,
		.gap_extension = 2,
	};
	for(i=0; i<sizeA; i++){
		for(j=0; j<sizeB; j++){
			mm_allocator_t* const mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);
			char* seq1 = strdup(seqsA[i]);
			char* seq2 = strdup(seqsB[j]);
			affine_wavefronts_t* affine_wavefronts = affine_wavefronts_new_complete(strlen(seq1),strlen(seq2),&affine_penalties,NULL,mm_allocator);
			affine_wavefronts_align(affine_wavefronts,seq1,strlen(seq1),seq2,strlen(seq2));
			//const int score = edit_cigar_score_gap_affine(&affine_wavefronts->edit_cigar,&affine_penalties);
			char* const pattern_alg = mm_allocator_calloc(mm_allocator,strlen(seq1)+strlen(seq2)+1,char,true);
			char* const text_alg = mm_allocator_calloc(mm_allocator,strlen(seq1)+strlen(seq2)+1,char,true);
			int k, alg_pos =0, pattern_pos= 0, text_pos =0;
			for (k=affine_wavefronts->edit_cigar.begin_offset;k<affine_wavefronts->edit_cigar.end_offset;++k) {
				switch (affine_wavefronts->edit_cigar.operations[k]) {
					case 'M':
						if (seq1[pattern_pos] != seq2[text_pos]) {
							pattern_alg[alg_pos] = seq1[pattern_pos];
							text_alg[alg_pos++] = seq2[text_pos];
						}else{
							pattern_alg[alg_pos] = seq1[pattern_pos];
							text_alg[alg_pos++] = seq2[text_pos];
						}
						pattern_pos++; text_pos++;
						break;
					case 'X':
						if (seq1[pattern_pos] != seq2[text_pos]) {
							pattern_alg[alg_pos] = seq1[pattern_pos++];
							text_alg[alg_pos++] = seq2[text_pos++];
						}else{
							pattern_alg[alg_pos] = seq1[pattern_pos++];
							text_alg[alg_pos++] = seq2[text_pos++];
						}
						break;
					case 'I':
						pattern_alg[alg_pos] = '-';
						text_alg[alg_pos++] = seq2[text_pos++];
						break;
					case 'D':
						pattern_alg[alg_pos] = seq1[pattern_pos++];
						text_alg[alg_pos++] = '-';
						break;
					default:
						break;
				}
			}
			k=0;
			while (pattern_pos < strlen(seq1)) {
				pattern_alg[alg_pos+k] = seq1[pattern_pos++];
				++k;
			}
			while (text_pos < strlen(seq2)) {
				text_alg[alg_pos+k] = seq2[text_pos++];
				++k;
			}
			int alignment_length = strlen(pattern_alg);
			int** DATA = (int **)malloc(2*sizeof(int *));
			for(k=0; k<2; k++){
				DATA[k] = (int *)malloc(alignment_length*sizeof(int));
			}
			int* mult = (int *)malloc(alignment_length*sizeof(int));
			alignment_length = populate_DATA(pattern_alg,text_alg,DATA,alignment_length,mult);
			Get_dist_JC(alignment_length,distMat,DATA,mult,i,j);
			free(mult);
			free(DATA[0]);
			free(DATA[1]);
			free(DATA);
			affine_wavefronts_delete(affine_wavefronts);
			mm_allocator_delete(mm_allocator);
			free(seq1);
			free(seq2);
		}
	}
	double totalDist=calculateAvg[0];
	double numberOfPairs=calculateAvg[1];
	for(i=0; i<sizeA; i++){
		for(j=0; j<sizeB; j++){
			totalDist = totalDist + distMat[i][j];
			numberOfPairs++;
		}
	}
	calculateAvg[0]=totalDist;
	calculateAvg[1]=numberOfPairs;
}
void calculateAverageDist(char** seqsA, char** seqsB, int sizeA, int sizeB, int longestSeq, double** distMat, double* calculateAvg){
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
			Get_dist_JC(alignment_length,distMat,DATA,mult,i,j);
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
			totalDist = totalDist + distMat[i][j];
			numberOfPairs++;
		}
	}
	calculateAvg[0]=totalDist;
	calculateAvg[1]=numberOfPairs;
	free(scoring);
	alignment_free(aln);
	needleman_wunsch_free(nw);
}
double findShortestDist_WFA(int index, char* seq, int clusterSize, int longestSeq, double** distMat, int** DATA, int* mult){
	int i,j;
	affine_penalties_t affine_penalties = {
		.match = 0,
		.mismatch =4,
		.gap_opening = 6,
		.gap_extension = 2,
	};
	for(i=0; i<clusterSize; i++){
		mm_allocator_t* const mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);
		char* seq1 = strdup(rootSeqs[index]);
		char* seq2 = strdup(seq);
		affine_wavefronts_t* affine_wavefronts = affine_wavefronts_new_complete(strlen(seq1),strlen(seq2),&affine_penalties,NULL,mm_allocator);
		affine_wavefronts_align(affine_wavefronts,seq1,strlen(seq1),seq2,strlen(seq2));
		const int score = edit_cigar_score_gap_affine(&affine_wavefronts->edit_cigar,&affine_penalties);
		char* const pattern_alg = mm_allocator_calloc(mm_allocator,strlen(seq1)+strlen(seq2)+1,char,true);
		char* const text_alg = mm_allocator_calloc(mm_allocator,strlen(seq1)+strlen(seq2)+1,char,true);
		int k, alg_pos =0, pattern_pos= 0, text_pos =0;
		for (k=affine_wavefronts->edit_cigar.begin_offset;k<affine_wavefronts->edit_cigar.end_offset;++k) {
			switch (affine_wavefronts->edit_cigar.operations[k]) {
				case 'M':
					if (seq1[pattern_pos] != seq2[text_pos]) {
						pattern_alg[alg_pos] = seq1[pattern_pos];
						text_alg[alg_pos++] = seq2[text_pos];
					}else{
						pattern_alg[alg_pos] = seq1[pattern_pos];
						text_alg[alg_pos++] = seq2[text_pos];
					}
					pattern_pos++; text_pos++;
					break;
				case 'X':
					if (seq1[pattern_pos] != seq2[text_pos]) {
						pattern_alg[alg_pos] = seq1[pattern_pos++];
						text_alg[alg_pos++] = seq2[text_pos++];
					}else{
						pattern_alg[alg_pos] = seq1[pattern_pos++];
						text_alg[alg_pos++] = seq2[text_pos++];
					}
					break;
				case 'I':
					pattern_alg[alg_pos] = '-';
					text_alg[alg_pos++] = seq2[text_pos++];
					break;
				case 'D':
					pattern_alg[alg_pos] = seq1[pattern_pos++];
					text_alg[alg_pos++] = '-';
					break;
				default:
					break;
			}
		}
		k=0;
		while (pattern_pos < strlen(seq1)) {
			pattern_alg[alg_pos+k] = seq1[pattern_pos++];
			++k;
		}
		while (text_pos < strlen(seq2)) {
			text_alg[alg_pos+k] = seq2[text_pos++];
			++k;
		}
		int alignment_length = strlen(pattern_alg);
		alignment_length = populate_DATA(pattern_alg,text_alg,DATA,alignment_length,mult);
		Get_dist_JC(alignment_length,distMat,DATA,mult,i,0);
		affine_wavefronts_delete(affine_wavefronts);
		mm_allocator_delete(mm_allocator);
		free(seq1);
		free(seq2);
	}
	double shortestDist = 1;
	for(i=0; i<clusterSize; i++){
		if (shortestDist > distMat[i][0]){
			shortestDist = distMat[i][0];
		}
	}
	return shortestDist;
}
double findShortestDist(int index, char* seq, int clusterSize, int longestSeq, nw_alignment* nw_struct, double** distMat, int** DATA, int* mult){
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
		//pthread_mutex_lock(&lock);
		needleman_wunsch_align(rootSeqs[index],seq,nw_struct->scoring,nw_struct->nw,nw_struct->aln);
		//pthread_mutex_unlock(&lock);
		int alignment_length = strlen(nw_struct->aln->result_a);
		//int** DATA = (int **)malloc(2*sizeof(int *));
		//int k;
		//for(k=0; k<2; k++){
		//	DATA[k] = (int *)malloc(alignment_length*sizeof(int));
		//}
		//int* mult = (int *)malloc(alignment_length*sizeof(int));
		alignment_length = populate_DATA(nw_struct->aln->result_a,nw_struct->aln->result_b,DATA,alignment_length,mult);
		//pthread_mutex_lock(&lock);
		Get_dist_JC(alignment_length,distMat,DATA,mult,i,0);
		//pthread_mutex_unlock(&lock);
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
void makeNewCluster(char*** clusterSeqs, int number_of_clusters, char* sequence_to_add,int kseqs){
	strcpy(clusterSeqs[number_of_clusters][0],sequence_to_add);
	strcpy(clusterSeqs[0][kseqs],sequence_to_add);
}
void addToCluster(char* sequence_to_add, int clusterNumber, resultsStruct *results, int sizeOfChunk){
	int i;
	for(i=sizeOfChunk-1; i>=0; i--){
		if (results->clusterNumber[i]!=-1){
			break;
		}
	}
	strcpy(results->accession[i+1],sequence_to_add);
	results->clusterNumber[i+1]=clusterNumber;
}
void printClusters(int start, int number_of_clusters, int number_of_total_seqs, Options opt, char** taxonomy, mystruct mstr, int numToPrint, struct hashmap taxMap, int hasTaxFile, char** seqNames, char** sequences){
	FILE *clusterFile, *clusterTaxFile;
	char fileName[FASTA_MAXLINE];
	char *directory = strdup(opt.output_directory);
	int i, j, k;
	char*** assignedSeqs = (char ***)malloc(number_of_clusters*sizeof(char **));
	for (i=0; i<number_of_clusters; i++){
		assignedSeqs[i] = (char **)malloc(MAXNUMINCLUSTER*sizeof(char *));
		for(j=0; j<MAXNUMINCLUSTER; j++){
			assignedSeqs[i][j] = (char *)malloc(MAXNAME*sizeof(char));
		}
	}
	int *clusterSizes = (int *)malloc(number_of_clusters*sizeof(int));
	for(i=0; i<number_of_clusters; i++){
		clusterSizes[i]=0;
	}
	for(i=0; i<numToPrint; i++){
		strcpy(assignedSeqs[mstr.str->clusterNumber[i]-1][clusterSizes[mstr.str->clusterNumber[i]-1]],mstr.str->accession[i]);
		clusterSizes[mstr.str->clusterNumber[i]-1]++;
	}
	for(i=1; i<number_of_clusters; i++){
		if (opt.slash==0 && opt.default_directory==0){
			snprintf(fileName,FASTA_MAXLINE,"%s/%d.fasta",directory,i+start-2);
		}else if (opt.slash==1){
			snprintf(fileName,FASTA_MAXLINE,"%s%d.fasta",directory,i+start-2);
		}else if (opt.default_directory==1){
			snprintf(fileName,FASTA_MAXLINE,"%d.fasta",i+start-2);
		}
		clusterFile = fopen(fileName, "a");
		if (clusterFile == NULL){ printf("Error opening cluster file!"); exit(1); }
		if (hasTaxFile==1){
		if (opt.slash==0 && opt.default_directory==0){
			snprintf(fileName,FASTA_MAXLINE, "%s/%d_taxonomy.txt",directory,i+start-2);
		}else if (opt.slash==1){
			snprintf(fileName,FASTA_MAXLINE, "%s%d_taxonomy.txt",directory,i+start-2);
		}else if (opt.default_directory==1){
			snprintf(fileName,FASTA_MAXLINE, "%d_taxonomy.txt",i+start-2);
		}
		clusterTaxFile = fopen(fileName, "a");
		if (clusterTaxFile == NULL){ printf("Error opening cluster taxonomy file!"); exit(1); }
		}
		/*for(j=0; j<clusterSizes[i]; j++){
			for(k=0; k<number_of_total_seqs; k++){
				if (strcmp(clusters[i][j],seqNames[k])==0){
					fprintf(clusterFile,">%s\n",seqNames[k]);
					fprintf(clusterFile,"%s\n",sequences[k]);
					fprintf(clusterTaxFile,"%s\t%s\n",seqNames[k],taxonomy[k]);
				}
			}
		}*/
		for(j=0; j<clusterSizes[i-1]; j++){
			fprintf(clusterFile,">%s\n",assignedSeqs[i-1][j]);
			//fprintf(clusterFile,"%s\n",(char *)hashmap_get(&map,assignedSeqs[i-1][j]));
			for(k=0; k<number_of_total_seqs; k++){
				if ( strcmp(seqNames[k],assignedSeqs[i-1][j])==0 ){
					break;
				}
			}
			fprintf(clusterFile,"%s\n",sequences[k]);
			if (hasTaxFile==1){
				fprintf(clusterTaxFile,"%s\t%s\n",assignedSeqs[i-1][j],(char *)hashmap_get(&taxMap,assignedSeqs[i-1][j]));
			}
		}
		fclose(clusterFile);
		if (hasTaxFile==1){
			fclose(clusterTaxFile);
		}
	}
	free(clusterSizes);
	for(i=0; i<number_of_clusters; i++){
		free(assignedSeqs[i]);
	}
	free(assignedSeqs);
}
void saveCLSTR(int start, int number_of_total_seqs, mystruct mstr, int numToPrint, char*** clstr, int** clstr_lengths, int max_length, char** seqNames, char** sequences){
	int i,j,k,l;
	int number_of_clusters;
	int index = -1;
	for(i=MAXNUMBEROFCLUSTERS-1; i>=0; i--){
		if(clstr[i][0][0]=='\0'){
			index = i;
		}
	}
	number_of_clusters=index;
	int *clusterSizes = (int *)malloc(number_of_clusters*sizeof(int));
	for(i=0; i<number_of_clusters; i++){
		clusterSizes[i]=0;
	}
	char*** assignedSeqs = (char ***)malloc(number_of_clusters*sizeof(char **));
	for (i=0; i<number_of_clusters; i++){
		assignedSeqs[i] = (char **)malloc(MAXNUMINCLUSTER*sizeof(char *));
		for(j=0; j<MAXNUMINCLUSTER; j++){
			assignedSeqs[i][j] = (char *)malloc(MAXNAME*sizeof(char));
		}
	}
	for(i=0; i<numToPrint; i++){
		strcpy(assignedSeqs[start-1+mstr.str->clusterNumber[i]-1][clusterSizes[start-1+mstr.str->clusterNumber[i]-1]],mstr.str->accession[i]);
		clusterSizes[start-1+mstr.str->clusterNumber[i]-1]++;
	}
	for(i=0; i<number_of_clusters; i++){
		for(j=0; j<clusterSizes[i]; j++){
			for(k=0; k<number_of_total_seqs; k++){
				if ( strcmp(seqNames[k],assignedSeqs[i][j])==0 ){
					break;
				}
			}
			int length;
			for(l=0; l<max_length; l++){
				if ( sequences[k][l]=='\0' ){
					length=l;
					break;
				}
			}
			int place=-1;
			for(l=number_of_total_seqs-1; l>=0; l--){
				if ( clstr[i][l][0]=='\0' ){
					place = l;
				}
			}
			strcpy(clstr[i][place],seqNames[k]);
			clstr_lengths[i][place]=length;
		}
	}
	free(clusterSizes);
	for(i=0; i<number_of_clusters; i++){
		free(assignedSeqs[i]);
	}
	free(assignedSeqs);
}
void printLessThanFour(int numberOfUnassigned, Options opt, struct hashmap taxMap, int starting_number_of_clusters,char** seqNames, char** sequences){
	int i;
	FILE *clusterFile, *clusterTaxFile;
	char fileName[FASTA_MAXLINE];
	char *directory = strdup(opt.output_directory);
	for(i=0; i<numberOfUnassigned; i++){
		if(opt.slash==0 && opt.default_directory==0){
			snprintf(fileName,FASTA_MAXLINE,"%s/%d.fasta",directory,i+starting_number_of_clusters-1);
		}else if (opt.slash==1){
			snprintf(fileName,FASTA_MAXLINE,"%s%d.fasta",directory,i+starting_number_of_clusters-1);
		}else if (opt.default_directory==1){
			snprintf(fileName,FASTA_MAXLINE,"%d.fasta",i+starting_number_of_clusters-1);
		}
		clusterFile = fopen(fileName,"w");
		if (clusterFile == NULL){ printf("Error opening fasta file!"); exit(1); }
		if ( opt.hasTaxFile == 1 ){
			if(opt.slash==0 && opt.default_directory==0){
				snprintf(fileName,FASTA_MAXLINE, "%s/%d_taxonomy.txt",directory,i+starting_number_of_clusters-1);
			}else if (opt.slash==1){
				snprintf(fileName,FASTA_MAXLINE,"%s%d_taxonomy.txt",directory,i+starting_number_of_clusters-1);
			}else if (opt.default_directory==1){
				snprintf(fileName,FASTA_MAXLINE,"%d_taxonomy.txt",i+starting_number_of_clusters-1);
			}
			clusterTaxFile = fopen(fileName, "w");
			if (clusterTaxFile==NULL){ printf("Error opening taxonomy file!"); exit(1); }
			fprintf(clusterTaxFile,"%s\t%s\n",seqNames[i],(char *)hashmap_get(&taxMap,seqNames[i]));
			fclose(clusterTaxFile);
		}
		fprintf(clusterFile,">%s\n",seqNames[i]);
		fprintf(clusterFile,"%s\n",sequences[i]);
		fclose(clusterFile);
	}	
}
void printLessThanFour_CLSTR(int numberOfUnAssigned, char*** clstr, int** clstr_lengths, int max_length, char** seqNames, char** sequences){
	int i,j;
	int index = -1;
	int length = -1;
	for(i=MAXNUMBEROFCLUSTERS-1; i>=0; i--){
		if(clstr[i][0][0]=='\0'){
			index = i;
		}
	}
	for(i=0; i<numberOfUnAssigned; i++){
		strcpy(clstr[index+i][0],seqNames[i]);
		for(j=0; j<max_length; j++){
			if (sequences[i][j] == '\0' ){
				length = j;
				break;
			}
		}
		clstr_lengths[index+i][0]=length;
	}
}
void readInTaxFile(FILE* taxonomyFile, char** taxonomy, struct hashmap taxMap){
	char buffer[FASTA_MAXLINE];
	char *accession, *lineTaxonomy;
	//char accession[MAXNAME];
	//char lineTaxonomy[FASTA_MAXLINE];
	int i=0;
	while( fgets(buffer,FASTA_MAXLINE,taxonomyFile) != NULL){
		accession = strtok(buffer,"\t");
		lineTaxonomy = strtok(NULL,"\n");
		assert(strlen(accession) <= MAXNAME);
		strcpy(taxonomy[i],lineTaxonomy);
		char *acc= strdup(accession);
		char *tax = strdup(lineTaxonomy);
		//strcpy(acc,accession);
		//strcpy(tax,lineTaxonomy);
		hashmap_put(&taxMap,acc,tax);
		i++;
	}
}
void print_taxonomy(int clusterNumber, int clusterSize, int* chooseK, int kseqs, char** seqNames, struct hashmap taxMap,int hasTaxFile){
	int i, j;
	printf("cluster %d:\n",clusterNumber);
	for(i=0; i<clusterSize; i++){
		for(j=0; j<kseqs; j++){
			if (strcmp(clusters[clusterNumber][i],seqNames[chooseK[j]])==0){
				if (hasTaxFile==1){
					//printf("%s\t%s\n",seqNames[chooseK[j]],taxonomy[chooseK[j]]);
					printf("%s\t%s\n",seqNames[chooseK[j]],(char *)hashmap_get(&taxMap,seqNames[chooseK[j]]));
				}else{
					printf("%s\n",seqNames[chooseK[j]]);
				}
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
void freeSeqsInClusterAvg(char*** seqsInClusterAvg, int number_of_kseqs, int number_of_clusters){
	int i,j;
	for(i=0; i<number_of_clusters-1; i++){
		for(j=0; j<number_of_kseqs; j++){
			free(seqsInClusterAvg[i][j]);
		}
		free(seqsInClusterAvg[i]);
	}
	free(seqsInClusterAvg);
}
void freeTreeArr(node** tree,int* clusterSizes){
	int i,j;
}
void freeClusters(int kseqs){
	int i,j;
	for(i=0; i<MAXNUMBEROFCLUSTERS; i++){
		for(j=0; j<kseqs; j++){
			free(clusters[i][j]);
		}
		free(clusters[i]);
	}
	free(clusters);
}
void freeSequences(int number_of_seqs,char** sequences,char** seqNames, char** taxonomy, char** sequencesToClusterLater){
	int i;
	for(i=0;i<number_of_seqs; i++){
		free(sequences[i]);
		free(seqNames[i]);
		free(taxonomy[i]);
		free(sequencesToClusterLater[i]);
		//free(assigned[i]);
	}
	free(sequences);
	free(seqNames);
	free(taxonomy);
	free(sequencesToClusterLater);
	//free(assigned);
}
void freeSeqsInCluster(char** seqsInCluster, int number_of_kseqs){
	int i;
	for(i=0; i<number_of_kseqs; i++){
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
void allocateMemForResults( resultsStruct *results, int sizeOfChunk, int num_threads, int number_of_clusters, double average, int* fasta_specs){
	int i,j;
	results->accession = (char **)malloc((sizeOfChunk)*sizeof(char *));
	results->savedForNewClusters = (char **)malloc((sizeOfChunk)*sizeof(char *));
	for(i=0; i<sizeOfChunk; i++){
		results->accession[i] = malloc((fasta_specs[2]+1)*sizeof(char));
		results->savedForNewClusters[i] = malloc((fasta_specs[2]+1)*sizeof(char));
		results->clusterNumber=-1;
		memset(results->accession[i],'\0',fasta_specs[2]+1);
		memset(results->savedForNewClusters[i],'\0',fasta_specs[2]+1);
	}
	results->clusterNumber = malloc(sizeOfChunk*sizeof(int));
	results->average = average;
	results->number_of_clusters = number_of_clusters;
	results->clusterSizes = (int *)malloc(MAXNUMBEROFCLUSTERS*sizeof(int));
}
void saveForLater(char* accToSave, resultsStruct *results, int sizeOfChunk){
	int i;
	for(i=sizeOfChunk-1; i>=0; i--){
		if (results->savedForNewClusters[i][0]!='\0'){
			break;
		}
	}
	strcpy(results->savedForNewClusters[i+1],accToSave);
}
void *runAssignToCluster(void *ptr){
	struct mystruct *mstr = (mystruct *) ptr;
	resultsStruct *results=mstr->str;
	int start=mstr->start;
	int end=mstr->end;
	int num_threads = mstr->num_threads;
	int sizeOfChunk=end-start;
	printf("sizeOfChunk is %d\n",sizeOfChunk);
	double average = mstr->average;
	//double average = 0.2;
	int number_of_clusters = mstr->number_of_clusters;
	int number_of_kseqs = mstr->number_of_kseqs;
	int largest_cluster = mstr->largest_cluster;
	int* fasta_specs = mstr->fasta_specs;
	char** sequences = mstr->sequences;
	char** seqNames = mstr->seqNames;
	char** taxonomy = mstr->taxonomy;
	//struct hashmap seqsToCompare = mstr->seqsToCompare;
	//struct hashmap assignedSeqs = mstr->assignedSeqs;
	//char** assignedSeqs = mstr->assignedSeqs;
	int* chooseK = mstr->chooseK;
	int i,j,k,l,m;
	int closestCluster;
	double distance=0;
	double shortest_distance=1;
	int* clusterSize = mstr->clusterSize;
	int next=0;
	int numAssigned=0;
	//int number_assigned_in_cluster=results->numassigned;
	//for(i=0; i<sizeOfChunk; i++){
	//	results->clusterNumber[i]=-1;
	//	memset(results->accession[i],'\0',fasta_specs[2]+1);
	//	memset(results->savedForNewClusters[i],'\0',fasta_specs[2]+1);
	//}
	//char** seqsInCluster = (char **)malloc(number_of_kseqs*sizeof(char *));
	//for(i=0; i<number_of_kseqs; i++){
	//	seqsInCluster[i] = (char *)malloc((fasta_specs[1]+1)*sizeof(char));
	//}
	//char*** clusterNames = mstr->clusterNames;
	nw_alignment* nw_struct = mstr->nw_struct;
	//nw_struct = initialize_nw(fasta_specs[3]);
	char*** clusterSeqs = mstr->clusterSeqs;
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
	//for(i=start; i<end; i++){
	for(i=0; i<sizeOfChunk; i++){
		shortest_distance=1;
		for(j=0; j<number_of_kseqs; j++){
			if (i+start==mstr->chooseK[j]){
				next=1;
			}
		}
		//for(j=0; j<number_of_kseqs; j++){
		//	if (strcmp(seqNames[i],clusters[0][j])==0){
		//		next=1;
		//	}
		//}
		//if ( 1==(int)hashmap_get(&seqsToCompare,seqNames[i]) ){
		//	next=1;
		//}
		//if ( 1==(int)hashmap_get(&assignedSeqs,seqNames[i]) ){
		//	next = 1;
		//}
		//for(j=0; j<end-start; j++){
		//	if (strcmp(assignedSeqs[j],seqNames[i])==0){
		//		next=1;
		//	}
		//}
		if (next != 1){
			for(j=1; j<number_of_clusters; j++){
				//for(k=0; k<clusterSize[j]; k++){
					//for(l=0; l<update_initial_cluster; l++){
						//if (strcmp(clusters[j][k],clusters[0][l])==0){
						//	strcpy(seqsInCluster[k],sequences[chooseK[l]]);
						//}
					//}
					//if ( 1==(int)hashmap_get(&seqsToCompare,clusters[j][k]) ){
					//	strcpy(seqsInCluster[k],(char *)hashmap_get(&map,clusters[j][k]));
					//}
				//}
				//pthread_mutex_lock(&lock);
				//distance=findShortestDist(clusterSeqs[j],sequences[i],clusterSize[j],fasta_specs[3],nw_struct,distMat2,DATA,mult);
				if (mstr->use_nw==0){
					distance=findShortestDist_WFA(j-1,sequences[i],1,fasta_specs[3],distMat2,DATA,mult);
				}else{
					distance=findShortestDist(j-1,sequences[i],1,fasta_specs[3],nw_struct,distMat2,DATA,mult);
				}
				//pthread_mutex_unlock(&lock);
				if (distance < shortest_distance){
					shortest_distance = distance;
					closestCluster=j;
				}
			}
			//printf("thread %d\t%s\t%d\t%lf\t%s\n",mstr->threadnumber,seqNames[i],closestCluster,shortest_distance,taxonomy[i]);
			if (shortest_distance > average){
				//printf("making new cluster for %s number of clusters now %d\n",seqNames[i],number_of_clusters);
				//if (num_threads==1){
				/*
				makeNewCluster(clusterSeqs,number_of_clusters,sequences[i],number_of_kseqs);
				addToCluster(seqNames[i],number_of_clusters,results,sizeOfChunk);
				clusterSize[number_of_clusters]=1;	
				number_of_clusters++;
				number_of_kseqs++;
				//newClusterSizes[numberOfNodesToCut]=1;
				//clusterSize[numberOfNodesToCut]=1;
				//update_initial_cluster++;
				//numberOfNodesToCut++;
				results->number_of_clusters++;
				numAssigned++;*/
				//}
				saveForLater(seqNames[i],results,sizeOfChunk);
			}else{
				//for(j=sizeOfChunk-1; j>=0; j--){
				//	if(results->assigned[j][0]=='\0'){
				//		break;
				//	}
				//}
				//strcpy(results->assigned[j],seqNames[i]);
				addToCluster(seqNames[i],closestCluster,results,sizeOfChunk);
				numAssigned++;
				results->numassigned++;
				//newClusterSizes[closestCluster]++;
			}
		}
		next=0;
	}
	free(DATA[0]);
	free(DATA[1]);
	free(DATA);
	free(mult);
	//free_nw(nw_struct);
	for(i=0; i<clusterSize[largest_cluster];i++){
		free(distMat2[i]);
	}
	free(distMat2);
	//freeSeqsInCluster(seqsInCluster,number_of_kseqs);
	mstr->numAssigned = numAssigned;
	pthread_exit(NULL);
}
int countNumUnassigned(char** sequencesToClusterLater, int* fasta_specs){
	int i;
	for(i=0; i<fasta_specs[0]; i++){
		if (sequencesToClusterLater[i][0]=='\0'){
			break;
		}
	}
	return i;
}
void clearSequencesToAssignLater(char** sequencesToClusterLater, int* fasta_specs){
	int i;
	for(i=0; i<fasta_specs[0]; i++){
		strcpy(sequencesToClusterLater[i],"\0");
	}
}
int findNewClusterSizes(int* newClusterSizes, int* fasta_specs){
	int i,j;
	int number_of_clusters=0;
	for(i=0; i<MAXNUMBEROFCLUSTERS; i++){
		if (clusters[i][0][0]=='\0'){
			number_of_clusters=i;
			break;
		}
		for(j=0; j<fasta_specs[0]; j++){
			if (clusters[i][j][0]=='\0'){
				newClusterSizes[i]=j;
				break;
			}
		}
	}
	return number_of_clusters;
}
void makeNewClusters(int numberOfUnAssigned, char** sequencesToClusterLater){
	int i,j;
	for(i=0; i<MAXNUMBEROFCLUSTERS; i++){
		if (clusters[i][0][0]=='\0'){
			break;
		}
	}
	for(j=0; j<numberOfUnAssigned; j++){
		strcpy(clusters[i][0],sequencesToClusterLater[j]);
		i++;
	}
}
/*int readInFastaClusters(Options opt, char** readInFiles, struct hashmap hash_for_files){
	struct dirent *de, *de2;
	DIR *dr = opendir(opt.output_directory);
	regex_t regex1;
	int reti;
	char msgbuf[100];
	int i, j;
	i=0;
	reti = regcomp(&regex1, "[0-9]*\.fasta$",0);
	if (reti){
		printf(stderr, "Could not compile regex\n");
		exit(1);
	}
	if (dr==NULL){
		printf("Could not open directory for reads");
		return 0;
	}
	while((de=readdir(dr))!=NULL){
		reti = regexec(&regex1, de->d_name, 0, NULL, 0);
		if (!reti){
			strcpy(readInFiles[i],opt.output_directory);
			strcat(readInFiles[i],"/");
			strcat(readInFiles[i],de->d_name);
			printf("file: %s\n",readInFiles[i]);
			hashmap_put(&hash_for_files,de->d_name,i);
			i++;
		}else if (reti == REG_NOMATCH){
		}else{
			regerror(reti, &regex1, msgbuf, sizeof(msgbuf));
			printf(stderr, "Regex match failed: %s\n", msgbuf);
			exit(1);
		}
	}
	closedir(dr);
	return i;
}*/
void printInitialClusters(int starting_number_of_clusters,int number_of_clusters, int* clusterSizes, Options opt, int total_number_of_sequences, struct hashmap taxMap, int hasTaxFile, char** seqNames, char** sequences){
	FILE *fasta;
	FILE *tax;
	char fileName[FASTA_MAXLINE];
	char *directory = strdup(opt.output_directory);
	int i,j,k;
	for(i=1; i<number_of_clusters; i++){
		if (opt.slash==0 && opt.default_directory==0){
			snprintf(fileName,FASTA_MAXLINE,"%s/%d.fasta",directory,i+starting_number_of_clusters-2);
		}else if (opt.slash==1){
			snprintf(fileName,FASTA_MAXLINE,"%s%d.fasta",directory,i+starting_number_of_clusters-2);
		}else if (opt.default_directory==1){
			snprintf(fileName,FASTA_MAXLINE,"%d.fasta",i+starting_number_of_clusters-2);
		}
		fasta = fopen(fileName, "w");
		if (fasta == NULL){ printf("Error opening fasta file!"); exit(1); }
		if (hasTaxFile==1){
		if (opt.slash==0 && opt.default_directory==0){
			snprintf(fileName,FASTA_MAXLINE, "%s/%d_taxonomy.txt",directory,i+starting_number_of_clusters-2);
		}else if (opt.slash==1){
			snprintf(fileName,FASTA_MAXLINE, "%s%d_taxonomy.txt",directory,i+starting_number_of_clusters-2);
		}else if (opt.default_directory==1){
			snprintf(fileName,FASTA_MAXLINE, "%d_taxonomy.txt",i+starting_number_of_clusters-2);
		}
		tax = fopen(fileName, "w");
		if (tax==NULL){ printf("Error opening taxonomy file!"); exit(1); }
		}
		for(j=0; j<clusterSizes[i]; j++){
			fprintf(fasta,">%s\n",clusters[i][j]);
			for(k=0; k<total_number_of_sequences; k++){
				if (strcmp(seqNames[k],clusters[i][j])==0 ){
					break;
				}
			}
			fprintf(fasta,"%s\n",sequences[k]);
			//fprintf(fasta,"%s\n",(char *)hashmap_get(&map,clusters[i][j]));
			if (hasTaxFile==1){
				fprintf(tax,"%s\t%s\n",clusters[i][j],(char *)hashmap_get(&taxMap,clusters[i][j]));
			}
		}
		fclose(fasta);
		if (hasTaxFile==1){
			fclose(tax);
		}
	}
}
void printInitialClusters_CLSTR( int starting_number_of_clusters, int number_of_clusters, int* clusterSizes, Options opt, int total_number_of_sequences, int max_length, char*** clstr, int** clstr_lengths, char** seqNames, char** sequences){
	//FILE *clstr;
	//char fileName[FASTA_MAXLINE];
	//char *directory = strdup(opt.output_directory);
	int i,j,k, l;
	//if (opt.slash==0 && opt.default_directory==0){
	//	snprintf(fileName,FASTA_MAXLINE,"%s/output.clstr",directory);
	//}else if (opt.slash==1){
	//	snprintf(fileName,FASTA_MAXLINE,"%soutput.clstr",directory);
	//}else if (opt.default_directory==1){
	//	snprintf(fileName,FASTA_MAXLINE,"output.clstr");
	//}
	//clstr = fopen(fileName, "a");
	//if (clstr == NULL ){ printf("Error opening cluster file!"); exit(1); }
	for(i=1; i<number_of_clusters; i++){
		//fprintf(clstr,">Cluster %d\n",i+starting_number_of_clusters-2);
		for(j=0; j<clusterSizes[i]; j++){
			int length=0;
			for(k=0; k<total_number_of_sequences; k++){
				if (strcmp(seqNames[k],clusters[i][j])==0 ){
					break;
				}
			}
			for(l=0; l<max_length; l++){
				if(sequences[k][l]=='\0'){
					length=l;
					break;
				}
			}
			strcpy(clstr[starting_number_of_clusters+i-2][j],seqNames[k]);
			clstr_lengths[starting_number_of_clusters+i-2][j] = length;
			//fprintf(clstr,"%d\t%dnt, >%s...\n",j,length,seqNames[k]);
		}
	}
}
void printCLSTR(Options opt, char*** clstr, int ** clstr_lengths, int total_number_of_seqs, int number_of_clusters){
	FILE *clstrFile;
	char fileName[FASTA_MAXLINE];
	char *directory = strdup(opt.output_directory);
	int i,j,k;
	if ( opt.output_file[0] == '\0' ){
		strcpy(opt.output_file,"output.clstr");
	}
	if (opt.slash==0 && opt.default_directory==0){
		snprintf(fileName,FASTA_MAXLINE,"%s/%s",directory,opt.output_file);
	}else if (opt.slash==1){
		snprintf(fileName,FASTA_MAXLINE,"%s%s",directory,opt.output_file);
	}else if (opt.default_directory==1){
		snprintf(fileName,FASTA_MAXLINE,"%s",opt.output_file);
	}
	clstrFile = fopen(fileName, "w");
	if (clstr == NULL ){ printf("Error opening cluster file!"); exit(1); }
	for(i=0; i<number_of_clusters; i++){
		if (clstr[i][0][0] != '\0'){
			fprintf(clstrFile,">Cluster %d\n",i);
		}
		for(j=0; j<total_number_of_seqs; j++){
			if ( clstr[i][j][0] != '\0'){
				fprintf(clstrFile,"%d\t%dnt, >%s...\n",j,clstr_lengths[i][j],clstr[i][j]);
			}
		}
	}
	fclose(clstrFile);
}
void makeClustersFromNodeCuts(node** tree, int node, int whichTree, int clusterNumber, int kseqs){
	int child1 = tree[whichTree][node].up[0];
	int child2 = tree[whichTree][node].up[1];
	int i;
	//if (node==-1){ return;}
	if ( tree[whichTree][node].nodeToCut==1 ){
		makeClustersFromNodeCuts(tree,child1,whichTree,clusterNumber+1,kseqs);
		makeClustersFromNodeCuts(tree,child2,whichTree,clusterNumber+2,kseqs);
		return;
	}
	if (tree[whichTree][node].up[0]==-1 && tree[whichTree][node].up[1]==-1){
		printf("clusterNumber: %d name: %s node: %d\n",clusterNumber,tree[whichTree][node].name,node);
		int count=0;
		for(i=kseqs-1; i>=0; i--){
			if (clusters[clusterNumber][i][0]=='\0'){
				count=i;
			}
		}
		strcpy(clusters[clusterNumber][count],tree[whichTree][node].name);
		return;
	}
	//if ( child1 == -1) { return; }
	//if ( child2 == -1) { return; }
	makeClustersFromNodeCuts(tree,child1,whichTree,clusterNumber,kseqs);
	makeClustersFromNodeCuts(tree,child2,whichTree,clusterNumber,kseqs);
	//return;
}
void findFirstNodeCut(node** tree, int node, int whichTree, int* answer){
	int child1 = tree[whichTree][node].up[0];
	int child2 = tree[whichTree][node].up[1];
	if ( tree[whichTree][node].nodeToCut==1 && node != -1){
		printf("returning node %d depth %d\n",node,tree[whichTree][node].depth);	
		if ( answer[0] == -1){
			answer[0]=node;
			answer[1]=tree[whichTree][node].depth;
		}else{
			//if ( answer[1] == tree[whichTree][node].depth ){
			//	answer[0]=-1;
			//}
			if ( answer[1] > tree[whichTree][node].depth ){
				answer[0]=node;
				answer[1]=tree[whichTree][node].depth;
			}
		}
	}
	if ( tree[whichTree][node].up[0]!=-1 && tree[whichTree][node].up[1]!=-1){
		findFirstNodeCut(tree,child1,whichTree,answer);
		findFirstNodeCut(tree,child2,whichTree,answer);
	}
}
int findMaxClade(node** tree, int number_of_leaves){
	int i,j;
	int numberOfNodesToCut = 0;
	for(i=0; i<number_of_leaves; i++){
		for(j=i+1; j<number_of_leaves; j++){
			int lca = findLCA(tree,i,j);
			printf("LCA is %d\n",lca);
			double* distance = (double *)malloc(1*sizeof(double));
			distance[0]=0;
			calculateDistance(tree,i,lca,distance);
			calculateDistance(tree,j,lca,distance);
			printf("distance between %d and %d is %lf\n",i,j,distance[0]);
			if (distance[0] > 0.4){
				tree[0][tree[0][lca].up[0]].nodeToCut = 1;
				tree[0][tree[0][lca].up[0]].nodeToCut = 1;
				//tree[0][tree[0][lca].down].nodeToCut=1;
				//tree[0][lca].nodeToCut = 1;
				numberOfNodesToCut++;
			}
		}
	}
	return numberOfNodesToCut;
}
void calculateDistance(node** tree, int leaf, int LCA,double* distance){
	if (leaf==LCA){ return;
	}else{
		if ( tree[0][leaf].bl > 0 ){
			distance[0] = tree[0][leaf].bl + distance[0];
		}
		calculateDistance(tree,tree[0][leaf].down,LCA,distance);
	}
}
int findLCA(node** tree, int nodeA, int nodeB){
	if (nodeA == nodeB){ return nodeA; }
	if (tree[0][nodeA].depth > tree[0][nodeB].depth){
		int tmp = nodeA;
		nodeA = nodeB;
		nodeB = tmp;
	}
	if ( tree[0][nodeB].down == -1 ){ return nodeB; }
	nodeB = tree[0][nodeB].down;
	return findLCA(tree,nodeA,nodeB);
}
void findLeaves(node** tree, int node, int whichTree, int* parentCuts, int** clusterarr, char*** cluster_seqs, int number_of_sequences,char** seqNames, char** sequences){
	int child1 = tree[whichTree][node].up[0];
	int child2 = tree[whichTree][node].up[1];
	if (tree[whichTree][node].up[0]==-1 && tree[whichTree][node].up[1]==-1){
		printf("finding node to cut parent of %d %s\n",node,tree[whichTree][node].name);
		int parentCut = findParentCut(tree,node,whichTree);
		if (parentCut==-1){ return; }
		//int *nodesInCluster = (int *)malloc(100*sizeof(int));
		//nodesInCluster = (int *)hashmap_get(&clusterhash,parentCut);
		int i=0;
		int index=-1;
		for(i=0; i<number_of_sequences; i++){
			if ( parentCuts[i]==parentCut){
				index=i;
			}
		}
		if (index==-1){
			for(i=number_of_sequences; i>=0; i--){
				if (parentCuts[i]==-1){
					index=i;
				}
			}
		}
		parentCuts[index]=parentCut;
		int placement=-1;
		for(i=number_of_sequences; i>=0; i--){
			if (clusterarr[index][i]==-1){
				placement=i;
			}
		}
		clusterarr[index][placement]=node;
		strcpy(clusters[index+1][placement],tree[whichTree][node].name);
		for(i=0; i<number_of_sequences; i++){
			if (strcmp(seqNames[i],tree[whichTree][node].name)==0){
				break;
			}
		}
		strcpy(cluster_seqs[index+1][placement],sequences[i]);
		//if (nodesInCluster==NULL){
		//	nodesInCluster[0] = node;
		//	int i=1;
		//	for(i=1; i<100; i++){
		//		nodesInCluster[i]=-1;
		//	}
		//}else{
		//	int count = 0;
		//	int i;
		//	for(i=100-1; i>=0; i--){
		//		if (nodesInCluster[i]==-1){
		//			count=i;
		//		}
		//	}
		//	nodesInCluster[count]=node;
		//}		
		//hashmap_put(&clusterhash,parentCut,nodesInCluster);
	}else{
		findLeaves(tree,child1,whichTree,parentCuts,clusterarr,cluster_seqs,number_of_sequences,seqNames,sequences);
		findLeaves(tree,child2,whichTree,parentCuts,clusterarr,cluster_seqs,number_of_sequences,seqNames,sequences);
	}
}
int findParentCut(node** tree, int node, int whichTree){
	if ( tree[whichTree][node].down == -1 ){ return -1; }
	int parent = tree[whichTree][node].down;
	if ( tree[whichTree][node].nodeToCut == 1){
		printf("We have node %d\n",node);
		return node;
	}else{
		findParentCut(tree,parent,whichTree);
	}
}
void addFirstCluster(node** tree, int node, int whichTree, int* parentCuts, int number_of_sequences){
	int i=0;
	int index = -1;
	for(i=number_of_sequences; i>=0; i--){
		if (parentCuts[i]==-1){
			index=i;
		}
	}
	int placement=-1;
	for(i=number_of_sequences; i>=0; i--){
		if ( clusters[index+1][i][0]=='\0' ){
			placement=i;
		}
	}
	strcpy(clusters[index+1][placement],tree[whichTree][node].name);
}
void findPathToRoot(node** tree, int node, int whichTree, int* foundPath){
	int parent = tree[whichTree][node].down;
	if (parent==-1){ return; }
	if (tree[whichTree][node].nodeToCut==1){
		foundPath[0]=1;
	}
	findPathToRoot(tree,parent,whichTree,foundPath);
}
void findLeavesOfNodeCut(node** tree, int node, int whichTree, int clusterNumber, int kseqs, int firstiter){
	int child1 = tree[whichTree][node].up[0];
	int child2 = tree[whichTree][node].up[1];
	int i;
	if (firstiter>0 && tree[whichTree][node].nodeToCut==1){
		return;
	}
	firstiter++;
	tree[whichTree][node].nodeToCut=1;
	if (tree[whichTree][node].up[0]==-1 && tree[whichTree][node].up[1]==-1){
		printf("clusterNumber: %d name: %s node: %d\n",clusterNumber,tree[whichTree][node].name,node);
		int count=0;
		for(i=kseqs-1; i>=0; i--){
			if (clusters[clusterNumber][i][0]=='\0'){
				count=i;
			}
		}
		strcpy(clusters[clusterNumber][count],tree[whichTree][node].name);
		return;
	}else{
		findLeavesOfNodeCut(tree,child1,whichTree,clusterNumber,kseqs,firstiter);
		findLeavesOfNodeCut(tree,child2,whichTree,clusterNumber,kseqs,firstiter);
	}
}
int findNumberOfClusters(){
	int i;
	for(i=0; i<MAXNUMBEROFCLUSTERS; i++){
		if (clusters[i][0][0]=='\0'){
			return i;
		}
	}
}
void shiftColumns(int kseqs){
	int i,j,k;
	char **tmp = (char**)malloc(kseqs*sizeof(char *));
	for(i=0; i<kseqs; i++){
		tmp[i] = (char *)malloc(MAXNAME*sizeof(char));
	}
	for(i=0; i<MAXNUMBEROFCLUSTERS; i++){
		if (clusters[i][0][0]=='\0'){
			for(j=i; j<MAXNUMBEROFCLUSTERS-1; j++){
				for(k=0; k<kseqs; k++){
					strcpy(tmp[k],clusters[j+1][k]);
					strcpy(clusters[j][k],tmp[k]);
					//strcpy(clusters[j+1][k],'\0');
					memset(clusters[j+1][k],'\0',MAXNAME);
				}
			}
		}
	}
	for(i=0; i<kseqs; i++){
		free(tmp[i]);
	}
	free(tmp);
}
void createTreesForClusters(node** treeArr, int number_of_clusters, int* clusterSize, char*** cluster_seqs, char*** clusters, int* fasta_specs, int* rootArr){
	//int rootArr[number_of_clusters];
	int i=0;
	int j=0;
	int k=0;
	for(i=1; i<number_of_clusters; i++){
		double** distMat = (double **)malloc((clusterSize[i]+1)*sizeof(double *));
		for(j=0; j<clusterSize[i]+1; j++){
			distMat[j] = (double *)malloc((clusterSize[i]+1)*sizeof(double));
		}
		for(j=0; j<clusterSize[i]+1; j++){
			for(k=0; k<clusterSize[i]+1; k++){
				distMat[j][k]=0;
			}
		}
		createDistMat(cluster_seqs[i],distMat,clusterSize[i],fasta_specs);
		if ( clusterSize[i] > 3 ){
			rootArr[i-1] = NJ(treeArr,distMat,clusterSize[i],i-1,i);
		}else{
			rootArr[i-1]=0;
		}
		for(j=0; j<clusterSize[i]+1; j++){
			free(distMat[j]);
		}
	}
}
void clearGlobals(){
	int i,j, k;
	for(i=0;i<STATESPACE;i++){
		RRVAL[i]=0;
		for(j=0; j<STATESPACE;j++){
			LRVEC[i][j]=0;
			RRVEC[i][j]=0;
			PMAT1[i][j]=0;
			PMAT2[i][j]=0;
		}
	}
	for(i=4;i<4;i++){
		RRVALnc[i]=0;
		for(i=4;i<4;i++){
			LRVECnc[i][j]=0;
			RRVECnc[i][j]=0;
		}
	}
	for (i=0;i<4;i++){
		PMATnc[0][i][4]=1.0;//missing data
		PMATnc[1][i][4]=1.0;//missing data
	}
	parameters[0]=0.0;
	for(i=1;i<10;i++){
		parameters[i]=1.0;
	}
}
void inittransitionmatrixnc(double pi[4]){
	int i, j;
	double sum, piT, RIVAL[4], RIVEC[4][4],  A[4][4], workspace[8];
	for (i=0; i<8; i++){
		workspace[i]=0;
	}
	A[0][1]=pi[1]*parameters[4];
	A[0][2]=pi[2]*parameters[5];
	A[0][3]=pi[3]*parameters[6];
	A[1][0]=pi[0]*parameters[4];
	A[1][2]=pi[2]*parameters[7];
	A[1][3]=pi[3]*parameters[8];
	A[2][0]=pi[0]*parameters[5];
	A[2][1]=pi[1]*parameters[7];
	A[2][3]=pi[3]; //unscaled rate of GT = 1.0
	A[3][0]=pi[0]*parameters[6];
	A[3][1]=pi[1]*parameters[8];
	A[3][2]=pi[2]; //unscaled rate of GT = 1.0
	for (i=0; i<4; i++){
		A[i][i]=0.0;
		sum=0.0;
		for (j=0; j<4; j++){
			sum = sum + A[i][j];
		}
		A[i][i] = -sum;
	}
	if (eigen(1, A[0], 4, RRVALnc, RIVAL, RRVECnc[0], RIVEC[0], workspace) != 0){
		printf("Transitions matrix did not converge or contained non-real values!\n");
		exit(-1);
	}
	for (i=0; i<4; i++){
		for (j=0; j<4; j++){
			LRVECnc[i][j] = RRVECnc[i][j];
		}
	}
	if (matinv(RRVECnc[0],4, 4, workspace) != 0){
		printf("Could not invert matrix!\nResults may not be reliable!\n");
	}
}
void maketransitionmatrixnc(int n, double t, int whichRoot){
	int i, j, k;
	double EXPOS[4];
	for (k=0; k<4; k++){
		EXPOS[k] = exp(t*RRVALnc[k]);
		if ( EXPOS[k] < 0.000000 ){
			printf("EXPOS[k]=%e\n",EXPOS[k]);
		}
	}
	for (i=0; i<4; i++){
		for (j=0; j<4; j++){
			PMATnc[n][i][j] = 0.0;
			for (k=0; k<4; k++){
				PMATnc[n][i][j] =  PMATnc[n][i][j] + RRVECnc[k][j]*LRVECnc[i][k]*EXPOS[k];
			}
		}
	}
}
void makeconnc(int node, double lambda, int whichRoot, int numbase, int numspec, int*** seqArr, int root){
	int i, j, child, seqn, site;
	double L, max;
	child = treeArr[whichRoot][node].up[0];
	if (treeArr[whichRoot][child].up[0]==-1){
		maketransitionmatrixnc(0, lambda*treeArr[whichRoot][child].bl,whichRoot);
		//seqn=child-numspec+1;
		seqn=child;
		for (site=0; site<numbase; site++){
			for (i=0; i<4; i++){
				treeArr[whichRoot][node].likenc[site][i] = PMATnc[0][i][seqArr[whichRoot][seqn][site]];
			}
			//if (site==0 && node ==root) printf("node 10 (b=%i): %lf %lf %lf %lf\n",seqArr[whichRoot][seqn][site],PMATnc[0][0][seqArr[whichRoot][seqn][site]],PMATnc[0][1][seqArr[whichRoot][seqn][site]],PMATnc[0][2][seqArr[whichRoot][seqn][site]],PMATnc[0][3][seqArr[whichRoot][seqn][site]]);
		}
	}else{
		makeconnc(child, lambda, whichRoot, numbase, numspec, seqArr, root);
		maketransitionmatrixnc(0, lambda*treeArr[whichRoot][child].bl,whichRoot);
		for (site=0; site<numbase; site++){
			for (i=0; i<4; i++){
				treeArr[whichRoot][node].likenc[site][i]=0.0;
				for (j=0; j<4; j++){
					treeArr[whichRoot][node].likenc[site][i] += PMATnc[0][i][j]*treeArr[whichRoot][child].likenc[site][j];
				}
			}
		}

	}
	child = treeArr[whichRoot][node].up[1];
	if (treeArr[whichRoot][child].up[1]==-1){
		maketransitionmatrixnc(0,lambda*treeArr[whichRoot][child].bl,whichRoot);
		//seqn=child-numspec+1;
		seqn=child;
		for (site=0; site<numbase; site++){
			for (i=0; i<4; i++){
				treeArr[whichRoot][node].likenc[site][i] = treeArr[whichRoot][node].likenc[site][i]*PMATnc[0][i][seqArr[whichRoot][seqn][site]];
			}
		}
	}else{
		makeconnc(child, lambda,whichRoot, numbase, numspec, seqArr, root);
		maketransitionmatrixnc(0,lambda*treeArr[whichRoot][child].bl,whichRoot);
		for (site=0; site<numbase; site++){
			max=0.0;
			for (i=0; i<4; i++){
				L=0.0;
				for (j=0; j<4; j++){
					L += PMATnc[0][i][j]*treeArr[whichRoot][child].likenc[site][j];
				}
				if ((treeArr[whichRoot][node].likenc[site][i] = treeArr[whichRoot][node].likenc[site][i]*L)>max){
					max = treeArr[whichRoot][node].likenc[site][i];
				}
			}
			if (max<0.00000000001) printf("Warning, max = %lf\n",max);
			for (i=0; i<4; i++){
				treeArr[whichRoot][node].likenc[site][i]=treeArr[whichRoot][node].likenc[site][i]/max;
			}
			UFCnc[site] = UFCnc[site] + log(max);
		}
	}
}
double getlike_gamma(double par[],int whichRoot, int numbase, int root, int numspec, int*** seqArr){
	double stand, L, loclike, **locloglike, max, pi[4], gampar[2], d, like = 0.0;
	int i, j, k;
	COUNT2++;
	stand = 1.0+par[1]+par[2]+par[3];
	pi[0]=par[1]/stand;
	pi[1]=par[2]/stand;
	pi[2]=par[3]/stand;
	pi[3]=1.0-pi[0]-pi[1]-pi[2];
	gampar[0]=gampar[1]=par[9]; //We are setting alpha=beta to keep a constant mean to avoid identifiability issues.  This is not the same as a standard gammma.
	UFCnc = malloc(numbase*(sizeof(double)));
	statevector = malloc(NUMCAT*(sizeof(double)));
	locloglike = malloc(numbase*(sizeof(double *)));
	for (i=0; i<numbase; i++){
		locloglike[i] = malloc(NUMCAT*(sizeof(double)));
	}
	definegammaquantiles(NUMCAT, gampar);
	statevector[0]=1.0;
	inittransitionmatrixnc(pi);
	for (j=0; j<NUMCAT; j++){
		for (i=0; i<numbase; i++){
			UFCnc[i]=0.0;
		}
		makeconnc(root, statevector[j],whichRoot,numbase,numspec,seqArr,root);
		for (i=0; i<numbase; i++){
			L=0.0;
			for (k=0;k<4;k++){
				L += treeArr[whichRoot][root].likenc[i][k]*pi[k];
			}
			if (L>0.0) locloglike[i][j] = log(L) + UFCnc[i];
		}
	}
	for (i=0; i<numbase; i++){
		loclike=0.0;
		max = -100000000000.0;
		for (j=0; j<NUMCAT; j++){
			if (locloglike[i][j]>max){ //underflow protection
				max=locloglike[i][j];
			}
		}
		for (j=0; j<NUMCAT; j++){
			d=locloglike[i][j]-max;
			if (d>-100){
				loclike += exp(d);
			}
		}
		like = like + log(loclike) + max;
	}
	free(statevector);
	for (i=0; i<numbase; i++){
		free(locloglike[i]);
	}
	free(locloglike);
	//printf("LIKE: %lf\n",like - (double)numbase*log((double)NUMCAT));
	//printf("\n");
	free(UFCnc);
	return -like + (double)numbase*log((double)NUMCAT);
}
double like_bl_Arr(double par[2], int whichRoot, int numbase, int root, int numspec, int*** seqArr){
	int i, j, s, base;
	double b, p, L=0.0;
	maketransitionmatrixnc(0,par[1],whichRoot);
	for (s=0; s<numbase; s++){
		p=0.0;
		if (treeArr[whichRoot][localnode].up[0]==-1){
			//base = seqArr[whichRoot][localnode-numspecArr[whichRoot]+1][s];
			base = seqArr[whichRoot][localnode][s];
			assert(base >= 0 && base <= 4);
			if (base<4){
				for (j=0; j<4; j++){
					p += localpi[base]*PMATnc[0][base][j]*templike_nc[s][j];
				}
			}else{
				p=1.0;
			}
		}else{
			for (i=0; i<4; i++){
				b=0;
				for (j=0; j<4; j++){
					b += PMATnc[0][i][j]*templike_nc[s][j];
				}
				p += b*localpi[i]*treeArr[whichRoot][localnode].likenc[s][i];
			}
		}
		L += log(p);
	} COUNT++;
	return -L;
}
void maxbl_nc(int node, int parent, double pi[4], int precision, int whichRoot, int numbase, int*** seqArr, int root, int numspec){
	double par[2], minpar[2], maxpar[2], L;
	par[1]=treeArr[whichRoot][node].bl; /*This stuff should probably be cleaned up*/
	minpar[1]=MINBL;
	maxpar[1]=MAXBL;
	localpi=pi;
	localnode=node;
	L = findmax_Arr(par, minpar, maxpar, 1, like_bl_Arr, precision, whichRoot, numbase, root, numspec, seqArr);
	treeArr[whichRoot][node].bl = par[1];
}
void recurse_estimatebranchlengths(int node, double pi[4], int precision, int whichRoot, int numbase, int*** seqArr, int root, int numspec){
	int i,j, s, parent, otherb, child1, child2;
	double max, bl;
	child1 = treeArr[whichRoot][node].up[0];
	parent = treeArr[whichRoot][node].down;
	bl = treeArr[whichRoot][node].bl;
	if ((otherb = treeArr[whichRoot][parent].up[0])==node){
		otherb = treeArr[whichRoot][parent].up[1];
	}
	maketransitionmatrixnc(1, treeArr[whichRoot][otherb].bl,whichRoot);
	for (s=0; s<numbase; s++){
		for (i=0; i<4; i++){
			max=0.0;
			if (treeArr[whichRoot][otherb].up[0]==-1){
				//templike_nc[s][i] = PMATnc[1][i][seqArr[whichRoot][otherb-numspec+1][s]]; /*if(s==0) printf("temp[s][%i]: %lf (b=%i, %lf), ",i,templike_nc[s][i],seq[otherb-numspec+1][s], PMATnc[1][i][seq[otherb-numspec+1][s]]);*/
				templike_nc[s][i] = PMATnc[1][i][seqArr[whichRoot][otherb][s]]; /*if(s==0) printf("temp[s][%i]: %lf (b=%i, %lf), ",i,templike_nc[s][i],seq[otherb-numspec+1][s], PMATnc[1][i][seq[otherb-numspec+1][s]]);*/
			}else{
				templike_nc[s][i]=0.0;
				for (j=0; j<4; j++){
					templike_nc[s][i]=templike_nc[s][i]+treeArr[whichRoot][otherb].likenc[s][j]*PMATnc[1][i][j];
				}
			}
			if ((templike_nc[s][i]=templike_nc[s][i]*treeArr[whichRoot][parent].posteriornc[s][i])>max){
				max=templike_nc[s][i];
			}
		}
	}
	maxbl_nc(node,  parent, pi,precision,whichRoot, numbase, seqArr, root, numspec);
	maketransitionmatrixnc(0, treeArr[whichRoot][node].bl,whichRoot);
	for (s=0; s<numbase; s++){
		max=0.0;
		for (i=0; i<4; i++){
			treeArr[whichRoot][node].posteriornc[s][i]=0.0;
			for (j=0; j<4; j++){
				treeArr[whichRoot][node].posteriornc[s][i] = treeArr[whichRoot][node].posteriornc[s][i] + PMATnc[0][i][j]*templike_nc[s][j];
			}
			if (treeArr[whichRoot][node].posteriornc[s][i]>max){//more hysterical underflow protection
				max=treeArr[whichRoot][node].posteriornc[s][i];
			}
		}
		for (i=0; i<4; i++){
			treeArr[whichRoot][node].posteriornc[s][i]=treeArr[whichRoot][node].posteriornc[s][i]/max;
		}
	}
	if (child1>-1){
		child2 = treeArr[whichRoot][node].up[1];
		recurse_estimatebranchlengths(child1, pi, precision, whichRoot, numbase, seqArr, root, numspec);
		recurse_estimatebranchlengths(child2, pi, precision, whichRoot, numbase, seqArr, root, numspec);
	}
}
void estimatebranchlengths(double par[10], int precision, int whichRoot, int numbase, int root, int numspec, int*** seqArr){
	int i, j, s, child1, child2;
	double stand, pi[4];
	for (i=0; i<10; i++){
		currentestimate[i]=par[i];
	}
	doNRinits(1);
	stand = 1.0+par[1]+par[2]+par[3];
	pi[0]=par[1]/stand;
	pi[1]=par[2]/stand;
	pi[2]=par[3]/stand;
	pi[3]=1.0-pi[0]-pi[1]-pi[2];
	child1 = treeArr[whichRoot][root].up[0];
	child2 = treeArr[whichRoot][root].up[1];
	templike_nc = malloc(numbase*(sizeof(double *)));
	for (i=0; i<numbase; i++){
		templike_nc[i]=malloc(4*(sizeof(double)));
	}
	if ((treeArr[whichRoot][child2].bl=treeArr[whichRoot][child2].bl+treeArr[whichRoot][child1].bl-MINBL)<MINBL){
		treeArr[whichRoot][child2].bl=MINBL;
	}
	treeArr[whichRoot][child1].bl=MINBL;
	getlike_gamma(par,whichRoot,numbase,root,numspec,seqArr); /*this is not necessary if the likelihood fucntion has already been called*/
	maketransitionmatrixnc(0, treeArr[whichRoot][child2].bl+treeArr[whichRoot][child1].bl,whichRoot);
	for (s=0; s<numbase; s++){
		for (i=0; i<4; i++){
			treeArr[whichRoot][root].posteriornc[s][i] = 1.0;/*tree[child1].likenc[s][i];*/
			if (treeArr[whichRoot][child2].up[0]>-1){
				treeArr[whichRoot][child1].posteriornc[s][i]=0.0;
				for (j=0; j<4; j++){
					treeArr[whichRoot][child1].posteriornc[s][i] += treeArr[whichRoot][child2].likenc[s][j]*PMATnc[0][i][j];
				}
			}else{
				//treeArr[whichRoot][child1].posteriornc[s][i] = PMATnc[0][i][seqArr[whichRoot][child2-numspec+1][s]];
				treeArr[whichRoot][child1].posteriornc[s][i] = PMATnc[0][i][seqArr[whichRoot][child2][s]];
			}
		}
		if (s==0){
			//printf("node %d initialization: %lf %lf %lf %lf\n",root,treeArr[whichRoot][root].posteriornc[s][0],treeArr[whichRoot][root].posteriornc[s][1],treeArr[whichRoot][root].posteriornc[s][2],treeArr[whichRoot][root].posteriornc[s][3]);
			//printf("Initial child 1 like: %lf %lf %lf %lf\n",treeArr[whichRoot][child1].likenc[s][0],treeArr[whichRoot][child1].likenc[s][1],treeArr[whichRoot][child1].likenc[s][2],treeArr[whichRoot][child1].likenc[s][3]);
		}
	}
	if (treeArr[whichRoot][child1].up[0]>-1){
		recurse_estimatebranchlengths(treeArr[whichRoot][child1].up[0], pi, precision, whichRoot, numbase, seqArr, root, numspec);
	}
	if (treeArr[whichRoot][child1].up[1]>-1){
		recurse_estimatebranchlengths(treeArr[whichRoot][child1].up[1], pi, precision, whichRoot, numbase, seqArr, root, numspec);
	}
	recurse_estimatebranchlengths(child2, pi, precision,whichRoot, numbase, seqArr, root, numspec);
	for (i=0; i<numbase; i++){
		free(templike_nc[i]);
	}
	free(templike_nc);
	freeNRinits(1);
}
double maximizelikelihoodnc_globals(double parameters[10], int precision, int whichRoot, int numbase, int*** seqArr, int root, int numspec){
	//optimization function starts counting at 1 so arrays have dimensionality n+1
	int i;
	double L, lowbound[10], upbound[10];
	doNRinits(10);
	for (i=1; i<10; i++){
		lowbound[i]=0.05;
		upbound[i]=20.0;
	}
	lowbound[9]=0.3;
	L=findmax_Arr(parameters, lowbound, upbound, 9, getlike_gamma, precision, whichRoot, numbase, root, numspec, seqArr);
	freeNRinits(10);
	return L;
}
void print_branch_lengths(node** treeArr, int node, int whichRoot){
	printf("NODE: %d BL: %lf\n",node,treeArr[whichRoot][node].bl);
	if ( treeArr[whichRoot][node].up[0] != -1){
		print_branch_lengths(treeArr,treeArr[whichRoot][node].up[0],whichRoot);
		print_branch_lengths(treeArr,treeArr[whichRoot][node].up[1],whichRoot);
	}
}
void estimatenucparameters(int whichRoot, int numbase, int root, int numspec, int*** seqArr){
	double L;
	COUNT=COUNT2=0;
	clearGlobals();
	//double precision[1];
	//precision[0]=0;
	//treeArr[whichRoot][root].bl = 0.00001;
	//print_branch_lengths(treeArr,root,whichRoot);
	printf("Initial likelihoodL value = %lf\n",-getlike_gamma(parameters,whichRoot,numbase,root,numspec,seqArr));
	//print_branch_lengths(treeArr,root,whichRoot);
	estimatebranchlengths(parameters,0, whichRoot, numbase, root, numspec,seqArr);
	L=maximizelikelihoodnc_globals(parameters,0,whichRoot,numbase,seqArr,root,numspec);
	//printf("Current ML value = %lf\n",L);
	estimatebranchlengths(parameters,0, whichRoot, numbase, root, numspec,seqArr);
	L=maximizelikelihoodnc_globals(parameters,0,whichRoot,numbase,seqArr,root,numspec);
	//print_branch_lengths(treeArr,root,whichRoot);
	//printf("Current ML value = %lf\n",L);
	estimatebranchlengths(parameters,1, whichRoot, numbase, root, numspec,seqArr);
	L=maximizelikelihoodnc_globals(parameters,0,whichRoot,numbase,seqArr,root,numspec);
	//printf("Current ML value = %lf\n",L);
	estimatebranchlengths(parameters,2, whichRoot, numbase, root, numspec,seqArr);
	L=maximizelikelihoodnc_globals(parameters,2,whichRoot,numbase,seqArr,root,numspec);
	//printf("Current ML value = %lf\n",L);
	estimatebranchlengths(parameters,2, whichRoot, numbase, root, numspec,seqArr);
	estimatebranchlengths(parameters,2, whichRoot, numbase, root, numspec,seqArr);
	printf("Current ML value= %lf\n",-getlike_gamma(parameters,whichRoot,numbase,root,numspec,seqArr));
}
void makeposterior_nc(int node, int whichRoot, int numbase, int*** seqArr){
	int i,j, s, parent, otherb, child1, child2, b;
	double bl, max;
	child1 = treeArr[whichRoot][node].up[0];
	child2 = treeArr[whichRoot][node].up[1];
	parent = treeArr[whichRoot][node].down;
	bl = treeArr[whichRoot][node].bl;
	maketransitionmatrixnc(0, bl,whichRoot);
	if ((otherb = treeArr[whichRoot][parent].up[0])==node){
		otherb = treeArr[whichRoot][parent].up[1];
	}
	maketransitionmatrixnc(1, treeArr[whichRoot][otherb].bl,whichRoot);
	for (s=0; s<numbase; s++){
		if (treeArr[whichRoot][otherb].up[0]>-1){
			for (i=0; i<4; i++){
				templike_nc[s][i]=0;
				for (j=0; j<4; j++){
					templike_nc[s][i] += treeArr[whichRoot][otherb].likenc[s][j]*PMATnc[1][i][j];
				}
				templike_nc[s][i]=templike_nc[s][i]*treeArr[whichRoot][parent].posteriornc[s][i];
			}
		}else{
			b=seqArr[whichRoot][otherb][s];
			for (i=0; i<4; i++){
				templike_nc[s][i] = PMATnc[1][i][b]*treeArr[whichRoot][parent].posteriornc[s][i];
			}
		}
		for (i=0; i<4; i++){
			treeArr[whichRoot][node].posteriornc[s][i]=0.0;
			max=0.0;
			for (j=0; j<4; j++){
				if ((treeArr[whichRoot][node].posteriornc[s][i] = treeArr[whichRoot][node].posteriornc[s][i] + PMATnc[0][i][j]*templike_nc[s][j])>max){
					max=treeArr[whichRoot][node].posteriornc[s][i];//more underflow protection
				}
			}
		}
		for (i=0; i<4; i++){
			treeArr[whichRoot][node].posteriornc[s][i]=treeArr[whichRoot][node].posteriornc[s][i]/max;
		}
	}
		if (treeArr[whichRoot][child1].up[0]>-1){
			makeposterior_nc(child1,whichRoot,numbase,seqArr);
		}
		if (treeArr[whichRoot][child2].up[0]>-1){
			makeposterior_nc(child2,whichRoot,numbase,seqArr);
		}
}
void getposterior_nc(int whichRoot, int numbase, int root, int numspec, int*** seqArr){
	int i, j, s, k, parent, b, notdonebefore;
	double p, sum, pi[4], stand, **templike;
	getlike_gamma(parameters,whichRoot,numbase,root,numspec,seqArr); /*need to call likelihood again*/
	templike_nc = malloc(numbase*(sizeof(double *)));
	for (i=0; i<numbase; i++){
		templike_nc[i]=malloc(4*(sizeof(double)));
	}
	stand = 1.0+parameters[1]+parameters[2]+parameters[3];
	pi[0]=parameters[1]/stand;
	pi[1]=parameters[2]/stand;
	pi[2]=parameters[3]/stand;
	pi[3]=1.0-pi[0]-pi[1]-pi[2];
	for (s=0; s<numbase; s++){
		for (i=0; i<4; i++){
			treeArr[whichRoot][root].posteriornc[s][i] = 1.0;
		}
	}
	if (treeArr[whichRoot][treeArr[whichRoot][root].up[0]].up[0]>-1){
		makeposterior_nc(treeArr[whichRoot][root].up[0],whichRoot,numbase,seqArr);
	}
	if (treeArr[whichRoot][treeArr[whichRoot][root].up[1]].up[0]>-1){
		makeposterior_nc(treeArr[whichRoot][root].up[1],whichRoot,numbase,seqArr);
	}
	for (j=0; j<2*numspec-1; j++){
		if (treeArr[whichRoot][j].up[0]>-1){
			for (s=0; s<numbase; s++) {
				sum = 0.0;
				for (i=0; i<4; i++){
					sum = sum + (treeArr[whichRoot][j].posteriornc[s][i]=treeArr[whichRoot][j].likenc[s][i]*treeArr[whichRoot][j].posteriornc[s][i]*pi[i]);
				}
				for (i=0; i<4; i++){
					treeArr[whichRoot][j].posteriornc[s][i] = treeArr[whichRoot][j].posteriornc[s][i]/sum;
				}
			}
		}else{
			for (s=0; s<numbase; s++) {
				notdonebefore=1;
				b=seqArr[whichRoot][j][s];
				if (b==4){
					if (notdonebefore==1) {
						maketransitionmatrixnc(0, treeArr[whichRoot][j].bl,whichRoot);
						notdonebefore=0;
					}
					parent=treeArr[whichRoot][j].down;
					sum = 0.0;
					for (i=0; i<4; i++){
						treeArr[whichRoot][j].posteriornc[s][i]=0.0;
						for (k=0; k<4; k++){
							treeArr[whichRoot][j].posteriornc[s][i] += treeArr[whichRoot][parent].posteriornc[s][k]*PMATnc[0][i][k];
						}
						sum = sum + (treeArr[whichRoot][j].posteriornc[s][i]=treeArr[whichRoot][j].posteriornc[s][i]*pi[i]);
					}
					for (i=0; i<4; i++){
						treeArr[whichRoot][j].posteriornc[s][i] = treeArr[whichRoot][j].posteriornc[s][i]/sum;
					}
				}else{
					for (i=0; i<4; i++){
						if (i==b){
							treeArr[whichRoot][j].posteriornc[s][i]=1.0;
						}else{
							treeArr[whichRoot][j].posteriornc[s][i]=0.0;
						}
					}
				}
			}
		}
	}
	for (i=0; i<numbase; i++){
		free(templike_nc[i]);
	}
	free(templike_nc);

}
void initialize_assignment_mem(type_of_PP**** PP, int numberOfRoots,int* numspec, int* numbase){
	int i, j, k;
	PP = (type_of_PP ****)malloc(numberOfRoots*sizeof(type_of_PP ***));
		for(i=0; i<numberOfRoots;i++){
			PP[i] = (type_of_PP ***)malloc((2*numspec[i+1]-1)*sizeof(type_of_PP **));
			for(j=0; j<2*numspec[i+1]-1; j++){
				PP[i][j] = (type_of_PP **)malloc((numbase[j])*sizeof(type_of_PP *));
				for(k=0; k<numbase[j];k++){
					PP[i][j][k] = (type_of_PP *)malloc(4*(sizeof(type_of_PP)));
				}
			}
		}
}
void store_PPs(type_of_PP**** PP, int numberOfRoots,int* numspec, int* numbase){
	int i, j, k, l;
	for(i=0; i<numberOfRoots; i++){
		if ( numspec[i+1] > 3){
		for (j=0; j<2*numspec[i+1]-1; j++){
			for (k=0; k<numbase[i]; k++){
				for (l=0; l<4; l++){
					//if ( treeArr[i][j].posteriornc[k][l] == -1 ){
					//	PP[i][j][k][l] = -1;
					//}else{
						PP[i][j][k][l] = (type_of_PP)(1.0-treeArr[i][j].posteriornc[k][l]);//This need to change if we change type for this variable
					//}
				}
			}
		}
		}
	}
}
void printRootSeqs(char** rootSeqs, type_of_PP**** PP, int numberOfRoots,int* numbase, int* rootArr, int *clusterSize, char*** cluster_seqs){
	type_of_PP minimum;
	int index,i,j,k;
	for(i=0; i<numberOfRoots;i++){
		printf("NUMBASE: %d\n",numbase[i]);
		if ( clusterSize[i+1] > 3){
		for(j=0;j<numbase[i];j++){
			minimum=PP[i][rootArr[i]][j][0];
			index=0;
			for(k=0;k<4;k++){
				if (minimum > PP[i][rootArr[i]][j][k]){
					minimum=PP[i][rootArr[i]][j][k];
					index=k;
				}
				//if( PP[i][rootArr[i]][j][k] == -1){
				//	index=-1;
				//}
			}
			char base;
			if (index==0){ base='A'; }
			else if (index ==1 ){ base='C'; }
			else if (index==2 ){ base='G'; }
			else if (index==3 ){base='T'; }
			rootSeqs[i][j]=base;
		}
		rootSeqs[i][numbase[i]]='\0';
		}else{
			strcpy(rootSeqs[i],cluster_seqs[i+1][0]);
		}
	}
	for(i=0;i<numberOfRoots;i++){
		printf("Root Sequence %d\t",i);
		for(j=0;j<numbase[i];j++){
			printf("%c",rootSeqs[i][j]);
		}
		printf("\n");
	}
}
int main(int argc, char **argv){
	Options opt;
	opt.number_of_clusters = 10;
	opt.number_of_kseqs=50;
	opt.slash=0;
	opt.default_directory=1;
	opt.numthreads=1;
	opt.hasTaxFile=0;
	opt.output_fasta=0;
	opt.clstr_format=1;
	opt.use_nw=0;
	strcpy(opt.output_directory,"");
	memset(opt.output_file,'\0',2000);
	parse_options(argc, argv, &opt);
	FILE* fasta_for_clustering;
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
	printf("Number of threads: %d\n",opt.numthreads);
	opt.numberOfLinesToRead=fasta_specs[0]/opt.numthreads;
	int i,j;
	char** seqNames = (char **)malloc(fasta_specs[0]*sizeof(char *));
	char** sequences = (char **)malloc(fasta_specs[0]*sizeof(char *));
	for(i=0; i<fasta_specs[0]; i++){
		seqNames[i]=(char *)malloc((fasta_specs[2]+1)*sizeof(char));
		sequences[i]=(char *)malloc((fasta_specs[1]+1)*sizeof(char));
	}
	//hashmap_init(&map, hashmap_hash_string, hashmap_compare_string, fasta_specs[0]);
	if (( fasta_for_clustering = fopen(opt.fasta,"r")) == (FILE *) NULL ) fprintf(stderr,"FASTA file could not be opened.\n");
	readInFasta(fasta_for_clustering,seqNames,sequences);
	fclose(fasta_for_clustering);
	char** taxonomy = (char **)malloc(fasta_specs[0]*sizeof(char *));
	for(i=0; i<fasta_specs[0]; i++){
		taxonomy[i]=(char *)malloc(FASTA_MAXLINE*sizeof(char));
	}
	FILE* taxonomyFile;
	struct hashmap taxMap;
	if (opt.hasTaxFile==1){
		hashmap_init(&taxMap, hashmap_hash_string, hashmap_compare_string, 4*fasta_specs[0]);
		if (( taxonomyFile = fopen(opt.taxonomy,"r")) == (FILE *) NULL ) fprintf(stderr,"TAXONOMY file could not be opened.\n");
		readInTaxFile(taxonomyFile,taxonomy,taxMap);
		fclose(taxonomyFile);
	}
	int* chooseK = (int *)malloc(MAXNUMBEROFKSEQS*sizeof(int));
	int numberOfUnAssigned=fasta_specs[0];
	int kseqs = opt.number_of_kseqs;
	char** sequencesToClusterLater = (char**)malloc(fasta_specs[0]*sizeof(char *));
	for(i=0; i<fasta_specs[0]; i++){
		sequencesToClusterLater[i] = (char *)malloc(MAXNAME*sizeof(char));
	}
	clusters = (char***)malloc(MAXNUMBEROFCLUSTERS*sizeof(char **));
	char*** cluster_seqs = (char***)malloc(MAXNUMBEROFCLUSTERS*sizeof(char **));
	for(i=0; i<MAXNUMBEROFCLUSTERS; i++){
		clusters[i]=(char **)malloc((fasta_specs[0])*sizeof(char *));
		cluster_seqs[i]=(char **)malloc((fasta_specs[0])*sizeof(char *));
		for(j=0; j<fasta_specs[0]; j++){
			clusters[i][j]=(char *)malloc(MAXNAME*sizeof(char));
			clusters[i][j][0]='\0';
			cluster_seqs[i][j]=(char *)malloc((fasta_specs[1]+1)*sizeof(char));
			cluster_seqs[i][j][0]='\0';
		}
	}
	int numberOfNodesToCut=0;
	int numberOfSequencesLeft=fasta_specs[0];
	int starting_number_of_clusters=1;
	char*** clstr;
	int** clstr_lengths;
	//struct hashmap assignedSeqs;
	//hashmap_init(&assignedSeqs, hashmap_hash_string, hashmap_compare_string, 4*fasta_specs[0]);
	char** assigned = (char **)malloc(fasta_specs[0]*sizeof(char **));
	for(i=0; i<fasta_specs[0]; i++){
		assigned[i] = (char *)malloc(MAXNAME*sizeof(char));
		assigned[i][0]='\0';
	}
	//START THE LOOP
	//numberOfUnAssigned=0;
	while(numberOfUnAssigned>3){
	printf("starting number of clusters: %d\n",starting_number_of_clusters);
	//struct hashmap seqsToCompare;
	//hashmap_init(&seqsToCompare, hashmap_hash_string, hashmap_compare_string, numberOfSequencesLeft);
	char*** cluster_seqs = (char***)malloc(MAXNUMBEROFCLUSTERS*sizeof(char **));
	for(i=0; i<MAXNUMBEROFCLUSTERS; i++){
		cluster_seqs[i]=(char **)malloc((fasta_specs[0])*sizeof(char *));
		for(j=0; j<fasta_specs[0]; j++){
			cluster_seqs[i][j]=(char *)malloc((fasta_specs[1]+1)*sizeof(char));
			cluster_seqs[i][j][0]='\0';
		}
	}
	//clearSequencesToAssignLater(sequencesToClusterLater,fasta_specs);
	numberOfNodesToCut=0;
	for(i=0; i<MAXNUMBEROFKSEQS; i++){
		chooseK[i]=-1;
	}
	if (numberOfUnAssigned < fasta_specs[0] ){
		//kseqs = (numberOfUnAssigned + (4-1))/4;
		//opt.number_of_clusters = (numberOfUnAssigned + (2-1))/2;
		if (kseqs > numberOfUnAssigned){
			kseqs=numberOfUnAssigned;
		}
		//opt.number_of_clusters = numberOfUnAssigned - 1;
		//opt.number_of_clusters = kseqs-1;
	}
	printf("kseq is %d\n",kseqs);
	printf("number_of_clusters is %d\n",opt.number_of_clusters);
	struct timespec tstart={0,0}, tend={0,0};
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	printf("Choosing %d sequences at random...\n",kseqs);
	int choosing=0;
		srand(time(0));
	while(choosing==0){
		int random_number = generateRandom(numberOfSequencesLeft-1);
		int index=0;
		for(i=kseqs; i>=0; i--){
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
			//hashmap_put(&seqsToCompare,seqNames[random_number],1);
			//hashmap_put(&assignedSeqs,seqNames[random_number],1);
			int place=-1;
			for(j=fasta_specs[0]-1; j>=0; j--){
				if ( assigned[j][0]=='\0' ){
					place=j;
				}
			}
			strcpy(assigned[place],seqNames[random_number]);
		}
		if (chooseK[kseqs-1]!=-1){
			choosing=1;
		}
	}
	for(i=0; i<kseqs; i++){
		strcpy(clusters[numberOfNodesToCut][i],seqNames[chooseK[i]]);
		strcpy(cluster_seqs[numberOfNodesToCut][i],sequences[chooseK[i]]);
	}
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %lf seconds\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	double** distMat = (double **)malloc((kseqs+1)*sizeof(double *));
	for(i=0; i<kseqs+1; i++){
		distMat[i] = (double *)malloc((kseqs+1)*sizeof(double));
	}
	for(i=0; i<kseqs+1; i++){
		for(j=0; j<kseqs+1; j++){
			distMat[i][j] = 0;
		}
	}
	char** seqsInCluster = (char **)malloc(MAXNUMBEROFKSEQS*sizeof(char *));
	for(i=0; i<MAXNUMBEROFKSEQS; i++){
		seqsInCluster[i] = (char *)malloc((fasta_specs[1]+1)*sizeof(char));
		//printf("seqsInCluster[%d]: %s\n",i,seqsInCluster[i]);
	}
	for(i=0; i<kseqs; i++){
		strcpy(seqsInCluster[i],sequences[chooseK[i]]);
	}
	printf("Creating distance matrix...\n");
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	if ( opt.use_nw ==0 ){
		createDistMat_WFA(seqsInCluster,distMat,kseqs,fasta_specs);
	}else{
		createDistMat(seqsInCluster,distMat,kseqs,fasta_specs);
	}
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %lf seconds\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	int* clusterSize = (int *)malloc(fasta_specs[3]*sizeof(int));
	clusterSize[0]=kseqs;
	//print_distance_matrix(distMat,clusters,0,opt.number_of_kseqs);
	node **tree = malloc((fasta_specs[3]+1)*sizeof(node *));
	for(i=0; i<fasta_specs[3]+1; i++){
		tree[i]=(node *)malloc((2*kseqs-1)*sizeof(node));
		for(j=0; j<2*kseqs-1; j++){
			tree[i][j].name=(char *)malloc(MAXNAME*sizeof(char));
			strcpy(tree[i][j].name,"internal");
			tree[i][j].nodeToCut=0;
		}
	}
	int root = NJ(tree,distMat,kseqs,0,0);
	get_number_descendants(tree,root,0);
	assignDepth(tree,tree[0][root].up[0],tree[0][root].up[1],1);
	//printtree(tree,0,opt.number_of_kseqs);
	double* branchLengths = (double *)malloc((2*kseqs-1)*sizeof(double));
	for(i=0; i<2*kseqs-1; i++){
		branchLengths[i] = tree[0][i].bl;
	}
	int* indexArray = (int *)malloc((2*kseqs-1)*sizeof(int));
	sortArray(branchLengths,indexArray,kseqs);
	printf("Longest Branch: %lf node: %d\n",branchLengths[0],indexArray[0]);
	//numberOfNodesToCut = findMaxClade(tree,kseqs);
	for(i=0; i<2*kseqs-1; i++){
		//if (tree[0][indexArray[i]].nd > 1 && numberOfNodesToCut < opt.number_of_clusters-1){
		if ( numberOfNodesToCut < opt.number_of_clusters-1 && tree[0][indexArray[i]].nd > 10){
		//if ( tree[0][indexArray[i]].bl > 0.03 && tree[0][indexArray[i]].nd > 1){
			printf("cutting at node %d\n",indexArray[i]);
			numberOfNodesToCut++;
			//printdescendants(tree,indexArray[i],numberOfNodesToCut,0,kseqs);
			tree[0][indexArray[i]].nodeToCut=1;
			//clearDescendants(tree,tree[0][indexArray[i]].up[0],0);
			//clearDescendants(tree,tree[0][indexArray[i]].up[1],0);
			//updateNumberOfDescendants(tree,indexArray[i],tree[0][indexArray[i]].nd,0);
			//findLeavesOfNodeCut(tree,indexArray[i],0,numberOfNodesToCut,kseqs,0);	
			printf("updating tree\n");
		}
	}
	printf("made %d cuts\n",numberOfNodesToCut);
	free(indexArray);
	free(branchLengths);
	//struct hashmap clusterhash;
	//hashmap_init(&clusterhash, hashmap_hash_string, hashmap_compare_string, fasta_specs[0]);
	int* parentcuts = (int *)malloc(fasta_specs[0]*sizeof(int));
	for(i=0; i<fasta_specs[0]; i++){
		parentcuts[i]=-1;
	}
	int** clusterarr = (int **)malloc(fasta_specs[0]*sizeof(int *));
	for(i=0; i<fasta_specs[0]; i++){
		clusterarr[i] = (int *)malloc(fasta_specs[0]*sizeof(int));
		for(j=0; j<fasta_specs[0]; j++){
			clusterarr[i][j]=-1;
		}
	}
	findLeaves(tree,root,0,parentcuts,clusterarr,cluster_seqs,fasta_specs[0],seqNames,sequences);
	//int* nodeArrays = (int *)malloc(100*sizeof(int));
	//int key;
	//hashmap_iter(&clusterhash);
		//int i=0;
		//for(i=0; i<100; i++){
		//	printf("key %d: %d\n",key,nodeArrays[i]);
		//}
	//}
	int* firstNodeCut=(int *)malloc(1*sizeof(int));
	firstNodeCut[0]=-1;
	//findFirstNodeCut(tree,root,0,firstNodeCut);
	//printf("first node cut is %d\n",firstNodeCut[0]);
	//if ( firstNodeCut[0] != -1){
		//addFirstCluster(tree,root,firstNodeCut[0],0,parentcuts);
	//}
	//makeClustersFromNodeCuts(tree,root,0,1,kseqs);
	int* foundPath = malloc(1*sizeof(int));
	foundPath[0]=0;
	for(i=0; i<kseqs; i++){
		foundPath[0]=0;
		findPathToRoot(tree,i,0,foundPath);
		printf("found path node %d: %d\n",i,foundPath[0]);
		if (foundPath[0]==0){
			addFirstCluster(tree,i,0,parentcuts,fasta_specs[0]-1);
		}
	}
	for(i=0; i<fasta_specs[0]; i++){
		free(clusterarr[i]);
	}
	free(clusterarr);
	free(parentcuts);
	free(foundPath);
	//shiftColumns(kseqs);
	//numberOfNodesToCut++;
	//if ( tree[0][root].nodeToCut != 1 ){
	//	printf("root is not 1\n");
	//	printdescendants(tree,root,numberOfNodesToCut,0,kseqs);
	//	numberOfNodesToCut++;
	//}
	numberOfNodesToCut=findNumberOfClusters();
	printf("number of nodes to cut: %d\n",numberOfNodesToCut);
	int largestCluster=0;
	int sum_kseq=0;
	for(i=1; i<numberOfNodesToCut; i++){
		clusterSize[i]=countNumInCluster(i,kseqs);
		printf("number in cluster%d: %d\n",i,clusterSize[i]);
		if (clusterSize[i] > largestCluster){
			largestCluster = clusterSize[i];
		}
		sum_kseq = sum_kseq + clusterSize[i];
	}
	printf("largest cluster is %d\n",largestCluster);
	printf("sum of clusters is %d\n",sum_kseq);
	treeArr = (node **)malloc((numberOfNodesToCut-1)*sizeof(node *));
	for(i=0; i<numberOfNodesToCut-1; i++){
		treeArr[i]=malloc((2*clusterSize[i+1]-1)*sizeof(node));
		for(j=0; j<2*clusterSize[i+1]-1; j++){
			treeArr[i][j].name = (char *)malloc(MAXNAME*sizeof(char));
			strcpy(treeArr[i][j].name,"internal");
			treeArr[i][j].likenc = malloc(fasta_specs[1]*sizeof(double *));
			treeArr[i][j].posteriornc = malloc(fasta_specs[1]*sizeof(double *));
			int n;
			for(n=0; n<fasta_specs[1]; n++){
				treeArr[i][j].likenc[n] = malloc(4*sizeof(double));
				treeArr[i][j].posteriornc[n] = malloc(4*sizeof(double));
			}
		}
	}
	for(i=0; i<4; i++){
		PMATnc[0][i][4]=1.0;
		PMATnc[0][i][4]=1.0;
	}
	//allocateMemForTreeArr(numberOfNodesToCut-1,clusterSize,treeArr,kseqs);
	for(i=0; i<clusterSize[numberOfNodesToCut-1]; i++){
		for(j=0; j<fasta_specs[0]; j++){
			if (strcmp(seqNames[j],clusters[numberOfNodesToCut-1][i])==0){
				break;
			}
		}
		strcpy(cluster_seqs[numberOfNodesToCut-1][i],sequences[j]);
	}
	int* rootArr = (int *)malloc(numberOfNodesToCut*sizeof(int));
	createTreesForClusters(treeArr, numberOfNodesToCut,clusterSize,cluster_seqs,clusters,fasta_specs,rootArr);
	char** kalign_args = (char **)malloc(2*sizeof(char *));
	kalign_args[0] = (char*)malloc(100*sizeof(char));
	kalign_args[1] = (char*)malloc(100*sizeof(char));
	strcpy(kalign_args[0],"test_CO1.1000.fasta");
	strcpy(kalign_args[1],"test_dbug/test_kalign_MSA.fasta");
	int*** seqArr = (int ***)malloc(numberOfNodesToCut*sizeof(int **));
	int k,l,m;
	int *numbase = (int *)malloc((numberOfNodesToCut-1)*sizeof(int));
	for(i=0; i<numberOfNodesToCut-1; i++){
		seqArr[i] = (int **)malloc(clusterSize[i+1]*sizeof(int *));
		for(j=0; j<clusterSize[i+1]; j++){
			seqArr[i][j] = (int *)malloc(5000*sizeof(int));
				for(k=0; k<5000; k++){
					seqArr[i][j][k]='\0';
				}
		}
		numbase[i] = 0;
		if (clusterSize[i+1] > 3){
			main_kalign(1,kalign_args,clusterSize[i+1],clusters[i+1],cluster_seqs[i+1],seqArr,numbase,i);
			estimatenucparameters(i,numbase[i],rootArr[i],clusterSize[i+1],seqArr);
			getposterior_nc(i,numbase[i],rootArr[i],clusterSize[i+1],seqArr);
		}else{
			numbase[i] = strlen(cluster_seqs[i+1][0]);
		}
	}
	type_of_PP**** PP = (type_of_PP ****)malloc((numberOfNodesToCut-1)*sizeof(type_of_PP ***));
	for(i=0; i<numberOfNodesToCut-1;i++){
		PP[i] = (type_of_PP ***)malloc((2*clusterSize[i+1]-1)*sizeof(type_of_PP **));
		for(j=0; j<2*clusterSize[i+1]-1; j++){
			PP[i][j] = (type_of_PP **)malloc((numbase[i])*sizeof(type_of_PP *));
			for(k=0; k<numbase[i];k++){
				PP[i][j][k] = (type_of_PP *)malloc(4*(sizeof(type_of_PP)));
			}
		}
	}
	//initialize_assignment_mem(PP,numberOfNodesToCut-1,clusterSize,numbase);
	store_PPs(PP,numberOfNodesToCut-1,clusterSize,numbase);
	rootSeqs = (char **)malloc((numberOfNodesToCut-1)*(sizeof(char *)));
	for(i=0; i<numberOfNodesToCut-1;i++){
		rootSeqs[i]=(char *)malloc((numbase[i]+1)*(sizeof(char)));
	}
	printRootSeqs(rootSeqs,PP,numberOfNodesToCut-1,numbase,rootArr,clusterSize,cluster_seqs);
	int count=0;
	double averageDist;
	//char*** seqsInClusterAvg = (char ***)malloc((opt.number_of_clusters-1)*sizeof(char **));
	//for(i=0; i<opt.number_of_clusters-1; i++){
	//	seqsInClusterAvg[i] = (char **)malloc(opt.number_of_kseqs*sizeof(char *));
	//	for(j=0; j<opt.number_of_kseqs; j++){
	//		seqsInClusterAvg[i][j] = (char *)malloc(fasta_specs[1]*sizeof(char));
	//	}
	//}
	double* calculateAvg = (double *)malloc(2*sizeof(double));
	calculateAvg[0]=0;
	calculateAvg[1]=0;
	//double*** distMatAvg = (double ***)malloc(numberOfNodesToCut*sizeof(double **));
	for(i=1; i<numberOfNodesToCut; i++){
		for(j=0; j<clusterSize[i]; j++){
			for(k=0; k<kseqs; k++){
				if (strcmp(clusters[i][j],clusters[0][k])==0){
					strcpy(seqsInCluster[j],sequences[chooseK[k]]);
					//strcpy(seqsInClusterAvg[i-1][j],sequences[chooseK[k]]);
				}
			}
		}
		for(l=0; l<kseqs+1; l++){
			for(m=0; m<kseqs+1; m++){
				distMat[l][m]=0;
			}
		}
		//if (i>1){
		//	distMatAvg[count] = (double **)malloc(clusterSize[i-1]*sizeof(double *));
		//	for(l=0; l<clusterSize[i-1]; l++){
		//		distMatAvg[count][l] = (double *)malloc(clusterSize[i]*sizeof(double));
		//	}
		//	averageDist = calculateAverageDist(seqsInClusterAvg[i-2],seqsInClusterAvg[i-1],clusterSize[i-1],clusterSize[i],fasta_specs[1],distMatAvg,count,calculateAvg);
		//	count++;
		//}
		//createDistMat(seqsInCluster,distMat,clusterSize[i],fasta_specs);
		//print_distance_matrix(distMat,clusters,i,clusterSize[i]);
		//root = NJ(tree,distMat,clusterSize[i],i);
		//get_number_descendants(tree,root,i);
		printf("cluster %d tree\n",i);
		//printtree(tree,i,clusterSize[i]);
		print_taxonomy(i,clusterSize[i],chooseK,kseqs,seqNames,taxMap,opt.hasTaxFile);
	}
	//freeSeqsInClusterAvg(seqsInClusterAvg,opt);
	double* averageDistances;
	averageDistances = malloc(1000*sizeof(double));
	char** seqsInClusterA = (char **)malloc(largestCluster*sizeof(char *));
	char** seqsInClusterB = (char **)malloc(largestCluster*sizeof(char *));
	for(i=0; i<largestCluster; i++){
		seqsInClusterA[i] = (char *)malloc((fasta_specs[1]+1)*sizeof(char));
		seqsInClusterB[i] = (char *)malloc((fasta_specs[1]+1)*sizeof(char));
	}
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	printf("Calculating pairwise average...\n");
	count=0;
	for(i=1; i<numberOfNodesToCut-1; i++){
		for(j=0; j<clusterSize[i]; j++){
			//strcpy(seqsInClusterA[j],(char *)hashmap_get(&map,clusters[i][j]));
			strcpy(seqsInClusterA[j],cluster_seqs[i][j]);
			//for(k=0; k<opt.number_of_kseqs; k++){
			//	if (strcmp(clusters[i][j],clusters[0][k])==0){
			//		strcpy(seqsInClusterA[j],sequences[chooseK[k]]);
			//	}
			//}
		}
		for(j=i+1; j<numberOfNodesToCut-1; j++){
			for(k=0; k<clusterSize[j]; k++){
				//strcpy(seqsInClusterB[k],(char *)hashmap_get(&map,clusters[j][k]));
				strcpy(seqsInClusterB[k],cluster_seqs[j][k]);
				//for(l=0; l<opt.number_of_kseqs; l++){
				//	if (strcmp(clusters[j][k],clusters[0][l])==0){
				//		strcpy(seqsInClusterB[k],sequences[chooseK[l]]);
				//	}
				//}
			}
			for(l=0; l<kseqs+1; l++){
				for(m=0; m<kseqs+1; m++){
					distMat[l][m]=0;
				}
			}
		//for(j=0; j<clusterSize[i]; j++){
			//for(k=0; k<clusterSize[i+1]; k++){
			if ( opt.use_nw==0 ){
				calculateAverageDist_WFA(seqsInClusterA,seqsInClusterB,clusterSize[i],clusterSize[j],2*fasta_specs[1],distMat,calculateAvg);
			}else{
				calculateAverageDist(seqsInClusterA,seqsInClusterB,clusterSize[i],clusterSize[j],2*fasta_specs[1],distMat,calculateAvg);
			}
			//}
		//}
			averageDistances[count]=calculateAvg[0]/calculateAvg[1];
			count++;
			calculateAvg[0]=0;
			calculateAvg[1]=0;
		}
	}
	for(i=0; i<largestCluster; i++){
		free(seqsInClusterA[i]);
		free(seqsInClusterB[i]);
	}
	free(seqsInClusterA);
	free(seqsInClusterB);
	for(i=0; i<MAXNUMBEROFCLUSTERS; i++){
		free(cluster_seqs[i]);
	}
	free(cluster_seqs);
	for(i=1; i< numberOfNodesToCut; i++){
		double** distMatCluster = (double **)malloc((clusterSize[i]+1)*sizeof(double *));
		for(j=0; j<clusterSize[i]+1; j++){
			distMatCluster[j] = (double *)malloc((clusterSize[i]+1)*sizeof(double));
		}
		for(j=0; j<clusterSize[i]; j++){
			for(k=0; k<clusterSize[i]; k++){
				distMatCluster[j][k]=0;
			}
		}
		char** cluster_sequences = (char **)malloc(clusterSize[i]*sizeof(char *));
		for(j=0; j<clusterSize[i]; j++){
			cluster_sequences[j] = (char *)malloc(fasta_specs[1]*sizeof(char));
		}
		int index=0;
		for(j=0; j<fasta_specs[0]; j++){
			for(k=0; k<clusterSize[i]; k++){
				if ( strcmp(seqNames[j],clusters[i][k])==0 ){
					strcpy(cluster_sequences[index],sequences[j]);
					index++;
				}
			}
		}
		//createDistMat(cluster_sequences,distMatCluster,clusterSize[i],fasta_specs);
		//root = NJ(tree,distMatCluster,clusterSize[i],i,i);
	}
	if ( opt.output_fasta==1 ){
		printInitialClusters(starting_number_of_clusters,numberOfNodesToCut,clusterSize,opt,fasta_specs[0],taxMap,opt.hasTaxFile,seqNames,sequences);
	}
	if ( opt.clstr_format==1 ){
		if (numberOfUnAssigned == fasta_specs[0]){
			clstr = (char ***)malloc(MAXNUMBEROFCLUSTERS*sizeof(char **));
			clstr_lengths = (int **)malloc(MAXNUMBEROFCLUSTERS*sizeof(int *));
			for(i=0; i<MAXNUMBEROFCLUSTERS; i++){
				clstr[i] = (char **)malloc(fasta_specs[0]*sizeof(char *));
				clstr_lengths[i] = (int *)malloc(fasta_specs[0]*sizeof(int));
				for(j=0; j<fasta_specs[0]; j++){
					clstr_lengths[i][j] = -1;
					clstr[i][j] = (char *)malloc(MAXNAME*sizeof(char));
					memset(clstr[i][j],'\0',MAXNAME);
				}
			}
		}
		printInitialClusters_CLSTR(starting_number_of_clusters,numberOfNodesToCut,clusterSize,opt,fasta_specs[0],fasta_specs[1],clstr,clstr_lengths,seqNames,sequences);
	}
	double averageAvg = 0;
	printf("count is %d\n",count);
	//double shortestDistBetweenClusters = averageDistances[0];
	for(i=0; i<count; i++){
		averageAvg = averageDistances[i] + averageAvg;
	}
	averageAvg = averageAvg/count;
	freeDistMat(distMat,kseqs+1);
	free(averageDistances);
	//freeDistMatAvg(distMatAvg,clusterSize,count);
	freeTreeArr(tree,clusterSize);
	double average = calculateAvg[0]/calculateAvg[1];
	free(calculateAvg);
	printf("average pairwise: %lf\n",averageAvg);
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %lf seconds\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	pthread_t threads[opt.numthreads];
	mystruct mstr[opt.numthreads];
	for(i=0; i<opt.numthreads; i++){
		mstr[i].str = malloc(sizeof(struct resultsStruct));
		mstr[i].average = averageAvg;
		mstr[i].clusterSize = clusterSize;
		mstr[i].number_of_kseqs = kseqs;
		mstr[i].numAssigned = 0;
		mstr[i].num_threads = opt.numthreads;
		mstr[i].number_of_clusters = numberOfNodesToCut;
		//mstr[i].clusterNames = (char ***)malloc(numberOfNodesToCut*sizeof(char **));
		mstr[i].fasta_specs = fasta_specs;
		mstr[i].nw_struct = initialize_nw(fasta_specs[3]);
		//mstr[i].seqNames = (char **)malloc(opt.numberOfLinesToRead*sizeof(char *));
		//mstr[i].sequences = (char **)malloc(opt.numberOfLinesToRead*sizeof(char *));
		//mstr[i].taxonomy = (char **)malloc(opt.numberOfLinesToRead*sizeof(char *));
		mstr[i].chooseK = (int *)malloc(kseqs*sizeof(int));
		for(j=0; j<kseqs; j++){
			mstr[i].chooseK[j] = chooseK[j];
		}
		//for(j=0; j<opt.numberOfLinesToRead; j++){
		//	mstr[i].seqNames[j] = (char *)malloc((fasta_specs[2]+1)*sizeof(char));
		//	mstr[i].sequences[j] = (char *)malloc((fasta_specs[1]+1)*sizeof(char));
		//	mstr[i].taxonomy[j] = (char *)malloc(FASTA_MAXLINE*sizeof(char));
		//}
		//for(j=0; j<numberOfNodesToCut; j++){
			//mstr[i].clusterNames[j]=(char **)malloc(clusterSize[j]*sizeof(char *));
			//for(k=0; k<clusterSize[j]; k++){
				//mstr[i].clusterNames[j][k]=(char *)malloc(MAXNAME*sizeof(char));
				//mstr[i].clusterSeqs[j][k]=(char *)malloc((fasta_specs[1]+1)*sizeof(char));
				//strcpy(mstr[i].clusterNames[j][k],clusters[j][k]);
				//strcpy(mstr[i].clusterSeqs[j][k],(char *)hashmap_get(&map,clusters[j][k]));
			//}
		//}
		mstr[i].clusterSeqs = (char ***)malloc((numberOfNodesToCut+PADDING)*sizeof(char **));
		for(j=0; j<numberOfNodesToCut+PADDING; j++){
			if (j<numberOfNodesToCut){
				if (j>0){
					mstr[i].clusterSeqs[j]=(char **)malloc(clusterSize[j]*sizeof(char *));
					for(k=0; k<clusterSize[j]; k++){
						mstr[i].clusterSeqs[j][k]=(char *)malloc((fasta_specs[1]+1)*sizeof(char));
						for(l=0; l<fasta_specs[0]; l++){
							if (strcmp(clusters[j][k],seqNames[l])==0){
								break;
							}
						}
						strcpy(mstr[i].clusterSeqs[j][k],sequences[l]);
						//strcpy(mstr[i].clusterSeqs[j][k],(char *)hashmap_get(&map,clusters[j][k]));
					}
				}else{
					mstr[i].clusterSeqs[j]=(char **)malloc(MAXNUMBEROFKSEQS*sizeof(char *));
					for(k=0; k<MAXNUMBEROFKSEQS; k++){
						mstr[i].clusterSeqs[j][k]=(char *)malloc((fasta_specs[1]+1)*sizeof(char));
						if (k<clusterSize[j]){
							for(l=0; l<fasta_specs[0]; l++){
								if (strcmp(seqNames[l],clusters[j][k])==0){
									break;
								}
							}
							strcpy(mstr[i].clusterSeqs[j][k],sequences[l]);
						}
					}
				}
			}else{
				mstr[i].clusterSeqs[j]=(char **)malloc(1*sizeof(char *));
				for(k=0; k<1; k++){
					mstr[i].clusterSeqs[j][k]=(char *)malloc((fasta_specs[1]+1)*sizeof(char));
				}
			}
		}
		//allocateMemForResults(mstr[i].str, opt.numberOfLinesToRead, opt.numthreads, numberOfNodesToCut, averageAvg);
		//mstr[i].str->accession = (char **)malloc((opt.numberOfLinesToRead)*sizeof(char *));
		//mstr[i].str->savedForNewClusters = (char **)malloc((opt.numberOfLinesToRead)*sizeof(char *));
		//mstr[i].str->clusterNumber = (int *)malloc(opt.numberOfLinesToRead*sizeof(int));
		//mstr[i].str->assigned = (char **)malloc((opt.numberOfLinesToRead*sizeof(char *)));
		mstr[i].str->numassigned=0;
		//for(j=0; j<opt.numberOfLinesToRead; j++){
			//mstr[i].str->accession[j] = malloc((fasta_specs[2]+1)*sizeof(char));
			//mstr[i].str->assigned[j] = malloc((fasta_specs[2]+1)*sizeof(char));
			//mstr[i].str->savedForNewClusters[j] = malloc((fasta_specs[2]+1)*sizeof(char));
			//mstr[i].str->clusterNumber[j]=-1;
			//memset(mstr[i].str->assigned[j],'\0',fasta_specs[2]+1);
			//memset(mstr[i].str->accession[j],'\0',fasta_specs[2]+1);
			//memset(mstr[i].str->savedForNewClusters[j],'\0',fasta_specs[2]+1);
		//}
		mstr[i].str->average = average;
		mstr[i].str->number_of_clusters = numberOfNodesToCut;
		mstr[i].str->clusterSizes = (int *)malloc(MAXNUMBEROFCLUSTERS*sizeof(int));
		for (j=0; j<numberOfNodesToCut; j++){
			mstr[i].str->clusterSizes[j] = clusterSize[j];
		}
		//mstr[i].seqsToCompare = seqsToCompare;
		//mstr[i].assignedSeqs = (char **)malloc(fasta_specs[0]*sizeof(char *));
		//for( j=0; j<fasta_specs[0]; j++){
		//	mstr[i].assignedSeqs[j] = (char *)malloc(MAXNAME*sizeof(char));
		//	mstr[i].assignedSeqs[j][0] = '\0';
		//}
		//int place;
		//for(j=fasta_specs[0]-1; j>=0; j--){
		//	if (assigned[j][0]=='\0'){
		//		place=j;
		//	}
		//}
		//for( j=0; j<place; j++){
		//	strcpy(mstr[i].assignedSeqs[j],assigned[j]);
		//}
		mstr[i].chooseK = (int *)malloc(kseqs*sizeof(int));
		for(j=0; j<kseqs; j++){
			mstr[i].chooseK[j] = chooseK[j];
		}
	}
	int closestCluster=0;
	int next=0;
	int update_initial_cluster = kseqs;
	int largest_cluster = findLargestCluster(clusterSize,numberOfNodesToCut);
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	printf("Starting to assign sequences to clusters\n");
	int start=0;
	int end=0;
	int returnLineNumber=0;
	int numberToAssign = opt.numberOfLinesToRead;
	int endOfAssign = fasta_specs[0];
	if (numberOfUnAssigned < fasta_specs[0]){
		numberToAssign = numberOfUnAssigned/opt.numthreads;
		endOfAssign = numberOfUnAssigned;
	}
	while (1){
		j=0;
		for(i=0; i<opt.numthreads; i++){
			start = j;
			end = j+numberToAssign;
			if (i==opt.numthreads-1)
				end=endOfAssign;
			mstr[i].start=start;
			mstr[i].end=end;
			mstr[i].largest_cluster=largest_cluster;
			mstr[i].threadnumber=i;
			mstr[i].use_nw = opt.use_nw;
			mstr[i].seqNames = (char **)malloc((end-start)*sizeof(char *));
			mstr[i].sequences = (char **)malloc((end-start)*sizeof(char *));
			mstr[i].taxonomy = (char **)malloc((end-start)*sizeof(char *));
			mstr[i].str->accession = (char **)malloc((end-start)*sizeof(char *));
			mstr[i].str->savedForNewClusters = (char **)malloc((end-start)*sizeof(char *));
			mstr[i].str->clusterNumber = (int *)malloc((end-start)*sizeof(char *));
			for(k=0; k<end-start; k++){
				mstr[i].seqNames[k] = (char *)malloc((fasta_specs[2]+1)*sizeof(char));
				mstr[i].sequences[k] = (char *)malloc((fasta_specs[1]+1)*sizeof(char));
				mstr[i].taxonomy[k] = (char *)malloc(FASTA_MAXLINE*sizeof(char));
				mstr[i].str->accession[k] = (char *)malloc((fasta_specs[2]+1)*sizeof(char));
				mstr[i].str->savedForNewClusters[k] = (char *)malloc((fasta_specs[2]+1)*sizeof(char));
			       	mstr[i].str->clusterNumber[k] = -1;
				memset(mstr[i].str->accession[k],'\0',fasta_specs[2]+1);
				memset(mstr[i].str->savedForNewClusters[k],'\0',fasta_specs[2]+1);
				strcpy(mstr[i].seqNames[k],seqNames[j+k]);
				strcpy(mstr[i].sequences[k],sequences[j+k]);
				strcpy(mstr[i].taxonomy[k],taxonomy[j+k]);
			}
			j=j+numberToAssign;
		}
		for(i=0; i<opt.numthreads; i++){
			size_t default_stack_size;
			pthread_attr_t stack_size_custom_attr;
			pthread_attr_init(&stack_size_custom_attr);
			pthread_attr_getstacksize(&stack_size_custom_attr,&default_stack_size);
			if (default_stack_size < MIN_REQ_SSIZE){
				pthread_attr_setstacksize(&stack_size_custom_attr,(size_t)MIN_REQ_SSIZE);
			}
			//pthread_create(&threads[i], NULL, runAssignToCluster, &mstr[i]);
			pthread_create(&threads[i], &stack_size_custom_attr, runAssignToCluster, &mstr[i]);
		}
		for(i=0; i<opt.numthreads; i++){
			pthread_join(threads[i], NULL);
		}
		//for(i=0; i<opt.numthreads; i++){
			//for( j=0; j<(mstr[i].end-mstr[i].start); j++){
				//fprintf(results);
			//}
		//}
		break;
	}
	k=0;
	char** actualSeqsToClusterLater = (char **)malloc(fasta_specs[0]*sizeof(char *));
	for(i=0; i<fasta_specs[0]; i++){
		actualSeqsToClusterLater[i] = (char *)malloc((fasta_specs[1]+1)*sizeof(char));
		memset(actualSeqsToClusterLater[i],'\0',fasta_specs[1]+1);
		memset(sequencesToClusterLater[i],'\0',MAXNAME);
	}
	for(i=0; i<opt.numthreads; i++){
		for(j=0; j<numberToAssign; j++){
			if (mstr[i].str->savedForNewClusters[j][0]!='\0'){
				int dont_add=0;
				for(m=0; m<fasta_specs[0]; m++){
					if (strcmp(sequencesToClusterLater[m],mstr[i].str->savedForNewClusters[j])==0){
						dont_add=1;
					}
				}
				if (dont_add==0){
					strcpy(sequencesToClusterLater[k],mstr[i].str->savedForNewClusters[j]);
					for(l=0; l<fasta_specs[0]; l++){
						if ( strcmp(seqNames[l],mstr[i].str->savedForNewClusters[j])==0 ){
							break;
						}
					}
					strcpy(actualSeqsToClusterLater[k],sequences[l]);
				//strcpy(actualSeqsToClusterLater[k],(char*)hashmap_get(&map,mstr[i].str->savedForNewClusters[j]));
				//hashmap_put(&map,sequencesToClusterLater[k],(char*)hashmap_get(&map,mstr[i].str->savedForNewClusters[j]));
					k++;
				}
			}
		}
	}
	//for(i=0; i<opt.numthreads; i++){
	//	for(j=0; j<mstr[i].str->numassigned; j++){
	//		int place;
	//		for(k=fasta_specs[0]-1; k>=0; k--){
	//			if (assigned[k]=='\0'){
	//				place=k;
	//			}
	//		}
	//		strcpy(assigned[place],mstr[i].str->assigned[j]);
			//hashmap_put(&assignedSeqs,mstr[i].str->assigned[j],1);
	//	}
	//}
	for(i=0; i<opt.numthreads; i++){
		if (opt.output_fasta==1){
			printClusters(starting_number_of_clusters,mstr[i].str->number_of_clusters,fasta_specs[0],opt,taxonomy,mstr[i],mstr[i].numAssigned,taxMap,opt.hasTaxFile,seqNames,sequences);
		}
		if(opt.clstr_format==1){
			saveCLSTR(starting_number_of_clusters,fasta_specs[0],mstr[i],mstr[i].numAssigned,clstr,clstr_lengths,fasta_specs[1],seqNames,sequences);
		}
	}
	starting_number_of_clusters = starting_number_of_clusters + mstr[0].str->number_of_clusters - 1;
	for(i=0; i<opt.numthreads; i++){
		start=l;
		end=l+numberToAssign;
		if(i==opt.numthreads-1)
			end=endOfAssign;
		free(mstr[i].chooseK);
		for(j=0; j<end-start; j++){
			free(mstr[i].seqNames[j]);
			free(mstr[i].sequences[j]);
			free(mstr[i].taxonomy[j]);
			free(mstr[i].str->accession[j]);
			//free(mstr[i].str->assigned[j]);
			free(mstr[i].str->savedForNewClusters[j]);
		}
		free(mstr[i].seqNames);
		free(mstr[i].sequences);
		free(mstr[i].taxonomy);
		free(mstr[i].nw_struct);
		free(mstr[i].str->accession);
		//free(mstr[i].str->assigned);
		free(mstr[i].str->savedForNewClusters);
		free(mstr[i].str->clusterNumber);
		free(mstr[i].str->clusterSizes);
		//for(j=0; j<numberOfNodesToCut; j++){
			//for(k=0; k<clusterSize[j]; k++){
			//	free(mstr[i].clusterNames[j][k]);
			//}
			//free(mstr[i].clusterNames[j]);
		//}
		//free(mstr[i].clusterNames);
		//for(j=0; j<fasta_specs[0]; j++){
		//	free(mstr[i].assignedSeqs[j]);
		//}
		//free(mstr[i].assignedSeqs);
		for(j=0; j<numberOfNodesToCut+PADDING; j++){
			if (j<numberOfNodesToCut){
				if (j>0){
					for(k=0; k<clusterSize[j]; k++){
						free(mstr[i].clusterSeqs[j][k]);
					}
				}else{
					for(k=0; k<MAXNUMBEROFKSEQS; k++){
						free(mstr[i].clusterSeqs[j][k]);
					}
				}
			}else{
				free(mstr[i].clusterSeqs[j][0]);
			}
			free(mstr[i].clusterSeqs[j]);
		}
		free(mstr[i].clusterSeqs);
		l=l+numberToAssign;
	}
	for(i=0; i<fasta_specs[0]; i++){
		memset(seqNames[i],'\0',fasta_specs[2]);
		memset(sequences[i],'\0',fasta_specs[2]);
	}
	numberOfUnAssigned=countNumUnassigned(sequencesToClusterLater,fasta_specs);
	numberOfSequencesLeft=numberOfUnAssigned;
	for(i=0; i<numberOfUnAssigned; i++){
		strcpy(seqNames[i],sequencesToClusterLater[i]);
		strcpy(sequences[i],actualSeqsToClusterLater[i]);
	}
	for(i=0; i<fasta_specs[0]; i++){
		free(actualSeqsToClusterLater[i]);
	}
	free(actualSeqsToClusterLater);
	free(clusterSize);
	clock_gettime(CLOCK_MONOTONIC, &tend);
	int* newClusterSizes = (int *)malloc(MAXNUMBEROFCLUSTERS*sizeof(int));
	numberOfNodesToCut=findNewClusterSizes(newClusterSizes,fasta_specs);
	printf("Finished! Took %lf seconds\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	free(newClusterSizes);
	//hashmap_destroy(&seqsToCompare);
	for(i=0; i<MAXNUMBEROFCLUSTERS; i++){
		for(j=0;j<fasta_specs[0];j++){
			memset(clusters[i][j],'\0',MAXNAME);
		}
	}
	if (numberOfUnAssigned<=3){
		if (opt.output_fasta==1){
			printLessThanFour(numberOfUnAssigned, opt, taxMap, starting_number_of_clusters,seqNames, sequences);
		}
		if (opt.clstr_format==1){
			printLessThanFour_CLSTR(numberOfUnAssigned, clstr, clstr_lengths,fasta_specs[1],seqNames, sequences);
		}
	}
	}
	printCLSTR(opt,clstr,clstr_lengths,fasta_specs[0],starting_number_of_clusters);
	//freeMemForAlign(DATA,fasta_specs[1],mult);
	//freeMemForDistMat(clusterSize,largest_cluster,distMat2);
	freeClusters(kseqs);
	freeSequences(fasta_specs[0],sequences,seqNames,taxonomy, sequencesToClusterLater);
	free(fasta_specs);
	free(chooseK);
	//hashmap_destroy(&map);
}
