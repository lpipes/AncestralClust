#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "needleman_wunsch.h"

#define FASTA_MAXLINE 5000
#define MAXNAME 30
#define KSEQ 100
#define DISTMAX 30.0
#define MINBL 0.00001
#define MAXNAME 30

typedef struct node{
	int down;
	int up[2];
	double frac[4];
	double bl;
	int nd;
	char* name;
	int nodeToCut;
}node;

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
			strcpy(sequences[i],buffer);
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
	if (numrealsites < 1) {printf("Wacky distance in function 'Get_dist_JC' between sequence %i and %i\n,i,j,"); exit(-1);}
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
	//for (i=0; i<KSEQ+1; i++){
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
void sortArray(double* branchLengths, int* indexArray){
	int i,j, tmp2;
	double tmp;
	//Bubble Sort
	for(i=0; i<2*KSEQ-1; i++){
		indexArray[i]=i;
	}
	for(i=0; i<2*KSEQ; i++){
		for(j=i+1;j<2*KSEQ-1;j++){
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
	/*for(i=0; i<2*KSEQ-1; i++){
		printf("branchLengths[%d]: %lf\n",i,branchLengths[i]);
	}*/
}
int get_number_descendants(node** tree, int node, int whichTree){
	if (tree[whichTree][node].up[0]==-1)  return (tree[whichTree][node].nd=1);
	else return (tree[whichTree][node].nd=(get_number_descendants(tree,tree[whichTree][node].up[0],whichTree)+get_number_descendants(tree,tree[whichTree][node].up[1],whichTree)));
}
void printdescendants(node** tree, int node, char*** clusters, int clusterNumber, int whichTree){
	int child1 = tree[whichTree][node].up[0];
	int child2 = tree[whichTree][node].up[1];
	int i;
	if (tree[whichTree][node].nodeToCut==1){
		return;
	}
	if (tree[whichTree][node].up[0]==-1 && tree[whichTree][node].up[1]==-1){
		int count=0;
		for(i=KSEQ-1; i>=0; i--){
			if (clusters[clusterNumber][i][0]=='\0'){
				count=i;
			}
		}
		printf("count %d node %d: %s\n",count,node,tree[whichTree][node].name);
		strcpy(clusters[clusterNumber][count],tree[whichTree][node].name);
	}else{
		printdescendants(tree,child1,clusters,clusterNumber,whichTree);
		printdescendants(tree,child2,clusters,clusterNumber,whichTree);
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
int countNumInCluster(char*** clusters, int index){
	int i,j;
	for(i=0; i<KSEQ; i++){
		if(clusters[index][i][0]!='\0'){
			j=i;
		}
	}
	return j+1;
}
void allocateMemForTreeArr(int numberOfClusters, int* clusterSize, node** treeArr){
	int i,j;
	treeArr=(node **)malloc(numberOfClusters*sizeof(node *));
	for(i=0; i<numberOfClusters; i++){
		treeArr[i]=malloc((2*KSEQ-1)*sizeof(node));
		for(j=0; j<2*KSEQ-1; j++){
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
}
void updateNumberOfDescendants(node** tree, int node, int descendants, int whichTree, int returnNode){
	int parent = tree[whichTree][node].down;
	if (node==returnNode){
		return;
	}
	if ( tree[whichTree][node].down == -1 ){
		return;
	}
	if (node != returnNode){
		tree[whichTree][node].nd=tree[whichTree][node].nd-descendants;
	}
	updateNumberOfDescendants(tree,parent,descendants,whichTree,returnNode);
}

int main(){
	FILE* fasta_for_clustering;
	if (( fasta_for_clustering = fopen("/space/s1/lenore/crux_db/crux_db2/new_blast/GazF1_GazR1/GazF1_GazR1_db_filtered/GazF1_GazR1_fasta_and_taxonomy/GazF1_GazR1_.rmAmbig.fasta","r")) == (FILE *) NULL ) fprintf(stderr,"FASTA file could not be opened.\n");
	int number_of_sequences = 0;
	int* fasta_specs = (int *)malloc(4*sizeof(int));
	setNumSeq(fasta_for_clustering,fasta_specs);
	fasta_specs[3] = 5; //NUMBER OF CLUSTERS
	printf("Number of sequences: %d\n",fasta_specs[0]);
	printf("Longest sequence: %d\n",fasta_specs[1]);
	printf("Longest name: %d\n",fasta_specs[2]);
	printf("Number of clusters: %d\n",fasta_specs[3]);
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
	if (( fasta_for_clustering = fopen("/space/s1/lenore/crux_db/crux_db2/new_blast/GazF1_GazR1/GazF1_GazR1_db_filtered/GazF1_GazR1_fasta_and_taxonomy/GazF1_GazR1_.rmAmbig.fasta","r")) == (FILE *) NULL ) fprintf(stderr,"FASTA file could not be opened.\n");
	readInFasta(fasta_for_clustering,seqNames,sequences);
	close(fasta_for_clustering);
	int* chooseK = (int *)malloc(KSEQ*sizeof(int));
	for(i=0; i<KSEQ; i++){
		chooseK[i]=-1;
	}
	int choosing=0;
		srand(time(0));
	while(choosing==0){
		int random_number = generateRandom(fasta_specs[0]);
		int index=0;
		for(i=KSEQ; i>=0; i--){
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
		}
		if (chooseK[KSEQ-1]!=-1){
			choosing=1;
		}
	}
	double** distMat = (double **)malloc((KSEQ+1)*sizeof(double *));
	for(i=0; i<KSEQ+1; i++){
		distMat[i] = (double *)malloc((KSEQ+1)*sizeof(double));
	}
	for(i=0; i<KSEQ+1; i++){
		for(j=0; j<KSEQ+1; j++){
			distMat[i][j] = 0;
		}
	}
	char** seqsInCluster = (char **)malloc(KSEQ*sizeof(char *));
	for(i=0; i<KSEQ; i++){
		seqsInCluster[i] = (char *)malloc(fasta_specs[1]*sizeof(char));
		strcpy(seqsInCluster[i],sequences[chooseK[i]]);
		//printf("seqsInCluster[%d]: %s\n",i,seqsInCluster[i]);
	}
	createDistMat(seqsInCluster,distMat,KSEQ,fasta_specs[1]);
	char*** clusters = (char***)malloc((fasta_specs[3]+1)*sizeof(char **));
	for(i=0; i<fasta_specs[3]+1; i++){
		clusters[i]=(char **)malloc((KSEQ)*sizeof(char *));
		for(j=0; j<KSEQ; j++){
			clusters[i][j]=(char *)malloc(MAXNAME*sizeof(char));
			clusters[i][j][0]='\0';
		}
	}
	for(i=0; i<KSEQ; i++){
		strcpy(clusters[0][i],seqNames[chooseK[i]]);
	}
	int* clusterSize = (int *)malloc(fasta_specs[3]*sizeof(int));
	clusterSize[0]=KSEQ;
	print_distance_matrix(distMat,clusters,0,KSEQ);
	node **tree = malloc((fasta_specs[3]+1)*sizeof(node *));
	for(i=0; i<fasta_specs[3]+1; i++){
		tree[i]=(node *)malloc((2*KSEQ-1)*sizeof(node));
		for(j=0; j<2*KSEQ-1; j++){
			tree[i][j].name=(char *)malloc(MAXNAME*sizeof(char));
			strcpy(tree[i][j].name,"internal");
			tree[i][j].nodeToCut=0;
		}
	}
	int root = NJ(tree,distMat,clusters,KSEQ,0);
	get_number_descendants(tree,root,0);
	printtree(tree,0,KSEQ);
	int* longestBranches = (int *)malloc(fasta_specs[3]*sizeof(int));
	double* branchLengths = (double *)malloc((2*KSEQ-1)*sizeof(double));
	for(i=0; i<2*KSEQ-1; i++){
		branchLengths[i] = tree[0][i].bl;
	}
	int* indexArray = (int *)malloc((2*KSEQ-1)*sizeof(int));
	sortArray(branchLengths,indexArray);
	printf("Longest Branch: %lf node: %d\n",branchLengths[0],indexArray[0]);
	int numberOfNodesToCut=0;
	for(i=0; i<2*KSEQ-1; i++){
		if (tree[0][indexArray[i]].nd > 4 && numberOfNodesToCut < 3){
			printf("cutting at node %d\n",indexArray[i]);
			numberOfNodesToCut++;
			printdescendants(tree,indexArray[i],clusters,numberOfNodesToCut,0);
			tree[0][indexArray[i]].nodeToCut=1;
			clearDescendants(tree,indexArray[i],0);
			updateNumberOfDescendants(tree,indexArray[i],tree[0][indexArray[i]].nd,0,tree[0][indexArray[i]].up[0]);
			printf("updating tree\n");
			printtree(tree,0,clusterSize[0]);
		}
	}
	printf("number of nodes to cut: %d\n",numberOfNodesToCut);
	printdescendants(tree,root,clusters,numberOfNodesToCut+1,0);
	for(i=0; i<numberOfNodesToCut+2; i++){
		clusterSize[i]=countNumInCluster(clusters,i);
		printf("number in cluster%d: %d\n",i,clusterSize[i]);
	}
	//node** treeArr;
	//allocateMemForTreeArr(numBranchesToFind,clusterSize,treeArr);
	int k,l,m;
	for(i=1; i<numberOfNodesToCut+2; i++){
		for(j=0; j<clusterSize[i]; j++){
			for(k=0; k<KSEQ; k++){
				if (strcmp(clusters[i][j],clusters[0][k])==0){
					strcpy(seqsInCluster[j],sequences[chooseK[k]]);
				}
			}
		}
		for(l=0; l<KSEQ+1; l++){
			for(m=0; m<KSEQ+1; m++){
				distMat[l][m]=0;
			}
		}
		createDistMat(seqsInCluster,distMat,clusterSize[i],fasta_specs[1]);
		print_distance_matrix(distMat,clusters,i,clusterSize[i]);
		root = NJ(tree,distMat,clusters,clusterSize[i],i);
		get_number_descendants(tree,root,i);
		printf("cluster %d tree\n",i);
		printtree(tree,i,clusterSize[i]);
	}
	//createDistMat(seqsInCluster,distMat,clusterSize[1],fasta_specs[1]);
	//print_distance_matrix(distMat,clusters,1,clusterSize[1]);
	//root = NJ(tree,distMat,clusters,clusterSize[1],1);
	//get_number_descendants(tree,root,1);
	//printf("cluster 1 tree\n");
	//printtree(tree,1,clusterSize[1]);
	//root = NJ(tree,distMat,clusters,clusterSize[2],2);
	//get_number_descendants(tree,root,2);
	//printf("cluster 2 tree\n");
	//printtree(tree,2,clusterSize[2]);
}
