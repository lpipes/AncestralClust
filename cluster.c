#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "needleman_wunsch.h"

#define FASTA_MAXLINE 5000
#define MAXNAME 30
#define KSEQ 10
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

void print_distance_matrix(double** distMat, char** seqNames, int* chooseK){
	int i,j;
	for(i=0; i<KSEQ; i++){
		printf("row %i %s\t",i,seqNames[chooseK[i]]);
		for(j=0; j<KSEQ; j++){
			printf("%lf\t",distMat[i][j]);
		}
		printf("\n");
	}
}
int populate_DATA(char* query1, char* query2, int** DATA, int alignment_length, int* mult){
	
	int i;
	printf("query1: %s\n",query1);
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
	printf("query2: %s\n",query2);
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
	for(i=0; i<alignment_length; i++){
		printf("%d",DATA[0][i]);
	}
	printf("\n");
	for(i=0; i<alignment_length; i++){
		printf("%d",DATA[1][i]);
	}
	printf("\n");
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
int NJ(node* tree, double** distMat,char** seqNames,int* chooseK){
	double *r, minval, D, u1;
	int child1, child2, i, j, n, min, pair[2], newnode, *nodenum;
	r = malloc(KSEQ*(sizeof(double)));
	nodenum = malloc(KSEQ*(sizeof(int)));
	n=KSEQ;
	for(i=0; i<KSEQ; i++){
		nodenum[i]=i;
		tree[i].up[0]=tree[i].up[1]=-1;
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
		newnode=2*KSEQ-n;
		child1 = nodenum[pair[0]];
		child2 = nodenum[pair[1]];
		tree[newnode].up[0]=child1;
		tree[newnode].up[1]=child2;
		tree[child1].down=newnode;
		tree[child2].down=newnode;
		if (tree[child1].up[0]==-1){
			strcpy(tree[child1].name,seqNames[chooseK[pair[0]]]);
		}
		if (tree[child2].up[1]==-1){
			strcpy(tree[child2].name,seqNames[chooseK[pair[1]]]);
		}
		// HERE WE DISALLOW NEGATIVE BRANCHLENGTHS!  IS THIS THE BEST THING TO DO?
		if ((u1 = (distMat[pair[0]][pair[1]] +r[pair[0]]-r[pair[1]])/2.0) < MINBL){
			tree[child1].bl = MINBL;
		}else{
			tree[child1].bl = u1;
		}
		if ((tree[child2].bl = distMat[pair[0]][pair[1]]-u1) < MINBL)/*(DM[pair[0]][pair[1]] +r[pair[1]]-r[pair[0]])/2.0;*/{
			tree[child2].bl = MINBL;
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
		updateChooseK(chooseK, pair[1]);
		updateChooseK(chooseK, pair[0]);
		n--;
	} while (n>2);
	newnode=2*KSEQ-2;
	child1 = nodenum[0];
	child2 = nodenum[1];
	tree[newnode].up[0]=child1;
	tree[newnode].up[1]=child2;
	tree[newnode].down=-1;
	tree[newnode].bl=-1.0;
	tree[child1].down=newnode;
	tree[child2].down=newnode;
	if (tree[child1].up[0]==-1 && strcmp(tree[child1].name,"internal")==0){
		strcpy(tree[child1].name,seqNames[chooseK[0]]);
	}
	if (tree[child2].up[0]==-1 && strcmp(tree[child2].name,"internal")==0){
		strcpy(tree[child1].name,seqNames[chooseK[0]]);
	}	
	if (distMat[0][1] > MINBL*2.0)  // HERE WE DISALLOW NEGATIVE BRANCHLENGTHS!  IS THIS THE BEST THING TO DO?
		tree[child1].bl = tree[child2].bl = distMat[0][1]/2.0;
	else tree[child1].bl = tree[child2].bl = MINBL;
	//for (i=0; i<KSEQ+1; i++){
	//	free(distMat[i]);
	//}
	//free(distMat);
	free(nodenum);
	free(r);
	return(newnode);
}
void updateChooseK(int* chooseK, int index){
	int i;
	int tmp[KSEQ];
	for(i=0; i<KSEQ; i++){
		tmp[i] = chooseK[i];
	}
	for(i=0; i<KSEQ-1; i++){
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
void printtree(node* tree){
	int i;
	for(i=0; i<2*KSEQ-1; i++){
		printf("%i: up: %i %i, down: %i, bl: %f, nd: %d, name: %s\n",i,tree[i].up[0],tree[i].up[1],tree[i].down,tree[i].bl/*totsites*/,tree[i].nd,tree[i].name);
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
	for(i=0; i<2*KSEQ-1; i++){
		printf("branchLengths[%d]: %lf\n",i,branchLengths[i]);
	}
}
int get_number_descendants(node* tree, int node){
	if (tree[node].up[0]==-1)  return (tree[node].nd=1);
	else return (tree[node].nd=(get_number_descendants(tree,tree[node].up[0])+get_number_descendants(tree,tree[node].up[1])));
}
void printdescendants(node* tree, int node, char*** clusters, int clusterNumber){
	int child1 = tree[node].up[0];
	int child2 = tree[node].up[1];
	int i;
	if (tree[node].nodeToCut==1){
		return;
	}
	if (tree[node].up[0]==-1 && tree[node].up[1]==-1){
		int count=0;
		for(i=KSEQ-1; i>=0; i--){
			if (clusters[clusterNumber][i][0]=='\0'){
				count=i;
			}
		}
		printf("count %d node %d: %s\n",count,node,tree[node].name);
		strcpy(clusters[clusterNumber][count],tree[node].name);
	}else{
		printdescendants(tree,child1,clusters,clusterNumber);
		printdescendants(tree,child2,clusters,clusterNumber);
	}
}
int countNumInCluster(char*** clusters, int index){
	int i,j;
	for(i=0; i<KSEQ; i++){
		if(clusters[index][i][0]!='\0'){
			j=i;
		}
	}
	return j;
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
		printf("on %d\n",i);
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
int main(){
	FILE* fasta_for_clustering;
	if (( fasta_for_clustering = fopen("/space/s1/lenore/crux_db/crux_db2/new_blast/GazF1_GazR1/GazF1_GazR1_db_filtered/GazF1_GazR1_fasta_and_taxonomy/GazF1_GazR1_.fasta","r")) == (FILE *) NULL ) fprintf(stderr,"FASTA file could not be opened.\n");
	int number_of_sequences = 0;
	int* fasta_specs = (int *)malloc(3*sizeof(int));
	setNumSeq(fasta_for_clustering,fasta_specs);
	printf("Number of sequences: %d\n",fasta_specs[0]);
	printf("Longest sequence: %d\n",fasta_specs[1]);
	printf("Longest name: %d\n",fasta_specs[2]);
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
	if (( fasta_for_clustering = fopen("/space/s1/lenore/crux_db/crux_db2/new_blast/GazF1_GazR1/GazF1_GazR1_db_filtered/GazF1_GazR1_fasta_and_taxonomy/GazF1_GazR1_.fasta","r")) == (FILE *) NULL ) fprintf(stderr,"FASTA file could not be opened.\n");
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
		printf("seqsInCluster[%d]: %s\n",i,seqsInCluster[i]);
	}
	createDistMat(seqsInCluster,distMat,KSEQ,fasta_specs[1]);
	print_distance_matrix(distMat,seqNames,chooseK);
	struct node *tree = malloc((2*KSEQ-1)*sizeof(node));
	for(i=0; i<2*KSEQ-1; i++){
		tree[i].name=(char *)malloc(MAXNAME*sizeof(char));
		strcpy(tree[i].name,"internal");
		tree[i].nodeToCut=0;
	}
	int root = NJ(tree,distMat,seqNames,chooseK);
	get_number_descendants(tree,root);
	printtree(tree);
	int numBranchesToFind = 2;
	int* longestBranches = (int *)malloc(numBranchesToFind*sizeof(int));
	double* branchLengths = (double *)malloc((2*KSEQ-1)*sizeof(double));
	for(i=0; i<2*KSEQ-1; i++){
		branchLengths[i] = tree[i].bl;
	}
	int* indexArray = (int *)malloc((2*KSEQ-1)*sizeof(int));
	sortArray(branchLengths,indexArray);
	printf("Longest Branch: %lf node: %d\n",branchLengths[0],indexArray[0]);
	char*** clusters = (char***)malloc(numBranchesToFind*sizeof(char **));
	for(i=0; i<numBranchesToFind; i++){
		clusters[i]=(char **)malloc((KSEQ)*sizeof(char *));
		for(j=0; j<KSEQ; j++){
			clusters[i][j]=(char *)malloc(MAXNAME*sizeof(char));
			clusters[i][j][0]='\0';
		}
	}
	int* clusterSize = (int *)malloc(numBranchesToFind*sizeof(int));
	for(i=0; i<2*KSEQ-1; i++){
		if (tree[indexArray[i]].nd > 4){
			printf("cutting at node %d\n",indexArray[i]);
			printdescendants(tree,indexArray[i],clusters,0);
			tree[indexArray[i]].nodeToCut=1;
			break;
		}
	}
	printf("cluster 2:\n");
	printdescendants(tree,root,clusters,1);
	clusterSize[0]=countNumInCluster(clusters,0);
	clusterSize[1]=countNumInCluster(clusters,1);
	printf("number in cluster1: %d\n",clusterSize[0]);
	printf("number in cluster2: %d\n",clusterSize[1]);
	node** treeArr;
	allocateMemForTreeArr(numBranchesToFind,clusterSize,treeArr);
	for(i=0; i<clusterSize[0]; i++){
		for(j=0; j<KSEQ; j++){
			if (strcmp(clusters[0][i],seqNames[chooseK[j]])==0){
				strcpy(seqsInCluster[i],sequences[chooseK[j]]);
			}
		}
	}
	for(i=0; i<KSEQ+1; i++){
		for(j=0; j<KSEQ+1; j++){
			distMat[i][j] = 0;
		}
	}
	createDistMat(seqsInCluster,distMat,clusterSize[0],fasta_specs[1]);
}
