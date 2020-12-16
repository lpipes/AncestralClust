#include <stdio.h>
#include "global.h"

int read_sequences_into_msa( char*** clusters, char*** cluster_seqs, struct msa** m, int* clusterSize, int whichCluster){
	struct msa* msa = NULL;
	struct msa_seq* seq_ptr = NULL;
	int i,j;
	for(i=0; i<clusterSize[whichCluster]; i++){
		strcpy(seq_ptr->name,cluster[whichCluster][i]);
		seq_ptr = msa->sequences[msa->numseq];
		msa->numseq++;
		strcpy(seq_ptr->seq[seq_ptr->len],cluster_seqs[whichCluster][i]);
		seq_ptr->len++;
	}
}
