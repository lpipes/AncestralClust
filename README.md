# AncestralClust
AncestralClust was developed to cluster divergent sequences. A neighbor-joining phylogenetic tree is built from a random number of user-specified initial sequences, the tree is split into initial clusters based on the longest branch lengths, each initial cluster is aligned with a multiple sequence alignment and the ancestral sequence from each cluster is inferred by the root node after being midpoint rooted. The rest of the sequences are assigned to the clusters based on their closest genetic distance to the ancestral sequence (or saved for a future iteration if the distance is greater than the average distance between clusters). The process reiterates until all sequences are assigned to a cluster.

<img src="https://github.com/lpipes/AncestralClust/blob/master/cluster_ancestral.png?raw=true">

ancestralclust [OPTIONS]
	
	-h, --help				usage: -i [Input FASTA] -r [Integer <= Number of Sequences] -k [Integer > 1]
	-i, --infile [REQUIRED]			fasta to cluster [expects sequence in 1 line (i.e., without line breaks)]
	-t, --infile_taxonomy [OPTIONAL]	taxonomy of fasta [Sorted in same order as FASTA, not used in clustering]
	-b, --number_of_clusters [REQUIRED]	number of initial clusters [default: 10]
	-r, --number_of_sequences [REQUIRED]	number of sequences in initial cluster [default: 100]
	-d, --directory				directory to print clusters [DIRECTORY MUST EXIST PRIOR TO RUNNING]
	-c, --threads				number of threads [default: 1]
	-o, --output_file			output file
	-f, --fasta_format			output fasta files for each cluster
	-u, --use_nw				use Needleman-Wunsch [default is WFA]
	-l, --number_of_lines_to_read		number of lines to read in from file [default: 10000]
	-p, --number_of_descendants		number of descendants to require to cut branch [default: 10]
	-q, --root_seqs				file to print root sequences
	-a, --set_average			set the average branch length [double > 0, default: calculates averge]
	

AncestralClust uses <a href="https://github.com/TimoLassmann/kalign">kalign3</a> to construct multiple sequence alignments, <a href="https://github.com/smarco/WFA">wavefront alignment algorithm</a> for pairwise alignments, and <a href="https://github.com/noporpoise/seq-align">needleman-wunsch alignment</a> for pairwise alignments if chosen by the user, and <a href="https://github.com/DavidLeeds/hashmap">David Leeds' hashmap</a> for taxonomy files if user chooses.

# Installation
To install

	git clone https://github.com/lpipes/AncestralClust.git
	cd AncestralClust
	make

On macOS, -fopenmp is required for installation.

If you get the following message after you `cd AncestralClust`:
```
make
make: Nothing to be done for 'all'.
```
Do `make clean` and `make` again.

# Usage
The only required arguments are the -i, --infile which is a FASTA file to be clustered, -b, --number_of_clusters which is the number of clusters desired, and -r, --number_of_sequences which is the number of initial sequences chosen at random for the initial clustering. The taxonomy file is optional and is not used to inform the clustering. The taxonomy file must be in the same order as the FASTA file with the following format (the FASTA header followed by a tab followed by the taxonomy):

	FASTA_header\tdomain;phylum;class;order;family;genus;species

AncestralClust also expects no line breaks in your FASTA file. You can remove linebreaks by (replace `test.fasta` with your infile FASTA):

	sed -i ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' test.fasta

AncestralClust also expects that your sequences only contain A, G, C, T, or N (no ambiguous nucleotides or gaps `-`). Some examples on how to run AncestralClust. Using 500 initial sequences with an approximate 10 clusters:

	ancestralclust -i SeqsToCluster.fasta -r 500 -b 10

Printing out the fasta files for each cluster:

	ancestralclust -i SeqsToCluster.fasta -r 500 -b 10 -d clusters -f

# Performance
AncestralClust tends to have more even clusters (as measured by the Coefficient of Variation) with higher relative NMI for every taxonomic level than leading clustering software UCLUST measured from three different metabarcode libraries (16S, 18S, and COI).
<img src="https://github.com/lpipes/AncestralClust/blob/master/RelativeNMI_species.png?raw=true">

# Limits
The maximum number of clusters is set to 100. To set this number higher change line 17 in global.h:

	#define MAXNUMBEROFCLUSTERS

The maximum number of sequences in a cluster is set to 10,000. To set this number higher change line 19 in global.h:

	#define MAXNUMINCLUSTER 10000

The maximum number of sequences randomly chosen in the initial clusters is set to 10,000. To set this number higher change line 18 in global.h:

	#define MAXNUMBEROFKSEQS 10000

The maximum length for sequences is set to 6000bp. To set this number higher change line 10 in global.h:

	#define FASTA_MAXLINE 6000

The minimum length for a sequence is set to 100bp. To set this number lower change line 11 in global.h:

	#define MIN_SEQ 100
 
After changing global.h, compile the program with
	
	make

# Citation
Pipes L, and Nielsen R (2021) AncestralClust: Clustering of Divergent Nucleotide Sequences by Ancestral Sequence Reconstruction using Phylogenetic Trees. biorxiv. 
<a href="https://www.biorxiv.org/content/10.1101/2021.01.08.426008v3">https://www.biorxiv.org/content/10.1101/2021.01.08.426008v3</a>

# Data from manuscript
The data from the manuscript is contained in this repository: <a href="https://github.com/lpipes/AncestralClust_manuscript">https://github.com/lpipes/AncestralClust_manuscript</a>
