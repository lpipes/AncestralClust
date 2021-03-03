# AncestralClust
AncestralClust was developed to cluster divergent sequences. A neighbor-joining phylogenetic tree is built from a random number of user-specified initial sequences, the tree is split into initial clusters based on the longest branch lengths, each initial cluster is aligned with a multiple sequence alignment and the ancestral sequence from each cluster is inferred by the root node after being midpoint rooted. The rest of the sequences are assigned to the clusters based on their closest genetic distance to the ancestral sequence (or saved for a future iteration if the distance is greater than the average distance between clusters). The process reiterates until all sequences are assigned to a cluster.

<img src="https://github.com/lpipes/AncestralClust/blob/master/cluster_ancestral.png?raw=true">

ancestralclust [OPTIONS]
	
	-h, --help			usage: -i file.fasta
	-i, --infile [REQUIRED]		fasta to cluster
	-t, --infile_taxonomy		taxonomy of fasta to cluster [sorted in same order as fasta]
	-n, --number_of_clusters	number of initial clusters [default: 10]
	-k, --number_of_sequences	number of sequences in initial cluster [default: 100]
	-d, --directory			directory to print clusters [directory must exist]
	-c, --threads			number of threads
	-o, --output_file		output file
	-f, --fasta_format		output fasta files for each cluster
	-u, --use_nw			use Needleman-Wunsch [default is WFA]
	-l, --number_of_lines_to_read	number of lines to read in from file
	-p, --number_of_descendants	number of descendants to require to cut branch [default: 10]	

AncestralClust uses <a href="https://github.com/TimoLassmann/kalign">kalign3</a> to construct multiple sequence alignments, <a href="https://github.com/smarco/WFA">wavefront alignment algorithm</a> for pairwise alignments, and <a href="https://github.com/noporpoise/seq-align">needleman-wunsch alignment</a> for pairwise alignments if chosen by the user, and <a href="https://github.com/DavidLeeds/hashmap">David Leeds' hashmap</a> for taxonomy files if user chooses.

# Installation
To install

	git clone https://github.com/lpipes/AncestralClust.git
	cd AncestralClust
	make

On macOS, -fopenmp is required for installation.

# Usage
The only required argument is the -i, --infile which is a FASTA file to be clustered. The taxonomy file is optional and is not used to inform the clustering. The taxonomy file must be in the same order as the FASTA file with the following format (the FASTA header followed by a tab followed by the taxonomy):

	FASTA_header\tdomain;phylum;class;order;family;genus;species

# Performance
AncestralClust tends to have more even clusters (as measured by the Coefficient of Variation) with higher relative NMI for every taxonomic level than leading clustering software UCLUST measured from three different metabarcode libraries (16S, 18S, and COI).
<img src="https://github.com/lpipes/AncestralClust/blob/master/RelativeNMI_species.png?raw=true">

# Citation
Pipes L, and Nielsen R (2021) AncestralClust: Clustering of Divergent Nucleotide Sequences by Ancestral Sequence Reconstruction using Phylogenetic Trees. biorxiv. 
<a href="https://www.biorxiv.org/content/10.1101/2021.01.08.426008v1">https://www.biorxiv.org/content/10.1101/2021.01.08.426008v1</a>
