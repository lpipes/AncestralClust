# AncestralClust
AncestralClust

ancestralclust [OPTIONS]
	
	-h, --help			usage: -i file.fasta -t file_taxonomy.txt -d output_directory
	-i, --infile			fasta to cluster
	-t, --infile_taxonomy		taxonomy of fasta to cluster (sorted)
	-n, --number_of_clusters	number of initial clusters
	-k, --number_of_sequences	number of sequences in initial cluster
	-d, --directory			directory to print clusters
	-c, --threads			number of threads
	-o, --output_file		output file
	-f, --fasta_format		output fasta files for each cluster
	-u, --use_nw			use Needleman-Wunsch (default is WFA)
	-l, --number_of_lines_to_read	number of lines to read in from file

AncestralClust uses <a href="https://github.com/TimoLassmann/kalign">kalign3</a> to construct multiple sequence alignments, <a href="https://github.com/smarco/WFA">wavefront alignment algorithm</a> for pairwise alignments, and <a href="https://github.com/noporpoise/seq-align">needleman-wunsch alignment</a> for pairwise alignments if chosen by the user, and <a href="https://github.com/DavidLeeds/hashmap">David Leeds' hashmap</a> for taxonomy files if user chooses.	
