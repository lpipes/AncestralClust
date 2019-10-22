#include "options.h"

static struct Options long_options[]=
{
	{"help", no_argument, 0, 'h'},
	{"infile", required_argument, 0, 'i'},
	{"infile_taxonomy", required_argument, 0, 't'},
	{"number_of_clusters", required_argument, 0, 'n'},
	{"number_of_sequences", required_argument, 0, 'k'},
	{"directory", required_argument, 0, 'd'}
};

char usage[] = "\nNJclust [OPTIONS]\n\
	\n\
	-h, --help			usage: -i file.fasta -t file_taxonomy.txt\n\
	-i, --infile			fasta to cluster\n\
	-t, --infile_taxonomy		taxonomy of fasta to cluster (sorted)\n\
	-n, --number_of_clusters	number of initial clusters\n\
	-k, --number_of_sequences	number of sequences in initial cluster\n\
	-d, --directory			directory to print clusters\n\
	\n";

void print_help_statement(){
	printf("%s", &usage[0]);
	return;
}

void parse_options(int argc, char **argv, Options *opt){
	int option_index, success;
	char c;
	if (argc==1){
		print_help_statement();
		exit(0);
	}
	while(1){
		c=getopt_long(argc,argv,"hn:k:d:i:t:",long_options, &option_index);
		if (c==-1) break;
		switch(c){
			case 'h':
				print_help_statement();
				exit(0);
				break;
			case 'i':
				success = sscanf(optarg, "%s", opt->fasta);
				if (!success)
					fprintf(stderr, "Invalid fasta file\n");
				break;
			case 't':
				success = sscanf(optarg, "%s", opt->taxonomy);
				if (!success)
					fprintf(stderr, "Invalid taxonomy file\n");
				break;
			case 'n':
				success = sscanf(optarg, "%d", &(opt->number_of_clusters));
				if (!success)
					fprintf(stderr, "Could not read number of clusters\n");
				break;
			case 'k':
				success = sscanf(optarg, "%d", &(opt->number_of_kseqs));
				if (!success)
					fprintf(stderr, "Could not read number of initial sequences to choose\n");
				break;
			case 'd':
				success = sscanf(optarg, "%s", opt->output_directory);
				if (!success)
					fprintf(stderr, "Invalid output directory\n");
				break;
		}
	}
}
