#include "options.h"

static struct option long_options[]=
{
	{"help", no_argument, 0, 'h'},
	{"infile", required_argument, 0, 'i'},
	{"infile_taxonomy", optional_argument, 0, 't'},
	{"number_of_clusters", optional_argument, 0, 'n'},
	{"number_of_sequences", optional_argument, 0, 'k'},
	{"directory", optional_argument, 0, 'd'},
	{"threads", optional_argument, 0, 'c'},
	{"output_file", optional_argument, 0, 'o'},
	{"clstr_format", optional_argument, 0, 's'},
	{"fasta_to_assign", optional_argument, 0, 'f'},
	{"use_nw", optional_argument, 0, 'u'},
	{"number_of_lines_to_read", optional_argument, 0, 'l'},
	{"number_of_descendants", optional_argument, 0, 'p'},
	{0,0,0,0}
};

char usage[] = "\nancestralclust [OPTIONS]\n\
	\n\
	-h, --help				usage: -i [Input FASTA] -r [Integer <= Number of Sequences] -k [Integer > 1]\n\
	-i, --infile [REQUIRED]			fasta to cluster\n\
	-t, --infile_taxonomy [OPTIONAL]	taxonomy of fasta [Sorted in same order as FASTA, not used in clustering]\n\
	-b, --number_of_clusters [REQUIRED]	number of initial clusters [default: 10]\n\
	-r, --number_of_sequences [REQUIRED]	number of sequences in initial cluster [default: 100]\n\
	-d, --directory				directory to print clusters [DIRECTORY MUST EXIST PRIOR TO RUNNING]\n\
	-c, --threads				number of threads [default: 1]\n\
	-o, --output_file			output file\n\
	-f, --fasta_format			output fasta files for each cluster\n\
	-u, --use_nw				use Needleman-Wunsch [default is WFA]\n\
	-l, --number_of_lines_to_read		number of lines to read in from file\n\
	-p, --number_of_descendants		number of descendants to require to cut branch [default: 10]\n\
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
		c=getopt_long(argc,argv,"hfun:b:d:i:t:c:o:l:p:r:",long_options, &option_index);
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
			case 'f':
				opt->output_fasta=1;
				break;
			case 'u':
				opt->use_nw=1;
				break;
			case 'l':
				success = sscanf(optarg, "%d", &(opt->numberOfLinesToRead));
				if (!success)
					fprintf(stderr, "Could not read number of lines to read\n");
				break;
			case 't':
				success = sscanf(optarg, "%s", opt->taxonomy);
				opt->hasTaxFile=1;
				if (!success)
					fprintf(stderr, "Invalid taxonomy file\n");
				break;
			case 'p':
				success = sscanf(optarg, "%d", &(opt->number_of_desc));
				if (!success)
					fprintf(stderr, "Could not read number of descendants\n");
				break;
			case 'b':
				success = sscanf(optarg, "%d", &(opt->number_of_clusters));
				if (!success)
					fprintf(stderr, "Could not read number of clusters\n");
				break;
			case 'r':
				success = sscanf(optarg, "%d", &(opt->number_of_kseqs));
				if (!success)
					fprintf(stderr, "Could not read number of initial sequences to choose\n");
				break;
			case 'd':
				success = sscanf(optarg, "%s", opt->output_directory);
				if (optarg[strlen(optarg)-1]=='/'){
					opt->slash=1;
				}
				opt->default_directory=0;
				if (!success)
					fprintf(stderr, "Invalid output directory\n");
				break;
			case 'o':
				success = sscanf(optarg, "%s", opt->output_file);
				if (!success)
					fprintf(stderr, "Invalid output file\n");
				break;
			case 'c':
				success = sscanf(optarg, "%d", &(opt->numthreads));
				if (!success)
					fprintf(stderr, "Could not read number of threads\n");
				break;
		}
	}
}
