#ifndef PARSE_ARGUMENTS_H

#define PARSE_ARGUMENTS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "parse_arguments.h"
#include "c_tools.h"


void print_usage(char *prog_loc) {
	fprintf(stderr, "usage:"
								"\t%s <mode> [options] file [file, ...]\n"
								"\t%s [-h | --help]\n\n"
								"\twhere <mode> is one of {hist, extract, both}\n\n"
					, prog_loc, prog_loc);
}


void print_help(char *prog_loc) {

	char *diagram;

	print_usage(prog_loc);

	fprintf(stderr, "modes:\n"
						"\thist : only count k-mers and print histogram\n"
						"\textract : extract reads with above 'cutoff' number of k-mers mapping to it\n"
						"\tboth : do both hist and extract\n\n"

					"options (default):\n"
						"\tapplicable in both functions:\n"
							"\t\t-i, --in : location of hash table file - optional: hash table will be computed if not provided\n"
							"\t\t-o, --out : file in which to store hash table - optional\n"
							"\t\t-q, --quiet : supress progress messages normally printed to stderr (false)\n"
							"\t\t-v, --verbose : print each k-mer as it is hashed (only really useful for debugging) (false)\n"
							"\t\t-c, --canonical : count canonical version of k-mers (i.e. the lowest scoring hash of the k-mer and its reverse complement) (false)\n"
							"\t\t-r, --region-size : number of bases in each region (15)\n"
							"\t\t-g, --interval-size : number of bases in gap between each region (0)\n\n"

						"\tonly applicable in extract function:\n"
							"\t\t-a, --min : minimum number of occurrences of k-mer for it to be masked on read (1)\n"
							"\t\t-b, --max : maximum number of occurrences of k-mer for it to be masked on read (999)\n"
							"\t\t-u, --cutoff : minimum number of k-mers mapped to read for read to be printed (see notes)\n"
							"\t\t-x, --max-difference : maximum difference between number of k-mer hits and number of possible k-mer hits for the read (see notes)\n"
							"\t\t-d, --disable-mask : leave bases not occurring in desired k-mer peaks unmasked when extracting reads (faslse)\n\n"

						"\tmisc:\n"
							"\t\t-h, --help : print this message\n\n"

					"notes:\n"
						"\t* If neither --cutoff nor --max-difference is specified but one is required, cutoff defaults to 50\n"
						"\t* A maximuim of one of --in and --out may be specified by the user\n"
						"\t* --quiet and --verbose are mutually exclusive\n"
						"\t* --min cannot be greater than --max\n"
						"\t* --region-size must be 1, 3, 5, or 15\n\n");

	diagram =	"diagram:\n"
				"\t---------------------------------------------------------------\n"
				"\t|                                                             |\n"
				"\t|    3     10     3     10     3     10     3     10     3    |\n"
				"\t|   |||----------|||----------|||----------|||----------|||   |\n"
                "\t|           ^-- interval       ^-- region                     |\n"
				"\t|   <---------------------- window ----------------------->   |\n"
				"\t|                                                             |\n"
				"\t---------------------------------------------------------------\n";

	fprintf(stderr, "%s\n", diagram);
	exit(EXIT_FAILURE);
}


argument_struct parse_arguments(int argc, char **argv) {

	argument_struct to_return;

	int arg_i; /* Argument parser for loop counter */
	bool argument_error = false; /* Set to true if we need to quit after all error checking has taken place */
	bool help_requested = false; /* Set to true if we need to print the help message after we have looked through the remaining arguments */

	to_return.print_hist = false;
	to_return.extract_reads = false;
	to_return.min_kmer_hits = -1;
	to_return.max_kmers_missed = -1;
	to_return.min_val = 0;
	to_return.max_val = 0;
	to_return.quiet = false;
	to_return.verbose = false;
	to_return.use_canonical = false;
	to_return.kmer_size = 0;
	to_return.mask = 2; /* 0 = no masking; 1 = strict mask; 2 = normal mask */
	to_return.where_to_save_hash_table = NULL;
	to_return.stored_hash_table_location = NULL;
	to_return.region_size = -1;
	to_return.interval_size = -1;
	to_return.index_first_file = argc - 1;

	if (argc <= 2) {
		if (argc == 2) {
			if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
				print_help(argv[0]);
			}
			else {
				print_usage(argv[0]);
				exit(EXIT_FAILURE);
			}
		}
		else {
			print_usage(argv[0]);
			exit(EXIT_FAILURE);
		}
	}

	if (!strcmp(argv[1], "both")) {
		to_return.print_hist = true;
		to_return.extract_reads = true;
	}

	else if (!strcmp(argv[1], "hist")) {
		to_return.print_hist = true;
	}

	else if (!strcmp(argv[1], "extract")) {
		to_return.extract_reads = true;
	}

	else {
		fprintf(stderr, "ERROR: Mode not recognised\n");
		print_usage(argv[0]);
		exit(EXIT_FAILURE);
	}

	for (arg_i = 2; arg_i < argc - 1; arg_i++) {

		if (!strcmp(argv[arg_i], "-h") || !strcmp(argv[arg_i], "--help")) {
			help_requested = true;
		}

		else if (!strcmp(argv[arg_i], "-u") || !strcmp(argv[arg_i], "--cutoff")) {
			if (!to_return.extract_reads) {
				fprintf(stderr, "ERROR: -u/--cutoff must not be specified in this mode\n");
				argument_error = true;
			}
			if (is_str_integer(argv[++arg_i])) {
				to_return.min_kmer_hits = atoi(argv[arg_i]);
				if (to_return.min_kmer_hits < 0) {
					fprintf(stderr, "ERROR: -u/--cutoff must be a non-negative integer\n");
					argument_error = true;
				}
			}
			else {
				fprintf(stderr, "ERROR: -u/--cutoff must be a non-negative integer\n");
				argument_error = true;
			}
		}

		else if (!strcmp(argv[arg_i], "-x") || !strcmp(argv[arg_i], "--max-difference")) {
			if (!to_return.extract_reads) {
				fprintf(stderr, "ERROR: -x/--max-difference must not be specified in this mode\n");
				argument_error = true;
			}
			if (is_str_integer(argv[++arg_i])) {
				to_return.max_kmers_missed = atoi(argv[arg_i]);
				if (to_return.max_kmers_missed < 0) {
					fprintf(stderr, "ERROR: -x/--max-difference must be a non-negative integer\n");
					argument_error = true;
				}
			}
			else {
				fprintf(stderr, "ERROR: -x/--max-difference must be a non-negative integer\n");
				argument_error = true;
			}
		}

		else if (!strcmp(argv[arg_i], "-a") || !strcmp(argv[arg_i], "--min")) {
			if (!to_return.extract_reads) {
				fprintf(stderr, "ERROR: -a/--min-val must not be specified in this mode\n");
				argument_error = true;
			}
			if (is_str_integer(argv[++arg_i])) {
				if (atoi(argv[arg_i]) < 0) {
					fprintf(stderr, "ERROR: -a/--min-val must be a non-negative integer\n");
					argument_error = true;
				}
				else {
					to_return.min_val = atoi(argv[arg_i]);
				}
			}
			else {
				fprintf(stderr, "ERROR: -a/--min-val must be a non-negative integer\n");
				argument_error = true;
			}
		}

		else if (!strcmp(argv[arg_i], "-b") || !strcmp(argv[arg_i], "--max")) {
			if (!to_return.extract_reads) {
				fprintf(stderr, "ERROR: -b/--max-val must not be specified in this mode\n");
				argument_error = true;
			}
			if (is_str_integer(argv[++arg_i])) {
				if (atoi(argv[arg_i]) < 0) {
					fprintf(stderr, "ERROR: -b/--max-val must be a non-negative integer\n");
					argument_error = true;
				}
				else {
					to_return.max_val = atoi(argv[arg_i]);
				}
			}
			else {
				fprintf(stderr, "ERROR: -b/--max-val must be a non-negative integer\n");
				argument_error = true;
			}
		}

		else if (!strcmp(argv[arg_i], "-q") || !strcmp(argv[arg_i], "--quiet")) {
			to_return.quiet = true;
		}

		else if (!strcmp(argv[arg_i], "-v") || !strcmp(argv[arg_i], "--verbose")) {
			to_return.verbose = true;
		}

		else if (!strcmp(argv[arg_i], "-c") || !strcmp(argv[arg_i], "--canonical")) {
			to_return.use_canonical = true;
		}

		else if (!strcmp(argv[arg_i], "-k") || !strcmp(argv[arg_i], "--kmer-size")) {
			if (to_return.kmer_size != 0) {
				fprintf(stderr, "ERROR: -k/--kmer-size specified more than once\n");
				argument_error = true;
			}
			if (is_str_integer(argv[++arg_i])) {
				to_return.kmer_size = atoi(argv[arg_i]);
				if (to_return.kmer_size != 13 && to_return.kmer_size != 15 && to_return.kmer_size != 17) {
					fprintf(stderr, "ERROR: -k/--kmer-size must be 13, 15, or 17\n");
					argument_error = true;
				}
			}
			else {
				fprintf(stderr, "ERROR: -k/--kmer-size must be 13, 15, or 17\n");
				argument_error = true;
			}
		}

		else if (!strcmp(argv[arg_i], "-d") || !strcmp(argv[arg_i], "--disable-mask")) {
			if (!to_return.extract_reads) {
				fprintf(stderr, "ERROR: Mask must not be specified in this mode\n");
				argument_error = true;
			}
			if (to_return.mask == 1) {
				fprintf(stderr, "ERROR: --disable-mask cannot be used with --strict-mask\n");
				argument_error = true;
			}
			else {
				to_return.mask = 0;
			}
		}
		else if (!strcmp(argv[arg_i], "-s") || !strcmp(argv[arg_i], "--strict-mask")) {
			if (!to_return.extract_reads) {
				fprintf(stderr, "ERROR: Mask must not be specified in this mode\n");
				argument_error = true;
			}
			if (to_return.mask == 0) {
				fprintf(stderr, "ERROR: --strict-mask cannot be used with --disable-mask\n");
				argument_error = true;
			}
			else {
				to_return.mask = 1;
			}
		}

		else if (!strcmp(argv[arg_i], "-i") || !strcmp(argv[arg_i], "--in")) {
			to_return.stored_hash_table_location = argv[++arg_i];
		}

		else if (!strcmp(argv[arg_i], "-o") || !strcmp(argv[arg_i], "--out")) {
			to_return.where_to_save_hash_table = argv[++arg_i];
		}

		else if (!strcmp(argv[arg_i], "-r") || !strcmp(argv[arg_i], "--region-size")) {
			if (is_str_integer(argv[++arg_i])) {
				to_return.region_size = atoi(argv[arg_i]);
				if (to_return.region_size < 0) {
					fprintf(stderr, "ERROR: -r/--region-size must be a non-negative integer\n");
					argument_error = true;
				}
			}
			else {
				fprintf(stderr, "ERROR: -r/--region-size must be a non-negative integer\n");
				argument_error = true;
			}
		}

		else if (!strcmp(argv[arg_i], "-g") || !strcmp(argv[arg_i], "--interval-size")) {
			if (is_str_integer(argv[++arg_i])) {
				to_return.interval_size = atoi(argv[arg_i]);
				if (to_return.interval_size < 0) {
					fprintf(stderr, "ERROR: -g/--interval-size must be a non-negative integer\n");
					argument_error = true;
				}
			}
			else {
				fprintf(stderr, "ERROR: -g/--interval-size must be a non-negative integer\n");
				argument_error = true;
			}
		}

		else {
			to_return.index_first_file = arg_i;
			break;
		}
	}


	/* ----- Error check user input ----- */


	if (to_return.kmer_size == 0) {
		fprintf(stderr, "ERROR: -k/--kmer-size must be specified\n");
		argument_error = true;
	}

	if (to_return.extract_reads) {
		if (to_return.min_val == 0 || to_return.max_val == 0) {
			fprintf(stderr, "ERROR: -a/--min-val and -b/--max-val must both be specified and non-zero when extracting reads\n");
			argument_error = true;
		}

		if (to_return.max_val < to_return.min_val) {
			fprintf(stderr, "ERROR: -a/--min-val must not be greater than -b/--max-val\n");
			argument_error = true;
		}
	}

	if (to_return.kmer_size == 15) {
		if (to_return.region_size != -1) {
			if (to_return.region_size != 1 && to_return.region_size != 3 && to_return.region_size != 5 && to_return.region_size != 15) {
				fprintf(stderr, "ERROR: Region size must be either 1, 3, 5, or 15\n");
				argument_error = true;
			}
		}
	}

	else {
		if (to_return.region_size != -1) {
			fprintf(stderr, "ERROR: -r/--region-size must not be specified if -k/--kmer-size is not 15\n");
			argument_error = true;
		}

		if (to_return.interval_size != -1) {
			fprintf(stderr, "ERROR: -g/--interval-size must not be specified if -k/--kmer-size is not 15\n");
			argument_error = true;
		}

		if (to_return.mask == 1) {
			fprintf(stderr, "ERROR: -s/--strict cannot be used if -k/--kmer-size is not 15\n");
			argument_error = true;
		}
	}

	if (to_return.stored_hash_table_location && to_return.where_to_save_hash_table) {
		fprintf(stderr, "ERROR: Cannot specify both -i/--in and -o/--out\n");
		argument_error = true;
	}

	if (to_return.quiet && to_return.verbose) {
		fprintf(stderr, "ERROR: Cannot enable both -q/--quiet and -v/--verbose modes\n");
		argument_error = true;
	}

	if (help_requested) {
		if (argument_error) {
			putchar('\n');
		}
		print_help(argv[0]);
	}

	if (argument_error) {
		putchar('\n');
		print_usage(argv[0]);
		exit(EXIT_FAILURE);
	}

	/* ----- End of user input error checking ----- */


	return to_return;
}

#endif /* PARSE_ARGUMENTS_H */
