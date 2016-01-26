#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>

#include "c_tools.h"
#include "fastlib.h"
#include "zkc2.h"

#define NUM_CELLS_HASH_TABLE 1073741824 /* = 4^15 */
#define HISTOGRAM_SIZE 10001


void read_hash_table_from_file(uint32_t *hash_table, char *hash_table_location, bool quiet) {

	FILE* input_file;

	if (!quiet) {
		fprintf(stderr, "Reading hash table from file\n");
	}
	input_file = fopen(hash_table_location, "rb");

	if (input_file == NULL) {
		fprintf(stderr, "WARNING: Failed to open hash table file\n");
		exit(EXIT_FAILURE);
	}
	
	if (fread(hash_table, sizeof(uint32_t), NUM_CELLS_HASH_TABLE, input_file) != NUM_CELLS_HASH_TABLE) {
		fprintf(stderr, "WARNING: Failed to load hash table from file\n");
		exit(EXIT_FAILURE);
	}

	if (!quiet) {
		fprintf(stderr, "Successfully read hash table from file\n");
	}

	return;

}


void write_hash_table_to_file(uint32_t *hash_table, char *hash_file_name, bool quiet) {

	FILE *out_file;

	if (!quiet) {
		fprintf(stderr, "Writing hash table to file\n");
	}
	out_file = fopen(hash_file_name, "wb");

	if (out_file == NULL) {
		fprintf(stderr, "WARNING: Failed to create hash table file - it has not been written\n");
		return;
	}

	if (fwrite(hash_table, sizeof(uint32_t), NUM_CELLS_HASH_TABLE, out_file) != NUM_CELLS_HASH_TABLE) {
		fprintf(stderr, "WARNING: Did not manage to write hash table to file\n");
		return;
	}

	if (!quiet) {
		fprintf(stderr, "Successfully wrote hash table to file\n");
	}

	return;
}


int hash_base (char base) {
	int IntBase;
	if     (base == 'A' || base == 'a') IntBase = 0;
	else if(base == 'C' || base == 'c') IntBase = 1;
	else if(base == 'G' || base == 'g') IntBase = 2;
	else if(base == 'T' || base == 't') IntBase = 3;
	else
	{
		IntBase = -1;
	}
	return IntBase;
}


seq_hash_return hash_sequence(char *seq, int region_size, int interval_size, int window_size) {

	uint32_t seq_hash = 0;
	int base_hash;
	int i;
	seq_hash_return to_return;

	to_return.found_n = false;

	for (i = 0; i < window_size; i++) {
		base_hash = hash_base(seq[i]);
		if (base_hash != -1) {
			if ((i % (region_size + interval_size)) < region_size) {
				seq_hash += base_hash;
				seq_hash <<= 2;
			}
		}
		else {
			to_return.found_n = true;
			seq_hash = 0;
			break;
		}
	}

	seq_hash >>= 2;
	to_return.hash = seq_hash;

	return to_return;
}


uint32_t hash_rc(uint32_t seq_hash) {

	uint32_t mask = -4; /* All bits should be set to 1 except least significant two */
	uint32_t rc_hash = 0;
	int i;

	for (i = 0; i < 14; i++) {
		rc_hash += ~((seq_hash & 3) | mask);
		rc_hash <<= 2;
		seq_hash >>= 2;
	}
	rc_hash += ~((seq_hash & 3) | mask);

	return rc_hash;
}


void decode_hash(uint32_t hash, unsigned int region_size, unsigned int window_size, unsigned int interval_size) {

	unsigned int j;

	hash <<= 2;
	 
	for (j = 0; j < window_size; j++) {
		if ((j % (region_size + interval_size)) < region_size) {
			if ((hash & (3 << 30)) == 0) {
				putchar('A');
			}
			else if ((hash & (3 << 30)) == (1UL << 30)) {
				putchar('C');
			}
			else if ((hash & (3 << 30)) == (2UL << 30)) {
				putchar('G');
			}
			else if ((hash & (3 << 30)) == (3UL << 30)) {
				putchar('T');
			}
			else {
				printf("ERROR\n");
			}
			hash <<= 2;
		}
		else {
			putchar('-');
		}
	}
}


void print_usage(char *prog_loc) {

	char *diagram;

	fprintf(stderr, "\nusage: %s mode [options] <file>\n\n"

					"modes:\n"
						"\thist : only count k-mers and print histogram\n"
						"\textract : extract reads with above 'cutoff' number of k-mers mapping to it\n"
						"\tboth : do both hist and extract\n\n"

					"options (default):\n"
						"\t--cutoff : minimum number of k-mers mapped to read for read to be printed (50)\n"
						"\t--min : minimum number of occurrences of k-mer for it to be masked on read (1)\n"
						"\t--max : maximum number of occurrences of k-mer for it to be masked on read (999)\n"
						"\t--out : file in which to store hash table - optional\n"
						"\t--in : location of hash table file - optional: hash table will be computed if not provided\n"
						"\t--quiet : supress progress messages normally printed to stderr (false)\n"
						"\t--verbose : print each k-mer as it is hashed (only really useful for debugging) (false)\n"
						"\t--rc : count canonical version of k-mers (i.e. the lowest scoring hash of the k-mer and its reverse complement) (false)\n"
						"\t--region-size : number of bases in each region (15)\n"
						"\t--interval : number of bases in gap between each region (0)\n\n"

					"notes:\n"
						"\t* A maximuim of one of --in and --out may be specified by the user\n"
						"\t* --quiet and --verbose are mutually exclusive\n"
						"\t* --min cannot be greater than --max\n"
						"\t* --region-size must be 1, 3, 5, or 15\n\n", prog_loc);

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


new_hashes shift_hash(uint32_t current_seq_hash, uint32_t current_rc_hash, int num_regions, int *base_hash_array) {

	new_hashes to_return;
	int iCount;
	uint32_t seq_mask; 
	uint32_t rc_mask;
	int jump = 30 / num_regions; /* Distance to next region */

	if (num_regions == 1) {
		seq_mask = 1073741820;	/* = 0011 1111 1111 1111 1111 1111 1111 1100 */
		rc_mask = 268435455;	/* = 0000 1111 1111 1111 1111 1111 1111 1111 */
	}

	else if (num_regions == 3) {
		seq_mask = 1070593020;	/* = 0011 1111 1100 1111 1111 0011 1111 1100 */
		rc_mask = 267648255;	/* = 0000 1111 1111 0011 1111 1100 1111 1111 */
	}

	else if (num_regions == 5) {
		seq_mask = 1022611260;	/* = 0011 1100 1111 0011 1100 1111 0011 1100 */
		rc_mask = 255652815;	/* = 0000 1111 0011 1100 1111 0011 1100 1111 */
	}

	else if (num_regions == 15) {
		seq_mask = 0;			/* = 0000 0000 0000 0000 0000 0000 0000 0000 */
		rc_mask = 0;			/* = 0000 0000 0000 0000 0000 0000 0000 0000 */
	}

	else {
		fprintf(stderr, "INTERNAL ERROR: Invalid number of regions\n");
		exit(EXIT_FAILURE);
	}

	current_seq_hash <<= 2;
	current_rc_hash >>= 2;

	/* Zero two least significant bits of each region, and two most significant bits overall, of 32-bit int */
	current_seq_hash &= seq_mask;
	/* Zero two most significant bits of each region, and two most significant bits overall, of 32-bit int */
	current_rc_hash &= rc_mask;

	to_return.new_hash = current_seq_hash;
	to_return.new_rc_hash = current_rc_hash;

	for (iCount = 0; iCount < num_regions; iCount++) {
		to_return.new_hash += base_hash_array[iCount] << (jump * (num_regions - iCount - 1));
		to_return.new_rc_hash += ((base_hash_array[num_regions - iCount - 1] ^ 3) << (28 - (iCount * jump)));
	}

	to_return.canonical_hash = (to_return.new_hash < to_return.new_rc_hash) ? to_return.new_hash : to_return.new_rc_hash;

	return to_return;
}


new_hashes hash_new_window(uint32_t current_seq_hash) {
	new_hashes to_return;

	to_return.new_hash = current_seq_hash;
	to_return.new_rc_hash = hash_rc(current_seq_hash);

	to_return.canonical_hash = (to_return.new_hash < to_return.new_rc_hash) ? to_return.new_hash : to_return.new_rc_hash;

	return to_return;
}


int main(int argc, char **argv) {

	FILE *input_file;
	seg_return ret; /* Returned data from get_next_seg (i.e. read info and bEOF) */
	int format;
	uint32_t *hash_table;
	uint32_t hash_val; 
	uint32_t rc_hash;
	uint32_t canonical_hash;
	uint32_t hash_to_use;
	int hash;
	seq_hash_return hash_seq;
	new_hashes new_hashes_triple;
	unsigned int kmer_hits;
	unsigned int min_val = 1;
	unsigned int max_val = 999;
	unsigned int cutoff = 50;
	unsigned int interval_size = 0;		/* Number of bases between regions      |||				|||            |||            |||            |||	*/
	unsigned int region_size = 15;		/* Number of bases in each region       ----------------------------------------------------------------	*/
	unsigned int window_size;			/* Number of bases in window             3      10       3       10     3      10      3      10      3		*/
	unsigned int num_regions;			/*										         ^-- interval           ^-- region							*/
										/*										<--------------------------- window --------------------------->	*/
	char *where_to_save_hash_table = "\0";
	char *stored_hash_table_location = "\0";
	unsigned long i, j; /* Counter */
	unsigned long new_base_loc;
	long hist[HISTOGRAM_SIZE];
	bool extract_reads = false;
	bool print_hist = false;
	bool quiet = false;
	bool verbose = false;
	bool use_canonical = false;
	unsigned long end_newest_kmer = 0; /* Index of the end of the most recently found k-mer word in the desired range. Set to -1 to avoid the first base being unmasked. */
	int arg_i;
	int new_base_hash_array[5];
	unsigned int iCount; 
	int phase; /* 0 => Pass 1; 1 => Pass 2 */

	phase = 999; /* Default phase, meaning that I have forgotten to set it to anything meaningful */

	if (argc <= 2) {
		print_usage(argv[0]);
	}

	if (!strcmp(argv[1], "both")) {
		print_hist = true;
		extract_reads = true;
	}

	else if (!strcmp(argv[1], "hist")) {
		print_hist = true;
	}

	else if (!strcmp(argv[1], "extract")) {
		extract_reads = true;
	}

	else {
		print_usage(argv[0]);
	}

	for (arg_i = 2; arg_i < argc - 1; arg_i++) {
		if ((!strcmp(argv[arg_i], "-c")) || (!strcmp(argv[arg_i], "--cutoff"))) {
			if (is_str_of_digits(argv[++arg_i])) {
				cutoff = atoi(argv[arg_i]);
			}
			else {
				print_usage(argv[0]);
			}
		}
		else if (!strcmp(argv[arg_i], "--min")) {
			if (is_str_of_digits(argv[++arg_i])) {
				min_val = atoi(argv[arg_i]);
			}
			else {
				print_usage(argv[0]);
			}
		}
		else if (!strcmp(argv[arg_i], "--max")) {
			if (is_str_of_digits(argv[++arg_i])) {
				max_val = atoi(argv[arg_i]);
			}
			else {
				print_usage(argv[0]);
			}
		}
		else if (!strcmp(argv[arg_i], "--quiet")) {
			quiet = true;
		}
		else if (!strcmp(argv[arg_i], "--verbose")) {
			verbose = true;
		}
		else if (!strcmp(argv[arg_i], "--rc")) {
			use_canonical = true;
		}
		else if (!strcmp(argv[arg_i], "--in")) {
			stored_hash_table_location = argv[++arg_i];
		}
		else if (!strcmp(argv[arg_i], "--out")) {
			where_to_save_hash_table = argv[++arg_i];
		}
		else if (!strcmp(argv[arg_i], "--region-size")) {
			if (is_str_of_digits(argv[++arg_i])) {
				region_size = atoi(argv[arg_i]);
			}
			else {
				print_usage(argv[0]);
			}
		}
		else if (!strcmp(argv[arg_i], "--interval")) {
			if (is_str_of_digits(argv[++arg_i])) {
				interval_size = atoi(argv[arg_i]);
			}
			else {
				print_usage(argv[0]);
			}
		}
		else {
			print_usage(argv[0]);
		}
	}

	if (region_size != 1 && region_size != 3 && region_size != 5 && region_size != 15) {
		fprintf(stderr, "ERROR: Region size must be either 1, 3, 5, or 15\n");
		exit(EXIT_FAILURE);
	}
	else {
		num_regions = (15 / region_size);
	}

	window_size = ((num_regions - 1) * interval_size) + 15;

	if ((strlen(stored_hash_table_location) != 0) && (strlen(where_to_save_hash_table) != 0)) {
		fprintf(stderr, "ERROR: Cannot specify both --in and --out\n");
		exit(EXIT_FAILURE);
	}

	if (quiet && verbose) {
		fprintf(stderr, "ERROR: Cannot enable both quiet and verbose modes\n");
		exit(EXIT_FAILURE);
	}

	if ((input_file = fopen(argv[argc - 1], "r")) == NULL) {
		fprintf(stderr, "ERROR: Could not open data file\n");
		exit(EXIT_FAILURE);
	}

	if (max_val < min_val) {
		fprintf(stderr, "ERROR: Minimum value must not be greater than maximum value\n");
		exit(EXIT_FAILURE);
	}

	/* Hash table is 4^15 cells which are guaranteed to be capable of holding the count of a k-mer 
	 * providing that the count does not exceed 2^32 - the minimum size of a long). 
	 */

	if ((hash_table = calloc(NUM_CELLS_HASH_TABLE, sizeof(uint32_t))) == NULL) {
		fprintf(stderr, "ERROR: Out of memory\n");
		exit(EXIT_FAILURE);
	}

	format = which_format(input_file);

	if (strlen(stored_hash_table_location) == 0) {
		if (!quiet) {
			fprintf(stderr, "Counting k-mers into hash table\n");
		}
		phase = 0;
	}
	else {

		read_hash_table_from_file(hash_table, stored_hash_table_location, quiet);

		if (print_hist) {
			phase = 1;
		}
		else if (extract_reads && !print_hist) {
			phase = 2;
		}
	}

	if (phase == 999) {
		fprintf(stderr, "INTERNAL ERROR: Phase has not been set correctly (is still at 999)\n");
		exit(EXIT_FAILURE);
	}

	/* 
	 * Phase 0 = construct hash table
	 * Phase 1 = print histogram
	 * Phase 2 = extract reads 
	 */

	while (true) {

		if (phase == 0 || phase == 2) {

			if (phase == 2) {
				if (!quiet) {
					fprintf(stderr, "Extracting reads with desired k-mer coverage\n");
				}

				rewind(input_file);
			}

			do {
				i = 0;

				if (phase == 2) {
					kmer_hits = 0;
					end_newest_kmer = 0;
				}

				ret = get_next_seg(input_file, format);

				if (ret.segment.length < window_size) {

					free(ret.segment.name);
					free(ret.segment.seq);
					if (format == 1) {
						free(ret.segment.qual);
					}
					continue;
				} 

				hash_seq = hash_sequence(ret.segment.seq, region_size, interval_size, window_size);

				while (hash_seq.found_n == true && i <= (ret.segment.length - window_size)) {
					i += 1;
					hash_seq = hash_sequence(ret.segment.seq + i, region_size, interval_size, window_size);
					if (phase == 2) {
						ret.segment.seq[i] = 'N';
					}
				}

				if (hash_seq.found_n == false) {

					new_hashes_triple = hash_new_window(hash_seq.hash);
					hash_val = new_hashes_triple.new_hash;
					rc_hash = new_hashes_triple.new_rc_hash;
					canonical_hash = new_hashes_triple.canonical_hash;

					hash_to_use = use_canonical ? canonical_hash : hash_val;

					i += window_size; 

					if (verbose) {
						decode_hash(hash_val, region_size, window_size, interval_size);
						putchar('\n');
					}

					if (phase == 0) {
						hash_table[hash_to_use] += 1;
					}

					else if (phase == 2) {
						if (hash_table[hash_to_use] >= min_val && hash_table[hash_to_use] <= max_val) {
							kmer_hits++;
							end_newest_kmer = i;
						}

						if ((end_newest_kmer == 0) || ((i - 14) > end_newest_kmer)) {
							ret.segment.seq[i - 14] = 'N';
						}
					}

					i = phase ? i + 1 : i; /* May not be necessary (different initial loop conditions were used previously) */
					for (; i < ret.segment.length; i++) {

						for (iCount = 0; iCount < num_regions - 1; iCount++) {
							/* Can guarantee that only the final new character hashed might be an 'N', as otherwise we would have already found it */
							new_base_loc = i - window_size + region_size + (iCount * (region_size + interval_size));
							hash = hash_base(ret.segment.seq[new_base_loc]);
							new_base_hash_array[iCount] = hash;
						}

						new_base_loc = i - window_size + region_size + (iCount * (region_size + interval_size));
						hash = hash_base(ret.segment.seq[new_base_loc]);

						if (hash != -1) {
							new_base_hash_array[iCount] = hash;
							new_hashes_triple = shift_hash(hash_val, rc_hash, num_regions, new_base_hash_array);
							hash_val = new_hashes_triple.new_hash;
							rc_hash = new_hashes_triple.new_rc_hash;
							canonical_hash = new_hashes_triple.canonical_hash;

							hash_to_use = use_canonical ? canonical_hash : hash_val;

							if (verbose) {
								decode_hash(hash_val, region_size, window_size, interval_size);
								putchar('\n');
							}

							if (phase == 0) {
								hash_table[hash_to_use] += 1;
							}

							else if (phase == 2) {
								if (hash_table[hash_to_use] >= min_val && hash_table[hash_to_use] <= max_val) {
									kmer_hits++;
									end_newest_kmer = i;
								}

								if ((end_newest_kmer == 0) || ((i - 14) > end_newest_kmer)) {
									ret.segment.seq[i - 14] = 'N';
								}
							}
						}

						else {

							if (phase == 2) {
								/* Before moving onto the next k-mer word, mask, if necessary, the remainder of the k-mer word which is going to be skipped */
								for (j = i - 14; j < i; j++) {
									if ((end_newest_kmer == 0) || (j > end_newest_kmer)) {
										ret.segment.seq[j] = 'N';
									}
								}
							}

							i += 1;
							hash_seq = hash_sequence(ret.segment.seq + i, region_size, interval_size, window_size);

							/* Keep hashing the sequence starting at the next base and moving along the window until we don't find any more 'N's */
							while (hash_seq.found_n == true && i < (ret.segment.length - window_size)) {
								if (phase == 2) {
									if ((end_newest_kmer == 0) || (i > end_newest_kmer)) {
										ret.segment.seq[i] = 'N';
									}
								}

								i += 1;
								hash_seq = hash_sequence(ret.segment.seq + i, region_size, interval_size, window_size);
							}

							if (hash_seq.found_n == true) {
								break;
							}

							else {
								new_hashes_triple = hash_new_window(hash_seq.hash);
								hash_val = new_hashes_triple.new_hash;
								rc_hash = new_hashes_triple.new_rc_hash;
								canonical_hash = new_hashes_triple.canonical_hash;

								hash_to_use = use_canonical ? canonical_hash : hash_val;

								/* Move to start of next 'N' free k-mer word */
								i += window_size;

								if (verbose) {
									decode_hash(hash_val, region_size, window_size, interval_size);
									putchar('\n');
								}

								if (phase == 0) {
									hash_table[hash_to_use] += 1;
								}

								else if (phase == 2) {
									if (hash_table[hash_to_use] >= min_val && hash_table[hash_to_use] <= max_val) {
										kmer_hits++;
										end_newest_kmer = i;
									}

									if ((end_newest_kmer == 0) || ((i - 14) > end_newest_kmer)) {
										ret.segment.seq[i - 14] = 'N';
									}
								}
							}
						}
					}

					if (phase == 2) {
						/* If necessary, mask the bases in the final k-mer word */
						for (j = i - 15; j < i; j++) {
							if ((end_newest_kmer == 0) || (j > end_newest_kmer)) {
								ret.segment.seq[j] = 'N';
							}
						}
					}
				}

				if (phase == 2) {
					if (kmer_hits >= cutoff) {
						printf(">%s %d\n%s\n", ret.segment.name, kmer_hits, ret.segment.seq);
					}
				}

				free(ret.segment.name);
				free(ret.segment.seq);
				if (format == 1) {
					free(ret.segment.qual);
				}

			} while (!ret.bEOF);

		}

		if (phase == 0) {
			if (strlen(where_to_save_hash_table) != 0) {
				write_hash_table_to_file(hash_table, where_to_save_hash_table, quiet);
			}

			if (print_hist) {
				phase = 1;
			}
			else if (extract_reads && !print_hist) {
				rewind(input_file);
				phase = 2;
			}
			else {
				break;
			}
		}

		else if (phase == 1) {

			if (!quiet) {
				fprintf(stderr, "Computing histogram\n");
			}

			for (i = 0; i < (HISTOGRAM_SIZE - 1); i++) {
				hist[i] = 0;
			}

			for (i = 0; i < NUM_CELLS_HASH_TABLE; i++) {
				if (hash_table[i] > 0) {
					if (hash_table[i] < HISTOGRAM_SIZE) {
						hist[hash_table[i] - 1]++;
					}
					else {
						hist[HISTOGRAM_SIZE - 1]++;
					}
				}
			}

			for (i = 0; i < (HISTOGRAM_SIZE - 1); i++) {
				if (hist[i] > 0) {
					printf("%lu %ld\n", i + 1, hist[i]);
				}
			}

			if (extract_reads) {
				rewind(input_file);
				phase = 2;
			}	
			else {
				break;
			}
		}

		else if (phase == 2) {
			break;
		}

		else {
			fprintf(stderr, "INTERNAL ERROR: Phase not correctly set\n");
			exit(EXIT_FAILURE);
		}
	}

	free(hash_table);
	fclose(input_file);

	return 0;
}
