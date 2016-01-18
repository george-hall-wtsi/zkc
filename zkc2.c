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
		fprintf(stderr, "Reading hash table from memory\n");
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
		fprintf(stderr, "Successfully read hash table from memory\n");
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


seq_hash_return hash_sequence (char *seq) {

	uint32_t seq_hash = 0;
	int base_hash;
	int i;
	seq_hash_return to_return;

	to_return.found_n = false;

	for (i = 0; i < 15; i++) {
		base_hash = hash_base(seq[i]);
		if (base_hash != -1) {
			seq_hash += base_hash;
			if (i != 14) {
				seq_hash <<= 2;
			}
		}
		else {
			to_return.found_n = true;
			seq_hash = 0;
			break;
		}
	}

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


void decode_hash(uint32_t hash) {

	char to_return[16];
	int i;
	to_return[15] = '\0';

	for (i = 14; i >= 0; i--) {
		if ((hash & 3) == 0) {
			to_return[i] = 'A';
		}
		else if ((hash & 3) == 1) {
			to_return[i] = 'C';
		}
		else if ((hash & 3) == 2) {
			to_return[i] = 'G';
		}
		else if ((hash & 3) == 3) {
			to_return[i] = 'T';
		}
		hash >>= 2;
	}

	printf("%s ", to_return); 
}


void print_usage(char *prog_loc) {
	fprintf(stderr, "\nusage: %s mode [options] <file>\n\n", prog_loc);
	fprintf(stderr, "modes:\n\thist : only count k-mers and print histogram\n\textract : extract reads with above 'cutoff' number of k-mers mapping to it\n\tboth : do both hist and extract\n\n");
	fprintf(stderr, "options - (default):\n\t--cutoff : minimum number of k-mers mapped to read for read to be printed (50)\n\t--min : minimum number of occurrences of k-mer for it to be masked on read (1)\n");
	fprintf(stderr, "\t--max : maximum number of occurrences of k-mer for it to be masked on read (999)\n");
	fprintf(stderr, "\t--out : file in which to store hash table  - optional\n\t--in : location of hash table file - optional: hash table will be computed if not provided\n");
	fprintf(stderr, "\tNote: A maximuim of one of --in and --out may be specified by the user\n\n");
	exit(EXIT_FAILURE);
}


new_hashes shift_hash(int new_base_hash, uint32_t current_seq_hash, uint32_t current_rc_hash) {

	new_hashes to_return;

	current_seq_hash <<= 2;
	/* Zero two most significant bits */
	current_seq_hash &= 1073741823;

	current_rc_hash >>= 2;

	to_return.new_hash = current_seq_hash;
	to_return.new_rc_hash = current_rc_hash;

	to_return.new_hash += new_base_hash;
	to_return.new_rc_hash += (new_base_hash & 3) << 28;
	to_return.new_rc_hash ^= (3UL << 28);

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
	int hash;
	seq_hash_return hash_seq;
	new_hashes new_hashes_triple;
	int kmer_hits;
	unsigned int min_val = 1;
	unsigned int max_val = 999;
	int cutoff = 50;
	char *where_to_save_hash_table = "\0";
	char *stored_hash_table_location = "\0";
	unsigned long i, j; /* Counter */
	long hist[HISTOGRAM_SIZE];
	bool extract_reads = false;
	bool print_hist = false;
	bool quiet = false;
	bool verbose = false;
	unsigned long end_newest_kmer = 0; /* Index of the end of the most recently found k-mer word in the desired range. Set to -1 to avoid the first base being unmasked. */
	int arg_i;

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
		else if (!strcmp(argv[arg_i], "--in")) {
			stored_hash_table_location = argv[++arg_i];
		}
		else if (!strcmp(argv[arg_i], "--out")) {
			where_to_save_hash_table = argv[++arg_i];
		}
		else {
			print_usage(argv[0]);
		}
	}

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

		/* First pass through file to count k-mers into hash_table */
		do {
			i = 0;
			ret = get_next_seg(input_file, format);

			hash_seq = hash_sequence(ret.segment.seq);

			while (hash_seq.found_n == true && i <= (ret.segment.length - 15)) {
				i += 1;
				hash_seq = hash_sequence(ret.segment.seq + i);
			}

			if (hash_seq.found_n == false) {

				new_hashes_triple = hash_new_window(hash_seq.hash);
				hash_val = new_hashes_triple.new_hash;
				rc_hash = new_hashes_triple.new_rc_hash;
				canonical_hash = new_hashes_triple.canonical_hash;

				hash_table[canonical_hash] += 1;
				i += 15; 

				for (; i < ret.segment.length; i++) {

					hash = hash_base(ret.segment.seq[i]);

					if (hash != -1) {
						new_hashes_triple = shift_hash(hash, hash_val, rc_hash);
						hash_val = new_hashes_triple.new_hash;
						rc_hash = new_hashes_triple.new_rc_hash;
						canonical_hash = new_hashes_triple.canonical_hash;

						hash_table[canonical_hash] += 1;

					}
					else {
						i += 1;
						hash_seq = hash_sequence(ret.segment.seq + i);
						while (hash_seq.found_n == true && i < (ret.segment.length - 15)) {
							i += 1;
							hash_seq = hash_sequence(ret.segment.seq + i);
						}

						if (hash_seq.found_n == true) {
							break;
						}

						else {
							new_hashes_triple = hash_new_window(hash_seq.hash);
							hash_val = new_hashes_triple.new_hash;
							rc_hash = new_hashes_triple.new_rc_hash;
							canonical_hash = new_hashes_triple.canonical_hash;

							hash_table[canonical_hash] += 1;
						}

						/* Move to start of next 'N' free k-mer word */
						i += 14;
					}
				}
			}

			free(ret.segment.name);
			free(ret.segment.seq);
			if (format == 1) {
				free(ret.segment.qual);
			}
		} while (!ret.bEOF);

		if (strlen(where_to_save_hash_table) != 0) {
			write_hash_table_to_file(hash_table, where_to_save_hash_table, quiet);
		}
	}

	else {
		read_hash_table_from_file(hash_table, stored_hash_table_location, quiet);
	}

	if (print_hist) {

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
	}

	if (extract_reads) {

		rewind(input_file);

		if (!quiet) {
			fprintf(stderr, "Extracting reads with desired k-mer coverage\n");
		}

		/* Second pass through file to extract reads with desired k-mer coverage */
		do {
			i = 0;
			kmer_hits = 0;
			end_newest_kmer = 0;
			ret = get_next_seg(input_file, format);

			if (verbose) {
				printf("Name: %s\nSeq: %s\n", ret.segment.name, ret.segment.seq);
			}


			hash_seq = hash_sequence(ret.segment.seq);

			while (hash_seq.found_n == true && i <= (ret.segment.length - 15)) {
				ret.segment.seq[i] = 'N';
				i += 1;
				hash_seq = hash_sequence(ret.segment.seq + i);
			}

			if (hash_seq.found_n == false) {

				new_hashes_triple = hash_new_window(hash_seq.hash);
				hash_val = new_hashes_triple.new_hash;
				rc_hash = new_hashes_triple.new_rc_hash;
				canonical_hash = new_hashes_triple.canonical_hash;

				i += 14; /* Move to end of current k-mer */

				if (hash_table[canonical_hash] >= min_val && hash_table[canonical_hash] <= max_val) {
					kmer_hits++;
					end_newest_kmer = i;
				}

				if (verbose) {
					decode_hash(hash_val);
					printf(" - %lu\n", hash_table[canonical_hash]);
				}

				if ((end_newest_kmer == 0) || ((i - 14) > end_newest_kmer)) {
					ret.segment.seq[i - 14] = 'N';
				}

				for (i = i + 1; i < ret.segment.length; i++) {

					hash = hash_base(ret.segment.seq[i]);

					if (hash != -1) {

						new_hashes_triple = shift_hash(hash, hash_val, rc_hash);
						hash_val = new_hashes_triple.new_hash;
						rc_hash = new_hashes_triple.new_rc_hash;
						canonical_hash = new_hashes_triple.canonical_hash;

						if (hash_table[canonical_hash] >= min_val && hash_table[canonical_hash] <= max_val) {
							kmer_hits++;
							end_newest_kmer = i;
						}

						if (verbose) {
							decode_hash(hash_val);
							printf(" - %lu\n", hash_table[canonical_hash]);
						}

						if ((end_newest_kmer == 0) || ((i - 14) > end_newest_kmer)) {
							ret.segment.seq[i - 14] = 'N';
						}
					}
					else {
						/* Before moving onto the next k-mer word, mask, if necessary, the remainder of the k-mer word which is going to be skipped */
						for (j = i - 14; j < i; j++) {
							if ((end_newest_kmer == 0) || (j > end_newest_kmer)) {
								ret.segment.seq[j] = 'N';
							}
						}

						/* Move to i to beginning of sequence after the 'N' we just found, and hash this sequence */
						i += 1;
						hash_seq = hash_sequence(ret.segment.seq + i);

						/* Keep hashing the sequence starting at the next base and moving along the window until we don't find any more 'N's */
						while (hash_seq.found_n == true && i <= (ret.segment.length - 15)) {
							if ((end_newest_kmer == 0) || (i > end_newest_kmer)) {
								ret.segment.seq[i] = 'N';
							}
							i += 1;
							hash_seq = hash_sequence(ret.segment.seq + i);
						}

						/* i pointing at beginning of new k-mer */

						if (hash_seq.found_n == true) {
							break;
						}

						else {
							new_hashes_triple = hash_new_window(hash_seq.hash);
							hash_val = new_hashes_triple.new_hash;
							rc_hash = new_hashes_triple.new_rc_hash;
							canonical_hash = new_hashes_triple.canonical_hash;

							/* Move to end of new k-mer word */
							i += 14;

							if (hash_table[canonical_hash] >= min_val && hash_table[canonical_hash] <= max_val) {
								kmer_hits++;
								end_newest_kmer = i;
							}

							if (verbose) {
								decode_hash(hash_val);
								printf(" - %lu\n", hash_table[canonical_hash]);
							}

							if ((end_newest_kmer == 0) || ((i - 14) > end_newest_kmer)) {
								ret.segment.seq[i - 14] = 'N';
							}
						}
					}
				}

				/* If necessary, mask the bases in the final k-mer word */
				for (j = i - 15; j < i; j++) {
					if ((end_newest_kmer == 0) || (j > end_newest_kmer)) {
						ret.segment.seq[j] = 'N';
					}
				}
			}

			if (kmer_hits >= cutoff) {
				printf(">%s %d\n%s\n", ret.segment.name, kmer_hits, ret.segment.seq);
			}

			free(ret.segment.name);
			free(ret.segment.seq);
			if (format == 1) {
				free(ret.segment.qual);
			}

		} while (!ret.bEOF);
	}

	free(hash_table);
	fclose(input_file);

	return 0;
}
