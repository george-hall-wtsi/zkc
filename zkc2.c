
/*******************************************************************************
 * Copyright (c) 2015-2016 Genome Research Ltd. 
 *  
 * Author: George Hall <gh10@sanger.ac.uk> 
 * 
 * This file is part of K-mer Toolkit. 
 * 
 * K-mer Toolkit is free software: you can redistribute it and/or modify it under 
 * the terms of the GNU General Public License as published by the Free Software 
 * Foundation; either version 3 of the License, or (at your option) any later 
 * version. 
 *  
 * This program is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
 * details. 
 *  
 * You should have received a copy of the GNU General Public License along with 
 * this program. If not, see <http://www.gnu.org/licenses/>. 
 ******************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <inttypes.h>

#include "c_tools.h"
#include "fastlib.h"
#include "zkc2.h"


void read_hash_table_from_file(uint32_t *hash_table, char *hash_table_location, bool quiet, uint64_t num_cells_hash_table) {

	FILE *input_file;

	if (!quiet) {
		fprintf(stderr, "Reading hash table from file\n");
	}
	input_file = fopen(hash_table_location, "rb");

	if (input_file == NULL) {
		fprintf(stderr, "ERROR: Failed to open hash table file\n");
		exit(EXIT_FAILURE);
	}
	
	if (fread(hash_table, sizeof(uint32_t), num_cells_hash_table, input_file) != num_cells_hash_table) {
		fprintf(stderr, "ERROR: Failed to load hash table from file\n");
		exit(EXIT_FAILURE);
	}

	if (fclose(input_file) != 0) {
		if (!quiet) {
			fprintf(stderr ,"WARNING: Failed to close input file - continuing anyway\n");
		}
	}

	else if (!quiet) {
		fprintf(stderr, "Successfully read hash table from file\n");
	}


	return;
}


void write_hash_table_to_file(uint32_t *hash_table, char *hash_file_name, bool quiet, uint64_t num_cells_hash_table) {

	FILE *out_file;

	if (!quiet) {
		fprintf(stderr, "Writing hash table to file\n");
	}
	out_file = fopen(hash_file_name, "wb");

	if (out_file == NULL) {
		fprintf(stderr, "WARNING: Failed to create hash table file - it has not been written\n");
		return;
	}

	if (fwrite(hash_table, sizeof(uint32_t), num_cells_hash_table, out_file) != num_cells_hash_table) {
		fprintf(stderr, "WARNING: Did not manage to write hash table to file\n");
		return;
	}

	if (fclose(out_file) != 0) {
		if (!quiet) {
			fprintf(stderr, "WARNING: Failed to close hash table file - continuing anyway\n");
		}
	}

	else if (!quiet) {
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


seq_hash_return hash_sequence(char *seq, unsigned int region_size, unsigned int interval_size, unsigned int window_size, int kmer_size) {

	uint64_t seq_hash = 0;
	int base_hash;
	seq_hash_return to_return;

	to_return.found_n = false;

	for (unsigned int base_index = 0; base_index < window_size; base_index++) {
		base_hash = hash_base(seq[base_index]);
		if (base_hash != -1) {
			if ((base_index % (region_size + interval_size)) < region_size) {
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


uint64_t hash_rc(uint64_t seq_hash, int kmer_size) {

	uint64_t mask = -4; /* All bits should be set to 1 except least significant two */
	uint64_t rc_hash = 0;

	for (int i = 0; i < (kmer_size - 1); i++) {
		rc_hash += ~((seq_hash & 3) | mask);
		rc_hash <<= 2;
		seq_hash >>= 2;
	}
	rc_hash += ~((seq_hash & 3) | mask);


	return rc_hash;
}


void decode_hash(uint64_t hash, unsigned int region_size, unsigned int window_size, unsigned int interval_size, int kmer_size) {

	int shift = 2 * (kmer_size - 1);

	for (unsigned int base_index = 0; base_index < window_size; base_index++) {
		if ((base_index % (region_size + interval_size)) < region_size) {
			if ((hash & (3ULL << shift)) == 0) {
				putc('A', stderr);
			}
			else if ((hash & (3ULL << shift)) == (1ULL << shift)) {
				putc('C', stderr);
			}
			else if ((hash & (3ULL << shift)) == (2ULL << shift)) {
				putc('G', stderr);
			}
			else if ((hash & (3ULL << shift)) == (3ULL << shift)) {
				putc('T', stderr);
			}
			else {
				fprintf(stderr, "INTERNAL ERROR: Something went wrong in decode_hash\n");
			}
			hash <<= 2;
		}
		else {
			putc('-', stderr);
		}
	}

	return;
}


void decode_all_hashes(uint64_t hash_val, uint64_t rc_hash, uint64_t canonical_hash, int region_size, int window_size, int interval_size, int kmer_size, uint64_t hash_to_use, uint32_t *hash_table) {
	fprintf(stderr, "Forward hash: ");
	decode_hash(hash_val, region_size, window_size, interval_size, kmer_size);
	fprintf(stderr, "\tReverse complement hash: ");
	decode_hash(rc_hash, region_size, window_size, interval_size, kmer_size);
	fprintf(stderr, "\tCanonical hash: ");
	decode_hash(canonical_hash, region_size, window_size, interval_size, kmer_size);
	fprintf(stderr, "\tValue of used hash in table: %" PRIu32, hash_table[hash_to_use]);


	return;
}


void print_usage(char *prog_loc) {

	char *diagram;

	fprintf(stderr, "\nusage: %s mode [options] <file>\n\n"

					"modes:\n"
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
							"\t\t-u, --cutoff : minimum number of k-mers mapped to read for read to be printed (50)\n"
							"\t\t-d, --disable-mask : leave bases not occurring in desired k-mer peaks unmasked when extracting reads (faslse)\n\n"

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


new_hashes shift_hash(uint64_t current_seq_hash, uint64_t current_rc_hash, unsigned int num_regions, int *base_hash_array, int kmer_size) {

	new_hashes to_return;
	uint64_t seq_mask; 
	uint64_t rc_mask;
	int jump = (2 * kmer_size) / num_regions; /* Distance to next region */

	if (kmer_size == 13) {
		seq_mask = 67108860;
		rc_mask = 16777215;
	}

	else if (kmer_size == 15) {
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
	}

	else if (kmer_size == 17) {
		seq_mask = 17179869180; /* 0011 1111 1111 1111 1111 1111 1111 1111 1100 */
		rc_mask = 4294967295; /* 0000 1111 1111 1111 1111 1111 1111 1111 1111 */
	}

	else {
		exit(1);
	}

	current_seq_hash <<= 2;
	current_rc_hash >>= 2;

	/* Zero two least significant bits of each region, and two most significant bits overall, of 32-bit int */
	current_seq_hash &= seq_mask;
	/* Zero two most significant bits of each region, and two most significant bits overall, of 32-bit int */
	current_rc_hash &= rc_mask;

	to_return.new_hash = current_seq_hash;
	to_return.new_rc_hash = current_rc_hash;

	for (int i = 0; i < num_regions; i++) {
		to_return.new_hash += base_hash_array[i] << (jump * (num_regions - i - 1));
		to_return.new_rc_hash += (((uint64_t) (base_hash_array[num_regions - i - 1] ^ 3)) << (((2 * kmer_size) - 2) - (i * jump)));
	}

	to_return.canonical_hash = (to_return.new_hash < to_return.new_rc_hash) ? to_return.new_hash : to_return.new_rc_hash;


	return to_return;
}


new_hashes hash_new_window(uint64_t current_seq_hash, int kmer_size) {

	new_hashes to_return;

	to_return.new_hash = current_seq_hash;
	to_return.new_rc_hash = hash_rc(current_seq_hash, kmer_size);

	to_return.canonical_hash = (to_return.new_hash < to_return.new_rc_hash) ? to_return.new_hash : to_return.new_rc_hash;


	return to_return;
}


int main(int argc, char **argv) {

	FILE *input_file;
	seg_return ret; /* Returned data from get_next_seg (i.e. read info and bEOF) */
	int format;
	uint32_t *hash_table;
	uint64_t hash_val; 
	uint64_t rc_hash;
	uint64_t canonical_hash;
	uint64_t hash_to_use;
	int hash;
	seq_hash_return hash_seq;
	new_hashes new_hashes_triple;
	unsigned int kmer_hits = 0;
	int min_val = 0;
	int max_val = 0;
	int cutoff = 50;
	int interval_size = -1;	/* Number of bases between regions      |||				|||            |||            |||            |||	*/
	int region_size = -1;	/* Number of bases in each region       ----------------------------------------------------------------	*/
	int window_size;		/* Number of bases in window             3      10       3       10     3      10      3      10      3		*/
	int num_regions;		/*										         ^-- interval           ^-- region							*/
							/*										<--------------------------- window --------------------------->	*/
	char *where_to_save_hash_table = NULL;
	char *stored_hash_table_location = NULL;
	uint64_t num_cells_hash_table;
	uint64_t base_index; 
	unsigned long new_base_loc;
	int histogram_size = 10001;
	long hist[histogram_size];
	bool extract_reads = false;
	bool print_hist = false;
	bool quiet = false;
	bool verbose = false;
	bool use_canonical = false;
	int mask = 2; /* 0 = no masking; 1 = strict mask; 2 = normal mask */
	unsigned long end_newest_kmer = 0; /* Index of the end of the most recently found k-mer word in the desired range. Set to 0 to avoid the first base being unmasked. */
	unsigned long *final_indices; /* Array holding the indices of the final base currently masked for each set of bases modulo (region_size + interval_size) */
	unsigned long final_index;
	int new_base_hash_array[5];
	unsigned int iCount; 
	int phase; /* 0 => Pass 1; 1 => Pass 2 */
	int kmer_size = 0; 
	long read_count = 0;
	long read_count_cutoff = 500000;

	phase = 999; /* Default phase = 999, meaning that I have forgotten to set it to anything meaningful */

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

	for (int arg_i = 2; arg_i < argc - 1; arg_i++) {

		if (!strcmp(argv[arg_i], "-u") || !strcmp(argv[arg_i], "--cutoff")) {
			if (!extract_reads) {
				fprintf(stderr, "ERROR: -u/--cutoff must not be specified in this mode\n");
				exit(EXIT_FAILURE);
			}
			if (is_str_integer(argv[++arg_i])) {
				cutoff = atoi(argv[arg_i]);
				if (cutoff < 0) {
					fprintf(stderr, "ERROR: -u/--cutoff must not be negative\n");
					exit(EXIT_FAILURE);
				}
			}
			else {
				print_usage(argv[0]);
			}
		}

		else if (!strcmp(argv[arg_i], "-a") || !strcmp(argv[arg_i], "--min")) {
			if (!extract_reads) {
				fprintf(stderr, "ERROR: -a/--min-val must not be specified in this mode\n");
				exit(EXIT_FAILURE);
			}
			if (is_str_integer(argv[++arg_i])) {
				min_val = atoi(argv[arg_i]);
				if (min_val < 0) {
					fprintf(stderr, "ERROR: -a/--min-val must not be negative\n");
					exit(EXIT_FAILURE);
				}
			}
			else {
				print_usage(argv[0]);
			}
		}

		else if (!strcmp(argv[arg_i], "-b") || !strcmp(argv[arg_i], "--max")) {
			if (!extract_reads) {
				fprintf(stderr, "ERROR: -b/--max-val must not be specified in this mode\n");
				exit(EXIT_FAILURE);
			}
			if (is_str_integer(argv[++arg_i])) {
				max_val = atoi(argv[arg_i]);
				if (max_val < 0) {
					fprintf(stderr, "ERROR: -b/--max-val must not be negative\n");
					exit(EXIT_FAILURE);
				}
			}
			else {
				print_usage(argv[0]);
			}
		}

		else if (!strcmp(argv[arg_i], "-q") || !strcmp(argv[arg_i], "--quiet")) {
			quiet = true;
		}

		else if (!strcmp(argv[arg_i], "-v") || !strcmp(argv[arg_i], "--verbose")) {
			verbose = true;
		}

		else if (!strcmp(argv[arg_i], "-c") || !strcmp(argv[arg_i], "--canonical")) {
			use_canonical = true;
		}

		else if (!strcmp(argv[arg_i], "-k") || !strcmp(argv[arg_i], "--kmer-size")) {
			if (kmer_size != 0) {
				fprintf(stderr, "ERROR: -k/--kmer-size specified more than once\n");
				exit(EXIT_FAILURE);
			}
			if (is_str_integer(argv[++arg_i])) {
				kmer_size = atoi(argv[arg_i]);
				if (kmer_size != 13 && kmer_size != 15 && kmer_size != 17) {
					fprintf(stderr, "ERROR: -k/--kmer-size must be 13, 15, or 17\n");
					exit(EXIT_FAILURE);
				}
			}
			else {
				print_usage(argv[0]);
			}
		}

		else if (!strcmp(argv[arg_i], "-d") || !strcmp(argv[arg_i], "--disable-mask")) {
			if (!extract_reads) {
				fprintf(stderr, "ERROR: Mask must not be specified in this mode\n");
				exit(EXIT_FAILURE);
			}
			if (mask == 1) {
				fprintf(stderr, "ERROR: --disable-mask cannot be used with --strict-mask\n");
				print_usage(argv[0]);
			}
			else {
				mask = 0;
			}
		}
		else if (!strcmp(argv[arg_i], "-s") || !strcmp(argv[arg_i], "--strict-mask")) {
			if (!extract_reads) {
				fprintf(stderr, "ERROR: Mask must not be specified in this mode\n");
				exit(EXIT_FAILURE);
			}
			if (mask == 0) {
				fprintf(stderr, "ERROR: --strict-mask cannot be used with --disable-mask\n");
				print_usage(argv[0]);
			}
			else {
				mask = 1;
			}
		}

		else if (!strcmp(argv[arg_i], "-i") || !strcmp(argv[arg_i], "--in")) {
			stored_hash_table_location = argv[++arg_i];
		}

		else if (!strcmp(argv[arg_i], "-o") || !strcmp(argv[arg_i], "--out")) {
			where_to_save_hash_table = argv[++arg_i];
		}

		else if (!strcmp(argv[arg_i], "-r") || !strcmp(argv[arg_i], "--region-size")) {
			if (is_str_integer(argv[++arg_i])) {
				region_size = atoi(argv[arg_i]);
				if (region_size < 0) {
					fprintf(stderr, "ERROR: -r/--region-size cannot be negative\n");
					exit(EXIT_FAILURE);
				}
			}
			else {
				print_usage(argv[0]);
			}
		}

		else if (!strcmp(argv[arg_i], "-g") || !strcmp(argv[arg_i], "--interval-size")) {
			if (is_str_integer(argv[++arg_i])) {
				interval_size = atoi(argv[arg_i]);
				if (interval_size < 0) {
					fprintf(stderr, "ERROR: -g/--interval-size cannot be negative\n");
					exit(EXIT_FAILURE);
				}
			}
			else {
				print_usage(argv[0]);
			}
		}

		else {
			print_usage(argv[0]);
		}
	}


	/* ----- Error check user input ----- */


	if (kmer_size == 0) {
		fprintf(stderr, "ERROR: -k/--kmer-size must be specified\n");
		exit(EXIT_FAILURE);
	}

	if (extract_reads) {
		if (min_val == 0 || max_val == 0) {
			fprintf(stderr, "ERROR: -a/--min-val and -b/--max-val must both be specified and non-zero when extracting reads\n");
			exit(EXIT_FAILURE);
		}

		if (max_val < min_val) {
			fprintf(stderr, "ERROR: -a/--min-val must not be greater than -b/--max-val\n");
			exit(EXIT_FAILURE);
		}
	}

	if (kmer_size == 15) {
		if (region_size != -1) {
			if (region_size != 1 && region_size != 3 && region_size != 5 && region_size != 15) {
				fprintf(stderr, "ERROR: Region size must be either 1, 3, 5, or 15\n");
				exit(EXIT_FAILURE);
			}
		}
	}

	else {
		if (region_size != -1) {
			fprintf(stderr, "ERROR: -r/--region-size must not be specified if -k/--kmer-size is not 15\n");
			exit(EXIT_FAILURE);
		}

		if (interval_size != -1) {
			fprintf(stderr, "ERROR: -g/--interval-size must not be specified if -k/--kmer-size is not 15\n");
			exit(EXIT_FAILURE);
		}

		if (mask == 1) {
			fprintf(stderr, "ERROR: -s/--strict cannot be used if -k/--kmer-size is not 15\n");
			exit(EXIT_FAILURE);
		}
	}

	if (stored_hash_table_location && where_to_save_hash_table) {
		fprintf(stderr, "ERROR: Cannot specify both -i/--in and -o/--out\n");
		exit(EXIT_FAILURE);
	}

	if (quiet && verbose) {
		fprintf(stderr, "ERROR: Cannot enable both -q/--quiet and -v/--verbose modes\n");
		exit(EXIT_FAILURE);
	}


	/* ----- End of user input error checking ----- */


	num_cells_hash_table = 1UL << (2 * kmer_size); /* = 4^kmer_size */

	if (region_size == -1) {
		region_size = kmer_size;
	}

	if (interval_size == -1) {
		interval_size = 0;
	}

	num_regions = (kmer_size / region_size);
	window_size = ((num_regions - 1) * interval_size) + kmer_size;
	
	if ((input_file = fopen(argv[argc - 1], "r")) == NULL) {
		fprintf(stderr, "ERROR: Could not open data file\n");
		exit(EXIT_FAILURE);
	}

	if (verbose) {
		fprintf(stderr, "Region size = %u Number regions = %d Interval size = %u Window size = %u K-mer size = %d\n", region_size, num_regions, interval_size, window_size, kmer_size);
	}

	/* Hash table is 4**kmer_size cells which are guaranteed to be capable of holding the count of a k-mer 
	 * providing that the count does not exceed 2^32 - the minimum size of a long). 
	 */

	if ((hash_table = calloc(num_cells_hash_table, sizeof(uint32_t))) == NULL) {
		fprintf(stderr, "ERROR: Out of memory\n"); exit(EXIT_FAILURE);
	}

	format = which_format(input_file);

	if (stored_hash_table_location == NULL) {
		if (!quiet) {
			fprintf(stderr, "Counting k-mers into hash table\n");
		}
		phase = 0;
	}
	else {
		read_hash_table_from_file(hash_table, stored_hash_table_location, quiet, num_cells_hash_table);

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

				if (mask == 1) {
					if ((final_indices = calloc(region_size + interval_size, sizeof(unsigned long))) == NULL) {
						fprintf(stderr, "ERROR: Ran out of memory\n");
						exit(EXIT_FAILURE);
					}
				}

				rewind(input_file);
			}

			if (!quiet) {
				fprintf(stderr, "One dot for each 500,000 reads processed\n");
			}

			read_count = 0;

			do {
				base_index = 0;

				if (phase == 2) {
					kmer_hits = 0;
					if (mask == 2) {
						end_newest_kmer = 0; 
					}
					else if (mask == 1) {
						for (int k = 0; k < region_size + interval_size; k++) {
							final_indices[k] = 0;
						}
					}
				}

				ret = get_next_seg(input_file, format);
				read_count++;

				if (ret.segment.length < window_size) {

					free(ret.segment.name);
					free(ret.segment.seq);
					if (format == 1) {
						free(ret.segment.qual);
					}
					if (!quiet) {
						if (read_count == read_count_cutoff) {
							read_count = 0;
							fprintf(stderr, ".");
						}
					}

					continue;
				} 

				if (verbose) {
					fprintf(stderr, "Read name: %s\n", ret.segment.name);
				}

				hash_seq = hash_sequence(ret.segment.seq, region_size, interval_size, window_size, kmer_size);

				while (hash_seq.found_n == true && base_index <= (ret.segment.length - window_size)) {
					base_index += 1;
					hash_seq = hash_sequence(ret.segment.seq + base_index, region_size, interval_size, window_size, kmer_size);
					if (phase == 2) {
						if (mask) {
							if (verbose) {
								fprintf(stderr, "(1) Masking at base_index = %lu\n", base_index);
							}
							ret.segment.seq[base_index] = 'N';
						}
					}
				}

				if (hash_seq.found_n == false) {

					new_hashes_triple = hash_new_window(hash_seq.hash, kmer_size);
					hash_val = new_hashes_triple.new_hash;
					rc_hash = new_hashes_triple.new_rc_hash;
					canonical_hash = new_hashes_triple.canonical_hash;

					hash_to_use = use_canonical ? canonical_hash : hash_val;

					base_index += window_size - 1; 

					if (verbose) {
						decode_all_hashes(hash_val, rc_hash, canonical_hash, region_size, window_size, interval_size, kmer_size, hash_to_use, hash_table);
						fprintf(stderr, " [1]\n");
					}

					if (phase == 0) {
						hash_table[hash_to_use] += 1;
					}

					else if (phase == 2) {
						if (hash_table[hash_to_use] >= min_val && hash_table[hash_to_use] <= max_val) {

							if (mask == 1) {
								for (int k = base_index - region_size + 1, l = 0; l < (region_size); l++) {
									if (verbose) {
										fprintf(stderr, "1: (k + l) %% (region_size + interval_size) = %" PRIu32 " k+l = %d\n", (k + l) % (region_size + interval_size), k+l);
									}
									final_indices[(k + l) % (region_size + interval_size)] = k + l;	
								}
							}

							else if (mask == 2) {
								end_newest_kmer = base_index;
							}

							kmer_hits++;
						}

						if (mask) {
							if (verbose) {
								fprintf(stderr, "(2) - Masking at base_index = %lu\n", base_index);
							}
							if (mask == 1) {
								final_index = final_indices[base_index - window_size + 1 % (region_size + interval_size)];
								if (verbose) {
									fprintf(stderr, "2: Final index = %lu\n", final_index);
								}
								if ((final_index == 0) || ((base_index - window_size + 1) > final_index)) {
									ret.segment.seq[base_index - window_size + 1] = 'N';
								}
							}
							else if (mask == 2) {
								if ((end_newest_kmer == 0) || ((base_index - window_size + 1) > end_newest_kmer)) {
									ret.segment.seq[base_index - window_size + 1] = 'N';
								}
							}
						}
					}

					for (base_index += 1; base_index < ret.segment.length; base_index++) {

						for (iCount = 0; iCount < num_regions - 1; iCount++) {
							/* Can guarantee that only the final new character hashed might be an 'N', as otherwise we would have already found it */
							new_base_loc = base_index - window_size + region_size + (iCount * (region_size + interval_size));
							hash = hash_base(ret.segment.seq[new_base_loc]);
							new_base_hash_array[iCount] = hash;
						}

						new_base_loc = base_index - window_size + region_size + (iCount * (region_size + interval_size));
						hash = hash_base(ret.segment.seq[new_base_loc]);

						if (hash != -1) {
							new_base_hash_array[iCount] = hash;
							new_hashes_triple = shift_hash(hash_val, rc_hash, num_regions, new_base_hash_array, kmer_size);
							hash_val = new_hashes_triple.new_hash;
							rc_hash = new_hashes_triple.new_rc_hash;
							canonical_hash = new_hashes_triple.canonical_hash;

							hash_to_use = use_canonical ? canonical_hash : hash_val;

							if (verbose) {
								decode_all_hashes(hash_val, rc_hash, canonical_hash, region_size, window_size, interval_size, kmer_size, hash_to_use, hash_table);
								fprintf(stderr, " [2]\n");
							}

							if (phase == 0) {
								hash_table[hash_to_use] += 1;
							}

							else if (phase == 2) {
								if (hash_table[hash_to_use] >= min_val && hash_table[hash_to_use] <= max_val) {
									if (mask == 1) {
										for (int k = base_index - region_size + 1, l = 0; l < (region_size); l++) {
											if (verbose) {
												fprintf(stderr, "2: (k + l) %% (region_size + interval_size) = %" PRIu32 " k+l = %d\n", (k + l) % (region_size + interval_size), k+l);
											}
											final_indices[(k + l) % (region_size + interval_size)] = k + l;	
										}
									}

									else if (mask == 2) {
										end_newest_kmer = base_index;
									}

									kmer_hits++;
								}

								if (mask) {
									if (verbose) {
										fprintf(stderr, "(3) Masking at base_index = %lu\n", base_index);
									}
									if (mask == 1) {
										final_index = final_indices[(base_index - window_size + 1) % (region_size + interval_size)];
										if (verbose) {
											fprintf(stderr, "3: base_index - window_size + 1 %% (region_size + interval_size) = %lu Final index = %lu\n", (base_index - window_size + 1) % (region_size + interval_size), final_index);
										}
										if ((final_index == 0) || ((base_index - window_size + 1) > final_index)) {
											ret.segment.seq[base_index - window_size + 1] = 'N';
										}
									}
									else if (mask == 2) {
										if ((end_newest_kmer == 0) || ((base_index - window_size + 1) > end_newest_kmer)) {
											ret.segment.seq[base_index - window_size + 1] = 'N';
										}
									}
								}
							}
						}

						else {

							if (phase == 2) {
								if (mask) {
									/* Before moving onto the next k-mer word, mask, if necessary, the remainder of the k-mer word which is going to be skipped */
									if (verbose) {
										fprintf(stderr, "(4) Masking at base_index = %lu\n", base_index);
									}
									if (mask == 1) {
										for (uint64_t i = base_index - window_size + 1; i < base_index; i++) {
											final_index = final_indices[i % (region_size + interval_size)];
											if (verbose) {
												fprintf(stderr, "4: Final index = %lu\n", final_index);
											}
											if ((final_index == 0) || (i > final_index)) {
												ret.segment.seq[i] = 'N';
											}
										}
									}
									else if (mask == 2) {
										for (uint64_t i = base_index - window_size + 1; i < base_index; i++) {
											if ((end_newest_kmer == 0) || (i > end_newest_kmer)) {
												ret.segment.seq[i] = 'N';
											}
										}
									}
								}
							}

							base_index += 1;
							hash_seq = hash_sequence(ret.segment.seq + base_index, region_size, interval_size, window_size, kmer_size);

							/* Keep hashing the sequence starting at the next base and moving along the window until we don't find any more 'N's */
							while (hash_seq.found_n == true && base_index < (ret.segment.length - window_size)) {
								if (phase == 2) {
									if (mask) {
										if (verbose) {
											fprintf(stderr, "(5) Masking at base_index = %lu\n", base_index);
										}
										if (mask == 1) {
											final_index = final_indices[base_index - window_size + 1 % (region_size + interval_size)];
											if (verbose) {
												fprintf(stderr, "5: Final index = %lu\n", final_index);
											}
											if ((final_index == 0) || (base_index > final_index)) {
												ret.segment.seq[base_index] = 'N';
											}
										}
										else if (mask == 2) {
											if ((end_newest_kmer == 0) || (base_index > end_newest_kmer)) {
												ret.segment.seq[base_index] = 'N';
											}
										}
									}
								}

								base_index += 1;
								hash_seq = hash_sequence(ret.segment.seq + base_index, region_size, interval_size, window_size, kmer_size);
							}

							if (hash_seq.found_n == true) {
								break;
							}

							else {
								new_hashes_triple = hash_new_window(hash_seq.hash, kmer_size);
								hash_val = new_hashes_triple.new_hash;
								rc_hash = new_hashes_triple.new_rc_hash;
								canonical_hash = new_hashes_triple.canonical_hash;

								hash_to_use = use_canonical ? canonical_hash : hash_val;

								/* Move to end of k-mer word */
								base_index += window_size - 1;

								if (verbose) {
									decode_all_hashes(hash_val, rc_hash, canonical_hash, region_size, window_size, interval_size, kmer_size, hash_to_use, hash_table);
									fprintf(stderr, " [3]\n");
								}

								if (phase == 0) {
									hash_table[hash_to_use] += 1;
								}

								else if (phase == 2) {
									if (hash_table[hash_to_use] >= min_val && hash_table[hash_to_use] <= max_val) {
										if (mask == 1) {
											for (int k = base_index - region_size + 1, l = 0; l < (region_size); l++) {
												if (verbose) {
													fprintf(stderr, "3: (k + l) %% (region_size + interval_size) = %" PRIu32 " k+l = %d\n", (k + l) % (region_size + interval_size), k+l);
												}
												final_indices[(k + l) % (region_size + interval_size)] = k + l;	
											}
										}

										else if (mask == 2) {
											end_newest_kmer = base_index;
										}

										kmer_hits++;
									}

									if (mask) {
										if (verbose) {
											fprintf(stderr, "(6) Masking at base_index = %lu\n", base_index);
										}
										if (mask == 1) {
											final_index = final_indices[base_index - window_size + 1 % (region_size + interval_size)];
											if (verbose) {
												fprintf(stderr, "6: Final index = %lu\n", final_index);
											}
											if ((final_index == 0) || ((base_index - window_size + 1) > final_index)) {
												ret.segment.seq[base_index - window_size + 1] = 'N';
											}
										}
										else if (mask == 2) {
											if ((end_newest_kmer == 0) || ((base_index - window_size + 1) > end_newest_kmer)) {
												ret.segment.seq[base_index - window_size + 1] = 'N';
											}
										}
									}
								}
							}
						}
					}

					if (phase == 2) {
						if (mask) {
							/* Mask the bases in the final k-mer word */
							if (verbose) {
								fprintf(stderr, "(7) Masking at base_index = %lu\n", base_index);
							}
							for (uint64_t i = base_index - window_size; i < base_index; i++) {
								if (mask == 1) {
									final_index = final_indices[i % (region_size + interval_size)];
									if (verbose) {
										fprintf(stderr, "7: Final index = %lu\n", final_index);
									}
									if ((final_index == 0) || (i > final_index)) {
										ret.segment.seq[i] = 'N';
									}
								}
								else if (mask == 2) {
									if ((end_newest_kmer == 0) || (i > end_newest_kmer)) {
										ret.segment.seq[i] = 'N';
									}
								}
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

				if (!quiet) {
					if (read_count == read_count_cutoff) {
						read_count = 0;
						fprintf(stderr, ".");
					}
				}

			} while (!ret.bEOF);

			if (!quiet) {
				fprintf(stderr, "\n");
			}

			if (phase == 2) {
				if (mask == 1) {
					free(final_indices);
				}
			}
		}

		if (phase == 0) {
			if (where_to_save_hash_table) {
				write_hash_table_to_file(hash_table, where_to_save_hash_table, quiet, num_cells_hash_table);
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

			for (uint64_t i = 0; i < histogram_size; i++) {
				hist[i] = 0;
			}

			for (uint64_t i = 0; i < num_cells_hash_table; i++) {
				if (hash_table[i] > 0) {
					if (hash_table[i] < histogram_size) {
						hist[hash_table[i] - 1]++;
					}
					else {
						hist[histogram_size - 1]++;
					}
				}
			}

			for (uint64_t i = 0; i < histogram_size; i++) {
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
