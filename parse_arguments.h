typedef struct argument_struct {
	bool print_hist;
	bool extract_reads;
	int min_kmer_hits;
	int max_kmers_missed;
	unsigned int min_val;
	unsigned int max_val;
	bool quiet;
	bool verbose;
	bool use_canonical;
	int kmer_size;
	int mask; /* 0 = no masking; 1 = strict mask; 2 = normal mask */
	char *where_to_save_hash_table;
	char *stored_hash_table_location;
	int region_size;
	int interval_size;
	int index_first_file;
} argument_struct;
argument_struct parse_arguments(int argc, char **argv);
