typedef struct {
	uint64_t hash;
	bool found_n; /* Set to 1 if we found an 'N' (i.e. this k-mer is not all 'A's, but rather had an 'N' in it) */
} seq_hash_return;

typedef struct {
	uint64_t new_hash;
	uint64_t new_rc_hash;
	uint64_t canonical_hash;
} new_hashes;

int hash_base(char base);
seq_hash_return hash_sequence(char *seq, unsigned int region_size, unsigned int interval_size, unsigned int window_size);
new_hashes hash_new_window(uint64_t current_seq_hash, int kmer_size);
new_hashes shift_hash(uint64_t current_seq_hash, uint64_t current_rc_hash, int num_regions, int *base_hash_array, int kmer_size);
uint64_t hash_rc(uint64_t seq_hash, int kmer_size);
void decode_hash(uint64_t hash, int region_size, int window_size, int interval_size, int kmer_size);
void read_hash_table_from_file(uint32_t *hash_table, char *hash_table_location, bool quiet, uint64_t num_cells_hash_table);
void write_hash_table_to_file(uint32_t *hash_table, char *hash_file_name, bool quiet, uint64_t num_cells_hash_table);
void compute_histogram(long *hist, bool quiet, unsigned int histogram_size, uint32_t *hash_table, uint64_t num_cells_hash_table);
