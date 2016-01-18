typedef struct {
	uint32_t hash;
	bool found_n; /* Set to 1 if we found an 'N' (i.e. this k-mer is not all 'A's, but rather had an 'N' in it) */
} seq_hash_return;

typedef struct {
	uint32_t new_hash;
	uint32_t new_rc_hash;
	uint32_t canonical_hash;
} new_hashes;

void print_usage(char *prog_loc);
int hash_base(char base);
seq_hash_return hash_sequence (char *seq);
new_hashes hash_new_window(uint32_t current_seq_hash);
new_hashes shift_hash(int new_base_hash, uint32_t current_seq_hash, uint32_t current_rc_hash);
uint32_t hash_rc(uint32_t seq_hash);
void decode_hash(uint32_t hash);
void read_hash_table_from_file(uint32_t *hash_table, char *hash_table_location, bool quiet);
void write_hash_table_to_file(uint32_t *hash_table, char *hash_file_name, bool quiet);
