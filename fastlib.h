#ifndef FASTLIB_H
#define FASTLIB_H

typedef struct {
	char *name;
	unsigned long length;
	char *seq;
	char *qual; /* Points to "\0" if fasta file */
} segment;

typedef struct {
	segment segment;
	bool bEOF;
} seg_return;

float calc_gc(char *seq);
void fastq_to_fasta(FILE *f);
void rename_reads(FILE *f, char *name, unsigned long min_length);
int which_format(FILE *f);
seg_return get_next_seg(FILE *f, int format);

#endif
