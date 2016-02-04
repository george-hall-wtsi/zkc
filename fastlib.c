/*******************************************************************************
 * Copyright (c) 2015 Genome Research Ltd. 
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

#include "c_tools.h"
#include "fastlib.h"


float calc_gc(char *seq) {

	size_t len = strlen(seq);

	int G_or_C = 0, 
		A_or_T = 0;

	if (*seq == '\0') {
		fprintf(stderr, "Length of sequence must be nonzero\n");
		exit(EXIT_FAILURE);
	}

	do {
		if (*seq == 'C' || *seq == 'G') {
			G_or_C++;
		}
		else if (*seq == 'A' || *seq == 'T') {
			A_or_T++;
		}
	} while (*++seq != '\0');

	return (100 * ((float) G_or_C / len));
}


void fastq_to_fasta(FILE *f) {

	seg_return new_segment;
	segment seg;

	do {
		new_segment = get_next_seg(f, 1);
		seg = new_segment.segment;

		printf(">%s\n%s\n", seg.name, seg.seq);

		free(seg.name);
		free(seg.seq);
		free(seg.qual);

	} while (!new_segment.bEOF);
}


void rename_reads(FILE *f, char *name, unsigned long min_length) {

	int format = which_format(f);
	int iCount = 0;
	bool bEOF = false;
	int name_len = strlen(name);
	char new_name[name_len + 12];
	char numbers[12];
	
	while (!bEOF) {
		seg_return new_segment = get_next_seg(f, format);
		segment seg = new_segment.segment;
		bEOF = new_segment.bEOF;

		if (seg.length >= min_length) {
			strcpy(new_name, name);
			sprintf(numbers, "%011d", iCount);
			strcat(new_name, numbers);
			new_name[name_len + 11] = '\0';

			if (format == 0) {
				printf(">%s\n%s\n", new_name, seg.seq); 
			}
			else if (format == 1) {
				printf("@%s\n%s\n+\n%s\n", new_name, seg.seq, seg.qual);
			}
			iCount++;
		}

		free(seg.name); 
		free(seg.seq);
		if (format == 1) {
			free(seg.qual);
		}
	}
}


int which_format(FILE *f) {

	/* Returns 0 if fasta, 1 if fastq */

	int format;

	line_return new_line;
	new_line = get_next_line(f);
	if (new_line.bEOF) {
		fprintf(stderr, "File too short!\n");
		exit(EXIT_FAILURE);
	}
	if (new_line.line[0] == '>') {format = 0;}
	else if (new_line.line[0] == '@') {format = 1;}
	else {
		fprintf(stderr, "Formatting error: File must be in fast(a/q) format (file provided does not begin with '>' or '@')\n");
		exit(EXIT_FAILURE);
	}

	free(new_line.line);
	rewind(f);	

	return format;
}

 
seg_return get_next_seg(FILE *f, int format) {

	seg_return to_return;
	segment new_segment;
	char *line, *tmp;
	bool bEOF = false; /* Set to true if EOF detected */
	unsigned long seq_len = 0;
	to_return.bEOF = false;
	unsigned int buffsize = 10000; /* Reads should rarely be longer than this => shouldn't need to realloc very often */

	if ((new_segment.seq = malloc(buffsize + 1)) == NULL) {
		fprintf(stderr, "Out of memory (malloc for sequence)\n");
		exit(EXIT_FAILURE);
	}


	/* Fasta file */
	if (format == 0) {

		line_return first_line = get_next_line(f);
		if (first_line.bEOF) {
			fprintf(stderr, "ERROR: Premature EOF (File ends without sequence)\n");
		}

		if ((new_segment.name = malloc(strlen(first_line.line))) == NULL) {
			fprintf(stderr, "ERROR: Out of memory (malloc for name)\n");
			exit(EXIT_FAILURE);
		}
		strcpy(new_segment.name, first_line.line + 1);
		free(first_line.line);

		new_segment.qual = "\0";

		/* Two state automaton:
		 *		State 0: Header
		 *		State 1: Bases
		 */

		while (true) {
			line_return new_line;
			new_line = get_next_line(f);

			line = new_line.line;
			bEOF = new_line.bEOF;

			if (bEOF) {
				to_return.bEOF = true;
				if (strlen(line) > 0) { /* Check we aren't dealing with a blank line at the end of a file */
					if ((seq_len + strlen(line) > buffsize)) {
						while ((seq_len + strlen(line) > buffsize)) {
							buffsize *= 2;
						}

						if ((tmp = realloc(new_segment.seq, buffsize + 1)) == NULL) {
							fprintf(stderr, "Out of memory (realloc for sequence)\n");
							exit(EXIT_FAILURE);
						}
						new_segment.seq = tmp;
					}

					seq_len += strlen(line);
					strcat(new_segment.seq, line);
				}

				free(line);
				break;
			}

			if (line[0] == '>') {
				/* This is the name of the next segment, so seek back to start of line */
				if (fseek(f, -strlen(line) - 1, SEEK_CUR) != 0) {
					fprintf(stderr, "ERROR: Failed to seek to correct position on file\n");
					exit(EXIT_FAILURE);
				}
				free(line);
				break;
			}

			else {

				if (seq_len == 0) {
					if ((seq_len + strlen(line) > buffsize)) {
						while ((seq_len + strlen(line) > buffsize)) {
							buffsize *= 2;
						}
						if ((tmp = realloc(new_segment.seq, buffsize + 1)) == NULL) {
							fprintf(stderr, "Out of memory (realloc for sequence)\n");
							exit(EXIT_FAILURE);
						}
					}

					strcpy(new_segment.seq, line);
				}

				else {
					if ((seq_len + strlen(line) > buffsize)) {
						while ((seq_len + strlen(line) > buffsize)) {
							buffsize *= 2;
						}

						if ((tmp = realloc(new_segment.seq, buffsize + 1)) == NULL) {
							fprintf(stderr, "Out of memory (realloc for sequence)\n");
							exit(EXIT_FAILURE);
						}
						new_segment.seq = tmp;
					}

					strcat(new_segment.seq, line);
				}

				seq_len += strlen(line);
				free(line);
			}
		}
	}
 
	/* Fastq file */
	else if (format == 1)  {

		/* Two state automaton: 
		 *		State 1: Bases
		 *		State 2: Quality scores
		 */

		int state = 1;
		bool bEOF = false; /* Set to true if EOF detected */
		char c;

		unsigned long qual_len = 0;
		
		line_return new_line;
		new_line = get_next_line(f);

		line = new_line.line;
		bEOF = new_line.bEOF;
		to_return.bEOF = bEOF;

		if ((new_segment.name = malloc(strlen(line) + 1)) == NULL) {
			fprintf(stderr, "ERROR: Out of memory\n");
			exit(EXIT_FAILURE);
		}
		strcpy(new_segment.name, line + 1);

		while (true) {

			free(line);
			new_line = get_next_line(f);
			line = new_line.line;
			bEOF = new_line.bEOF;
			to_return.bEOF = bEOF;


			/* Bases */
			if (state == 1) {
				if (line[0] == '+') {
					/* i.e. this is in fact the header for quality scores */
					state = 2;
					qual_len = 0;
				}

				else {
					if (seq_len == 0) {
						if ((seq_len + strlen(line) > buffsize)) {
							while ((seq_len + strlen(line) > buffsize)) {
								buffsize *= 2;
							}

							if ((tmp = realloc(new_segment.seq, buffsize + 1)) == NULL) {
								fprintf(stderr, "Out of memory (realloc for sequence)\n");
								exit(EXIT_FAILURE);
							}
							new_segment.seq = tmp;
						}

						strcpy(new_segment.seq, line); /* This line is causing a seqfault */
					}
					else {
						if ((seq_len + strlen(line) > buffsize)) {
							while ((seq_len + strlen(line) > buffsize)) {
								buffsize *= 2;
							}

							if ((tmp = realloc(new_segment.seq, buffsize + 1)) == NULL) {
								fprintf(stderr, "Out of memory (realloc for sequence)\n");
								exit(EXIT_FAILURE);
							}
							new_segment.seq = tmp;
						}

						strcat(new_segment.seq, line);
					}

					seq_len += strlen(line);
				}
			} 

			/* Quality Values */ 
			else if (state == 2) { 
				if (qual_len == 0) {
					if ((new_segment.qual = malloc(buffsize + 1)) == NULL) {
						fprintf(stderr, "Out of memory (malloc for qualities)\n");
						exit(EXIT_FAILURE);
					}

					if (strlen(line) > buffsize) {
						while (strlen(line) > buffsize) {
							buffsize *= 2;
						}

						if ((tmp = realloc(new_segment.qual, buffsize + 1)) == NULL) {
							fprintf(stderr, "Out of memory (realloc for qualities)\n");
							exit(EXIT_FAILURE);
						}
						new_segment.qual = tmp;
					}

					strcpy(new_segment.qual, line);
				}

				else if (qual_len < seq_len) {
					if ((qual_len + strlen(line) > buffsize)) {
						while ((qual_len + strlen(line) > buffsize)) {
							buffsize *= 2;
						}

						if ((tmp = realloc(new_segment.qual, buffsize + 1)) == NULL) {
							fprintf(stderr, "Out of memory (realloc for qualities)\n");
							exit(EXIT_FAILURE);
						}
						new_segment.qual = tmp;
					}

					strcat(new_segment.seq, line);
				}	

				qual_len += strlen(line);

				if (qual_len == seq_len) {
					/* End of quality values */
					free(line);

					/* Check for EOF */
					if ((c = getc(f)) != '@') {
						to_return.bEOF = true;
					}
					else {
						if (ungetc(c, f) == EOF) {
							fprintf(stderr, "ERROR: ungetc in get_next_seg failed\n");
							exit(EXIT_FAILURE);
						}
					}

					break;
				}

				else if (qual_len > seq_len) {
					fprintf(stderr, "ERROR: Quality string too long\n");
					exit(EXIT_FAILURE);
				}
			}

			else {
				fprintf(stderr, "ERROR: Automaton should never get here!!\n");
				exit(EXIT_FAILURE);
			}
		}
	}

	else {
		/* Should never reach here - Format needs to equal 0 or 1 */
		fprintf(stderr, "ERROR: Format variable needs to equal 0 or 1, currently equals %d\n", format);
		fprintf(stderr, "This is almost certainly an error in the program, not with the data\n");
		exit(EXIT_FAILURE);
	}

	new_segment.length = seq_len;
	to_return.segment = new_segment;

	return to_return;
}

