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
#include <stdbool.h>
#include <ctype.h>

#include "c_tools.h"


line_return get_next_line(FILE* f) {

	int c;
	char* line;
	size_t buffsize = 5000;
	unsigned long i = 0;

	line_return to_return;
	to_return.bEOF = false; /* Set to true if EOF detected	*/

	if ((line = malloc(buffsize)) == NULL) {
		fprintf(stderr, "ERROR: Out of memory\n");
		exit(EXIT_FAILURE);
	}

	while ((c = getc(f)) != '\n') {

		if (c == EOF) {
			to_return.bEOF = true;
			break;
		}

		if ((i + 1) == buffsize) {
			char* tmp = realloc(line, buffsize *= 2);

			if (tmp == NULL) {
				fprintf(stderr, "ERROR: Out of memory\n");
				free(line); 
				exit(EXIT_FAILURE);
			}

			line = tmp;
		}

		line[i++] = (char) c;
	}

	if ((c = getc(f)) == EOF) {
		to_return.bEOF = true;
	}
	else {
		if ((ungetc(c, f)) == EOF) {
			fprintf(stderr, "ERROR: ungetc failed in c_tools.c get_next_line function\n");
			exit(EXIT_FAILURE);
		}
	}

	line[i] = '\0';

	to_return.line = line;

	return to_return;
}


bool is_str_integer(char *str) {
	
	/* Return true if char array 'str' represents an integer */

	if (*str == '-') {
		if (*(++str) == '\0') {
			return false;
		}
	}

	while (*str != '\0') {
		if (!isdigit(*(str++))) {
			return false;
		}
	}

	return true;
}
