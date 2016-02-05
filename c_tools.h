#ifndef C_TOOLS_H

#define C_TOOLS_H

typedef struct {
	char* line;
	bool bEOF;
} line_return;

line_return get_next_line(FILE* f);
bool is_str_integer(char* str);

#endif
