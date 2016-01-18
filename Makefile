CFLAGS = -O3 -ansi -Wall -Wextra -pedantic 
CC = cc

SRCS = zkc2.c c_tools.c fastlib.c
OBJS = $(SRCS:.c=.o)

.PHONY: clean

all: $(OBJS) 
	$(CC) $(CFLAGS) -o zkc2 $(OBJS)

# Options to be used when building with Valgrind
val : CFLAGS = -O0 -g
val: all
	
clean:
	rm -f *.o
