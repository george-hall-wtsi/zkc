CFLAGS = -Wall -Wextra -O3 
CC = cc

SRCS = zkc2.c c_tools.c fastlib.c
OBJS = $(SRCS:.c=.o)
	
zkc2-test: $(OBJS)
	$(CC) $(CFLAGS) -o zkc2-test $(OBJS)

# Production version
zkc2: $(OBJS)
	$(CC) $(CFLAGS) -o zkc2 $(OBJS)

# Aliases
val: CFLAGS = -O0 -g -std=c99
val: zkc2-test
prod: zkc2
	
clean:
	rm -f *.o

%.o : %.c 
	$(CC) $(CFLAGS) -c $< -o $@
