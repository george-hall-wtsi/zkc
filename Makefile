CFLAGS = -Wall -Wextra -O3
CC = cc

SRCS = zkc2.c c_tools.c fastlib.c parse_arguments.c
OBJS = $(SRCS:.c=.o)
	
zkc2-test: $(OBJS)
	$(CC) $(CFLAGS) -o zkc2-test $(OBJS)

# Production version
zkc2: $(OBJS)
	$(CC) $(CFLAGS) -o zkc2 $(OBJS)

# Aliases
debug: CFLAGS = -Wall -Wextra -O0 -g
debug: zkc2-test
prod: zkc2
	
clean:
	rm -f *.o

%.o : %.c 
	$(CC) $(CFLAGS) -c $< -o $@
