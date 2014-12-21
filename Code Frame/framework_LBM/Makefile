# Include files
SOURCES=initLB.c visualLB.c boundary.c collision.c streaming.c computeCellValues.c main.c helper.c

# Compiler
# --------
CC=gcc

CFLAGS=-Werror -pedantic -Wall

# Linker flags
# ------------
LDFLAGS= 

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=lbsim

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ 

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)


$(OBJECTS): %.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@
