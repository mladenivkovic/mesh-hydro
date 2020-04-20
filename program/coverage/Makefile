include defines.mk
include ../bin/processing.mk

# executable name
EXEC = hydro  





#==================================
# COMPILER OPTIONS AND DIRECTIVES
#==================================


#Compiler to use
CC=gcc -I$(IDIR)


# flags for gcc compiler
CFLAGS=  -Wall -fprofile-arcs -ftest-coverage -fPIC -O0

CFLAGS += $(DEFINES)

LDFLAGS= -lm -lgcov








#==================================
# RECEPIES
#==================================

default: $(EXEC)

$(EXEC): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $<






#==================================
#==================================

.PHONY: clean # keeps makefile from doing something with a file named clean

clean:
	rm -f $(EXEC) *.o
