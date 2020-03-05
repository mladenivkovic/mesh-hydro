include defines.mk
include ../bin/processing.mk

# executable name
EXEC = hydro  


OBJECTS = main.o gas.o params.o io.o utils.o cell.o $(HYDROOBJ) $(LIMITEROBJ) $(RIEMANNOBJ)




#==================================
# COMPILER OPTIONS AND DIRECTIVES
#==================================


#Compiler to use
CC=gcc -I$(IDIR)


# flags for gcc compiler
CFLAGS= -g -O0 -pedantic -std=c11 \
		-fsignaling-nans -ffinite-math-only -fno-signed-zeros \
		-Wextra -Werror -Warray-bounds -Wno-unused-parameter
# CFLAGS= -O3

CFLAGS += $(DEFINES)

LDFLAGS= -lm








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
