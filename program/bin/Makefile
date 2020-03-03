#==================================
# COMPILATION TIME PARAMETERS
#==================================

# executable name
EXEC = hydro  

# how many dimensions to compile for
# Choices: 1, 2
NDIM = 1
# NDIM = 2


# set which Riemann solver to use.
# Choices: EXACT, TRRS, TSRS, HLLC
# RIEMANN = NONE
RIEMANN = EXACT
# RIEMANN = TRRS
# RIEMANN = TSRS
# RIEMANN = HLLC


# set which hydro solver to use.
# Choices: GODUNOV, ADVECTION_PWCONST, ADVECTION_PWLIN
SOLVER = GODUNOV
# SOLVER = ADVECTION_PWCONST
# SOLVER = ADVECTION_PWLIN

# set which flux limiter to use.
# Choices: NONE, MINMOD, SUPERBEE, VANLEER, MC
LIMITER = NONE
# LIMITER = MINMOD
# LIMITER = SUPERBEE
# LIMITER = MC
# LIMITER = VANLEER



# if there is a defines.mk present, overwrite the
# definitions above with its contents
ifneq ("$(wildcard defines.mk)","")
include defines.mk
endif


# process the definitions
include processing.mk






#==================================
# COMPILER OPTIONS AND DIRECTIVES
#==================================


#Compiler to use
CC=gcc -I$(IDIR) -static


# flags for gcc compiler
CFLAGS= -g -O0 -pedantic -std=c11 \
		-fsignaling-nans -ffinite-math-only -fno-signed-zeros \
		-Wextra -Werror -Warray-bounds #-Wno-unused-parameter
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
	rm -f *.o
