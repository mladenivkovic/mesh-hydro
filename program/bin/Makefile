#==================================
# COMPILATION TIME PARAMETERS
#==================================

# executable name
EXEC = hydro  

# how many dimensions to compile for
# Choices: 1, 2
# NDIM = 1
# NDIM = 2


# set which Riemann solver to use.
# Choices: EXACT, TRRS, TSRS, HLL, HLLC
# RIEMANN = EXACT
# RIEMANN = NONE


# set which hydro solver to use.
# Choices: GODUNOV_UPWIND, ADVECTION_PWCONST, ADVECTION_PWLIN
# SOLVER=GODUNOV_UPWIND
# SOLVER = ADVECTION_PWCONST
SOLVER = ADVECTION_PWLIN

# set which flux limiter to use.
# Choices: NONE, MINMOD, SUPERBEE, VANLEER, MC
# LIMITER = NONE
# LIMITER = MINMOD
# LIMITER = SUPERBEE
LIMITER = MC
# LIMITER = VANLEER




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

riemann: $(RIEMANN_OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) -DNDIM=1

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $<






#==================================
#==================================

.PHONY: clean # keeps makefile from doing something with a file named clean

clean:
	rm -f $(EXEC) *.o
