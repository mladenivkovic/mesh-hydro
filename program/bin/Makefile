#==================================
# COMPILATION TIME PARAMETERS
#==================================

# executable name
EXEC = hydro  

# how many dimensions to compile for
# Choices: 1, 2
# NDIM = 1
NDIM = 2

# set which hydro solver to use.
# Choices: GODUNOV, WAF, MUSCL, ADVECTION_PWCONST, ADVECTION_PWLIN, ADVECTION_WAF
SOLVER = GODUNOV
# SOLVER = WAF
# SOLVER = MUSCL
# SOLVER = ADVECTION_PWCONST
# SOLVER = ADVECTION_PWLIN
# SOLVER = ADVECTION_WAF


# set which Riemann solver to use.
# Choices: EXACT, TRRS, TSRS, HLLC
# RIEMANN = NONE
# RIEMANN = EXACT
# RIEMANN = TRRS
# RIEMANN = TSRS
RIEMANN = HLLC

# set which flux limiter to use.
# Choices: NONE, MINMOD, SUPERBEE, VANLEER, MC
LIMITER = NONE
# LIMITER = MINMOD
# LIMITER = SUPERBEE
# LIMITER = MC
# LIMITER = VANLEER


# set whether and which source terms to use
# Choices: NONE, CONSTANT, RADIAL
SOURCES = NONE
# SOURCES = CONSTANT
# SOURCES = RADIAL


# set which integrator to use (only for sources)
# Choices: RK2, RK4
INTEGRATOR = NONE
# INTEGRATOR = RK2
# INTEGRATOR = RK4






# if there is a defines.mk present, overwrite the
# definitions above with its contents
ifneq ("$(wildcard defines.mk)","")
$(warning WARNING: ================================================)
$(warning WARNING: READING IN DEFINITIONS FROM defines.mk)
$(warning WARNING: ================================================)
include defines.mk
endif


# process the definitions
include processing.mk






#==================================
# COMPILER OPTIONS AND DIRECTIVES
#==================================

# Compiler to use

# ---------------------------------------------------------------------
#  USING GCC

CC=gcc -I$(IDIR) #-static

# flags for gcc compiler
CFLAGS= -g -O0 -pedantic -std=c11 \
		-fsignaling-nans -ffinite-math-only -fno-signed-zeros \
		-Wall -Wextra -Werror -Warray-bounds -Wno-unused-parameter
# CFLAGS= -O3


# ---------------------------------------------------------------------
#  USING INTEL ONEAPI

# CC=icx -I$(IDIR) #-static

# flags for gcc compiler
# CFLAGS= -g -O0 -pedantic -std=c11 \
#         -Wall -Wextra -Werror -Wno-unused-parameter
# CFLAGS= -O3


# ---------------------------------------------------------------------


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
