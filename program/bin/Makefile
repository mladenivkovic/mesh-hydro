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
# Choices: EXACT, TRRS, TSRS, HLL, HLLC
# RIEMANN = EXACT
RIEMANN = NONE


# set which hydro solver to use.
# Choices: GODUNOV_UPWIND, ADVECTION_PWCONST, ADVECTION_PWLIN
# SOLVER=GODUNOV_UPWIND
SOLVER = ADVECTION_PWCONST
# SOLVER = ADVECTION_PWLIN

# set which flux limiter to use.
# Choices: NONE, MINMOD, SUPERBEE, VANLEER, MC
LIMITER = NONE
# LIMITER = MINMOD
# LIMITER = SUPERBEE
# LIMITER = MC
# LIMITER = VANLEER






#==============================
# Preprocessing defines
#==============================

# set default values if not defined otherwise
ifndef RIEMANN
RIEMANN = EXACT
endif
ifndef SOLVER
SOLVER = GODUNOV_UPWIND
endif
ifndef LIMITER
LIMITER = NONE
endif
ifndef NDIM
NDIM = 2
endif


# set advection flag if solver is advection
ifeq ($(strip $(SOLVER)), ADVECTION_PWLIN)
ADVECTION = true
RIEMANN = NONE
endif
ifeq ($(strip $(SOLVER)), ADVECTION_PWCONST)
ADVECTION = true
RIEMANN = NONE
LIMITER = NONE
endif


# transform defines into integers where needed
ifeq ($(strip $(SOLVER)), ADVECTION_PWCONST)
# ifeq ($(strip $(SOLVER)), ADVECTION_PWCONST)
SOLVERINT = 1
endif
ifeq ($(strip $(SOLVER)), ADVECTION_PWLIN)
SOLVERINT = 2
endif
ifeq ($(strip $(SOLVER)), GODUNOV_UPWIND)
SOLVERINT = 3
endif


ifeq ($(strip $(RIEMANN)), NONE)
RIEMANNINT = 0
endif
ifeq ($(strip $(RIEMANN)), EXACT)
RIEMANNINT = 1
endif
ifeq ($(strip $(RIEMANN)), HLLC)
RIEMANNINT = 2
endif
ifeq ($(strip $(RIEMANN)), HLL)
RIEMANNINT = 3
endif
ifeq ($(strip $(RIEMANN)), TRRS)
RIEMANNINT = 4
endif
ifeq ($(strip $(RIEMANN)), TSRS)
RIEMANNINT = 5
endif


ifeq ($(strip $(LIMITER)), NONE)
LIMITERINT = 0
endif
ifeq ($(strip $(LIMITER)), MINMOD)
LIMITERINT = 1
endif
ifeq ($(strip $(LIMITER)), SUPERBEE)
LIMITERINT = 2
endif
ifeq ($(strip $(LIMITER)), VANLEER)
LIMITERINT = 3
endif
ifeq ($(strip $(LIMITER)), MC)
LIMITERINT = 4
endif


COMPILEDATE=$(shell date "+%F %T")





DEFINES= -DNDIM=$(NDIM) -DSOLVER=$(SOLVERINT) -DRIEMANN=$(RIEMANNINT) -DLIMITER=$(LIMITERINT) -DCOMPDATE="$(COMPILEDATE)" 


ifdef ADVECTION
DEFINES += -DADVECTION
endif




#==================================
# FILE LISTS
#==================================

SRCDIR=../src

#include paths. Will be followed in that order.
VPATH=$(SRCDIR):$(SRCDIR)/limiter:$(SRCDIR)/solver:$(SRCDIR)/riemann

#include directories for headers
# IDIR=$(VPATH)
IDIR=$(SRCDIR)


ifeq ($(strip $(SOLVER)), ADVECTION_PWCONST)
	HYDROOBJ=advection_pwconst.o
endif
ifeq ($(strip $(SOLVER)), ADVECTION_PWLIN)
	HYDROOBJ=advection_pwlin.o
endif
ifeq ($(strip $(SOLVER)), GODUNOV_UPWIND)
	HYDROOBJ=godunov_upwind.o
endif


ifeq ($(strip $(LIMITER)), NONE)
	LIMITEROBJ=no_limiter.o
endif
ifeq ($(strip $(LIMITER)), MINMOD)
	LIMITEROBJ=minmod.o
endif
ifeq ($(strip $(LIMITER)), SUPERBEE)
	LIMITEROBJ=superbee.o
endif
ifeq ($(strip $(LIMITER)), VANLEER)
	LIMITEROBJ=van_leer.o
endif
ifeq ($(strip $(LIMITER)), MC)
	LIMITEROBJ=monotonized_central_difference.o
endif


# OBJECTS = $(notdir $(patsubst %.c, %.o, $(wildcard $(SRCDIR)/*.c)))
# OBJECTS = main.o gas.o godunov.o io.o params.o  utils.o $(RIEMANNFILE)
OBJECTS = main.o gas.o params.o io.o utils.o cell.o $(HYDROOBJ) $(LIMITEROBJ)








#==================================
# COMPILER OPTIONS AND DIRECTIVES
#==================================


#Compiler to use
CC=gcc -I$(IDIR) -static


# flags for gcc compiler
# CFLAGS= -g -O0 -Wall -pedantic -std=c11 -Wextra -Werror -fsignaling-nans  -ffinite-math-only -fno-signed-zeros -Warray-bounds # debug
CFLAGS= -O3  		 	   # performance

CFLAGS += $(DEFINES)

LDFLAGS= -lm -fopenmp








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
	rm -f $(OBJECTS) $(EXEC) riemann*.o
