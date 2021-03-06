# Makefile to be included by other makefiles where the definitions and variables are processed.
# Also has some global definitions, like srcdir and vpath.
# Not intended to be used by itself!
# Mainly used so all the different makefiles can just import it and I don't have to keep track
# of what I updated where.


#==============================
# Preprocessing defines
#==============================

#----------------------------------------------------
# set default values if not defined otherwise
#----------------------------------------------------

ifndef SOLVER
SOLVER = GODUNOV
endif

ifndef RIEMANN
RIEMANN = EXACT
endif

ifndef LIMITER
LIMITER = NONE
endif

ifndef NDIM
NDIM = 1
endif

ifndef SOURCES
SOURCES = NONE
endif









#----------------------------------------------------
# Do some checks
#----------------------------------------------------

# make sure Riemann solvers are selected if 
# hydro is being solved

ifeq ($(strip $(SOLVER)), GODUNOV)
ifeq ($(strip $(RIEMANN)), NONE)
RIEMANN = EXACT
endif
endif
ifeq ($(strip $(SOLVER)), WAF)
ifeq ($(strip $(RIEMANN)), NONE)
RIEMANN = EXACT
endif
endif
ifeq ($(strip $(SOLVER)), MUSCL)
ifeq ($(strip $(RIEMANN)), NONE)
RIEMANN = EXACT
endif
endif


# make sure an integrator is selected if
# sources are selected
ifneq ($(strip $(SOURCES)), NONE)
ifndef INTEGRATOR
INTEGRATOR = RK2
endif
ifeq ($(strip $(INTEGRATOR)), NONE)
INTEGRATOR = RK2
endif
else # if SOURCES = NONE
INTEGRATOR = NONE
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
ifeq ($(strip $(SOLVER)), ADVECTION_WAF)
ADVECTION = true
RIEMANN = NONE
endif







# ------------------------------------------------
# transform defines into integers where needed
# ------------------------------------------------

ifeq ($(strip $(SOLVER)), ADVECTION_PWCONST)
SOLVERINT = 11
endif
ifeq ($(strip $(SOLVER)), ADVECTION_PWLIN)
SOLVERINT = 12
endif
ifeq ($(strip $(SOLVER)), ADVECTION_WAF)
SOLVERINT = 13
endif
ifeq ($(strip $(SOLVER)), GODUNOV)
SOLVERINT = 21
endif
ifeq ($(strip $(SOLVER)), WAF)
SOLVERINT = 22
endif
ifeq ($(strip $(SOLVER)), MUSCL)
SOLVERINT = 23
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
ifeq ($(strip $(RIEMANN)), TRRS)
RIEMANNINT = 3
endif
ifeq ($(strip $(RIEMANN)), TSRS)
RIEMANNINT = 4
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



ifeq ($(strip $(SOURCES)), NONE)
SOURCESINT = 0
endif
ifeq ($(strip $(SOURCES)), CONSTANT)
SOURCESINT = 1
endif
ifeq ($(strip $(SOURCES)), RADIAL)
SOURCESINT = 2
endif





#-----------------------------------
# Set up Definition Flags
#-----------------------------------

COMPILEDATE:=$(shell date "+%F %T")


DEFINES= -DNDIM=$(NDIM) -DSOLVER=$(SOLVERINT) -DRIEMANN=$(RIEMANNINT) -DLIMITER=$(LIMITERINT) -DSOURCE=$(SOURCESINT) -DCOMPDATE="$(COMPILEDATE)" 


ifdef ADVECTION
DEFINES += -DADVECTION
endif

ifneq ($(strip $(SOURCES)), NONE)
DEFINES += -DWITH_SOURCES
endif





#==================================
# FILE LISTS
#==================================

SRCDIR=../src

#include paths. Will be followed in that order.
VPATH=$(SRCDIR):$(SRCDIR)/limiter:$(SRCDIR)/solver:$(SRCDIR)/riemann:$(SRCDIR)/sources:$(SRCDIR)/integrate

#include directories for headers
IDIR=$(SRCDIR)


ifeq ($(strip $(SOLVER)), ADVECTION_PWCONST)
	HYDROOBJ=advection_pwconst.o
endif
ifeq ($(strip $(SOLVER)), ADVECTION_PWLIN)
	HYDROOBJ=advection_pwlin.o
endif
ifeq ($(strip $(SOLVER)), GODUNOV)
	HYDROOBJ=godunov.o
endif
ifeq ($(strip $(SOLVER)), ADVECTION_WAF)
	HYDROOBJ=advection_waf.o
endif
ifeq ($(strip $(SOLVER)), WAF)
	HYDROOBJ=waf.o
endif
ifeq ($(strip $(SOLVER)), MUSCL)
	HYDROOBJ=muscl.o
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


ifeq ($(strip $(RIEMANN)), NONE)
	RIEMANNOBJ=
endif
ifeq ($(strip $(RIEMANN)), EXACT)
	RIEMANNOBJ=riemann-exact.o riemann.o
endif
ifeq ($(strip $(RIEMANN)), TRRS)
	RIEMANNOBJ=riemann-trrs.o riemann.o
endif
ifeq ($(strip $(RIEMANN)), TSRS)
	RIEMANNOBJ=riemann-tsrs.o riemann.o
endif
ifeq ($(strip $(RIEMANN)), HLL)
	RIEMANNOBJ=riemann-hll.o riemann.o
endif
ifeq ($(strip $(RIEMANN)), HLLC)
	RIEMANNOBJ=riemann-hllc.o riemann.o
endif



ifeq ($(strip $(SOURCES)), NONE)
	SRCOBJ=
endif
ifeq ($(strip $(SOURCES)), CONSTANT)
	SRCOBJ=sources.o sources-constant.o
endif
ifeq ($(strip $(SOURCES)), RADIAL)
	SRCOBJ=sources.o sources-radial.o
endif



ifeq ($(strip $(INTEGRATOR)), NONE)
	INTOBJ=
endif
ifeq ($(strip $(INTEGRATOR)), RK2)
	INTOBJ=integrate-runge-kutta-2.o
endif
ifeq ($(strip $(INTEGRATOR)), RK4)
	INTOBJ=integrate-runge-kutta-4.o
endif






OBJECTS = main.o gas.o params.o io.o utils.o cell.o solver.o limiter.o $(HYDROOBJ) $(LIMITEROBJ) $(RIEMANNOBJ) $(SRCOBJ) $(INTOBJ)
RIEMANN_OBJECTS = main-riemann.o gas.o params.o io.o utils.o cell.o $(RIEMANNOBJ)
