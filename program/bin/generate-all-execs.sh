#!/bin/bash


#-----------------------------------------------
# This script generates all possible
# version of the hydro code and names
# the different executables appropriately.
#-----------------------------------------------

genmakefile(){
    # generate makefile "header" which 
    # will be included into the actual 
    # Makefile.
    # $1: ndim
    # $2: solver
    # $3: riemann solver
    # $3: limiter

    f=defines.mk
    echo "# generated with generate-all-execs.sh" > $f
    echo ""             >> $f
    echo "NDIM = $1"    >> $f
    echo "SOLVER = $2"  >> $f
    echo "RIEMANN = $3" >> $f
    echo "LIMITER = $4" >> $f
    echo "EXEC = $5"    >> $f

}


make -f Makefile-individual-execnames clean


# Vanilla advection
for ndim in 1 2; do
    for solver in ADVECTION_PWLIN ADVECTION_PWCONST; do

        make -f Makefile-individual-execnames clean
        genmakefile $ndim $solver NONE NONE hydro-$solver-"$ndim"D
        make -f Makefile-individual-execnames
    done
done



# slope limiters
for ndim in 1 2; do
    for LIMITER in MINMOD SUPERBEE MC VANLEER; do

        make -f Makefile-individual-execnames clean
        genmakefile $ndim ADVECTION_PWLIN NONE $LIMITER hydro-advection-"$LIMITER"-"$ndim"D
        make -f Makefile-individual-execnames
    done
done

