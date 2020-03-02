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


make clean
make -f Makefile-Riemann clean


function advection(){

    #---------------------------
    # Vanilla advection
    #---------------------------
    for ndim in 1 2; do
        for solver in ADVECTION_PWLIN ADVECTION_PWCONST; do
            make clean
            genmakefile $ndim $solver NONE NONE hydro-$solver-"$ndim"D
            make
        done
    done



    #---------------------------
    # slope limiters
    #---------------------------
    for ndim in 1 2; do
        for LIMITER in MINMOD SUPERBEE MC VANLEER; do
            make clean
            genmakefile $ndim ADVECTION_PWLIN NONE $LIMITER hydro-advection-"$LIMITER"-"$ndim"D
            make
        done
    done
    rm -f defines.mk
}



function godunov(){

    #---------------------------
    # Godunov
    #---------------------------
    for ndim in 1 2; do
        for RIEMANN in EXACT TRRS TSRS HLLC; do
            make clean
            genmakefile $ndim GODUNOV $RIEMANN NONE hydro-godunov-"$RIEMANN"-"$ndim"D
            make
        done
    done
    rm -f defines.mk
}



function riemann(){

    #---------------------------
    # Code as a Riemann Solver
    #---------------------------

    for RIEMANN in EXACT TRRS TSRS; do
        make -f Makefile-Riemann clean
        genmakefile 1 NONE $RIEMANN NONE riemann-"$RIEMANN"
        make -f Makefile-Riemann
    done
    rm -f defines.mk
}



advection
godunov
riemann
