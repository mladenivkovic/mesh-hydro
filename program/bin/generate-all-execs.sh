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


function errexit() {
    # usage: errexit $? "optional message string"
    if [[ "$1" -ne 0 ]]; then
        myecho "ERROR OCCURED. ERROR CODE $1"
        if [[ $# > 1 ]]; then
            myecho "$2"
        fi
        exit $1
    else
        return 0
    fi
}



make clean
make -f Makefile-Riemann clean


function advection(){

    #---------------------------
    # Vanilla advection
    #---------------------------
    for ndim in 1 2; do
        for solver in ADVECTION_PWLIN ADVECTION_PWCONST ADVECTION_WAF; do
            make clean
            genmakefile $ndim $solver NONE NONE hydro-$solver-"$ndim"D
            make
            errexit $?
        done
    done



    #---------------------------
    # slope limiters
    #---------------------------
    for ndim in 1 2; do
        for solver in ADVECTION_PWLIN ADVECTION_WAF; do
            for LIMITER in MINMOD SUPERBEE MC VANLEER; do
                make clean
                genmakefile $ndim $solver NONE $LIMITER hydro-"$solver"-"$LIMITER"-"$ndim"D
                make
                errexit $?
            done
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
            errexit $?
        done
    done
    rm -f defines.mk
}



function waf(){

    #---------------------------
    # Godunov
    #---------------------------
    for ndim in 1 2; do
        for RIEMANN in EXACT TRRS TSRS HLLC; do
            for LIMITER in MINMOD SUPERBEE MC VANLEER; do
                make clean
                genmakefile $ndim WAF $RIEMANN NONE hydro-WAF-"$RIEMANN"-"$LIMITER"-"$ndim"D
                make
                errexit $?
            done
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
        errexit $?
    done
    rm -f defines.mk
}



advection
godunov
riemann
waf
