#!/bin/bash

# compare available slope limiter's effects on 1D advection.
# you need to set following variables:
# 
# EXECDIR: the directory containing executables of the program created by the /program/bin/generate-all-execs.sh script
# ICDIR:   the directory containing the IC file you require
# EVALDIR: the directory containing python evaluation and plotting scripts, i.e. /py/evaluations
#
# you also need to prepare your own IC files. Scripts to do that are in /py/IC

EXECDIR=
ICDIR=
EVALDIR=


for limiter in MC MINMOD SUPERBEE VANLEER; do
    $EXECDIR/hydro-advection-$limiter-1D paramfile.txt $ICDIR/advection-1D-four-shapes.dat | tee -a output.log
done

$EXECDIR/hydro-ADVECTION_PWLIN-1D paramfile.txt $ICDIR/advection-1D-four-shapes.dat | tee -a output.log

$EVALDIR/advection-1D-compare-limiters.py

