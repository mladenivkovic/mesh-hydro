#!/bin/bash

#-------------------------------------------
# Automatically compile and run the hydro
# program with all solvers and various IC
# files.
# Then plot the solutions, and create a TeX
# file where the images are compared.
#-------------------------------------------



#===================================
genparamfile() {
#===================================
    # generate parameter file.
    # $1 = nsteps
    # $2 = tmax
    # $3 = output frequency: after how many steps to write
    # $4 = dt_out
    # $5 = basename
    # $6 = Ccfl

    f=paramfile.txt
    echo "// parameter file for hydro program" > $f
    echo ""                 >> $f
    echo "verbose = 2"      >> $f
    echo "nx = 100"         >> $f
    echo "nsteps = $1"      >> $f
    echo "tmax = $2"        >> $f
    echo "foutput = $3"     >> $f
    echo "dt_out = $4"      >> $f
    echo "basename = $5"    >> $f
    echo "ccfl = $6"        >> $f
}


#==========================
genmakefile(){
#==========================
    # generate makefile "header" which 
    # will be included into the actual 
    # Makefile.
    # $1: ndim
    # $2: solver
    # $3: riemann solver
    # $3: limiter

    f=defines.mk
    echo "# generated with run.sh" > $f
    echo ""             >> $f
    echo "NDIM = $1"    >> $f
    echo "SOLVER = $2"  >> $f
    echo "RIEMANN = $3" >> $f
    echo "LIMITER = $4" >> $f

}



function errexit() {
    # usage: errexit $? "optional message string"
    if [[ "$1" -ne 0 ]]; then
        myecho "ERROR OCCURED. ERROR CODE $1"
        if [[ $# > 1 ]]; then
            myecho "$2"
        fi
        traceback 1
        exit $1
    else
        return 0
    fi
}



plotdir=../../py/plotting


# cleanup first
rm *.out *.log
# rm *.png


# create global log file
touch output.log


#==========================
# ADVECTION
#==========================


#-------------------
# 1D PWCONST
#-------------------

# genmakefile ndim solver riemann limiter
genmakefile 1 ADVECTION_PWCONST NONE NONE
make clean && make | tee -a output.log
errexit $?

# positive velocity
# genparamfile nsteps tmax foutput dt_out basename ccfl
genparamfile 0 10 0 1.0 advection-1D-pwconst 0.5

./hydro paramfile.txt ./IC/advection-1D.dat | tee -a output.log
errexit $?
$plotdir/plot_all_density.py advection-1D-pwconst-0000.out

# negative velocity
genparamfile 0 10 0 1.0 advection-1D-pwconst-negvel 0.5
./hydro paramfile.txt ./IC/advection-1D-negvel.dat | tee -a output.log
errexit $?
$plotdir/plot_all_density.py advection-1D-pwconst-negvel-0000.out


#-------------------
# 2D PWCONST
#-------------------

# genmakefile ndim solver riemann limiter
genmakefile 2 ADVECTION_PWCONST NONE NONE
make clean && make | tee -a output.log
errexit $?

# positive velocity
# genparamfile nsteps tmax foutput dt_out basename ccfl
genparamfile 0 1 0 0 advection-2D-pwconst 0.1

./hydro paramfile.txt ./IC/advection-2D.dat | tee -a output.log
errexit $?
$plotdir/plot_density.py advection-2D-pwconst-0001.out

# negative velocity
genparamfile 0 1 0 0 advection-2D-pwconst-negvel 0.1
./hydro paramfile.txt ./IC/advection-2D-negvel.dat | tee -a output.log
errexit $?
$plotdir/plot_density.py advection-2D-pwconst-negvel-0001.out



#-------------------
# 1D PWLIN
#-------------------

# genmakefile ndim solver riemann limiter
genmakefile 1 ADVECTION_PWLIN NONE NONE
make clean && make | tee -a output.log
errexit $?

# positive velocity
# genparamfile nsteps tmax foutput dt_out basename ccfl
genparamfile 0 10 0 1.0 advection-1D-pwlin 0.5

./hydro paramfile.txt ./IC/advection-1D.dat | tee -a output.log
errexit $?
$plotdir/plot_all_density.py advection-1D-pwlin-0000.out

# negative velocity
genparamfile 0 10 0 1.0 advection-1D-pwlin-negvel 0.5
./hydro paramfile.txt ./IC/advection-1D-negvel.dat | tee -a output.log
errexit $?
$plotdir/plot_all_density.py advection-1D-pwlin-negvel-0000.out


#-------------------
# 2D pwlin
#-------------------

# genmakefile ndim solver riemann limiter
genmakefile 2 ADVECTION_PWLIN NONE NONE
make clean && make | tee -a output.log
errexit $?

# positive velocity
# genparamfile nsteps tmax foutput dt_out basename
genparamfile 0 1 0 0 advection-2D-pwlin 0.1

./hydro paramfile.txt ./IC/advection-2D.dat | tee -a output.log
errexit $?
$plotdir/plot_density.py advection-2D-pwlin-0001.out

# negative velocity
genparamfile 0 1 0 0 advection-2D-pwlin-negvel 0.1
./hydro paramfile.txt ./IC/advection-2D-negvel.dat | tee -a output.log
errexit $?
$plotdir/plot_density.py advection-2D-pwlin-negvel-0001.out





#----------------------------
# 1D pwlin with limiters
#----------------------------

for LIMITER in MINMOD SUPERBEE MC VANLEER; do

    # genmakefile ndim solver riemann limiter
    genmakefile 1 ADVECTION_PWLIN NONE $LIMITER
    make clean && make | tee -a output.log
    errexit $?

    # positive velocity
    # genparamfile nsteps tmax foutput dt_out basename ccfl
    genparamfile 0 10 0 1.0 advection-1D-$LIMITER 0.5

    ./hydro paramfile.txt ./IC/advection-1D.dat | tee -a output.log
    errexit $?
    $plotdir/plot_all_density.py advection-1D-$LIMITER-0000.out

    # negative velocity
    genparamfile 0 10 0 1.0 advection-1D-pwlin-negvel 0.5
    ./hydro paramfile.txt ./IC/advection-1D-negvel.dat | tee -a output.log
    errexit $?
    $plotdir/plot_all_density.py advection-1D-$LIMITER-negvel-0000.out

done





#----------------------------
# 2D pwlin with limiters
#----------------------------

for LIMITER in MINMOD SUPERBEE MC VANLEER; do

    # genmakefile ndim solver riemann limiter
    genmakefile 2 ADVECTION_PWLIN NONE $LIMITER
    make clean && make | tee -a output.log
    errexit $?

    # positive velocity
    # genparamfile nsteps tmax foutput dt_out basename ccfl
    genparamfile 0 1 0 0 advection-2D-$LIMITER 0.1

    ./hydro paramfile.txt ./IC/advection-2D.dat | tee -a output.log
    errexit $?
    $plotdir/plot_density.py advection-2D-$LIMITER-0001.out

done






#---------------
# create TeX
#---------------
cd TeX
./run.sh
cd .
