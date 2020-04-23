#!/bin/bash

#--------------------------
# Run coverage tests
# made for gcc and gcovr
#--------------------------

# cleanup first
rm *html *json *gcda *gcno

# output file number
out=0

genmakefile(){
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



genmakefile_sources(){
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
    echo "SOURCES = $2"  >> $f
    echo "INTEGRATOR = $3"  >> $f
    echo "SOLVER = $4"  >> $f
    echo "RIEMANN = HLLC" >> $f
    echo "LIMITER = SUPERBEE" >> $f
}



function run_gcovr(){
    # run gcovr for a test case, create json intermediate output
    # then raise the output index
    # usage: run_gcovr
    fname=gcov_run"$out".json
    out=$((out + 1))
    gcovr . --root ../src/ --json -o "$fname"
}

function errexit() {
    # usage: errexit $? "optional message string"
    if [[ "$1" -ne 0 ]]; then
        echo "ERROR OCCURED. ERROR CODE $1"
        if [[ $# > 1 ]]; then
            echo "$2"
        fi
        exit $1
    else
        return 0
    fi
}



hydro_solvers="GODUNOV WAF MUSCL"
advection_solvers="ADVECTION_WAF ADVECTION_PWLIN ADVECTION_PWCONST"
ndims="1 2"
limiters="NONE MC SUPERBEE MINMOD VANLEER"
# limiters="NONE SUPERBEE MINMOD VANLEER"
riemann_solvers="EXACT HLLC TRRS TSRS"

paramfiles="paramfile-boundary-0  paramfile-boundary-2  paramfile-force-dt  paramfile-tmax paramfile-basename  paramfile-boundary-1  paramfile-dtout paramfile-foutput paramfile-dtoutlist"
icfiles_hydro="IC/riemann-left-vacuum.dat  IC/riemann-right-vacuum.dat  IC/riemann-sod-shock.dat  IC/riemann-vacuum-generating.dat"






#----------------------
# Run tests
#----------------------

# methods, solvers, limiters, dimensions...

echo "Running deault use coverage"

for ndim in $ndims; do
    for limiter in $limiters; do

        # advection
        for solver in $advection_solvers; do
            # advection tests
            make clean > /dev/null
            genmakefile $ndim $solver NONE $limiter #RIEMANN=NONE
            make > /dev/null
            ./hydro paramfile-default ./IC/advection-"$ndim"D.dat > /dev/null
            # ./hydro paramfile-default ./IC/advection-"$ndim"D.dat
            # errexit $?
            run_gcovr
        done

        # hydro
        for solver in $hydro_solvers; do
            for riemann in $riemann_solvers; do
                make clean > /dev/null
                genmakefile $ndim $solver $riemann $limiter
                make > /dev/null
                for icfile in $icfiles_hydro ./IC/advection-"$ndim"D.dat; do
                    ./hydro paramfile-default $icfile > /dev/null
                    # ./hydro paramfile-default $icfile
                    # errexit $?
                    run_gcovr
                done
            done
        done

    done
done


# Riemann solver tests

echo "Running Riemann solver coverage"

for riemann in $riemann_solvers; do
    make -f Makefile-Riemann clean > /dev/null
    genmakefile 1 NONE $riemann NONE
    make -f Makefile-Riemann > /dev/null
    for icfile in $icfiles_hydro; do
        for paramfile in paramfile-tmax paramfile-basename; do
            ./riemann $paramfile $icfile > /dev/null
            run_gcovr
        done
    done
done


# parameter IO tests

echo "running parameter I/O coverage"

for ndim in $ndims; do
    make clean > /dev/null
    genmakefile $ndim GODUNOV EXACT MINMOD;
    make > /dev/null
    for paramfile in $paramfiles; do
        ./hydro $paramfile ./IC/advection-"$ndim"D.dat > /dev/null
        run_gcovr
    done
done



# sources and integrators tests

echo "running sources and integrators coverage"

for ndim in $ndims; do
    for SOURCE in CONSTANT RADIAL; do
        for INTEGR in RK2 RK4; do
            for SOLVER in GODUNOV WAF MUSCL; do
                make clean > /dev/null
                genmakefile_sources $ndim $SOURCE $INTEGR $SOLVER;
                make > /dev/null
                ./hydro paramfile-sources ./IC/advection-"$ndim"D.dat > /dev/null
                run_gcovr
            done
        done
    done
done






#-------------------------
# generate html doc
#-------------------------

# generate cmd line args for all tracefiles
tracefileline=""
for i in `seq 0 $((out - 1))`; do 
    tracefileline="$tracefileline"" --add-tracefile ""gcov_run""$i".json
done


gcovr . --root ../src/ $tracefileline --html-details -o coverage.html
echo "finished writing coverage.html"


# clean up after yourself
rm *gcda *gcno *json *.out
