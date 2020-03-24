#!/bin/bash
#
# A script to automatize running stuff.
# use with -d to call program with gdb.
# Otherwise, uncomment the right lines that you want to use, or write your own.

PROGNAME=hydro
make clean
make

# PROGNAME=riemann
# make -f Makefile-Riemann clean
# make -f Makefile-Riemann

if [ $? -eq 0 ]; then # only continue if no error happened

    CMD='./hydro'

    if [ $# -gt 0 ]; then
        case "$1" in

            -d | --debug | deb | d | --d)
                CMD='gdb --args ./hydro'
            ;;

            *)
                echo "Unrecognized argument '"$1"'"
            ;;

        esac
    fi


    # RUNNING THE PROGRAM
    #------------------------

    # HYDRO
    #-----------------

    # $CMD paramfile.txt ../../IC/advection/advection-1D-step-100.dat | tee output.log
    # $CMD paramfile.txt ../../IC/advection/advection-1D-step-500.dat
    # $CMD paramfile.txt ../../IC/advection/advection-1D-step-100-negative-velocity.dat | tee output.log
    # $CMD paramfile.txt ../../IC/advection/advection-1D-step-1000.dat | tee output.log
    # $CMD paramfile.txt ../../IC/advection/advection-2D-step-100-ux-uy-negative-velocity.dat | tee output.log
    # $CMD paramfile.txt ../../IC/advection/advection-2D-step-100-ux-uy.dat | tee output.log
    # $CMD paramfile.txt ../../IC/advection/advection-2D-step-100-ux-uy-small-velocity.dat | tee output.log
    # $CMD paramfile.txt ../../IC/advection/advection-2D-step-100-ux.dat | tee output.log
    # $CMD paramfile.txt ../../IC/advection/advection-2D-step-100-uy.dat | tee output.log
    # $CMD paramfile.txt ../../IC/advection/advection-1D-four-shapes.dat | tee output.log
    # $CMD paramfile.txt ../../IC/advection/advection-2D-four-shapes.dat | tee output.log
    # $CMD paramfile.txt ../../IC/two-state/sod_test.dat
    # $CMD paramfile.txt ../../IC/two-state/sod_test_modified.dat
    # $CMD paramfile.txt ../../IC/two-state/vacuum_generating.dat
    $CMD paramfile.txt ../../IC/2D/kelvin-helmholtz-256.dat
    # $CMD paramfile.txt ../../IC/2D/corner-explosion-100.dat



    # RIEMANN
    #-----------------
    # $CMD paramfile.txt ../../IC/two-state/vacuum_generating.dat
    # $CMD paramfile.txt ../../IC/two-state/two_shocks.dat
    # $CMD paramfile.txt ../../IC/two-state/sod_test.dat
    # $CMD paramfile.txt ../../IC/two-state/left_vacuum.dat



    # PLOTTING
    #----------------------

    # python3 ../../py/plotting/plot_density.py *0001.out
    # python3 ../../py/plotting/plot_result.py *0001.out
    # python3 ../../py/plotting/plot_all_density.py *.out
    # python3 ../../py/plotting/plot_all_density_individually.py *.out
    python3 ../../py/plotting/plot_all_results_individually.py *.out

    # python3 ../../py/plotting/plot_riemann_result.py tsrs-solver-0001.out ../../IC/two-state/sod_test.dat


fi
