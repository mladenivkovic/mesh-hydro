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
    # $CMD paramfile.txt ../../IC/2D/kelvin-helmholtz-256.dat
    # $CMD paramfile.txt ../../IC/2D/rayleigh-taylor-256.dat
    # $CMD paramfile.txt ../../IC/2D/rayleigh-taylor-512.dat
    # $CMD paramfile.txt ../../IC/2D/kelvin-helmholtz-128.dat
    # $CMD paramfile.txt ../../IC/2D/corner-explosion-100.dat
    # $CMD paramfile.txt ../../IC/2D/center-explosion-100.dat
    # $CMD paramfile.txt ../../IC/1D/uniform-1D-100.dat
    # $CMD paramfile.txt ../../IC/2D/uniform-2D-100.dat
    $CMD paramfile.txt testing-IC-1D.dat
    python ../../py/plotting/plot_all_results_individually.py *out



    # RIEMANN
    #-----------------
    # icfile=../../IC/two-state/123problem.dat
    # icfile=../../IC/two-state/left_vacuum.dat
    # icfile=../../IC/two-state/right_vacuum.dat
    # icfile=../../IC/two-state/vacuum_generating.dat
    # icfile=../../IC/two-state/sod_test.dat
    # icfile=../../IC/two-state/sod_test_modified.dat
    # $CMD paramfile.txt $icfile | tee output.log
    # python3 ../../py/plotting/plot_all_riemann_results.py *-0001.out $icfile
    # python3 ../../py/plotting/plot_riemann_result.py *-0001.out $icfile




    # PLOTTING
    #----------------------

    # python3 ../../py/plotting/plot_density.py *0001.out
    # python3 ../../py/plotting/plot_result.py *0001.out
    # python3 ../../py/plotting/plot_all_results_individually.py *.out
    # python3 ../../py/plotting/plot_all_density.py *.out
    # python3 ../../py/plotting/plot_all_density_individually.py *.out


fi
