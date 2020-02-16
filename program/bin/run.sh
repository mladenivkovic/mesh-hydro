#!/bin/bash

make clean
make
if [ $? -eq 0 ]; then 

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

    # $CMD paramfile.txt ../../IC/advection/advection-1D-step-10.dat | tee output.log
    # $CMD paramfile.txt ../../IC/advection/advection-1D-step-100.dat | tee output.log
    # $CMD paramfile.txt ../../IC/advection/advection-1D-step-10-negative-velocity.dat | tee output.log
    # $CMD paramfile.txt ../../IC/advection/advection-1D-step-100-negative-velocity.dat | tee output.log
    $CMD paramfile.txt ../../IC/advection/advection-1D-step-1000.dat | tee output.log
    # $CMD paramfile.txt ../../IC/advection/advection-2D-step-100-ux-uy-negative-velocity.dat | tee output.log
    # $CMD paramfile.txt ../../IC/advection/advection-2D-step-100-ux-uy.dat | tee output.log
    # $CMD paramfile.txt ../../IC/advection/advection-2D-step-100-ux-uy-small-velocity.dat | tee output.log
    # $CMD paramfile.txt ../../IC/advection/advection-2D-step-100-ux.dat | tee output.log
    # $CMD paramfile.txt ../../IC/advection/advection-2D-step-100-uy.dat | tee output.log
    # $CMD paramfile.txt ../../IC/advection/advection-2D-step-10.dat | tee output.log
    # $CMD paramfile.txt ./restart.dat | tee output.log
    # $CMD paramfile.txt ../../IC/advection/advection-1D-four-shapes.dat | tee output.log
    # $CMD paramfile.txt ../../IC/advection/advection-2D-four-shapes.dat | tee output.log


    # ../../py/plotting/plot_all_density.py *.out
    # ../../py/plotting/plot_all_density_individually.py *.out
    # ../../py/plotting/plot_result.py *0001.out
    # eog *.png

fi
