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
    # $5 = basename. if == "NO_BASENAME" it will be skipped
    # $6 = Ccfl

    f=paramfile.txt
    echo "// parameter file for hydro program" > $f
    echo "// generated by run.sh" >> $f
    echo ""                 >> $f
    echo "verbose = 1"      >> $f
    echo "nx = 100"         >> $f
    echo "nstep_log = 50"   >> $f
    echo "nsteps = $1"      >> $f
    echo "tmax = $2"        >> $f
    echo "foutput = $3"     >> $f
    echo "dt_out = $4"      >> $f
    if [ "$5" != "NO_BASENAME" ]; then
        echo "basename = $5"    >> $f
    fi
    echo "ccfl = $6"        >> $f
}




#===================================
genparamfile_transmissive() {
#===================================
    # generate parameter file.
    # $1 = nsteps
    # $2 = tmax
    # $3 = output frequency: after how many steps to write
    # $4 = dt_out
    # $5 = basename. if == "NO_BASENAME" it will be skipped
    # $6 = Ccfl

    f=paramfile.txt
    echo "// parameter file for hydro program" > $f
    echo "// generated by run.sh" >> $f
    echo ""                 >> $f
    echo "verbose = 1"      >> $f
    echo "nx = 100"         >> $f
    echo "nstep_log =  50"  >> $f
    echo "nsteps = $1"      >> $f
    echo "tmax = $2"        >> $f
    echo "foutput = $3"     >> $f
    echo "dt_out = $4"      >> $f
    if [ "$5" != "NO_BASENAME" ]; then
        echo "basename = $5"    >> $f
    fi
    echo "ccfl = $6"        >> $f
    echo "boundary = 2"     >> $f
}




#===================================
genparamfile_sources() {
#===================================
    # generate parameter file.
    # $1 = nsteps
    # $2 = tmax
    # $3 = basename

    f=paramfile.txt
    echo "// parameter file for hydro program" > $f
    echo "// generated by run.sh" >> $f
    echo ""                 >> $f
    echo "verbose = 1"      >> $f
    echo "nx = 100"         >> $f
    echo "nstep_log =  50"  >> $f
    echo "nsteps = $1"      >> $f
    echo "tmax = $2"        >> $f
    if [ "$3" != "NO_BASENAME" ]; then
        echo "basename = $3"    >> $f
    fi
    echo "ccfl = 0.9"        >> $f
    echo "boundary = 1"     >> $f
    echo "src_const_acc_x = 1.0" >> $f
    echo "src_const_acc_y = 2.0" >> $f
    echo "src_const_acc_r = -1.0" >> $f
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






#==========================
genmakefile_sources(){
#==========================
    # generate makefile "header" which
    # will be included into the actual
    # Makefile.
    # $1: ndim
    # $2: solver
    # $3: source
    # $3: integrator

    f=defines.mk
    echo "# generated with run.sh" > $f
    echo ""             >> $f
    echo "NDIM = $1"    >> $f
    echo "SOLVER = $2"  >> $f
    echo "RIEMANN = HLLC" >> $f
    echo "LIMITER = SUPERBEE" >> $f
    echo "SOURCES = $3" >> $f
    echo "INTEGRATOR = $4" >> $f

}





function myecho(){
    echo "================ $1"
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



generate_tex(){
    cd TeX
    ./run.sh
    cd .
}


plotdir=../../python_module/scripts/plotting


# cleanup first
rm -f *.out *.log
# rm *.png




#==========================
# FUNCTIONS FOR TESTS
#==========================

# so they'll be easier to call for specific purposes/tests in groups


advection_pwconst_1D(){
    #-------------------
    # 1D PWCONST
    #-------------------

    # genmakefile ndim solver riemann limiter
    genmakefile 1 ADVECTION_PWCONST NONE NONE
    make clean && make
    errexit $?

    # positive velocity
    # genparamfile nsteps tmax foutput dt_out basename ccfl
    genparamfile 0 10 0 1.0 advection-1D-pwconst 0.8

    ./hydro paramfile.txt ./IC/advection-1D.dat
    errexit $?
    $plotdir/plot_all_density.py advection-1D-pwconst-0000.out
    errexit $?

    # negative velocity
    genparamfile 0 10 0 1.0 advection-1D-pwconst-negvel 0.8
    ./hydro paramfile.txt ./IC/advection-1D-negvel.dat
    errexit $?
    $plotdir/plot_all_density.py advection-1D-pwconst-negvel-0000.out
    errexit $?
}




advection_pwconst_2D(){
    #-------------------
    # 2D PWCONST
    #-------------------

    # genmakefile ndim solver riemann limiter
    genmakefile 2 ADVECTION_PWCONST NONE NONE
    make clean && make
    errexit $?

    # positive velocity
    # genparamfile nsteps tmax foutput dt_out basename ccfl
    genparamfile 0 1 0 0 advection-2D-pwconst 0.8
    ./hydro paramfile.txt ./IC/advection-2D.dat
    errexit $?
    $plotdir/plot_density.py advection-2D-pwconst-0001.out
    errexit $?

    # negative velocity
    genparamfile 0 1 0 0 advection-2D-pwconst-negvel 0.8
    ./hydro paramfile.txt ./IC/advection-2D-negvel.dat
    errexit $?
    $plotdir/plot_density.py advection-2D-pwconst-negvel-0001.out
    errexit $?

    # x only
    genparamfile 0 1 0 0 advection-2D-pwconst-x 0.8
    ./hydro paramfile.txt ./IC/advection-2D-x.dat
    errexit $?
    $plotdir/plot_density.py advection-2D-pwconst-x-0001.out
    errexit $?

    # y only
    genparamfile 0 1 0 0 advection-2D-pwconst-y 0.8
    ./hydro paramfile.txt ./IC/advection-2D-y.dat
    errexit $?
    $plotdir/plot_density.py advection-2D-pwconst-y-0001.out
    errexit $?
}




advection_pwlin_1D(){

    #-------------------
    # 1D PWLIN
    #-------------------

    # genmakefile ndim solver riemann limiter
    genmakefile 1 ADVECTION_PWLIN NONE NONE
    make clean && make
    errexit $?

    # positive velocity
    # genparamfile nsteps tmax foutput dt_out basename ccfl
    genparamfile 0 10 0 1.0 advection-1D-pwlin 0.8

    ./hydro paramfile.txt ./IC/advection-1D.dat
    errexit $?
    $plotdir/plot_all_density.py advection-1D-pwlin-0000.out
    errexit $?

    # negative velocity
    genparamfile 0 10 0 1.0 advection-1D-pwlin-negvel 0.8
    ./hydro paramfile.txt ./IC/advection-1D-negvel.dat
    errexit $?
    $plotdir/plot_all_density.py advection-1D-pwlin-negvel-0000.out
    errexit $?

}



advection_pwlin_2D(){
    #-------------------
    # 2D pwlin
    #-------------------

    # genmakefile ndim solver riemann limiter
    genmakefile 2 ADVECTION_PWLIN NONE NONE
    make clean && make
    errexit $?

    # positive velocity
    # genparamfile nsteps tmax foutput dt_out basename
    genparamfile 0 1 0 0 advection-2D-pwlin 0.8

    ./hydro paramfile.txt ./IC/advection-2D.dat
    errexit $?
    $plotdir/plot_density.py advection-2D-pwlin-0001.out
    errexit $?

    # negative velocity
    genparamfile 0 1 0 0 advection-2D-pwlin-negvel 0.8
    ./hydro paramfile.txt ./IC/advection-2D-negvel.dat
    errexit $?
    $plotdir/plot_density.py advection-2D-pwlin-negvel-0001.out
    errexit $?

    # x only
    genparamfile 0 1 0 0 advection-2D-pwlin-x 0.8
    ./hydro paramfile.txt ./IC/advection-2D-x.dat
    errexit $?
    $plotdir/plot_density.py advection-2D-pwlin-x-0001.out
    errexit $?

    # y only
    genparamfile 0 1 0 0 advection-2D-pwlin-y 0.8
    ./hydro paramfile.txt ./IC/advection-2D-y.dat
    errexit $?
    $plotdir/plot_density.py advection-2D-pwlin-y-0001.out
    errexit $?
}




advection_pwlin_limiters_1D(){
    #----------------------------
    # 1D pwlin with limiters
    #----------------------------

    for LIMITER in MINMOD SUPERBEE MC VANLEER; do

        # genmakefile ndim solver riemann limiter
        genmakefile 1 ADVECTION_PWLIN NONE $LIMITER
        make clean && make
        errexit $?

        # positive velocity
        # genparamfile nsteps tmax foutput dt_out basename ccfl
        genparamfile 0 10 0 1.0 advection-1D-$LIMITER 0.8

        ./hydro paramfile.txt ./IC/advection-1D.dat
        errexit $?
        $plotdir/plot_all_density.py advection-1D-$LIMITER-0000.out
        errexit $?
    done
}




advection_pwlin_limiters_2D(){
    #----------------------------
    # 2D pwlin with limiters
    #----------------------------

    for LIMITER in MINMOD SUPERBEE MC VANLEER; do

        # genmakefile ndim solver riemann limiter
        genmakefile 2 ADVECTION_PWLIN NONE $LIMITER
        make clean && make
        errexit $?

        # positive velocity
        # genparamfile nsteps tmax foutput dt_out basename ccfl
        genparamfile 0 1 0 0 advection-2D-$LIMITER 0.8

        ./hydro paramfile.txt ./IC/advection-2D.dat
        errexit $?
        $plotdir/plot_density.py advection-2D-$LIMITER-0001.out
        errexit $?

    done

}




advection_waf_1D(){

    #-------------------
    # 1D PWLIN
    #-------------------

    # genmakefile ndim solver riemann limiter
    genmakefile 1 ADVECTION_WAF NONE NONE
    make clean && make
    errexit $?

    # positive velocity
    # genparamfile nsteps tmax foutput dt_out basename ccfl
    genparamfile 0 10 0 1.0 advection-1D-waf 0.8

    ./hydro paramfile.txt ./IC/advection-1D.dat
    errexit $?
    $plotdir/plot_all_density.py advection-1D-waf-0000.out
    errexit $?

    # negative velocity
    genparamfile 0 10 0 1.0 advection-1D-waf-negvel 0.8
    ./hydro paramfile.txt ./IC/advection-1D-negvel.dat
    errexit $?
    $plotdir/plot_all_density.py advection-1D-waf-negvel-0000.out
    errexit $?

}



advection_waf_2D(){
    #-------------------
    # 2D WAF advection
    #-------------------

    # genmakefile ndim solver riemann limiter
    genmakefile 2 ADVECTION_WAF NONE NONE
    make clean && make
    errexit $?

    # positive velocity
    # genparamfile nsteps tmax foutput dt_out basename
    genparamfile 0 1 0 0 advection-2D-waf 0.8

    ./hydro paramfile.txt ./IC/advection-2D.dat
    errexit $?
    $plotdir/plot_density.py advection-2D-waf-0001.out
    errexit $?

    # negative velocity
    genparamfile 0 1 0 0 advection-2D-waf-negvel 0.8
    ./hydro paramfile.txt ./IC/advection-2D-negvel.dat
    errexit $?
    $plotdir/plot_density.py advection-2D-waf-negvel-0001.out
    errexit $?

    # x only
    genparamfile 0 1 0 0 advection-2D-waf-x 0.8
    ./hydro paramfile.txt ./IC/advection-2D-x.dat
    errexit $?
    $plotdir/plot_density.py advection-2D-waf-x-0001.out
    errexit $?

    # y only
    genparamfile 0 1 0 0 advection-2D-waf-y 0.8
    ./hydro paramfile.txt ./IC/advection-2D-y.dat
    errexit $?
    $plotdir/plot_density.py advection-2D-waf-y-0001.out
    errexit $?
}





advection_waf_limiters_1D(){
    #----------------------------------
    # 1D WAF advection with limiters
    #----------------------------------

    for LIMITER in MINMOD SUPERBEE MC VANLEER; do

        # genmakefile ndim solver riemann limiter
        genmakefile 1 ADVECTION_WAF NONE $LIMITER
        make clean && make
        errexit $?

        # positive velocity
        # genparamfile nsteps tmax foutput dt_out basename ccfl
        genparamfile 0 10 0 1.0 advection-waf-1D-$LIMITER 0.8

        ./hydro paramfile.txt ./IC/advection-1D.dat
        errexit $?
        $plotdir/plot_all_density.py advection-waf-1D-$LIMITER-0000.out
        errexit $?


        # negative velocity
        # genparamfile nsteps tmax foutput dt_out basename ccfl
        # genparamfile 0 10 0 1.0 advection-waf-1D-$LIMITER-negvel 0.8
        #
        # ./hydro paramfile.txt ./IC/advection-1D-negvel.dat
        # errexit $?
        # $plotdir/plot_all_density.py advection-waf-1D-$LIMITER-negvel-0000.out
        # errexit $?
    done
}





advection_waf_limiters_2D(){
    #-----------------------------------
    # 2D WAF advection with limiters
    #-----------------------------------

    for LIMITER in MINMOD SUPERBEE MC VANLEER; do

        # genmakefile ndim solver riemann limiter
        genmakefile 2 ADVECTION_WAF NONE $LIMITER
        make clean && make
        errexit $?

        # positive velocity
        # genparamfile nsteps tmax foutput dt_out basename ccfl
        genparamfile 0 1 0 0 advection-2D-waf-$LIMITER 0.8

        ./hydro paramfile.txt ./IC/advection-2D.dat
        errexit $?
        $plotdir/plot_density.py advection-2D-waf-$LIMITER-0001.out
        errexit $?

    done

}










riemann_vacuum(){
    #----------------------------
    # Test the vacuum solver
    #----------------------------

    for RIEMANN in EXACT; do

        # genmakefile ndim solver riemann limiter
        genmakefile 1 NONE $RIEMANN NONE
        make -f Makefile-Riemann clean && make -f Makefile-Riemann
        errexit $?

        # genparamfile nsteps tmax foutput dt_out basename ccfl
        genparamfile 0 0.01 0 0 "NO_BASENAME" 1
        ./riemann paramfile.txt ./IC/riemann-left-vacuum.dat
        errexit $?
        $plotdir/plot_riemann_result.py riemann-left-vacuum-*0001.out ./IC/riemann-left-vacuum.dat
        errexit $?

        genparamfile 0 0.01 0 0 "NO_BASENAME" 1
        ./riemann paramfile.txt ./IC/riemann-right-vacuum.dat
        errexit $?
        $plotdir/plot_riemann_result.py riemann-right-vacuum-*0001.out ./IC/riemann-right-vacuum.dat
        errexit $?

        genparamfile 0 0.01 0 0 "NO_BASENAME" 1
        ./riemann paramfile.txt ./IC/riemann-vacuum-generating.dat
        errexit $?
        $plotdir/plot_riemann_result.py riemann-vacuum-generating*0001.out ./IC/riemann-vacuum-generating.dat
        errexit $?

    done

}






riemann_solver(){
    #-----------------------------
    # Test various Riemann solvers
    #-----------------------------

    for RIEMANN in EXACT; do

        # genmakefile ndim solver riemann limiter
        genmakefile 1 NONE $RIEMANN NONE
        make -f Makefile-Riemann clean && make -f Makefile-Riemann
        errexit $?

        for icprefix in riemann-sod-shock riemann-sod-shock-reverse; do
            # genparamfile nsteps tmax foutput dt_out basename ccfl
            genparamfile 0 0.25 0 0 "NO_BASENAME" 1
            ./riemann paramfile.txt ./IC/"$icprefix".dat
            errexit $?

            $plotdir/plot_riemann_result.py "$icprefix"-RIEMANN-EXACT-0001.out ./IC/"$icprefix".dat
            errexit $?
        done

    done
}












function godunov_1D() {
    #-------------------------------------------
    # 1D godunov stuff
    #-------------------------------------------


    # for RIEMANN in EXACT TSRS; do
    for RIEMANN in EXACT TRRS TSRS HLLC; do

        # genmakefile ndim solver riemann limiter
        genmakefile 1 GODUNOV $RIEMANN NONE
        make clean && make
        errexit $?

        # non-vacuum
        for icprefix in riemann-sod-shock riemann-sod-shock-reverse; do
            # genparamfile_transmissive nsteps tmax foutput dt_out basename ccfl
            genparamfile_transmissive 0 0.2 0 0 "$icprefix"-GODUNOV-1D-$RIEMANN 0.9

            ./hydro paramfile.txt ./IC/"$icprefix".dat
            errexit $?
        done

        # vacuum
        for icprefix in riemann-left-vacuum riemann-right-vacuum riemann-vacuum-generating; do
            genparamfile_transmissive 0 0.01 0 0 "$icprefix"-GODUNOV-1D-$RIEMANN 0.9

            ./hydro paramfile.txt ./IC/"$icprefix".dat
            errexit $?
        done
    done

    for icprefix in riemann-left-vacuum riemann-right-vacuum riemann-vacuum-generating riemann-sod-shock riemann-sod-shock-reverse; do
        ./overplot_riemann_solvers.py "$icprefix" GODUNOV-1D junk
        errexit $?
    done
}








function godunov_2D() {
    #-------------------------------------------
    # 2D godunov stuff
    #-------------------------------------------


    for RIEMANN in EXACT TRRS TSRS HLLC; do

        # genmakefile ndim solver riemann limiter
        genmakefile 2 GODUNOV $RIEMANN NONE
        make clean && make
        errexit $?

        # non vacuum
        for icprefix in riemann-sod-shock riemann-sod-shock-reverse; do
            # genparamfile_transmissive nsteps tmax foutput dt_out basename ccfl
            genparamfile_transmissive 0 0.2 0 0 "$icprefix"-GODUNOV-2D-$RIEMANN 0.9

            ./hydro paramfile.txt ./IC/"$icprefix".dat
            errexit $?
        done

        # vacuum
        for icprefix in riemann-left-vacuum riemann-right-vacuum riemann-vacuum-generating; do
            genparamfile_transmissive 0 0.01 0 0 "$icprefix"-GODUNOV-2D-$RIEMANN 0.9

            ./hydro paramfile.txt ./IC/"$icprefix".dat
            errexit $?
        done
    done

    for icprefix in riemann-left-vacuum riemann-right-vacuum riemann-vacuum-generating riemann-sod-shock riemann-sod-shock-reverse; do
        ./overplot_riemann_solvers.py "$icprefix" GODUNOV-2D junk
        errexit $?
    done


    # do kelvin helmholtz only for one
    genmakefile 2 GODUNOV HLLC NONE
    make clean && make
    errexit $?

    # genparamfile nsteps tmax foutput dt_out basename ccfl
    genparamfile 0 3 0 0 "GODUNOV-kelvin-helmholtz" 0.8
    ./hydro paramfile.txt ./IC/kelvin-helmholtz-128.dat
    $plotdir/plot_result.py GODUNOV-kelvin-helmholtz-0001.out
    errexit $?

}






function waf_1D() {
    #-------------------------------------------
    # 1D WAF stuff without limiters
    #-------------------------------------------

    # this one needs a bit of special attention because the oscillations can crash the code

    for RIEMANN in EXACT TRRS TSRS HLLC; do
        # genmakefile ndim solver riemann limiter
        genmakefile 1 WAF $RIEMANN NONE
        make clean && make
        errexit $?

        # non-vacuum
        for icprefix in riemann-sod-shock riemann-sod-shock-reverse; do
            # genparamfile_transmissive nsteps tmax foutput dt_out basename ccfl
            genparamfile_transmissive 0 0.1 0 0 "$icprefix"-WAF-1D-$RIEMANN-NO_LIMITER 0.5

            ./hydro paramfile.txt ./IC/"$icprefix".dat
            errexit $?
        done

        # vacuum

        # !!!!!!!!!!!! WAF without limiters introduces strong oscillations, and the code can't finish.

        # for icprefix in riemann-left-vacuum riemann-right-vacuum riemann-vacuum-generating; do
        #     genparamfile_transmissive 0 0.005 0 0 "$icprefix"-WAF-1D-$RIEMANN-NO_LIMITER 0.2
        #
        #     ./hydro paramfile.txt ./IC/"$icprefix".dat
        #     errexit $?
        # done

    done

    for icprefix in riemann-sod-shock riemann-sod-shock-reverse; do
        ./overplot_riemann_solvers.py "$icprefix" WAF-1D NO_LIMITER
        errexit $?
    done

}






function waf_1D_limiters() {
    #-------------------------------------------
    # 1D WAF stuff with limiters
    #-------------------------------------------

    for LIMITER in MINMOD SUPERBEE MC VANLEER; do

        for RIEMANN in EXACT TRRS TSRS HLLC; do
            # genmakefile ndim solver riemann limiter
            genmakefile 1 WAF $RIEMANN $LIMITER
            make clean && make
            errexit $?

            # non-vacuum
            for icprefix in riemann-sod-shock riemann-sod-shock-reverse; do
                # genparamfile_transmissive nsteps tmax foutput dt_out basename ccfl
                genparamfile_transmissive 0 0.2 0 0 "$icprefix"-WAF-1D-$RIEMANN-$LIMITER 0.9

                ./hydro paramfile.txt ./IC/"$icprefix".dat
                errexit $?
            done

            # vacuum
            for icprefix in riemann-left-vacuum riemann-right-vacuum riemann-vacuum-generating; do
                genparamfile_transmissive 0 0.01 0 0 "$icprefix"-WAF-1D-$RIEMANN-$LIMITER 0.9

                ./hydro paramfile.txt ./IC/"$icprefix".dat
                errexit $?
            done

        done

        # plotting
        for icprefix in riemann-left-vacuum riemann-right-vacuum riemann-vacuum-generating riemann-sod-shock riemann-sod-shock-reverse; do
            ./overplot_riemann_solvers.py "$icprefix" WAF-1D $LIMITER
            errexit $?
        done

    done

}







function waf_2D_limiters() {
    #-------------------------------------------
    # 2D WAF stuff with limiters
    #-------------------------------------------

    for LIMITER in MINMOD SUPERBEE MC VANLEER; do

        for RIEMANN in EXACT TRRS TSRS HLLC; do
            # genmakefile ndim solver riemann limiter
            genmakefile 2 WAF $RIEMANN $LIMITER
            make clean && make
            errexit $?

            # non-vacuum
            for icprefix in riemann-sod-shock riemann-sod-shock-reverse; do
                # genparamfile_transmissive nsteps tmax foutput dt_out basename ccfl
                genparamfile_transmissive 0 0.2 0 0 "$icprefix"-WAF-2D-$RIEMANN-$LIMITER 0.9

                ./hydro paramfile.txt ./IC/"$icprefix".dat
                errexit $?
            done
        done

        for icprefix in riemann-sod-shock riemann-sod-shock-reverse; do
            ./overplot_riemann_solvers.py "$icprefix" WAF-2D $LIMITER
            errexit $?
        done

    done

    # do kelvin helmholtz only for one
    genmakefile 2 WAF HLLC VANLEER
    make clean && make
    errexit $?

    # genparamfile nsteps tmax foutput dt_out basename ccfl
    genparamfile 0 3 0 0 "WAF-kelvin-helmholtz" 0.9
    ./hydro paramfile.txt ./IC/kelvin-helmholtz-128.dat
    $plotdir/plot_result.py WAF-kelvin-helmholtz-0001.out
    errexit $?

}







function muscl_1D(){
    #-------------------------------------------
    # 1D MUSCL stuff without limiters
    #-------------------------------------------

    # this one needs a bit of special attention because the oscillations can crash the code

    for RIEMANN in EXACT TSRS HLLC; do
        # genmakefile ndim solver riemann limiter
        genmakefile 1 MUSCL $RIEMANN NONE
        make clean && make
        errexit $?

        # non-vacuum
        for icprefix in riemann-sod-shock riemann-sod-shock-reverse; do
            # genparamfile_transmissive nsteps tmax foutput dt_out basename ccfl
            genparamfile_transmissive 0 0.05 0 0 "$icprefix"-MUSCL-1D-$RIEMANN-NO_LIMITER 0.7

            ./hydro paramfile.txt ./IC/"$icprefix".dat
            errexit $?
        done

        # vacuum

        # !!!!!!!!!!!! MUSCL without limiters introduces strong oscillations, and the code can't finish.

        # for icprefix in riemann-left-vacuum riemann-right-vacuum riemann-vacuum-generating; do
        #     genparamfile_transmissive 0 0.005 0 0 "$icprefix"-MUSCL-1D-$RIEMANN-NO_LIMITER 0.2
        #
        #     ./hydro paramfile.txt ./IC/"$icprefix".dat
        #     errexit $?
        # done

    done

    for icprefix in riemann-sod-shock riemann-sod-shock-reverse; do
        ./overplot_riemann_solvers.py "$icprefix" MUSCL-1D NO_LIMITER
        errexit $?
    done


}





function muscl_1D_limiters() {
    #-------------------------------------------
    # 1D MUSCL stuff with limiters
    #-------------------------------------------

    for LIMITER in MINMOD SUPERBEE VANLEER; do

        for RIEMANN in EXACT TSRS HLLC; do
            # genmakefile ndim solver riemann limiter
            genmakefile 1 MUSCL $RIEMANN $LIMITER
            make clean && make
            errexit $?

            # non-vacuum
            for icprefix in riemann-sod-shock riemann-sod-shock-reverse; do
                # genparamfile_transmissive nsteps tmax foutput dt_out basename ccfl
                genparamfile_transmissive 0 0.2 0 0 "$icprefix"-MUSCL-1D-$RIEMANN-$LIMITER 0.8

                ./hydro paramfile.txt ./IC/"$icprefix".dat
                errexit $?
            done

            # vacuum

            # vacuum isn't well handled with MUSCL. Skip the tests.

            # for icprefix in riemann-left-vacuum riemann-right-vacuum riemann-vacuum-generating; do
            #     genparamfile_transmissive 0 0.01 0 0 "$icprefix"-MUSCL-1D-$RIEMANN-$LIMITER 0.5
            #
            #     ./hydro paramfile.txt ./IC/"$icprefix".dat
            #     errexit $?
            # done

        done

        # plotting
        for icprefix in riemann-sod-shock riemann-sod-shock-reverse; do
        # for icprefix in riemann-left-vacuum riemann-right-vacuum riemann-vacuum-generating riemann-sod-shock riemann-sod-shock-reverse; do
            ./overplot_riemann_solvers.py "$icprefix" MUSCL-1D $LIMITER
            errexit $?
        done

    done

}







function muscl_2D_limiters() {
    #-------------------------------------------
    # 2D MUSCL stuff with limiters
    #-------------------------------------------

    for LIMITER in MINMOD SUPERBEE VANLEER; do

        for RIEMANN in EXACT TSRS HLLC; do
            # genmakefile ndim solver riemann limiter
            genmakefile 2 MUSCL $RIEMANN $LIMITER
            make clean && make
            errexit $?

            # non-vacuum
            for icprefix in riemann-sod-shock riemann-sod-shock-reverse; do
                # genparamfile_transmissive nsteps tmax foutput dt_out basename ccfl
                genparamfile_transmissive 0 0.2 0 0 "$icprefix"-MUSCL-2D-$RIEMANN-$LIMITER 0.9

                ./hydro paramfile.txt ./IC/"$icprefix".dat
                errexit $?
            done
        done

        for icprefix in riemann-sod-shock riemann-sod-shock-reverse; do
            ./overplot_riemann_solvers.py "$icprefix" MUSCL-2D $LIMITER
            errexit $?
        done

    done

    # do kelvin helmholtz only for one
    genmakefile 2 MUSCL HLLC VANLEER
    make clean && make
    errexit $?

    # genparamfile nsteps tmax foutput dt_out basename ccfl
    genparamfile 0 3 0 0 "MUSCL-kelvin-helmholtz" 0.9
    ./hydro paramfile.txt ./IC/kelvin-helmholtz-128.dat
    $plotdir/plot_result.py MUSCL-kelvin-helmholtz-0001.out
    errexit $?

}






function sources() {
    #-------------------------------------------
    # tests for sources and integrators
    #-------------------------------------------


    for SOURCE in RADIAL CONSTANT; do
        for INTEGRATOR in RK2 RK4; do

            for NDIM in 1 2; do

                for SOLVER in WAF GODUNOV MUSCL; do

                    genmakefile_sources $NDIM $SOLVER $SOURCE $INTEGRATOR
                    make clean && make
                    errexit $?

                    basename=SOURCES-$SOLVER-$SOURCE-$INTEGRATOR-"$NDIM"D

                    genparamfile_sources 0 0.5 $basename
                    ./hydro paramfile.txt ./IC/uniform-"$NDIM"D.dat
                    errexit $?

                    if [ $NDIM == 2 ]; then
                        $plotdir/plot_result.py $basename-0001.out
                    fi
                done
            done
        done
        ./overplot_fname_label.py SOURCES-*"$SOURCE"*-1D-0001.out SOURCES-$SOURCE-1D.png
    done


}










#=====================================
# Now actually run the tests
#=====================================

# first check whether we can even plot stuff

advection_pwconst_1D
advection_pwconst_2D

advection_pwlin_1D
advection_pwlin_2D
advection_pwlin_limiters_1D
advection_pwlin_limiters_2D

advection_waf_1D
advection_waf_2D
advection_waf_limiters_1D
advection_waf_limiters_2D

riemann_vacuum
riemann_solver

godunov_1D
godunov_2D

waf_1D
waf_1D_limiters
waf_2D_limiters

muscl_1D
muscl_1D_limiters
muscl_2D_limiters

sources


#=====================================
# create TeX
#=====================================
generate_tex
