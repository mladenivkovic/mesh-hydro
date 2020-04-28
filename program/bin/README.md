bin directory
====================================


Contents:

    - ``arbitrary-1D.dat``: An example IC file of the arbitrary file format in 1D.
    - ``arbitrary-2D.dat``: An example IC file of the arbitrary file format in 2D.
    - ``example-paramfile.txt``: An example parameter file, contains all possible parameters and explaination of them.
    - ``example_output_times_file.txt``: An example output times file if the ``toutlist`` parameter is used.
    - ``generate-all-execs.sh``: A script to compile all possible combinations of solvers, limiters, dimensions, etc. Every executable file will be named recognizably and differenlty.
    - ``Makefile``: The Makefile used to compile with ``run.sh``, or just your default ``make clean && make``. You can change definitions at the top of the file.
    - ``Makefile-Riemann``: The Makefile used to compile the code to use it as a Riemann solver.
    - ``paramfile.txt``: A ready-to-go parameter file. Set up such that when running vanilla ``./run.sh``, it produces something nice.
    - ``processing.mk``: Some processing rules employed in all makefiles, written to be included in them.
    - ``run.sh``: A script to automatically compile, run, and plot results. Some possibilities are commented out. Run with ``-d`` to run the program with ``gdb``
    - ``twostate.dat``: An example IC file of the two-state format.
