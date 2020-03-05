#!/usr/bin/env python3

#------------------------------------------------------
# Check whether the module files can be imported.
# If not, tell me which directory I should add to my
# pythonpath.
#------------------------------------------------------



def try_to_import():
    """
    Try to import the hydro modules.
    Tell me where the module path is if import fails.
    """
    try:
        import hydro_io
        import hydro_plotting
        import hydro_riemann
        import hydro_utils

    except ImportError:

        import os
        cwd = os.getcwd()
        pre, mesh, post = cwd.partition("mesh-hydro")
        meshdir = os.path.join(pre, mesh)
        moddir = os.path.join(meshdir, "py", "module")

        print("ERROR: Couldn't find the hydro plotting module files.")
        print("ERROR: To use the python plotting, evaluation, or IC generating scripts,")
        print("ERROR: You first need to add the directory to your PYTHONPATH.")
        print("ERROR: The module is located at:")
        print("ERROR:      ", moddir)
        print("ERROR:")
        print("ERROR: If you are using bash, you could for example add")
        print('ERROR: export PYTHONPATH="$PYTHONPATH":{0}'.format(moddir))
        print("ERROR: to your .bashrc file.")

        quit(2)


if __name__ == "__main__":
    try_to_import()
