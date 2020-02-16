#!/usr/bin/env python3

#----------------------------
# Various small utilities
#----------------------------


from sys import argv
import os

def get_only_cmdlinearg():
    """
    Get the only expected cmdline argument. If there is none
    or too many, throw error.
    """

    argc = len(argv)
    if argc != 2:
        print(argc, argv)
        raise ValueError("I expect exactly one argument: The filename.")
    else:
        return argv[1]



def get_all_files_with_same_basename(fname):
    """
    Get a list of all files with the same basename as given file <fname>.
    Basename in this case means everything before _XXXX.out, which is the
    format of the hydro output files.
    """

    # first generate basename
    basename = fname[:-9]
    
    basedir = os.path.dirname(basename)
    if basedir == '':
        basedir = os.getcwd()

    allfiles = os.listdir(basedir)

    filelist = []

    for f in allfiles:
        if f.startswith(basename) and f.endswith(".out"):
            filelist.append(f)

    filelist.sort()

    return filelist
