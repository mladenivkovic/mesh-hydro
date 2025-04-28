import os

DIR_OF_THIS_SCRIPT = os.path.abspath( os.path.dirname( __file__ ) )


def _add_include(path, include="include"):
    """
    Add the "include" subdir to a path. E.g.
    /some/path -> /some/path/include
    """
    return os.path.join(os.path.abspath(path), include)


def _parse_defines_mk():
    """
    Parse the defines from the defines.mk file, if it exists.
    """

    defines = []

    makefile = os.path.join(DIR_OF_THIS_SCRIPT, "program/bin/defines.mk")
    if not os.path.exists(makefile):
        return defines

    f = open(makefile, "r")
    makelines = f.readlines()
    for line in makelines:
        if not "=" in line:
            continue
        stripped = line.strip()
        defline = stripped.replace(" ", "")
        print(defline)

        defines.append("-D"+defline)

    f.close()

    return defines




flags = [
    '-Wall',
    '-Wextra',
    '-Wno-unused-parameter'
    '-pedantic',
    '-x', 'c',
]

defines = [
    #  "-DWITH_MPI",
    ]

include = [
    #  "-I", ".",
    #  "-I", "..",                                         # for project builds
    #  "-I", DIR_OF_THIS_SCRIPT,                           # add config.h
    "-I", os.path.join(DIR_OF_THIS_SCRIPT, "program/src"),
        ]

libs = [
    #  MPI_ROOT,
    #  HDF5_ROOT,
        ]



def Settings( **kwargs ):

    if kwargs[ 'language' ] == 'cfamily':

        all_includes = include

        fname = kwargs[ "filename" ]
        absfname = os.path.abspath(fname)
        fdir = os.path.dirname(absfname)

        for lib in libs:
            all_includes.append("-I")
            all_includes.append(_add_include(lib))

        all_defines = defines + _parse_defines_mk()

        final_flags = flags + all_defines + all_includes

        return {
            'flags': final_flags,
          }
