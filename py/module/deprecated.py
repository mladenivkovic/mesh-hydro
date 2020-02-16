#!/usr/bin/env python3

#---------------------------------------
# Stuff that shouldn't be used any more.
#---------------------------------------




def read_arbitrary_ic_metadata(fname):
    """
    Read nx, ndim for arbitrary file format of IC files.
    """

    f = open(fname)

    nx = None
    ndim = None

    linecount = 0
    while True:
        # safety measure
        linecount += 1
        if linecount > 1000:
            print("================ Got to linecount = 1000 without having filetype found. Wtf? Exiting now.")
            quit(1)

        line = f.readline()
        clean = remove_C_style_comments(line)
        print("line:", line)
        print("clean:", clean)
        if line_is_empty(clean):
            #  print("clean is empty")
            continue

        else:
            # filetype MUST be first non-comment non-empty line
            # close file here, since we exit in any case
            name, eq, value = clean.partition("=")
            nstr = name.strip()
            if nstr == "filetype":
                continue
            elif nstr == "nx":
                nx = int(value)
            elif nstr == "ndim":
                ndim = int(value)
            else:
                raise ValueError("Unknown name: '{0}'".format(name))

        if nx is not None and ndim is not None:
            break

    f.close()


    return nx, ndim
