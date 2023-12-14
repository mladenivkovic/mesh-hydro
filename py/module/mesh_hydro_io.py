#!/usr/bin/env python3

# -------------------------------------------------------
# Module that contains io-routines for the hydro outputs
# and IC files.
# -------------------------------------------------------


import copy
import numpy as np


def write_ic(fname, ndim, rho, u, p):
    """
    Write an (arbitrary type) IC file.
    fname:  filename to be written
    ndim:   number of dimensions
    rho:    numpy array for density
    u:      numpy array for velocity
    p:      numpy array for pressure

    returns:
        Nothing
    """

    f = open(fname, "w")
    nx = rho.shape[0]

    f.write("filetype = arbitrary\n")
    f.write("ndim = {0:d}\n".format(ndim))
    f.write("nx = {0:d}\n".format(nx))
    f.write("\n")

    if ndim == 1:
        for i in range(nx):
            f.write("{0:12.6f} {1:12.6f} {2:12.6f}\n".format(rho[i], u[i], p[i]))

    elif ndim == 2:
        for j in range(nx):
            for i in range(nx):
                f.write(
                    "{0:12.6f} {1:12.6f} {2:12.6f} {3:12.6f}\n".format(
                        rho[i, j], u[i, j, 0], u[i, j, 1], p[i, j]
                    )
                )

    return


def read_output(fname):
    """
    Read the given output file.
    returns:
        ndim:       integer of how many dimensions we have
        rho:        numpy array for density
        u:          numpy array for velocity. In 1D: is 1D array. In 2D: is 2D array containing
                    both ux and uy
        p:          numpy array for pressure
        t:          time of the output
        step:       current step of the simulation
    """

    check_file_exists(fname)

    f = open(fname)

    nx = None
    ndim = None
    t = None
    step = None

    linecount = 0
    while True:
        # safety measure
        linecount += 1
        if linecount > 1000:
            print(
                "================ Got to linecount = 1000 without having found all metadata. Wtf? Exiting now."
            )
            print("got nx:", nx)
            print("got ndim", ndim)
            print("got t:", t)
            print("got step:", step)
            quit(1)

        line = f.readline()
        clean = remove_python_style_comments(line)
        if line_is_empty(clean):
            continue

        else:
            if clean.strip().startswith("x"):  # header name descriptions
                continue
            name, eq, value = clean.partition("=")
            nstr = name.strip()
            if nstr == "nx":
                nx = int(value)
            elif nstr == "ndim":
                ndim = int(value)
            elif nstr == "t":
                t = float(value)
            elif nstr == "nsteps":
                step = int(value)
            else:
                raise ValueError("Unknown name: '{0}'".format(name))

        if nx is not None and ndim is not None and t is not None and step is not None:
            break

    f.close()

    if ndim == 1:
        rho, u, p = np.loadtxt(
            fname, usecols=[1, 2, 3], dtype=np.float, unpack=True, skiprows=linecount
        )

    elif ndim == 2:
        rho, ux, uy, p = np.loadtxt(
            fname, usecols=[2, 3, 4, 5], dtype=np.float, unpack=True, skiprows=linecount
        )

        rho = rho.reshape((nx, nx))
        p = p.reshape((nx, nx))
        ux = ux.reshape((nx, nx))
        uy = uy.reshape((nx, nx))
        u = np.stack((ux, uy), axis=2)

    return ndim, rho, u, p, t, step


def read_ic(fname, nx=100):
    """
    Top-level function to read in the given IC file. File is passed as string fname.
    It figures out the dimensions etc by itself. Returns the relevant data as numpy
    arrays, which can be either 1D or 2D, depending on IC.

    If the file type is two-state type, it will generate a 1D numpy array with nx cells.

    returns:
        ndim:       integer of how many dimensions we have
        twostate:   if this is a two-state type IC file
        rho:        numpy array for density
        u:          numpy array for velocity. In 1D: is 1D array. In 2D: is 2D array containing
                    both ux and uy
        p:          numpy array for pressure
    """

    check_file_exists(fname)
    twostate = get_ic_filetype(fname)

    if twostate:
        ndim = 1
        rho, u, p = read_twostate_ic(fname, nx)

    else:
        ndim, rho, u, p = read_arbitrary_ic(fname)

    return ndim, twostate, rho, u, p


def read_arbitrary_ic(fname):
    """
    read the complete arbitrary format IC file
    """

    f = open(fname)
    data = f.readlines()
    f.close()

    got_ftype = False
    got_nx = False
    got_ndim = False
    got_header = False

    i = 0
    j = 0

    for line in data:
        clean = remove_C_style_comments(line)
        clean = remove_newline(clean)
        if line_is_empty(clean):
            continue

        if not got_header:
            name, eq, value = clean.partition("=")
            name = name.strip()
            if name == "filetype":
                got_ftype = True
                continue
            elif name == "ndim":
                ndim = int(value)
                got_ndim = True
            elif name == "nx":
                nx = int(value)
                got_nx = True
            else:
                print("Unrecognized value name:", name)

            got_header = got_ftype and got_nx and got_ndim
            if got_header:
                if ndim == 1:
                    rho = np.empty((nx), dtype=np.float)
                    u = np.empty((nx), dtype=np.float)
                    p = np.empty((nx), dtype=np.float)

                elif ndim == 2:
                    rho = np.empty((nx, nx), dtype=np.float)
                    u = np.empty((nx, nx, 2), dtype=np.float)
                    p = np.empty((nx, nx), dtype=np.float)

        else:
            vals = split_columns(clean)

            if ndim == 1:
                if len(vals) != 3:
                    print("Got wrong number of values in the line.")
                    print("Line was: ", line, end="")
                    print("Cleaned line is: '{0}'".format(clean))
                    print("I expect 3 values")
                    print("I got:", len(vals), vals)
                    quit(1)

                rho[i] = float(vals[0])
                u[i] = float(vals[1])
                p[i] = float(vals[2])
                i += 1

            elif ndim == 2:
                if len(vals) != 4:
                    print("Got wrong number of values in the line.")
                    print("Line was: ", line, end="")
                    print("Cleaned line is: '{0}'".format(clean))
                    print("I expect 4 values")
                    print("I got:", len(vals), vals)
                    quit(1)

                rho[j, i] = float(vals[0])
                u[j, i, 0] = float(vals[1])
                u[j, i, 1] = float(vals[2])
                p[j, i] = float(vals[3])
                i += 1
                if i == nx:
                    i = 0
                    j += 1

    # checks
    if ndim == 1:
        if i != nx:
            print("Got too few values in x direction. Got i=", i, "should be", nx)
            quit(1)
    if ndim == 2:
        if j != nx:
            print("Got too few values in y direction. Got i=", j, "should be", nx)
            quit(1)
        if i != 0:
            print("Got too few values in y direction. Got i=", i, "should be", 0)
            quit(1)

    return ndim, rho, u, p


def read_twostate_ic(fname, nx):
    """
    Read complete two-state format style.
    Return rho, u, p arrays with nx elements.
    """

    rhoL = None
    uL = None
    pL = None
    rhoR = None
    uR = None
    pR = None

    f = open(fname)
    data = f.readlines()
    f.close()

    for line in data:
        clean = remove_C_style_comments(line)
        clean = remove_newline(clean)
        if line_is_empty(clean):
            continue
        name, eq, value = clean.partition("=")
        name = name.strip()
        if name == "rho_L":
            rhoL = float(value)
        elif name == "u_L":
            uL = float(value)
        elif name == "p_L":
            pL = float(value)
        elif name == "rho_R":
            rhoR = float(value)
        elif name == "u_R":
            uR = float(value)
        elif name == "p_R":
            pR = float(value)
        elif name == "filetype":
            continue
        else:
            print("Unrecognized value name:", name)

    # check that we got everything
    for val, name in [
        (rhoL, "rhoL"),
        (uL, "uL"),
        (pL, "pL"),
        (rhoR, "rhoR"),
        (uR, "uR"),
        (pR, "pR"),
    ]:
        if val is None:
            print("Finished reading file, I got no value for", name)

    # now allocate rho, u, p arrays

    rho = np.empty((nx), dtype=np.float)
    u = np.empty((nx), dtype=np.float)
    p = np.empty((nx), dtype=np.float)

    nxhalf = nx // 2
    rho[:nxhalf] = rhoL
    u[:nxhalf] = uL
    p[:nxhalf] = pL
    rho[nxhalf:] = rhoR
    u[nxhalf:] = uR
    p[nxhalf:] = pR

    return rho, u, p


def get_ic_filetype(fname):
    """
    Get the filetype of the IC file. Returns True if two-state style file, False if arbitrary.
    """

    f = open(fname)

    linecount = 0
    while True:
        # safety measure
        linecount += 1
        if linecount > 1000:
            print(
                "================ Got to linecount = 1000 without having filetype found. Wtf? Exiting now."
            )
            quit(1)

        line = f.readline()
        clean = remove_C_style_comments(line)
        clean = remove_newline(clean)
        if line_is_empty(clean):
            continue

        else:
            # filetype MUST be first non-comment non-empty line
            # close file here, since we exit in any case
            f.close()
            ftstr, eq, ftype = clean.partition("=")
            if ftstr.strip() != "filetype":
                raise ValueError(
                    "First non-empty non-comment line is '{0}', but is should contain the filetype.".format(
                        clean
                    )
                )
            ftclean = ftype.strip()
            if ftclean == "two-state":
                return True
            elif ftclean == "arbitrary":
                return False
            else:
                raise ValueError("Unknown file type '{0}'".format(ftclean))


def remove_python_style_comments(line):
    """
    Remove # if the line starts with it
    """

    if line[0] == "#":
        return line[1:].strip()

    return line


def remove_C_style_comments(line):
    """
    Remove all comments (//, /* .. */) from the line
    actually just look for a slash. It has no other business in there.
    """

    clean = line

    for i, c in enumerate(line):
        if c == "/":
            clean = line[:i]
            break

    return clean


def line_is_empty(line):
    """
    Check whether line contains only spaces and/or newline chars
    """
    clean = line.strip()
    if len(clean) == 0:
        return True
    else:
        return False


def check_file_exists(fname):
    """
    Check that file exists, throw error if not.
    """
    import os

    if not os.path.isfile(fname):
        raise ValueError("Given file {0} doesn't exist.".format(fname))
    return


def remove_newline(line):
    """
    If newline character is the last character of the line, remove it.
    """
    if len(line) > 0:
        if line[-1] == "\n" or line[-1] == "\r":
            return (line[:-1]).strip()

    return line


def split_columns(line, delim=" "):
    """
    Split given line (string) into columns defined by the delimiter delim.
    If the delimiter occurs multiple times back - to - back, treat it as
    only one column delimiter.

    returns:
        splits: List of strings, representing columns

    """

    splits = []

    start = 0
    stop = 0

    while stop < len(line) and start < len(line):
        if line[start] == delim:
            start += 1
            continue

        if line[stop] != delim:
            if line[stop] == "\n" or line[stop] == "\r":
                # add value without newline
                splits.append(line[start:stop])
                start = stop
                break
            else:
                stop += 1
                continue

        else:  # if line[stop] == delim
            splits.append(line[start:stop])
            while stop < len(line):
                if line[stop] == delim:
                    stop += 1
                else:
                    if line[stop] == "\n" or line[stop] == "\r":
                        # add value without newline
                        splits.append(line[start:stop])
                        start = stop
                    break
            start = stop

    # add final column
    if start != stop:
        splits.append(line[start:])

    return splits
