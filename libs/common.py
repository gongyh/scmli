import sys
import os

def existing_file(fname):
    fname = os.path.realpath(fname)
    if os.path.isfile(fname):
        return fname
    else:
        raise TypeError('not a valid file "%s"' %fname)

def existing_dir(dname):
    dname = os.path.realpath(dname)
    if os.path.isdir(dname):
        return dname
    else:
        raise TypeError('not a valid directory "%s"' %dname)
