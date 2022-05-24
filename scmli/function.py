#!/usr/bin/python3


def print_c (a):
    print(a)

def pcr_pipline(r1,r2):
    import os
    os.system('fastqc -o ./ ' + r1 + ' ' + r2)
