#!/usr/bin/env python3
# encoding: utf-8
# pae3.py
# juanfc 2019-07-08
# builds gnuplot command line to plot the contents of a directory
# pae3.py directory [-s [0;1;4]] [-f "[30:60;2:20;50]"]

import sys
import os
from pathlib import Path
import argparse

BASEDIR = "results"

def buildList(pyExpr):
    r = []
    tmp = pyExpr.split(";")
    for item in tmp:
        if item.count(':') == 1:
            first, last = map(int,item.split(':'))
            r +=  list(map(str,list(range(first, last))))
        elif item.count(':') == 2:
            first, last, step = map(int,item.split(':'))
            r +=  list(map(str,list(range(first, last, step))))
        else:
            r.append(item)

    return r

def fileNumbering(n):
    return "%04d" % n

def plotall(n):
    s = "gnuplot -persist -e 'plot "
    for i in range(1,n):
        fname = fileNumbering(i)
        s += '"results/'+inDir+"/"+fname+'.txt"'+ ' using 11 with lines title "' + fname +'",'
    s = s.rstrip(',') + "'"
    print(s)
    os.system(s)
    # os.system("gnuplot -persist -e 'plot \"results/%s\" using 11 with lines title \"%s\"'" % (ffname, fname))


# #### #
# MAIN #
# #### #

if "-h" == sys.argv[1]:
    print("""Examples:
          pae3.py directory [-s [0;1;4]] [-f "[30:60;2:20;50]"]


          You must provided the input directory
          if it it starts in / it'd be an absolute path
          in other case, it will be searched inside 'results' directory

          IMPORTANT:
              Observe in the examples:
              To express a list of different values or ranges, separate them with ;
              but to prevent the Terminal shell from an undesired interpretation, surround it in ';' or ";"

          FIRST SPECIES IS NUMBER 0, LAST: THE TOTAL - 1
    """)
    exit(0)




inDir = None
testMode = False
species = None
files = None

i = 1
while i in range(1, len(sys.argv)):

    if "-t" == sys.argv[i]:
        testMode = True

    elif "-s" == sys.argv[i]:
        species = buildList(sys.argv[i+1])

    elif "-f" == sys.argv[i]:
        files = buildList(sys.argv[i+1])

    else inDir = sys.argv[i]

    i += 1

commandList = build("ae3.py ", ranges, " ")

#
# DIRECTORY LISTING OF FILES pfiles
#

if not inDir:
    print("You must give an input directory in parameters.")
    print("If the dir starts in / it'd be an absolute path")
    print("in other case, it will be created inside 'results' directory")
    print("pae3.py directory [-s [0;1;4]] [-f \"[30:60;2:20;50]\"]")
    exit(1)

if inDir.startswith('/'):
    inDir = Path(inDir)
else:
    inDir = Path(BASEDIR) / inDir

if not inDir.is_dir():
    print("You must an existing input directory in parameters.")
    print("The one provided: '", inDir, "'", sep="")
    print("does not exists.")
    exit(1)

prevDir = Path()
os.chdir(inDir)
inDir = Path()

if files:
    pfiles = buildFileList(inDir, files)
else:
    pfiles = inDir.glob("*.txt")

non_existing = list(filter(lambda x: not x.exists(), pfiles))
if non_existing:
    print("ERROR: NEXT FILES DO NOT EXIST:", non_existing)
    exit(1)

#
# species
#
actualNumberOfSpeciesInFiles = countNOfSpecies(pfiles)
if species:
    if list(filter(lambda x: x >= actualNumberOfSpeciesInFiles, species)):
        print("ERROR: Number of species in those files is: ", actualNumberOfSpeciesInFiles)
        print("You are asking to print a specie not there")
        exit(1)
else:
    species = list(range(actualNumberOfSpeciesInFiles))

plotall(pfiles, species)

def plotall(pfiles, species):
    s = "gnuplot -persist -e 'plot "
    sspecies = # 11+13+15+17, +14
    for file in pfiles:
        fname = fileNumbering(i)
        for nspecies in species:
            s += str(file) + stringSp(nspecies) + " title " + str(file) + ","
    s = s.rstrip(',') + "'"
    print(s)
    os.system(s)

# plot 'results/0001.txt' using 11 with lines,12 with lines title 'SP Sharksucker'
# plot "results/PR.txt" using ($11)+($12) with lines title "SP Una", "results/TR.txt" using ($11)+($12) with lines title "SP dos"
# def plotall(n):
    # s = "gnuplot -persist -e 'plot "
    # for i in range(1,n):
    #     fname = fileNumbering(i)
    #     s += '"results/'+outDir+"/"+fname+'.txt"'+ ' using 11 with lines title "' + fname +'",'
    # s = s.rstrip(',') + "'"
    # print(s)
    # os.system(s)
    # # os.system("gnuplot -persist -e 'plot \"results/%s\" using 11 with lines title \"%s\"'" % (ffname, fname))

