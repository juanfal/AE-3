#!/usr/bin/env python3
# encoding: utf-8
# rae3.py
# juanfc 2019-07-08
# Builds the ranges for each parameter and make the corresponding
# ae3.py calls

# TODO: --species can have more than one range… so it has to be expanded by itself
#       naming the output folder and each file

# SOLUTION:
#       --species has been treated especially.  As an special parameter, sorry
#       its ranges do not admit , yet
#       rae3.py --algo=[1:3] --uno=[a] lacagaste --none=nada --species=1,a=[2:4],b=[2:6]

# FOR LISTS INSIDE RANGES
#
#                               USE ;
#
#       rae3.py --algo=[1:3] --uno=[a] lacagaste --none=nada --species="1,2,a=[2:4;7],b=[2:6]"

# WORKING!:
#       rae3.py --outDir=kkkkkkkk 0619Amen2 --species="-1,IndirectOffspring=[-8:-4]"

import sys
import os
from pathlib import Path
import argparse

def frange(x, y, jump):
  while x < y:
    yield x
    x += jump

def buildList(pyExpr):
    r = []
    tmp = pyExpr.split(";")
    for item in tmp:
        if item.count(':') == 1:
            first, last = item.split(':')
            if '.' in first or '.' in last:
                first = float(first)
                last = float(last)
            else:
                first = int(first)
                last = int(last)
            r +=  list(map(str,list(frange(first, last))))
        elif item.count(':') == 2:
            first, last, step = item.split(':')
            if '.' in first or '.' in last or '.' in step:
                first = float(first)
                last = float(last)
                step = float(step)
            else:
                first = int(first)
                last = int(last)
                step = int(step)
            r +=  list(map(str,list(frange(first, last, step))))
        else:
            r.append(item)

    return r

def expandArg(anArg):
    r = []
    if "[" not in anArg:
        r = ["", [anArg]]
    else:
        r = anArg.strip()[:-1].split("[")
        # print(r)
        r = [r[0], buildList(r[1])]

    return r

def build(base, ranges, sep):
    tsep = ""
    commandList = [base]
    for theName, theRange in ranges:
        temp = []
        for i in theRange:
            for prev in commandList:
                temp += [prev + tsep + theName + str(i)]
        commandList = temp
        tsep = sep
    return commandList


def fileNumbering(n):
    return "%04d" % n

if "-h" == sys.argv[1]:
    print("""Examples:
          rae3.py --outDir=kkkkkkkk 0619Amen2 --species="-1,IndirectOffspring=[-8:-4;2],0,DirectOffspring=[1;2]"
          rae3.py --outDir=kkkkkkkk 0619Amen2 --NumberOfCells=[100:10000:100] --species="-1,IndirectOffspring=[-8:-4;2],0,DirectOffspring=[1;2]"
          rae3.py --outDir=kkkkkkkk -t 0619Amen2 --NumberOfCells=[100:103] --species="-1,IndirectOffspring=[-8:-4;2],0,DirectOffspring=[1;2]"
          rae3.py --outDir=kkkkkkkk -t 0619Amen2 --NumberOfCells="[100:110:3;8;10:14:2]" --species="-1,IndirectOffspring=[-8:-4:2;2],0,DirectOffspring=[1;2]"


          You must express an output directory in the first parameter:
                    --outDir=nameOfTheDirectory
          if it it starts in / it'd be an absolute path
          in other case, it will be created inside 'results' directory

          IMPORTANT:
              Observe in the last examples:
              To express a list of different values or ranges, separate them with ;
              but to prevent the Terminal shell from an undesired interpretation, surround it in ';' or ";"
    """)
    exit(0)




outDir = ""
testMode = False
ranges = []
for arg in sys.argv[1:]:

    if "-t" == arg:
        testMode = True
        continue

    if arg.startswith("--outDir="):
        _, outDir = arg.split("=")
        continue

    if not arg.startswith("--species"):
        ranges.append(expandArg(arg))
    else:
        ranges.append(["", build("--species=", map(expandArg, arg[10:].split(",")), "," )])

commandList = build("ae3.py --setRandomSeed=1 ", ranges, " ")


if not outDir:
    print("You must give an output directory in parameters:")
    print("          --outDir=nameOfTheDirectory")
    print("if the dir starts in / it'd be an absolute path")
    print("in other case, it will be created inside 'results' directory")
    exit(1)

outDir = Path(outDir)
if not outDir.root:
    finalOutDir = Path("results") / outDir
if not testMode:
    finalOutDir.mkdir(parents=True, exist_ok=True)

outListFName = finalOutDir / "_list.txt"

n = 1
print(outListFName)
if not testMode:
    listing = open(outListFName, "w")
for command in commandList:
    fname = fileNumbering(n)
    command += " --outDir=" + str(outDir)
    command += " --outFName=" + fname
    print(command)
    if not testMode:
        print(command, file=listing)
        os.system(command)
    n += 1


# --numGen int
# --verbose
# --saveExcel

# i - initFile
# r - --setRandomSeed int

# v - --varia
# n - --NumberOfCells int
# s - --NumberOfRsrcsInEachCell float
# d - --Distribution 'str'
#     --species 'str'
#     "NumberOfItems": 100,
#     "DirectOffspring": 5,

#     "GroupPartner":   "A",
#     "PhenotypicFlexibility": 0.0,

#     "AssociatedSpecies":   "",
#     "IndirectOffspring": 2,
#     "FitnessVariationLimit": 0





# --outDir 'str'
# --outFName 'str'


# print(sys.argv)
# status = os.system('ae3.py 0619Amen2 --species="-1,IndirectOffspring=-8"')
# print(status)

# print(buildList("1,2,4:8,100,1000:1010"))
# print(expandArg("--algo=[1,2,4:8,100,1000:1010]"))
# print(expandArg("--algo=[1,2,3]"))
# print(expandArg("--algo=[a]"))
# print(expandArg("--algo=a=[1:6],b=[-1,8:10]"))
# print(expandArg("--algo=nada"))

# print(ranges)
# ranges = [expandArg("--algo=[1:3]"), expandArg("--uno=[a]"), expandArg("lacagaste"), expandArg("--none=nada"), expandArg("--vars=[1,2,4:8]")]
# --algo=[1:3] --uno=[a] lacagaste --none=nada --vars=[1,2,4:8]

