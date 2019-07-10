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
# WORKING:
#       rae3.py --outDir=kkkkkkkk 0619Amen2 --species="-1,IndirectOffspring=[-8:-4]"

import sys
import os
from pathlib import Path
import argparse


def buildList(pyExpr):
    r = []
    tmp = pyExpr.split(";")
    for item in tmp:
        if ":" in item:
            first, last = map(int,item.split(':'))
            r +=  map(str,list(range(first, last)))
        else:
            r.append(item)

    return r

def expandArg(anArg):
    r = []
    if "[" not in anArg:
        r = ["", [anArg]]
    else:
        r = anArg.strip()[:-1].split("[")
        print(r)
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
    print("""Example:
          rae3.py --outDir=kkkkkkkk 0619Amen2 --species="-1,IndirectOffspring=[-8:-4;2],0,DirectOffspring=[1;2]"

          You must express an output directory in the first parameter:
                    --outDir=nameOfTheDirectory
          if it it starts in / it'd be an absolute path
          in other case, it will be created inside 'results' directory

          Inside --species=" .. "
          lists of values have to be separated by ;
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
        if not outDir.startswith("/"):
            outDir = os.path.join("results", outDir)
        continue

    if not arg.startswith("--species"):
        ranges.append(expandArg(arg))
    else:
        ranges.append(["", build("--species=", map(expandArg, arg[10:].split(",")), "," )])

commandList = build("ae3.py ", ranges, " ")


if not outDir:
    print("You must express an output directory in parameters:")
    print("          --outDir=nameOfTheDirectory")
    print("if it it starts in / it'd be an absolute path")
    print("in other case, it will be created inside 'results' directory")
    exit(1)

Path(outDir).mkdir(parents=True, exist_ok=True)



n = 1
print(outDir+"/"+"_list.txt")
listing = open(outDir+"/_list.txt", "a")
for command in commandList:
    fname = fileNumbering(n)
    command += " --saveExcel"
    command += " --outFName=" + fname
    print(command, file=listing)
    print(command)
    if not testMode:
        os.system(command)
    n += 1


# --numGen int
# --verbose
# --saveExcel

# i - initFile
# r - --setRandomSeed int

# c - --consume 'str'
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

