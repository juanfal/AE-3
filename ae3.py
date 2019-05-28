#!/usr/bin/env python3
# encoding: utf-8
# ae3.py
# Carlos Villagrasa, Javier Falgueras
# juanfc 2019-02-16

__version__ = 0.057 # 2019-05-28

import os
import sys
import json
import argparse
import textwrap
from pprint import pprint

from numpy import *
from numpy.random import seed, randint, shuffle, sample, permutation, normal

import xlsxwriter
from datetime import datetime



# Directories 'data' and 'results' are supposed being there
# finalForms goes to results/
if not os.path.isdir("results"): os.mkdir("results")

INDIVIDUAL = 0   # they consume df
ACTOR      = 1   # they consume df+if of recipient
RECIPIENT  = 2   # they consume df (it doesn't consume directly)
RECIPROCAL = 3   # they consume (intra 2x df+if) (inter df1+if1+df2+if2)
gNumberOfForms = 4

# Types of distribution
NEIGHBOURS_DISTRIBUTION = 0 # randomly distribute among neighbours around
RANDOM_GLOBAL_AVG   = 1 # random change around global average


# ###############################################################################
#                                MAIN SUBPROGRAMS
# ###############################################################################

def newWorld():
    """Build our world.
     First index:  species index
     Second index: cell index
     Third index:  form of association index
                    4 is the number of possible forms of association one species can have
                    (independent, actor, receptor, reciprocal).
                    Each of these numbers tells the number of members of that species
                    behaving that way
    So we have:
        gWorld[iSpecies, iCell, iForm]

    """
    return (zeros((gNumberOfSpecies, gNumberOfCells, gNumberOfForms), dtype=int),
            zeros((gNumberOfSpecies, gNumberOfForms), dtype=int),
            zeros((gNumberOfSpecies, gNumberOfForms), dtype=int)
            )

# def averagedDist(iSpecies, pos, length):
#     """Returns the average of elements around pos"""
#     otherSide = pos - length
#     if otherSide < 0:
#         s = gWorld[iSpecies, otherSide:, INDIVIDUAL].sum() + \
#             gWorld[iSpecies, 0:pos+length+1, INDIVIDUAL].sum()
#     else:
#         s = gWorld[iSpecies, pos-length:pos+length+1, INDIVIDUAL].sum()
#     mean      = s // (2*length+1)
#     return mean

def randomDist(nitems):
    """Returns a list of nitems numbers randomly distributed in gNumberOfCells"""

    r = sort(randint(0, nitems+1, gNumberOfCells-1))
    return concatenate((r,[nitems])) - concatenate(([0],r))

# For faster traversing of the world, vectorize the randomDist() function
# https://hackernoon.com/speeding-up-your-code-2-vectorizing-the-loops-with-numpy-e380e939bed3
vRandomDist = vectorize(randomDist, signature='()->(m)')

def doInitialSpreading():
    """In fact this is unneccesary as first step in the generation process
                will do it (again)"""
    if gConf["Distribution"].endswith("%"):
        gWorld[:, :, 0] = vRandomDist([gConf["species"][i]["NumberOfItems"] for i in range(gNumberOfSpecies)])
    else:
        # distribute all equally
        for iSpecies in range(gNumberOfSpecies):
            n = gConf["species"][iSpecies]["NumberOfItems"]
            meanPerCell = n // gNumberOfCells
            remainder   = n %  gNumberOfCells
            gWorld[iSpecies, :, 0] = [meanPerCell for _ in range(gNumberOfCells)]
            # the excess of items are equally given out to the first cells
            for iCell in range(remainder):
                gWorld[iSpecies, iCell, 0] += 1

def doDistribute():
    """Spread each each species in each cell"""
    # We suppose here we have only INDIVIDUALs (form index 0)
    global gDistribution
    totalOfEachSpecies = gWorld[:,:,INDIVIDUAL].sum(axis=1)

    for iSpecies in range(gNumberOfSpecies):
        distType, distVal = getDist(iSpecies)

        tempDist = zeros((gNumberOfCells), dtype=int)
        #             NEIGHBOURS_DISTRIBUTION
        if distType == NEIGHBOURS_DISTRIBUTION:
            distLen = int(gNumberOfCells*distVal/100/2)
            for iCell in range(gNumberOfCells):
                nitems = gWorld[iSpecies, iCell, INDIVIDUAL]
                # TODO: if the number of items nitems
                # is large, we could consider directly writing the
                # average on each cell around
                for _ in range(nitems):
                    dist = randint(-distLen, distLen+1)
                    p = iCell+dist
                    if p >= gNumberOfCells:
                        p %= gNumberOfCells
                    tempDist[p] += 1
            gWorld[iSpecies,:,INDIVIDUAL] = tempDist


        else:       # RANDOM_GLOBAL_AVG: x_a (1-\sigma) + xm \sigma
            nitems = gWorld[iSpecies,:,INDIVIDUAL].sum()
            average  = nitems / gNumberOfCells
            sigma = distVal/100
            wildDist = vRandomDist(totalOfEachSpecies)
            wildDist = wildDist * (1-sigma) + sigma * average
            gWorld[:, :, INDIVIDUAL] = wildDist


def doGrouping():
    """Form the groups following the group partner"""
    for iCell in range(gNumberOfCells):
        stillPossibleGroupInCell = True
        n = 0
        while stillPossibleGroupInCell:
            printv("Cell grouping %d" % n)
            n += 1
            stillPossibleGroupInCell = False
            permutedListOrigGroups = permutation(gWithPartnerList)
            printv("Grouping permutation list:", permutedListOrigGroups)
            for iSpecies in permutedListOrigGroups:
                phenFlex = gConf["species"][iSpecies]["PhenotypicFlexibility"]
                i_partner = iGetPartner(iSpecies)          # i partner for the group
                i_grouped = iGetGroupStartingIn(iSpecies)  # i already formed group that iSpecies starts
                # number of iSpecies items in that iCell
                ni = gWorld[iSpecies,iCell,INDIVIDUAL]
                ni_alreadyGrouped = gWorld[i_grouped,iCell,INDIVIDUAL]

                # item to group with
                if i_partner != iSpecies:
                    ni_partner = gWorld[i_partner,iCell,INDIVIDUAL]
                else:
                    # if itself, only half available
                    ni_alreadyGrouped *= 2
                    ni_partner = ni//2

                # Phenotypic complexity tells the % (0-1) from the total
                # existing particular species, that should be grouped, so we
                # compare the current amount of items (grouped + ungrouped)
                # and check if there is room for more grouping
                itemsAloneAndGrouped = ni_alreadyGrouped + ni
                ni_toGroup = noNeg(int(itemsAloneAndGrouped * phenFlex) - ni_alreadyGrouped)
                ni_feasible = min(ni_toGroup, ni_partner)
                #printv("ni_toGroup: %d with %d, ni_feasible: %d" % (ni_toGroup, i_grouped, ni_feasible))
                if ni_feasible > 0:
                    gWorld[ iSpecies,iCell,INDIVIDUAL] -= ni_feasible
                    gWorld[i_partner,iCell,INDIVIDUAL] -= ni_feasible

                    gWorld[i_grouped,iCell,INDIVIDUAL] += ni_feasible

                    stillPossibleGroupInCell = True

def doAssociation():
    """
    performs the association following their AssociatedSpecies
    """
    orderOfAssoc = permutation(gListOfAssociationActors)
    for iCell in range(gNumberOfCells):
        for iActor in orderOfAssoc:
            noIndiv = gWorld[iActor, iCell, INDIVIDUAL]
            iAssoc = iAssociatedTarget(iActor) # must exist
            # Kind of association
            # IRAR: Individual, Recipient, Actor, Reciprocal
            if iActor == iAssoc: # is Reciprocal intraspecific  A><A
                nWidows = noIndiv % 2
                gWorld[iActor, iCell, RECIPROCAL] = noIndiv - nWidows
                gWorld[iActor, iCell, INDIVIDUAL] = nWidows
            else:

                # either A><B or A>B>C, A>B
                # Here you need the minimum of the two numbers, so
                # the efficiency of this association gets cut down to that minimum
                ni_feasible = min(noIndiv, gWorld[iAssoc, iCell, INDIVIDUAL])
                gWorld[iActor, iCell, INDIVIDUAL] -= ni_feasible
                gWorld[iAssoc, iCell, INDIVIDUAL] -= ni_feasible

                if iActor == iAssociatedTarget(iAssoc): # cannot exist
                    #  A><B reciprocal interspecific
                    gWorld[iActor, iCell, RECIPROCAL] += ni_feasible
                    gWorld[iAssoc, iCell, RECIPROCAL] += ni_feasible
                else:  # A>B>C or A>B  simple association
                    gWorld[iActor, iCell,     ACTOR] += ni_feasible
                    gWorld[iAssoc, iCell, RECIPIENT] += ni_feasible
    for iSpecies in range(gNumberOfSpecies):
        gStatsAnt[iSpecies,INDIVIDUAL] = gWorld[iSpecies,:,INDIVIDUAL].sum()
        gStatsAnt[iSpecies,     ACTOR] = gWorld[iSpecies,:,     ACTOR].sum()
        gStatsAnt[iSpecies, RECIPIENT] = gWorld[iSpecies,:, RECIPIENT].sum()
        gStatsAnt[iSpecies,RECIPROCAL] = gWorld[iSpecies,:,RECIPROCAL].sum()

def doConsumeAndOffspring():
    """
    Set up a queue with every form of every species, for each cell
    since recipients receive from actors, both consume at the same time
    A><A    they eat as couple.  Gauss phenot. var. does "NOT" affect.
    A><B    they eat as couple.  Gauss phenot. var. does affect
    We have to queue
        INDIVIDUAL  always
        RECIPIENT   never
        ACTOR       always
        RECIPROCAL  A><A (intras) : half of them, and take out the other half
                    A><B (inters) : all A, but take out all B

    """
    # if gaussian phenotypic variability is considered, we can
    # wait until all the cells are evaluated and then, take the averages
    # for each species and alive item to reset direct/indirect fitness
    if gArgs["varia"]:
        newDirFit  = zeros((gNumberOfSpecies), dtype=[("sum", "int"), ("N", "int")])
    gStatsPost.fill(0)
    for iCell in range(gNumberOfCells):
        """Enqueue, consume resources and have offspring"""

        # ENQUEUEING
        queue = doEnqueueing(iCell)
        # queue = doMultilevelSelection(queue)


        # Lets go eating and have offspring:
        #   - Considers the gaussian phenotypic variability
        #   - Computes an sets offspring as individuals from the other forms
        rsrc = gConf["NumberOfRsrcsInEachCell"]
        i = 0
        while rsrc > 0 and i < queue.size:
            iSpecies = queue[i]["iSp"]
            iForm = queue[i]["iForm"]

            dirFit   = gConf["species"][iSpecies]["DirectOffspring"]
            indirFit = gConf["species"][iSpecies]["IndirectOffspring"]
            stdDev   = gConf["species"][iSpecies]["StandardDeviation"]
            toEat = 0.0
            if iForm == INDIVIDUAL:
                toEat = dirFit
                if rsrc >= toEat:
                    gStatsPost[iSpecies, INDIVIDUAL] += 1

                    gWorld[iSpecies, iCell, INDIVIDUAL] += dirFit
                    if gArgs["varia"]:
                        newDirFit[iSpecies]["sum"] += dirFit * dirFit
                        newDirFit[iSpecies]["N"] += dirFit

            else:
                iAssoc = iAssociatedTarget(iSpecies) # must exist
                assocDirFit = gConf["species"][iAssoc]["DirectOffspring"]
                if iForm == ACTOR:
                    toEat = dirFit + abs(indirFit) \
                             + assocDirFit
                    if rsrc >= toEat:

                        gStatsPost[iSpecies,ACTOR] += 1
                        gWorld[iSpecies, iCell, INDIVIDUAL] += dirFit

                        if gArgs["varia"]:
                            dirFit, indirFit = gauss(dirFit, indirFit, stdDev)
                            newDirFit[iSpecies]["sum"] += dirFit * dirFit
                            newDirFit[iSpecies]["N"] += dirFit

                        assocOffspring = noNeg(indirFit + assocDirFit)

                        gStatsPost[iAssoc,RECIPIENT] += 1
                        gWorld[iAssoc, iCell, INDIVIDUAL] += assocOffspring

                        if gArgs["varia"]:
                            newDirFit[iAssoc]["sum"] += assocOffspring * assocDirFit
                            newDirFit[iAssoc]["N"] += assocOffspring


                else: # RECIPROCAL
                    assocIndirFit = gConf["species"][iAssoc]["IndirectOffspring"]
                    toEat = dirFit    + abs(indirFit) + \
                          assocDirFit + abs(assocIndirFit)
                    if rsrc >= toEat:

                        assocStdDev = gConf["species"][iAssoc]["StandardDeviation"]

                        if gArgs["varia"]:
                            dirFit, indirFit = gauss(dirFit, indirFit, stdDev)
                            assocDirFit, assocIndirFit = gauss(assocDirFit, assocIndirFit, assocStdDev)

                        # dir nuevo (dirFit+assocIndirFit) * dirFit
                        # n   dirFit+assocIndirFit
                        offspring = noNeg(dirFit + assocIndirFit)

                        gStatsPost[iSpecies,RECIPROCAL] += 1
                        gWorld[iSpecies, iCell, INDIVIDUAL] += offspring

                        if gArgs["varia"]:
                            newDirFit[iSpecies]["sum"] += offspring * dirFit
                            newDirFit[iSpecies]["N"] += offspring

                        # dir nuevo (assocDirFit+indirFit) * assocDirFit
                        # n   assocDirFit+indirFit
                        offspring = noNeg(assocDirFit + indirFit)

                        gStatsPost[iAssoc,RECIPROCAL] += 1
                        gWorld[iAssoc, iCell, INDIVIDUAL] += offspring

                        if gArgs["varia"]:
                            newDirFit[iAssoc]["sum"] += offspring * assocDirFit
                            newDirFit[iAssoc]["N"] += offspring

            rsrc -= toEat
            i += 1

    if gArgs["varia"]:
        for iSpecies in range(gNumberOfSpecies):
            # CHANGE THE GLOBAL DirectOffspring gConf parameter
            if newDirFit[iSpecies]["N"]:
                incFit = gConf["species"][iSpecies]["DirectOffspring"] + gConf["species"][iSpecies]["IndirectOffspring"]
                dirFit = round(newDirFit[iSpecies]["sum"] / newDirFit[iSpecies]["N"])
                gConf["species"][iSpecies]["DirectOffspring"] = dirFit
                gConf["species"][iSpecies]["IndirectOffspring"] = incFit - dirFit

def doEnqueueing(iCell):
    # Setting up the queue to feed out
    nTotalInCell = gWorld[:,iCell,:].sum()   # size of the queue, larger than necessary
    # some forms wont enqueue

    queue = zeros((nTotalInCell), dtype=[("iSp", "int"), ("iForm", "int")])
    iQueueTop = 0
    for iSpecies in range(gNumberOfSpecies):
        iAssoc = iAssociatedTarget(iSpecies) # can not exist

        n = gWorld[iSpecies, iCell, INDIVIDUAL]
        if n > 0:
            queue[iQueueTop:iQueueTop+n] = n * [(iSpecies, INDIVIDUAL)]
            iQueueTop += n
            gWorld[iSpecies, iCell, INDIVIDUAL] = 0

        n = gWorld[iSpecies, iCell, ACTOR]
        if n > 0:
            queue[iQueueTop:iQueueTop+n] = n * [(iSpecies, ACTOR)]
            iQueueTop += n
            gWorld[iSpecies, iCell, ACTOR] = 0
            gWorld[  iAssoc, iCell, RECIPIENT] -= n # iAssoc must exist

        n = gWorld[iSpecies, iCell, RECIPROCAL]
        if n > 0:
            # iAssoc must exist
            if iAssoc == iSpecies: # intraspecific
                # n must be even
                queue[iQueueTop:iQueueTop+n//2] = (n//2)*[(iSpecies, RECIPROCAL)]
                iQueueTop += n//2
                gWorld[iSpecies, iCell, RECIPROCAL] = 0
            else:                  # interspecific
                queue[iQueueTop:iQueueTop+n] = n * [(iSpecies, RECIPROCAL)]
                iQueueTop += n
                gWorld[iSpecies, iCell, RECIPROCAL] = 0
                gWorld[  iAssoc, iCell, RECIPROCAL] -= n
    # Shrink the queue down to the real number of occupants
    queue.resize((iQueueTop))
    shuffle(queue)
    # print(queue)
    return queue

def doUngroup():
    for iCell in range(gNumberOfCells):
        for iOrigGroup in gWithPartnerList:
            printv("Ungrouping iOrigGroup:", iOrigGroup)
            iGroup = iGetGroupStartingIn(iOrigGroup)
            phenFlex = gConf["species"][iGroup]["PhenotypicFlexibility"]
            #JAVIER: la desagrupación es por la FlexPheno del patner principal,
            #debería ser por la FlexPhen del grupo A|B
            # Juan: Está puesto que la phenPlex es la del iGroup
            # Date cuenta que gOrigGroupedList es la lista de los que tienen | en el id
            # así que se toma la phenFlex de esos
            #Javier: Por ejemplo, la PhenotypicFlexibility de A es 0.8, y la PhenotypicFlexibility
            # de A|B es 0.1, la phenFlex de la linea 533 debe de ser 0.1

            # number of iGroup items in that iCell
            ni = gWorld[iGroup,iCell,INDIVIDUAL] # form index is 0 always

            ni_unGroup = int(round(ni * (1.0 - phenFlex)))

            if ni_unGroup > 0:
                gWorld[iOrigGroup,iCell,INDIVIDUAL] += ni_unGroup
                iPartner = iGetPartner(iOrigGroup)
                gWorld[iPartner,iCell,INDIVIDUAL] += ni_unGroup

                gWorld[iOrigGroup,iCell,INDIVIDUAL] -= ni_unGroup
                # printv("phenFlex de [%d] %s (%s): %.2f. To ungroup: %d -> %d -> %d [%d, %d]" %
                #        (iOrigGroup, gConf["species"][iOrigGroup]["id"],
                #        partnersList, phenFlex, ni, ni_unGroup,
                #        gWorld[iOrigGroup,iCell,INDIVIDUAL], gWorld[0,iCell,INDIVIDUAL],
                #        gStatsPost[iOrigGroup, INDIVIDUAL]))

def doMultilevelSelection(q):
    # perc = gConf["MultilevelDeath1Percent"]

    # i = 0
    # while i < q.size:
    #     iSpecies = q[i]["iSp"]
    #     iForm = q[i]["iForm"]
    #     # toExtinct =
    #     i += 1

    return q

def calcEgoism():
    """Computes the global variable gEgoism.
    It doesn't need to know what items are there for the estimation is based
    only on the corresponding direct fitnesses and association indirect
    fitnesses of each species/form
    """
    global gEgoism # because at start it is not a global structure/object
    if not gArgs["egoism"]:
        return
    if not gEgoism:
        gEgoism = zeros((gNumberOfSpecies, gNumberOfForms), dtype=float)
    else:
        gEgoism.fill(0)

    for iSpecies in range(gNumberOfSpecies):
        dfit = gConf["species"][iSpecies]["DirectOffspring"]
        gEgoism[iSpecies, INDIVIDUAL] = dfit

        if gConf["species"][iSpecies]["AssociatedSpecies"]:
            ifit = gConf["species"][iSpecies]["IndirectOffspring"]
            iAssoc = iAssociatedTarget(iSpecies)
            assoc_ifit = gConf["species"][iAssoc]["IndirectOffspring"]

            # rancour makes associated ind. fit. change sign
            gEgoism[iSpecies,      ACTOR] = dfit - abs(ifit)
            gEgoism[iSpecies,  RECIPIENT] = dfit + assoc_ifit
            gEgoism[iSpecies, RECIPROCAL] = dfit - abs(ifit) + assoc_ifit


    # Do scaling 0.0 - 1.0
    minEgo = amin(gEgoism) ; printv("minEgo: ", minEgo)
    maxEgo = amax(gEgoism) ; printv("maxEgo: ", maxEgo)
    if maxEgo > 0:
        gEgoism -= minEgo
        gEgoism /= maxEgo-minEgo
    printv("scaled egoism: ", gEgoism)
    # printv("rank:", gEgoism.argsort(axis=1).argsort(axis=0))
    #
    # printv("rank:", reshape(gEgoism.flatten().unique().argsort().argsort(), (gNumberOfSpecies, gNumberOfForms) ))
    # printv("rank:", reshape(gEgoism.flatten().unique().argsort().argsort(), (gNumberOfSpecies, gNumberOfForms) ))
    # TODO: NOS INTERESA LA JERARQUÍA, rank, MÁS QUE EL VALOR ESCALADO DEL EGOÍSMO
    # https://stackoverflow.com/questions/5284646/rank-items-in-an-array-using-python-numpy

# ###############################################################################
#                                     TOOLS
# ###############################################################################

def iFrom_id(id):
    """Returns index of species from its id, or -1"""
    i = gNumberOfSpecies - 1
    while i >= 0 and gConf["species"][i]["id"] != id:
        i -= 1
    return i

def iGetPartner(iSpecies):
    "returns index of the partner grouping"
    toFindId = gConf["species"][iSpecies]["GroupPartner"]
    return iFrom_id(toFindId)

def iGetGroupStartingIn(iSpecies):
    toFindId = gConf["species"][iSpecies]["id"] + '|' + gConf["species"][iSpecies]["GroupPartner"]
    return iFrom_id(toFindId)

def getListOfOrigGroups():
    theList = []
    for iSpecies in range(gNumberOfSpecies):
        partner = gConf["species"][iSpecies]["GroupPartner"]
        if partner != "":
            theList += [iSpecies]

    # Puts complex groups first
    theList.sort(key = lambda i:
            gConf["species"][i]["id"].count('|') +
            gConf["species"][i]["GroupPartner"].count('|'), reverse=True)
    return theList

def collectAssociationActors():
    "Return a list of indexes of actors which are associated with any other"
    return [i for i in range(gNumberOfSpecies) if gConf["species"][i]["AssociatedSpecies"] != ""]

def iAssociatedTarget(iSpecies):
    """returns either the index of associated partner or
    if not associated is there, -1"""
    return iFrom_id(gConf["species"][iSpecies]["AssociatedSpecies"])

def getDist(iSpecies):
    """from 45n or 50n returns the integer and the type
    n: NEIGHBOURS_DISTRIBUTION
    r: RANDOM_GLOBAL_AVG
    """
    if "Distribution" in gConf["species"][iSpecies]:
        dist = gConf["species"][iSpecies]["Distribution"]
    else:
        dist = gConf["Distribution"]

    if dist.endswith("n"): # average neighbours distribution
        distType = NEIGHBOURS_DISTRIBUTION
    else:
        distType = RANDOM_GLOBAL_AVG

    distSpan = int(dist[:-1])
    return distType, distSpan

def randomNormal(mean, stddev):
    """Draw random samples from a normal (Gaussian) distribution, with
    mean and standard deviation"""
    return normal(mean, stddev)

def gauss(direct, indirect, stddev):
    inclusive = direct + indirect
    gaussDirect = noNeg(randomNormal(direct, stddev))
    gaussIndirect = inclusive - gaussDirect
    return gaussDirect, gaussIndirect

def noNeg(n):
    if n < 0: n = 0
    return n

def ranking(l, value):
    rank = unique(l.flatten())
    return 1+where(rank==value)[0][0]

def initExcel():
    global gExcelCellHeader, gExcelCellID, gThedatetime
    gThedatetime = datetime.now().strftime("%Y%m%d-%H%M%S")
    excelOut = os.path.join("results", gInitConfFile + '_' + gThedatetime + ".xlsx")

    workbook = xlsxwriter.Workbook(excelOut)
    worksheet = workbook.add_worksheet()
    gExcelCellHeader = workbook.add_format({
                                           'align':'center', 'bg_color': '#CCCCFF', 'bold': True})
    gExcelCellID = workbook.add_format({'align':'center', 'bg_color': '#33FF99', 'bold': True})

    workbook.set_properties({
        'title':    'Evolutionary Automata',
        'subject':  gInitConfFile + '_' + gThedatetime,
        'author':   'Javier Falgueras, Juan Falgueras, Santiago Elena',
        'manager':  'Javier Falgueras, Juan Falgueras, Santiago Elena',
        'company':  'Univ. Valencia + Univ. Málaga',
        'category': 'Research Thesis',
        'keywords': 'Evolution Automatata',
        'comments': " ".join(sys.argv)})
    return workbook, worksheet

def saveExcel(sh, numGen):
    """Saving in Excel"""
    # id, NumberOfItems, DirectOffspring, IndirectOffspring,
    # AssociatedSpecies, StandardDeviation, INDIVIDUAL, ACTOR, RECIPIENT,
    # RECIPROCAL
    txtOutName = os.path.join("results", gInitConfFile + '_' + gThedatetime + ".txt")
    txtOut = open(txtOutName, "a")
    globalsW = 3
    if numGen == 1:
        globalsHeader =  ["NCel", "RsCel", "Dist"]
        sh.write_row(0,0, globalsHeader, gExcelCellHeader)
    sh.write_row(numGen,0, [gConf["NumberOfCells"], gConf["NumberOfRsrcsInEachCell"], gConf["Distribution"]])
    print(gConf["NumberOfCells"], "\t", gConf["NumberOfRsrcsInEachCell"],"\t", gConf["Distribution"], "\t", file=txtOut, end="")


    totNumCol = 15
    if numGen == 1:
        header =  ["ID", "D", "I", ">", "σ", "Gr", "Ph", "IND", "2", "ACT", "2", "RNT", "2", "RCL", "2"]
        for iSpecies in range(gNumberOfSpecies):
            sh.write_row(0,globalsW+iSpecies*totNumCol, header, gExcelCellHeader)


    for iSpecies in range(gNumberOfSpecies):
        sh.write(numGen, globalsW+iSpecies*totNumCol,
                 gConf["species"][iSpecies]["id"], gExcelCellID)
        print(gConf["species"][iSpecies]["id"] + "\t", file=txtOut, end="")

        toWrite =[
             gConf["species"][iSpecies]["DirectOffspring"],
             gConf["species"][iSpecies]["IndirectOffspring"],
             gConf["species"][iSpecies]["AssociatedSpecies"],
             gConf["species"][iSpecies]["StandardDeviation"],
             gConf["species"][iSpecies]["GroupPartner"],
             gConf["species"][iSpecies]["PhenotypicFlexibility"],
             gStatsAnt[ iSpecies,INDIVIDUAL],
             gStatsPost[iSpecies,INDIVIDUAL],
             gStatsAnt[ iSpecies,ACTOR],
             gStatsPost[iSpecies,ACTOR],
             gStatsAnt[ iSpecies,RECIPIENT],
             gStatsPost[iSpecies,RECIPIENT],
             gStatsAnt[ iSpecies,RECIPROCAL],
             gStatsPost[iSpecies,RECIPROCAL]
            ]
        sh.write_row(numGen, globalsW+iSpecies*totNumCol+1, toWrite)
        print("\t".join(map(str,toWrite)), file=txtOut, end="")
        if iSpecies < gNumberOfSpecies-1:
            print("\t", file=txtOut, end="")

        if numGen == 1:
            ori = iSpecies*totNumCol + 1 + globalsW
            end = ori + totNumCol - 2
            worksheet.set_column(ori, end, None, None, {'level': 1, 'hidden': True})

    print(file=txtOut)

def saveConf():
    """Save conf in a new _cont.json file
    if not there, or rewrite previous _cont.json file if was the initial conf loaded"""

    # gThedatetime = datetime.now().strftime("%Y%m%d-%H%M%S.%f")
    newConfFile = gInitConfFile
    if not newConfFile.endswith(gContExt):
        newConfFile += gContExt
    newConfFile = os.path.join("data", newConfFile + ".json")
    # print("(Re)writing", newConfFile)
    with open(newConfFile, 'w') as outfile:
        json.dump(gConf, outfile, sort_keys = True, indent = 4,
                   ensure_ascii = False)

def checkConf(conf):
    """Verifies conf returning inconsistencies or an empty string"""
    # TO DO
    # verify gConf so
    #   - groups agree
    #
    problems = ""

    items = {"NumberOfCells", "NumberOfRsrcsInEachCell",
        "MultilevelDeath1Percent", "LambdaForEgoism", "species"}

    itemsspecies = {"id", "NumberOfItems", "DirectOffspring",
        "GroupPartner", "PhenotypicFlexibility", "AssociatedSpecies",
        "IndirectOffspring", "StandardDeviation"}

    if not items.issubset(conf.keys()):
        problems += "Some essential item(s) is not defined\n\t Check the next items are all there:\n\t%s\n" % str(list(items))

    previousIds = set([])
    longListSpecies = len(conf["species"])
    if type(conf["species"]) != list or longListSpecies == 0:
        problems += "Species must be a list with species inside\n"
    for i in range(longListSpecies):
        theId = conf["species"][i]["id"]
        partnerId = conf["species"][i]["GroupPartner"]
        if theId in previousIds:
            problems += "Id '%s' REPEATED\n" % theId
        else:
            previousIds.add(theId)
        if partnerId != "":
            if iGetPartner(i) == -1:
                problems += "Species %d, id: '%s' has not the partner species '%s' in the conf\n" % \
                    (i, theId, partnerId)
            if iGetGroupStartingIn(i) == -1:
                problems += "Group id: '%s' is not configured\n" % \
                    (theId + '|' + partnerId)
        for aId in theId.split('|'):
            if iFrom_id(aId) == -1:
                problems += "Component id: '%s' from group '%s' is not configured\n" % \
                    (aId, theId)


        if not itemsspecies.issubset(conf["species"][i].keys()):
            problems += "Some essential item(s) of species %s is not defined\n\t Check the next items for all species are all there:\n\t%s\n" % (theId, str(list(itemsspecies)))

    if problems:
        problems = "ERROR. The configuration file '%s' has inconsistencies:\n\n%s" % (gInitConfFile, problems)
    return problems

def printv(*args, **kwargs):
    if gArgs["verbose"]:
        if type(args[0]) == str:
            print(*args)
        else:
            pprint(args)

# ARGS AND CONF PROCESSING

def getCommandLineArgs():
    """Parses and returns command line arguments"""

    theArgParser = argparse.ArgumentParser(description="* Evolutionary Automata *",
                    formatter_class=argparse.RawTextHelpFormatter)


    # initial configuration to load from a file
    theArgParser.add_argument(
        "initFile",
        nargs="?",
        default="defaultInit",
        help=textwrap.dedent("""\
        When provided it sets the file with the initial configuration.
        All the configuration files are inside the './data' directory.
        Do not add the .extension to the file name.
        If not provided, it defaults to 'defaultInit'""")
    )

    # Number of generations to run
    theArgParser.add_argument(
        "--numGen", type=int,
        default=10,
        metavar="int",
        help="Sets the number of generations to run, default: 10")

    # We want gaussian phenotypic variability
    theArgParser.add_argument(
        "--varia", help=textwrap.dedent("""\
        Changes offspring number following a random normal distribution
          around previous offspring with the standard deviation provided for
          each species.
        False if not provided"""),
        action='store_true')

    # Setting random seed to 0 to repeat pseudorandom values
    theArgParser.add_argument(
        "--setRandomSeed", type=int,
        default=-1,
        metavar="int",
        help=textwrap.dedent("""\
        Sets the random seed to a fixed initial value so
          to repeat same random sequences. If not, each
          running will start with different seeds"""),
        )

    # We want gaussian phenotypic variability
    theArgParser.add_argument(
        "--verbose", help=textwrap.dedent("""\
        Gives as much detailed information as it can"""),
        default=False,
        action='store_true')

    # # We want egoism multilevel selection
    # theArgParser.add_argument(
    #     "--egoism", help=textwrap.dedent("""\
    #     Considers egoism of each item in multilevel selection"""),
    #     action='store_true')

    # We want to save final status
    theArgParser.add_argument(
        "--saveExcel", help=textwrap.dedent("""\
        Save stats in 'Excel' file and in a txt file.
        See Excel header for meaning of txt columns"""),
        action='store_true')

    theArgParser.add_argument(
        "--NumberOfCells", type=int,
        default=argparse.SUPPRESS, metavar="int",
        help="Sets the total number of cells in our world")

    theArgParser.add_argument(
        "--NumberOfRsrcsInEachCell", type=float,
        default=argparse.SUPPRESS, metavar="float",
        help="Sets the items of resources in each cell for each generation")
    # theArgParser.add_argument(
    #     "--MultilevelDeath1Percent", type=int,
    #     default=argparse.SUPPRESS, metavar="int",
    #     help="Sets MultilevelDeath1Percent to. Range 0.0 to 1.0")
    # theArgParser.add_argument(
    #     "--LambdaForEgoism", type=float,
    #     default=argparse.SUPPRESS, metavar="float",
    #     help="Sets the coef. for egoism application")


    theArgParser.add_argument(
        "--Distribution", type=str,
        metavar="'str'",
        default="100%",
        help=textwrap.dedent("""\
        It allows specify the desired type of distribution between:
          a) Random           from 0%% to 100%% (with suffix %%)
          b) Length to cover the average
                           integer from 0 to number of cells N//2
                           if N//2 then is like 0%% in (a)

        When random, items are taken from cells forgetting
          their previous position and then distributed randomly
          0%%   means null random, it is equivalent to N//2
                averaged in case (b)
          100%% means _totally_ random numbers en each cell

        When averaged through length, the number of cells given
          is taken from the left, from the right and the
          [i] cell to average and this mean is the
          new value for [i]
             [i-n][i-(n-1)][…][i][i+1][…][i+n]
          The list is considered circular
          0 means no cells around are considered,
            so cells are isolated
          3 for example, means 3x2+1 cells
            are averaged: the central one plus the 3 at
            the left plus the 3 at the right""")
    )


    theArgParser.add_argument(
        "--species", type=str,
        default=argparse.SUPPRESS, metavar="'str'",
        help=textwrap.dedent("""\
    Sets one or more species parameters
    Surround everything with " "
    No spaces

    Put first the 0-index of the species:
        --species="0,DirectOffspring=8,1,IndirectOffspring=1"
    Changes the DirectOffspring of the first species and the
    IndirectOffspring of the second.

    A series of indexes separated by spaces is valid:
        --species="0,3,1,DirectOffspring=8"
    changes DirectOffspring of species 0, 3 and 1

    A Python standard range (no 'steps' here) is valid:
        --species="0,8:10,DirectOffspring=8"
    as with Python 8:10 means 8,9 (10 is not included)

    When an index is negative, is suppose from the end:
        --species="0,3,1,8:-2,DirectOffspring=8"
    where 8:-2  is  8,9,10,11,12,13  if there were 15 species.

    An empty index in ranges means a extreme:
        --species=":,DirectOffspring=8"
    means change all species
    """))


    return vars(theArgParser.parse_args())

def readInitConfFile(fileName):
    """Load config from init file"""
    with open(fileName) as f:
        conf = json.load(f)

    return conf

def replaceArgsInConf(conf, args):
    """Replace the conf args with the args provided in the command line"""
    def parseSpeciesArgs(s):
        def isInt(anyNumberOrString):
            try:
                int(anyNumberOrString)
                return True
            except ValueError :
                return False

        lenSpecies = len(conf["species"])
        l = s.split(',')
        indxsRead = True
        for arg in l:
            if isInt(arg) or ':' in arg:   # int (+-) or range
                if indxsRead:
                    ns = []
                    indxsRead = False
                if isInt(arg):
                    ns.append(int(arg))
                else:
                    if arg == ':':
                        first = 0
                        last = lenSpecies
                    elif arg.startswith(':'):
                        first = 0
                        last = int(arg[1:])
                    elif arg.endswith(':'):
                        first = int(arg[1:])
                        last = lenSpecies
                    else:
                        first, last = map(int,arg.split(':'))
                    if first < 0: first += lenSpecies
                    if last  < 0: last  += lenSpecies
                    ns += list(range(first, last))
                    printv(ns)
            else:
                indxsRead = True
                k, val=arg.split('=')
                for n in ns:
                    t = type(conf["species"][n][k])
                    conf["species"][n][k] = t(val)

    # substitute init file conf specific parameters
    # provided through the args in the command line
    for param in args:
        if param in conf:
            if param == "species":
                parseSpeciesArgs(args["species"])
            else:
                t = type(conf[param])
                conf[param] = t(args[param])
    return conf


# ###############################################################################
#                                    MAIN
# ###############################################################################

# GLOBALS

gArgs = getCommandLineArgs()
gInitConfFile = gArgs["initFile"]  # base file name
gInitConfCompName = os.path.join("data", gInitConfFile + ".json")
printv(gArgs)

gConf = replaceArgsInConf(readInitConfFile(gInitConfCompName), gArgs)

gNumberOfSpecies = len(gConf["species"])

s = checkConf(gConf)
if s:
    print(s)
    sys.exit(1)

printv(gConf)

gNumberOfCells   = gConf["NumberOfCells"]


gEgoism          = []
gExcelSaved      = []
gContExt         = "_cont"
gListOfAssociationActors = [] # list of species starting association
gListOfAssociationActors = collectAssociationActors() # list of species starting association
# if gConf["Distribution"].endswith("%"):
#     gDistribution   = float(gConf["Distribution"][:-1])
# else:
#     gDistribution   = int(gConf["Distribution"][:-1])


gWithPartnerList = getListOfOrigGroups()
printv("gWithPartnerList:", gWithPartnerList)


#
# LETS DO IT
#

if gArgs["setRandomSeed"] != -1:
    seed(int(gArgs["setRandomSeed"]))

gWorld, gStatsAnt, gStatsPost = newWorld()

if "egoism" in gArgs:
    calcEgoism()

doInitialSpreading() # pprint(gWorld)



if gArgs["saveExcel"]:
    workbook, worksheet = initExcel()

set_printoptions(formatter={'int': '{: 7d}'.format})
for genNumber in range(1, gArgs["numGen"]+1):

    doDistribute()             # ; pprint(gWorld)
    #print("%3d: %s Tot Spread:   %d" % (genNumber, gWorld[:,:,0].sum(axis=1),  gWorld[:,:,0].sum(axis=1).sum()))
    doGrouping()                # ; pprint(gWorld)
    print("%3d: %s Tot Grouping: %d" % (genNumber, gWorld[:,:,0].sum(axis=1),  gWorld[:,:,0].sum(axis=1).sum()))
    doAssociation()             # ; pprint(gWorld)
    #print("%3d: %s Tot Associac: %d" % (genNumber, gWorld[:,:,0].sum(axis=1),  gWorld[:,:,0].sum(axis=1).sum()))

    doConsumeAndOffspring()     # ; pprint(gWorld)
    #print("%3d: %s Tot Consume:  %d" % (genNumber, gWorld[:,:,0].sum(axis=1),  gWorld[:,:,0].sum(axis=1).sum()))
    doUngroup()
    #print("%3d: %s Tot Ungroup:  %d" % (genNumber, gWorld[:,:,0].sum(axis=1),  gWorld[:,:,0].sum(axis=1).sum()))
    if gArgs["saveExcel"]:
        saveExcel(worksheet, genNumber)
    #print("%3d: %s Tot: %d" % (genNumber, gWorld[:,:,0].sum(axis=1),  gWorld[:,:,0].sum(axis=1).sum()))


saveConf()

if gArgs["saveExcel"]:
    workbook.close()

