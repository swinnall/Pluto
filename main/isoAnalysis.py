" Analysis module within IAP "

import csv
import os
import sys
import numpy as np
from numpy import diff
from numpy import linalg
from scipy.signal import savgol_filter
import warnings
import re
import ast
import sys
import math as mt
import config
import genPlot


def importSampleData(info):

    # number of header rows
    nHeaders = 2

    # number of isotherms to plot
    nFiles  = len(info) - nHeaders

    # variable initialisations
    fileNames  = {init: 0 for init in range(nFiles)}
    equip      = {init: 0 for init in range(nFiles)}
    nLipids    = {init: 0 for init in range(nFiles)}
    lipidType  = {init: 0 for init in range(nFiles)}
    lipidRatio = {init: 0 for init in range(nFiles)}
    conc       = {init: 0 for init in range(nFiles)}
    volAdded   = {init: 0 for init in range(nFiles)}
    l          = {init: 0 for init in range(nFiles)}

    # assign data
    for i in range(nHeaders,len(info)):

        # correct for indexing
        j = i - nHeaders

        fileNames[j]  = info[i][0]
        equip[j]      = info[i][1]
        nLipids[j]    = info[i][2]
        lipidType[j]  = info[i][3]
        lipidRatio[j] = info[i][4]
        conc[j]       = info[i][5]
        volAdded[j]   = info[i][6]
        l[j]          = info[i][7]

    return nFiles, fileNames, equip, nLipids, lipidType, lipidRatio, conc, volAdded, l



def getFile(fileDIR,equipParams):

    with open(fileDIR, newline = '') as f:
        reader = csv.reader(f, delimiter=equipParams[2])
        data = list(reader)

    return data



def shiftIsotherm(P):

    # take absolute value to add postive number
    deltaP = abs(min(P))

    # update values in memory
    for i in range(0,len(P)):
        P[i] += deltaP

    return P, deltaP



def importData(equipParams, fname, plotDIR):

    # initialise variables
    t = []
    A = []
    P = []

    # read file into memory data
    fileDIR = "../input/00/" + fname + equipParams[1]
    data = getFile(fileDIR,equipParams)

    # Kibron
    if equipParams[1] == ".NTB":

        # extracting relevant data; key: [row][column]
        for j in range(equipParams[0],len(data)):
            t.append(float(data[j][8]))
            A.append(float(data[j][1])*(10**-6))
            P.append(float(data[j][5]))

    # Nima
    if equipParams[1] == ".txt":

        for j in range(equipParams[0],len(data)):
            t.append(float(data[j][0]))
            A.append(float(data[j][1])*(10**-4))
            P.append(float(data[j][5]))



    # set tStart to 0
    t0 = t[0]
    for i in range(len(t)):
        t[i] = t[i] - t0


    # shift data up so minimum = 0
    if config.shiftP == True:

        P, deltaP = shiftIsotherm(P)

        # define file save path name
        logDIR = plotDIR + "/dataChangeLog.txt"

        # write file
        header = ['filename', 'deltaP']
        data   = [fname, deltaP]

        if not os.path.exists(logDIR):

            with open(logDIR, 'w') as f:
                writer = csv.writer(f)
                writer.writerow(header)
                writer.writerow(data)

        else:
            with open(logDIR, 'a') as f:
                writer = csv.writer(f)
                writer.writerow(data)

    return t, A, P



def calcAreaPerMolecule(i, A_list, lipidMW, lipidType, nLipids, lipidRatio, conc, volAdded):

    # Aoagadro's Number
    NA = 6.0221409E23

    # initialisers
    Mw       = 0
    ratioTot = 0
    Am       = []

    # extract lipid types and ratios from string, returns list of one tuple per var
    if int(nLipids.get(i)[0]) > 1:
        types  = re.split(':',lipidType.get(i))
        ratios = re.split(':',lipidRatio.get(i))

        # sum all ratios for creating fractional components
        for j in range(len(ratios)):
            ratioTot += int(ratios[j])

        # call Mw of each lipid type, multiply by fraction, sum weights for total
        for j in range(int(nLipids.get(i)[0])):
            Mw += lipidMW.get(types[j]) * int(ratios[j]) / ratioTot

    # different condition if there is only one input lipid
    elif int(nLipids.get(i)[0]) == 1:
        types = lipidType.get(i)

        # lipids are singular
        Mw = lipidMW.get(types)

    else:
        print("Error: Undefined number of lipids in isotherm %d, row %d in txt input file." %(i,i+1))
        sys.exit()


    # calculate area per molecule
    for j in range(0,len(A_list)):
        Am.append(   ( A_list[j] / (ast.literal_eval(conc.get(i)) * ast.literal_eval(volAdded.get(i))) ) * (Mw/NA) )


    return Am



def calcElasticity(i, Am_list, P_list):

    # initialise lists
    E_list    = []
    Espl_list = []

    # calculates the derivative for each isotherm
    dPdAm = np.diff(P_list)/np.diff(Am_list)


    # calculate Cs^-1; start from 1st value to avoid spikes
    for j in range(2,len(dPdAm)):
        E_list.append((-Am_list[j])*dPdAm[j])


    # if length is even, truncate as savgol must have odd window length
    if len(E_list) % 2 == 0:
        E_list.pop()


    # remove noise introduced via derivative perturbations
    try:
        Espl_list = savgol_filter(E_list, len(E_list), config.nPoly, mode='interp')  # variable, window size (length), polynomial order, mode

    # remove infinity and nan values due to Nima having repeated area values
    except linalg.LinAlgError:
        print("LinAlgError")
        for j in range(len(E_list)):
            if E_list[j] == float("inf") or E_list[j] == float("-inf") or mt.isnan(E_list[j]) == True:
                E_list[j] = 0


    # high pass filter: remove values less than 0.01 (final smoothing of function)
    deleteMe = []
    for j in range(1,len(Espl_list)-1):

        if Espl_list[j] < 0.01:

            # store indices to be deleted
            deleteMe.append(j)

    # delete elements from the lists
    Espl_list = np.delete(Espl_list, deleteMe)


    return E_list, Espl_list



def calcPercArea(i, A_list):

    # initialise list
    percA_list = []

    # initial (max) area
    A0 = A_list[0]

    for j in range(len(A_list)):

        # sq and sqrt to give positive values
        diffSq = (A_list[j] - A0)**2
        diff   = mt.sqrt(diffSq)

        percA_list.append( (diff/A0) * 100 )

    return percA_list



def calcNormP(t_list, P_list):

    # initialise temporary lists
    tempNormT_list = []
    tempNormP_list = []

    # for every nth value, store in temp lists
    nth = 1000
    for j in range(0,len(P_list),nth):
        tempNormT_list.append(t_list[j])
        tempNormP_list.append(P_list[j])

    P0 = P_list[0]
    for j in range(len(tempNormP_list)):
        tempNormP_list[j] = tempNormP_list[j]/P0

    return tempNormT_list, tempNormP_list



def reduceNpoints(x, y):

    # initialise lists
    redX = []
    redY = []

    for j in range(0,len(x),config.nth):
        redX.append(x[j])
        redY.append(y[j])

    return redX, redY



def checkGrad(x, y):

    # generate list of numbers
    n = [j for j in range(len(x))]

    # differentiate x (the molecular area in this case)
    dxdn = np.diff(x)/np.diff(n)

    # initalise lists, 0 is row number of initial maximum
    minRowNums = []
    maxRowNums = [0]

    # iterate across length of dxdn except final value to stop for loop
    for j in range(0,len(dxdn)-1):

        # store row numbers where gradient changes -ve to +ve (compression to expansion)
        if dxdn[j] < 0 and dxdn[j+1] > 0:
            minRowNums.append(j)

        # store row numbers where gradient changes +ve to -ve (expansion to compression)
        if (dxdn[j] > 0 and dxdn[j+1] < 0) or j == (len(dxdn)-2):
            maxRowNums.append(j)

    return minRowNums, maxRowNums



def findCycles(x, y, l):

    # reduce number of x and y points to smooth data
    redX, redY = reduceNpoints(x, y)

    # differentiates and scans for dxdn=0,
    minRowNums, maxRowNums = checkGrad(redX, redY)

    # every cycle has a minimum compression point, number of mins = nCycles
    nCycles = len(minRowNums)

    # intialise new dict for storing compression values
    cycleX = {init: [] for init in range(nCycles)}
    cycleY = {init: [] for init in range(nCycles)}
    cycleL = {init: [] for init in range(nCycles)}

    # store pressure for each type in each cycle
    for cycle in range(nCycles):

        # cycles go maxA to maxA of next cycle
        for j in range(maxRowNums[cycle], maxRowNums[cycle+1]):

            # store values
            cycleX[cycle].append(redX[j])
            cycleY[cycle].append(redY[j])

            cycleL[cycle] = "Cycle " + str(cycle+1)

            # stop inner for loop during final cycle at point of max compression
            if cycle == nCycles and j == maxRowNums[cycle+1]:
                break

    return cycleX, cycleY, cycleL



def splitCycles(type, cycleX, cycleY, l):

    # extract number of cycles
    nCycles = len(cycleX)

    # intialise new dict for storing compression values
    halfCycleX = {init: [] for init in range(nCycles)}
    halfCycleY = {init: [] for init in range(nCycles)}
    halfCycleL = {init: [] for init in range(nCycles)}

    # store pressure for each type in each cycle
    for cycle in range(nCycles):

        Pmax = max(cycleY.get(cycle))

        for i in range(len(cycleX.get(cycle))):
            if cycleY.get(cycle)[i] == Pmax:
                jMinArea = i


        jMaxArea = len(cycleY.get(cycle))

        if type == " - compressions":

            # compressions go maxA to minA
            for j in range(0, jMinArea):

                # store values
                halfCycleX[cycle].append(cycleX.get(cycle)[j])
                halfCycleY[cycle].append(cycleY.get(cycle)[j])
                halfCycleL[cycle] = "Compression " + str(cycle+1)

                # stop inner for loop during final cycle at point of min compression
                if j == jMinArea:
                    print("call")
                    break


        if type == " - expansions":

            # expansions go minA to maxA
            for j in range(jMinArea, jMaxArea):

                # store values
                halfCycleX[cycle].append(cycleX.get(cycle)[j])
                halfCycleY[cycle].append(cycleY.get(cycle)[j])
                halfCycleL[cycle] = "Expansion " + str(cycle+1)

                # stop inner for loop during final cycle at point of min compression
                if j == jMaxArea:
                    break

    #print(halfCycleY)
    return halfCycleX, halfCycleY, halfCycleL



def main(info, title, plotDIR):

    # filter warnings
    warnings.filterwarnings("ignore")


    # import molecular weight database
    lipidMW = config.lipidMw


    # user input & sample instructions information
    nFiles, fileNames, equip, nLipids, lipidType, lipidRatio, conc, volAdded, l = importSampleData(info)


    # initialise dicts, store 1 file in each list
    t     = {new_list: [] for new_list in range(nFiles)}
    A     = {new_list: [] for new_list in range(nFiles)}
    P     = {new_list: [] for new_list in range(nFiles)}
    Am    = {new_list: [] for new_list in range(nFiles)}
    E     = {new_list: [] for new_list in range(nFiles)}
    Espl  = {new_list: [] for new_list in range(nFiles)}
    percA = {new_list: [] for new_list in range(nFiles)}
    normT = {new_list: [] for new_list in range(nFiles)}
    normP = {new_list: [] for new_list in range(nFiles)}


    # initialise master dicts, store dict of cycles in each element
    master_P  = {new_dict: {0: []} for new_dict in range(nFiles)}
    master_Am = {new_dict: {0: []} for new_dict in range(nFiles)}
    master_L = {new_dict: {0: []} for new_dict in range(nFiles)}


    ## perform analysis on each file and store in dict for genPlot.py
    for i in range(nFiles):

        # set filename
        fname = fileNames.get(i)

        # organise equipment
        if equip.get(i) in ["Kibron", "kibron"]:
            equipParams = (41,".NTB",";")
        elif equip.get(i) in ["Nima", "nima"]:
            equipParams = (1,".txt","\t")


        # import data from given file
        t_list, A_list, P_list = importData(equipParams, fname, plotDIR)
        t[i] = t_list
        A[i] = A_list
        P[i] = P_list


        # calculate area per molecule
        Am_list = calcAreaPerMolecule(i, A_list, lipidMW, lipidType, nLipids, lipidRatio, conc, volAdded)

        # convert units to Angstroms
        for j in range(len(Am_list)):
            Am_list[j] = Am_list[j] * (10**20)
        Am[i] = Am_list


        # pass list, get a dict of cycles, store each dict in a masterDict
        if config.plotPressure == False and config.plotIsotherm == True:
            cycleAm, cycleP, cycleL = findCycles(Am.get(i), P.get(i), l.get(i))
            master_Am[i] = cycleAm
            master_P[i]  = cycleP
            master_L[i]  = cycleL


        if config.plotNormInjection == True:
            tempNormT_list, tempNormP_list = calcNormP(t_list, P_list)
            normT[i] = tempNormT_list
            normP[i] = tempNormP_list


        # calculate percentage area
        if config.plotArea == True:
            percA = calcPercArea(A_list)
            percA[i] = percA_list


        # calculate elasticity
        if config.plotElasticity == True:
            suffix = " - compressions"

            # pass dict of cycles to splitCycles
            halfCycleAm, halfCycleP, halfCycleL = splitCycles(suffix, cycleAm, cycleP, l.get(i))

            # pass the first half of the files cycles as list
            E_list, Espl_list = calcElasticity(i, halfCycleAm.get(0), halfCycleP.get(0))
            E[i]    = E_list
            Espl[i] = Espl_list


    ## plot instructions
    if config.plotElasticity == True:
        axLabels = {"x": "$\pi$ ($mNm^{-1}$)", "x1": "Molecular Area ($\AA$$^2$ / Molecule)", "y": "$C_s^{-1} (mNm^{-1})$"}
        suffix   = " - elasticity"
        key = 2; vars = (nFiles, equip, l, axLabels, suffix, title, plotDIR, (master_P.get(i),master_Am.get(i)), Espl)
        genPlot.main(key,vars)

    # set default settings
    axLabels = {"x": "Molecular Area ($\AA$$^2$ / Molecule)", "y": "$\pi$ ($mNm^{-1}$)"}
    key      = 1

    if config.plotIsotherm == True:
        suffix   = " - isotherm"
        vars     = (nFiles, equip, l, axLabels, suffix, title, plotDIR, (Am,0), P)
        genPlot.main(key,vars)

    if config.plotPressure == True:
        axLabels = {"x": "auto calculated in genPlot", "y": "$\pi$ ($mNm^{-1}$)"}
        suffix   = " - pressure"
        vars     = (nFiles, equip, l, axLabels, suffix, title, plotDIR, (t,0), P)
        genPlot.main(key,vars)

    if config.plotNormInjection == True:
        axLabels = {"x": "auto calculated in genPlot", "y": "$\Delta\pi$"}
        suffix   = " - normInjPressure"
        vars     = (nFiles, equip, l, axLabels, suffix, title, plotDIR, (normT,0), normP)
        genPlot.main(key,vars)

    if config.plotArea == True:
        axLabels = {"x": "auto calculated in genPlot", "y": "$\Delta$A (%)"}
        suffix   = " - area"
        vars     = (nFiles, equip, l, axLabels, suffix, title, plotDIR, (t,0), percA)
        genPlot.main(key,vars)

    for i in range(nFiles):
        newTitle = title + " - " + l.get(i)

        if config.plotCompressions == True:
            suffix = " - compressions"
            halfCycleAm, halfCycleP, halfCycleL = splitCycles(suffix, master_Am.get(i), master_P.get(i), l.get(i))
            vars = (len(master_Am.get(i)), equip, halfCycleL, axLabels, suffix, newTitle, plotDIR, (halfCycleAm,0), halfCycleP)
            genPlot.main(key,vars)

        if config.plotExpansions == True:
            suffix = " - expansions"
            halfCycleAm, halfCycleP, halfCycleL = splitCycles(suffix, master_Am.get(i), master_P.get(i), l.get(i))
            vars = (len(master_Am.get(i)), equip, halfCycleL, axLabels, suffix, newTitle, plotDIR, (halfCycleAm,0), halfCycleP)
            genPlot.main(key,vars)

        if config.plotCycles == True:
            suffix = " - cycles"
            vars = (len(master_Am.get(i)), equip, master_L.get(i), axLabels, suffix, newTitle, plotDIR, (master_Am.get(i),0), master_P.get(i))
            genPlot.main(key,vars)


    # program executed
    print('\nAnalysis Complete! Have a nice day :)')

    return



if __name__ == '__main__':
    print("~Running isoAnalysis.py~\n")
    main()
