" Process surface pressure data "

# import packages
import csv, os, sys, re, ast
import warnings
import math as mt
import numpy as np
from numpy import diff
from numpy import linalg
from scipy.signal import savgol_filter

# import Pluto modeuls
import config
import database
from genFunc import modSelection, getFile


## Import Functions
def importSampleInfo(instructionsFile, nFiles):

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
    for i in range(nFiles):

        fileNames[i]  = instructionsFile["filename"][i]
        equip[i]      = instructionsFile["equipment"][i]
        nLipids[i]    = instructionsFile["nLipids"][i]
        lipidType[i]  = instructionsFile["lipidType"][i]
        lipidRatio[i] = instructionsFile["lipidRatio"][i]
        conc[i]       = instructionsFile["conc"][i]
        volAdded[i]   = instructionsFile["volAdded"][i]
        l[i]          = instructionsFile["label"][i]

    return fileNames, equip, nLipids, lipidType, lipidRatio, conc, volAdded, l



def importData(instrumentName, equipParams, fname, inputDIR, plotDIR):

    # read file into memory data
    fileDIR = inputDIR + '/' + fname + equipParams[1]

    # get data file as pandas dataframe
    data = getFile(path=fileDIR,nSkip=equipParams[0],delim=equipParams[2])

    if instrumentName == "KIBRON":
        t = data[data.columns.values[8]]
        A = data[data.columns.values[1]] * (10**-6)
        P = data[data.columns.values[5]]

    if instrumentName == "NIMA":
        t = data[data.columns.values[0]]
        A = data[data.columns.values[1]] * (10**-4)
        P = data[data.columns.values[5]]

    if instrumentName == "CUSTOM":
        t = data[data.columns.values[0]]
        A = data[data.columns.values[2]]  * (10**-6)
        P = data[data.columns.values[1]]


    # set tStart to 0
    if config.checkT0 == True:
        t0 = t[0]
        for i in range(len(t)):
            t[i] = t[i] - t0


    # shift data up so minimum = 0
    if config.shiftP == True and min(P) < 0:

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



## Data Manipulation Functions
def shiftIsotherm(P):

    # take absolute value to add postive number
    deltaP = abs(min(P))

    # update values in memory
    for i in range(0,len(P)):
        P[i] += deltaP

    return P, deltaP



## Calculation Functions
def calcAreaPerMolecule(i, A_list, lipidType, nLipids, lipidRatio, conc, volAdded):

    # import molecular weight database
    lipidMW = database.lipidMw

    # Aoagadro's Number
    NA = 6.0221409E23

    # initialisers
    Mw       = 0
    ratioTot = 0
    Am       = []

    # extract lipid types and ratios from string, returns list of one tuple per var
    if nLipids.get(i) > 1:
        types  = re.split(':',lipidType.get(i))
        ratios = re.split(':',lipidRatio.get(i))

        # sum all ratios for creating fractional components
        for j in range(len(ratios)):
            ratioTot += int(ratios[j])

        # call Mw of each lipid type, multiply by fraction, sum weights for total
        for j in range(int(nLipids.get(i))):
            Mw += lipidMW.get(types[j]) * int(ratios[j]) / ratioTot

    # different condition if there is only one input lipid
    elif nLipids.get(i) == 1:
        types = lipidType.get(i)

        # lipids are singular
        Mw = lipidMW.get(types)

    else:
        print("Error: Undefined number of lipids in isotherm %d, row %d in txt input file." %(i,i+1))
        sys.exit()


    # calculate area per molecule
    for j in range(0,len(A_list)):
        Am.append(   ( A_list[j] / (conc.get(i) * volAdded.get(i)) ) * (Mw/NA) * (10**20))

    return Am



def calcElasticity(Am_list, P_list):

    # initialise list
    E_list    = []

    # calculates the derivative for each isotherm
    dPdAm = np.diff(P_list)/np.diff(Am_list)

    # calculate Cs^-1; start from 1st value to avoid spikes
    for j in range(2,len(dPdAm)):
        E_list.append((-Am_list[j])*dPdAm[j])

    return E_list



def calcPercArea(A_list):

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



## Cycle Extraction Functions
def reduceNpoints(x, y):

    # initialise lists
    redX = []
    redY = []

    # data reduction for gradient calc.; take every 'Nth' point
    gradNth   = 20

    for j in range(0,len(x),gradNth):
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
    if len(x) > 2E3:
        redX, redY = reduceNpoints(x, y)
    else: redX, redY = x, y

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

    return halfCycleX, halfCycleY, halfCycleL



# Main
def main(instructionsFile, title, inputDIR, plotDIR):

    # filter warnings
    warnings.filterwarnings("ignore")

    # number of isotherms to plot
    nFiles  = len(instructionsFile)

    # user input & sample instructions information
    fileNames, equip, nLipids, lipidType, lipidRatio, conc, volAdded, l = importSampleInfo(instructionsFile, nFiles)

    # give analysis choice to user
    analysisOptions = ['plotIsotherm','plotCompressions','plotExpansions','plotCycles','plotPressure','plotArea','plotElasticity']

    analysisType, analysisRunning = modSelection(analysisOptions)

    # initialise dicts, store 1 file in each list
    t     = {new_list: [] for new_list in range(nFiles)}
    A     = {new_list: [] for new_list in range(nFiles)}
    P     = {new_list: [] for new_list in range(nFiles)}
    Am    = {new_list: [] for new_list in range(nFiles)}
    E     = {new_list: [] for new_list in range(nFiles)}
    percA = {new_list: [] for new_list in range(nFiles)}

    stitched_t  = {new_list: [] for new_list in range(nFiles)}
    stitched_P  = {new_list: [] for new_list in range(nFiles)}
    stitched_Am = {new_list: [] for new_list in range(nFiles)}

    elasticityCompression_Am = {new_list: [] for new_list in range(nFiles)}
    elasticityCompression_P  = {new_list: [] for new_list in range(nFiles)}
    elasticityCompression_L  = {new_list: [] for new_list in range(nFiles)}

    # initialise master dicts, store dict of cycles in each element
    master_P  = {new_dict: {0: []} for new_dict in range(nFiles)}
    master_Am = {new_dict: {0: []} for new_dict in range(nFiles)}
    master_L  = {new_dict: {0: []} for new_dict in range(nFiles)}


    ## perform analysis on each file and store in dict for genPlot.py
    for i in range(nFiles):
        fname = fileNames.get(i)

        # organise equipment
        if equip.get(i).upper() == "KIBRON":
            equipParams = (41,".NTB",";")
        elif equip.get(i).upper() == "NIMA":
            equipParams = (1,".txt","\t")
        elif equip.get(i).upper() == "CUSTOM":
            equipParams = (1,".txt","\t")

        # import data from given file i, immediate store t (need A & P lists as well)
        t[i], A[i], P[i] = importData(equip.get(i).upper(), equipParams, fname, inputDIR, plotDIR)

        # calculate percentage area
        if analysisType == 'plotArea':
            percA_list = calcPercArea(A.get(i))
            percA[i] = percA_list

        if analysisType == 'plotPressure' and config.subtractRef == True:

            refData = getFile(path=config.pressureRef,nSkip=equipParams[0],delim=equipParams[2])

            if equip.get(i).upper() == "KIBRON":
                t_ref = refData[refData.columns.values[8]]
                P_ref = refData[refData.columns.values[5]]

            if equip.get(i).upper() == "NIMA":
                t_ref = refData[refData.columns.values[0]]
                P_ref = refData[refData.columns.values[5]]

            if equip.get(i).upper() == "CUSTOM":
                t_ref = refData[refData.columns.values[0]]
                P_ref = refData[refData.columns.values[1]]

            ## isolate region of interest, a function that truncates P according to config parameters
            injROI = config.injROI[i]
            P[i]   = P.get(i)[injROI[0]:injROI[1]]
            t[i]   = np.linspace(0,injROI[1],len(P.get(i)),endpoint=True)

        # pass list, get a dict of cycles, store each dict in a masterDict
        if analysisType in ['plotIsotherm','plotCompressions','plotExpansions','plotCycles','plotElasticity']:

            # calculates area per molecule for given file i
            Am[i] = calcAreaPerMolecule(i, A.get(i), lipidType, nLipids, lipidRatio, conc, volAdded)

            cycleAm, cycleP, cycleL = findCycles(Am.get(i), P.get(i), l.get(i))
            master_Am[i] = cycleAm
            master_P[i]  = cycleP
            master_L[i]  = cycleL

            # selects cycles based on config input
            # iterate along each chosen cycle and list of values within each cycle

            for cycleNum in config.useCycles:

                # append every value to file i of stichedX
                for value in range(len(cycleAm.get(cycleNum))):
                    stitched_Am[i].append( cycleAm.get(cycleNum)[value] )
                    stitched_P[i].append( cycleP.get(cycleNum)[value] )

                for k in range(len(stitched_Am)):
                    stitched_t[i].append(t[k])


        # calculate elasticity
        if analysisType == 'plotElasticity':

            # pass dict of cycles to splitCycles, returns dicts containing all half cycles in file
            suffix = " - compressions"
            halfCycleAm, halfCycleP, halfCycleL = splitCycles(suffix, cycleAm, cycleP, l.get(i))
            elasticityCompression_Am[i] = halfCycleAm.get(0)
            elasticityCompression_P[i]  = halfCycleP.get(0)
            elasticityCompression_L[i]  = halfCycleL.get(0)

            # pass the compression of cycle 0 as list; later adapt to selected cycle in config
            E[i] = calcElasticity(halfCycleAm.get(0), halfCycleP.get(0))


 ## Plot instructions
    if analysisType == 'plotIsotherm':
        key      = (1,1)
        axLabels = {"x": "Molecular Area ($\AA$$^2$ / Molecule)", "y": "$\pi$ ($mNm^{-1}$)"}
        suffix   = " - isotherm"
        vars     = [nFiles, equip, l, axLabels, title, plotDIR, [stitched_Am,0], stitched_P]

    if analysisType == 'plotElasticity':
        axLabels = {"x": "$\pi$ ($mNm^{-1}$)", "x1": "Mol. Area ($\AA$$^2$/Molecule)", "y": "$C_s^{-1} (mNm^{-1})$"}
        suffix   = " - elasticity"
        key = (1,2); vars = [nFiles, equip, l, axLabels, title, plotDIR, [elasticityCompression_P,elasticityCompression_Am], E]

    if analysisType  == 'plotPressure':
        key      = (1,1)
        axLabels = {"x": "auto calculated in genPlot", "y": "$\pi$ ($mNm^{-1}$)"}
        suffix   = " - pressure"
        vars     = [nFiles, equip, l, axLabels, title, plotDIR, [t,0], P]

    if analysisType == 'plotArea':
        key      = (1,1)
        axLabels = {"x": "auto calculated in genPlot", "y": "$\Delta$A (%)"}
        suffix   = " - area"
        vars     = [nFiles, equip, l, axLabels, title, plotDIR, [t,0], percA]


    # set default settings for subsequent similar fits
    key      = (1,1)
    axLabels = {"x": "Molecular Area ($\AA$$^2$ / Molecule)", "y": "$\pi$ ($mNm^{-1}$)"}

    for i in range(nFiles):
        newTitle = title + " - " + l.get(i)

        if analysisType == 'plotCompressions':
            suffix = " - compressions"
            halfCycleAm, halfCycleP, halfCycleL = splitCycles(suffix, master_Am.get(i), master_P.get(i), l.get(i))
            vars = [len(master_Am.get(i)), equip, halfCycleL, axLabels, newTitle, plotDIR, [halfCycleAm,0], halfCycleP]

        if analysisType == 'plotExpansions':
            suffix = " - expansions"
            halfCycleAm, halfCycleP, halfCycleL = splitCycles(suffix, master_Am.get(i), master_P.get(i), l.get(i))
            vars = [len(master_Am.get(i)), equip, halfCycleL, axLabels, newTitle, plotDIR, [halfCycleAm,0], halfCycleP]

        if analysisType == 'plotCycles':
            suffix = " - cycles"
            vars = [len(master_Am.get(i)), equip, master_L.get(i), axLabels, newTitle, plotDIR, [stitched_Am,0], stitched_P]

    return key, vars, suffix



if __name__ == '__main__':
    print("~Running isoAnalysis.py~\n")
    main()
