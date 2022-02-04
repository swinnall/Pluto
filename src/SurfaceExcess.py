" Module for calculating surface excess "

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


def modSelection(analysisOptions):

    # ask user to pick one of the analysisOptions
    print("\n~~~\nAnalysis Options:")
    for i,option in enumerate(analysisOptions):
        print("%d: %s" %(i+1,option))
    print("~~~\n")
    
    analysisChoice = input("Which analysis would you like to do? Pick the associated number (1-%d) or 'q' to reutrn to landing page:\n  " %len(analysisOptions) )

    if analysisChoice == 'q':
        print("Returning to Pluto landing page.\n\n")
        analysisType    = 'n/a'
        analysisRunning = False

    elif analysisChoice == 'Q':
        print("Session closed.")
        sys.exit()

    elif analysisChoice in [str(i) for i in range(len(analysisOptions)+1)]:
        analysisType = analysisOptions[int(analysisChoice)-1]
        print("You picked %s.py\n" %analysisType)
        analysisRunning = True

    else:
        print("Not a valid response. Returning to Pluto landing page.\n\n")
        analysisType    = 'n/a'
        analysisRunning = False

    return analysisType, analysisRunning


## Import Functions
def importSampleData(instructionsFile, nFiles):

    # variable initialisations
    fileNames  = {init: 0 for init in range(nFiles)}
    l          = {init: 0 for init in range(nFiles)}

    # assign data
    for i in range(nFiles):

        fileNames[i]  = instructionsFile["fname"][i]
        l[i]          = instructionsFile["label"][i]

    return fileNames, l



def getFile(fileDIR,equipParams):

    with open(fileDIR, newline = '') as f:
        reader = csv.reader(f, delimiter=equipParams[2])
        data = list(reader)

    # filter out empty lines (seemingly randomly introduced in Nima files)
    data = [x for x in data if x != []]

    return data



def importData(equipParams, fname, inputDIR, plotDIR):

    # initialise variables
    t      = []
    gammaL = []
    gammaP = []

    # read file into memory data
    fileDIR = inputDIR + '/' + fname + equipParams[1]
    data = getFile(fileDIR,equipParams)

    for j in range(equipParams[0],len(data)):
        t.append(float(data[j][0]))
        gammaL.append(float(data[j][1]))
        gammaP.append(float(data[j][2]))

    return t, gammaL, gammaP


def reduceNpoints(x, y):

    # initialise lists
    redX = []
    redY = []

    for j in range(0,len(x),config.gammaNth):
        redX.append(x[j])
        redY.append(y[j])

    return redX, redY


def smoothData(genList):

    # initialise list for func output
    smoothList = []


    # if length is even, truncate as savgol must have odd window length
    if len(genList) % 2 == 0:
        genList.pop()


    # remove noise introduced via derivative perturbations
    try:
        smoothList = savgol_filter(genList, len(genList), config.gamma_nPoly, mode='interp')  # variable, window size (length), polynomial order, mode

    # remove infinity and nan values due to Nima having repeated area values
    except linalg.LinAlgError:
        print("LinAlgError")
        for j in range(len(genList)):
            if genList[j] == float("inf") or genList[j] == float("-inf") or mt.isnan(genList[j]) == True:
                genList[j] = 0

    except ValueError:
        print("ValueError: Polyorder must be less than window length. Attempting adjustment at 0.5*nPoly = %d." %int(config.nPoly/2))
        smoothList = savgol_filter(genList, len(genList), int(config.nPoly/2), mode='interp')  # variable, window size (length), polynomial order, mode


    # high pass filter: remove values less than 0.01 (final smoothing of function)
    deleteMe = []
    for j in range(1,len(smoothList)-1):

        if smoothList[j] < 0.01:

            # store indices to be deleted
            deleteMe.append(j)

    # delete elements from the lists
    smoothList = np.delete(smoothList, deleteMe)

    return smoothList




# Main
def main(instructionsFile, title, inputDIR, plotDIR):

    # filter warnings
    warnings.filterwarnings("ignore")


    # give analysis choice to user
    analysisOptions = ['plotGammaL','plotGammaP']

    analysisRunning = True
    while analysisRunning:
        analysisType, analysisRunning = modSelection(analysisOptions)

        if analysisRunning == False:
            break

        # number of files to plot
        nFiles  = len(instructionsFile)

        # user input & sample instructions information
        fileNames, l = importSampleData(instructionsFile, nFiles)

        # initialise dicts, store 1 file in each list
        t      = {new_list: [] for new_list in range(nFiles)}
        gammaL = {new_list: [] for new_list in range(nFiles)}
        gammaP = {new_list: [] for new_list in range(nFiles)}

        ## perform analysis on each file and store in dict for genPlot.py
        for i in range(nFiles):

            # get filename
            fname = fileNames.get(i)

            # organise equipment
            equipParams = (1,".txt","\t")

            # import data from given file i
            t_list, gammaL_list, gammaP_list = importData(equipParams, fname, inputDIR, plotDIR)


            # smooth data
            if analysisType == 'plotGammaL':

                # reduce number of points
                reduced_t, reduced_L = reduceNpoints(t_list, gammaL_list)
                t[i]      = reduced_t
                gammaL[i] = reduced_L

                # calculate smoothed function for plotting
                #temp = smoothData(gammaL_list)
                #gammaL[2*i+1] = temp


            if analysisType == 'plotGammaP':
                reduced_t, reduced_P = reduceNpoints(t_list, gammaP_list)
                t[i]      = reduced_t
                gammaP[i] = reduced_P

                #temp = smoothData(gammaP_list)
                #gammaP[2*i+1] = temp



	    ## Plot instructions
        if analysisType == 'plotGammaL':
            key      = (1,1)
            axLabels = {"x": "Time (min)", "y": "$\Gamma_{Lipid}$"}
            suffix   = " - gammaL"
            equip    = "N/A"
            vars     = (nFiles, equip, l, axLabels, title, plotDIR, (t,0), gammaL)

            # currently only allows 1x1, 2x1, 2x2 subplot types
            if config.plotMultiPanel == True:
                nRow = len(config.key)
                nCol = len(config.key[0])
                key = (nRow,nCol)
            genPlot.main(key,vars,suffix)


        if analysisType == 'plotGammaP':
            key      = (1,1)
            axLabels = {"x": "Time (min)", "y": "$\Gamma_{PolyA}$ (molecule $m^{-2}$)"}
            suffix   = " - gammaP"
            equip    = "N/A"
            vars     = (nFiles, equip, l, axLabels, title, plotDIR, (t,0), gammaP)

            # currently only allows 1x1, 2x1, 2x2 subplot types
            if config.plotMultiPanel == True:
                nRow = len(config.key)
                nCol = len(config.key[0])
                key = (nRow,nCol)
            genPlot.main(key,vars,suffix)


        # program executed
        print('\nAnalysis Complete!\n')

    return



if __name__ == '__main__':
    print("~Running gammaAnalysis.py~\n")
    main()