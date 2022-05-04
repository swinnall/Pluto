" Module for plotting surface excess "

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
from genFunc import modSelection, getFile


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



def importData(fname, inputDIR, plotDIR):

    # read file into memory data
    fileDIR = inputDIR + '/' + fname + '.txt'

    # get sample data as pandas df
    data = getFile(path=fileDIR, nSkip=0, delim='\t')

    # parse variables
    t      = data[data.columns.values[0]]
    gammaL = data[data.columns.values[1]]
    gammaP = data[data.columns.values[2]]

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


            # import data from given file i
            t_list, gammaL_list, gammaP_list = importData(fname, inputDIR, plotDIR)


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
            axLabels = {"x": "Time (min)", "y": "$\Gamma_{Lipid}$ (molecule $m^{-2}$)"}
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
    sys.exit()
    return



if __name__ == '__main__':
    print("~Running gammaAnalysis.py~\n")
    main()
