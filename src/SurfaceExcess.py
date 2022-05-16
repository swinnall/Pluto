" Process surface excess data "

# import packages
import csv, os, sys, re, ast, sys
import warnings
import numpy as np
import math as mt

# import Pluto modules
import config
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



def importData(fname, inputDIR, analysisType):

    # read file into memory data
    fileDIR = inputDIR + '/' + fname + '.txt'

    # get sample data as pandas df
    data = getFile(path=fileDIR, nSkip=0, delim='\t')

    # parse variables
    t = data[data.columns.values[0]]

    if analysisType == "plotGammaL":
        gamma = data[data.columns.values[1]]
    if analysisType == "plotGammaP":
        gamma = data[data.columns.values[2]]

    return t, gamma



# Main
def main(instructionsFile, title, inputDIR, plotDIR):

    # filter warnings
    warnings.filterwarnings("ignore")

    # give analysis choice to user
    analysisOptions = ['plotGammaL','plotGammaP']

    # select analysis type
    analysisType, analysisRunning = modSelection(analysisOptions)

    # number of files to plot
    nFiles  = len(instructionsFile)

    # user input & sample instructions information
    fileNames, l = importSampleData(instructionsFile, nFiles)

    # initialise dicts, store 1 file in each list
    t     = {new_list: [] for new_list in range(nFiles)}
    gamma = {new_list: [] for new_list in range(nFiles)}

    ## perform analysis on each file and store in dict for genPlot.py
    for i in range(nFiles):

        # get filename
        fname = fileNames.get(i)

        # import data from given file i
        t[i], gamma[i] = importData(fname, inputDIR, analysisType)

    ## Plot instructions
    if analysisType == 'plotGammaL':
        key      = (1,1)
        axLabels = {"x": "Time (min)", "y": "$\Gamma_{Lipid}$ (molecule $m^{-2}$)"}
        suffix   = " - gammaL"
        equip    = "N/A"
        vars     = [nFiles, equip, l, axLabels, title, plotDIR, [t,0], gamma]

    if analysisType == 'plotGammaP':
        key      = (1,1)
        axLabels = {"x": "Time (min)", "y": "$\Gamma_{PolyA}$ (molecule $m^{-2}$)"}
        suffix   = " - gammaP"
        equip    = "N/A"
        vars     = [nFiles, equip, l, axLabels, title, plotDIR, [t,0], gamma]

    return key, vars, suffix



if __name__ == '__main__':
    main()
