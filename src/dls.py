" Process DLS data "

# import packages
import sys, os
import warnings
import numpy as np

# import Pluto modules
import config
from genFunc import modSelection, getFile


def importSampleInfo(instructionsFile):

    # number of isotherms to plot
    nSamples  = len(instructionsFile)

    # separate input data into relevant ID lists of dicts
    fileInfoList  = []
    for i in range(nSamples):
       fileInfoList.append({'fname': instructionsFile["filename"][i], 'standard': instructionsFile["standard"][i], 'solvent': instructionsFile["solvent"][i], 'label': instructionsFile["label"][i]})

    return nSamples, fileInfoList



def importCountData(nSamples, fileInfoList):

    ## hardcoding paths
    inputRoot = '../../UoM-Data-Repository/input/TR-DLS/'

    inputPaths = []
    for sampleNum in range(nSamples):
        inputPaths.append(inputRoot + fileInfoList[sampleNum].get("fname"))


    # initialise variables for time plots
    x        = {new_list: [] for new_list in range(nSamples)}
    y        = {new_list: [] for new_list in range(nSamples)}
    label    = {new_list: 0 for new_list in range(nSamples)}


    for sampleNum in range(nSamples):

        time   = []
        meanCR = []
        fileCount = 0
        for filename in os.listdir(inputPaths[sampleNum]):

            if filename.endswith(".ASC"):

                filepath = inputPaths[sampleNum] + '/' + filename
                with open(filepath) as f:

                    for line in f:
                        l = line.strip("\n").split("\t")

                        if l[0] == "Time :":
                            date_time = l[1].strip().strip('""')
                            sec_time = sum(x * int(t) for x, t in zip([3600, 60, 1], date_time.split(":")))
                            if fileCount == 0:
                                startTime = sec_time
                            time.append(sec_time - startTime)

                        if l[0] == "MeanCR0 [kHz]   :":
                            meanCR0 = float(l[1].strip())
                        if l[0] == "MeanCR1 [kHz]   :":
                            meanCR1 = float(l[1].strip())

                f.close()
                meanCR.append(meanCR0+meanCR1)
                fileCount += 1

        x[sampleNum] = time
        y[sampleNum] = meanCR
        label[sampleNum] = fileInfoList[sampleNum].get("label")

    return x, y, label



def main(instructionsFile, title, inputDIR, outputPath):

    # filter warnings
    warnings.filterwarnings("ignore")

    # give analysis choice to user
    analysisOptions = ['plotCount']

    analysisType, analysisRunning = modSelection(analysisOptions)

    # user input & sample instructions information
    nSamples, fileInfoList = importSampleInfo(instructionsFile)

    if nSamples < 1:
        print("Fatal Error: No file input.")
        sys.exit()

    # import count rate data
    x, y, label = importCountData(nSamples, fileInfoList)

    # plot data
    equip = 'null'
    key   = (1,1)
    if analysisType == 'plotCount':
        axLabels = {"x": "T", "y": "Mean Count Rate"}
        suffix   = " - TR DLS"
        vars     = [len(x), equip, label, axLabels, title, outputPath, [x,0], y]

    return key, vars, suffix



if __name__ == '__main__':
    main()
