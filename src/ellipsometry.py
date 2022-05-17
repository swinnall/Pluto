" Process ellipsometry data "

# import packages
import sys
import warnings
import numpy as np
from statistics import mean

# import Pluto modules
import config
import genPlot
from genFunc import modSelection, getFile


def importSampleInfo(instructionsFile):

    # number of isotherms to plot
    nFiles  = len(instructionsFile)

    # store measurement info in file list
    fileInfoList  = []

    # separate input data into relevant ID lists of dicts
    for i in range(nFiles):
        fileInfoList.append({'fname': instructionsFile["filename"][i], 'refname': instructionsFile["refname"][i], 'label': instructionsFile["label"][i]})

    return nFiles, fileInfoList



def importData(inputDIR, nFiles, fileInfoList):

    # initialise variables for time plots
    refPsi   = []
    refDelta = []
    x        = {new_list: [] for new_list in range(nFiles)}
    yPsi     = {new_list: [] for new_list in range(nFiles)}
    yDelta   = {new_list: [] for new_list in range(nFiles)}
    label    = {new_list: 0 for new_list in range(nFiles)}


    for i in range(nFiles):

        ## Extract time series data
        fileDIR = inputDIR + '/' + fileInfoList[i].get('fname') + ".ds.dat"#".txt"
        file_df  = getFile(fileDIR,0,"\t")

        # store label information
        label[i] = fileInfoList[i].get('label')

        # store x-axis data (time or AOI)
        x[i] = file_df[file_df.columns.values[11]]

        # store corresponding y-axis data
        yPsi[i]   = file_df[file_df.columns.values[4]]
        yDelta[i] = file_df[file_df.columns.values[3]]

        # if delta < 10 it must be in radians therefore convert to degrees
        if yDelta.get(i)[0] < 10:
            for j in range(len(yDelta.get(i))):
                    yDelta[i][j] = yDelta.get(i)[j] * 180 / np.pi

        #print(yDelta.get(i))
        # for instructions where reference data was provided: subtract from measurement
        if fileInfoList[i].get('refname').upper() != "NULL":

            # get reference data file as dataframe
            fileDIR = inputDIR + '/' + fileInfoList[i].get('refname') +".ds.dat"#".txt"
            ref_df  = getFile(fileDIR,0,"\t")

            # extract raw ref data from dataframe
            refPsi   = ref_df[ref_df.columns.values[4]]
            refDelta = ref_df[ref_df.columns.values[3]]


            if config.ellipsBufferRef == True:

                # average reference values
                avRefPsi = mean(refPsi)
                avRefDel = mean(refDelta)

                # convert to degrees if necessary
                if yDelta.get(i)[0] < 10:
                    avRefPsi = avRefPsi * 180 / np.pi
                    avRefDel = avRefDel * 180 / np.pi

                # subtract reference from measurement
                yPsi[i]   = file_df[file_df.columns.values[4]] - avRefPsi
                yDelta[i] = file_df[file_df.columns.values[3]] - avRefDel

            elif config.ellipsLipidRef == True:

                # subtract reference from measurement
                yPsi[i]   = file_df[file_df.columns.values[4]] - refPsi
                yDelta[i] = file_df[file_df.columns.values[3]] - refDelta

            else:
                print("Fatal Error: Expected Null reference but found file.\n [Buffer & lipid ref == False]")
                sys.exit()
                
    return x, yPsi, yDelta, label



def main(instructionsFile, title, inputDIR, outputPath):

    # filter warnings
    warnings.filterwarnings("ignore")

    # give analysis choice to user
    analysisOptions = ['plotAOIpsi','plotAOIdelta','plotTimePsi','plotTimeDelta']

    analysisType, analysisRunning = modSelection(analysisOptions)

    # user input & sample instructions information
    nFiles, fileInfoList = importSampleInfo(instructionsFile)

    if nFiles < 1:
        print("Fatal Error: No file input.")
        sys.exit()

    # import data
    x, yPsi, yDelta, label = importData(inputDIR, nFiles, fileInfoList)

    # plot data
    equip = 'null'
    key   = (1,1)
    if analysisType == 'plotAOIpsi':
        axLabels = {"x": "AOI (\N{DEGREE SIGN})", "y": "$\Psi$ (\N{DEGREE SIGN})"}
        suffix   = " - psi AOI"
        vars     = [nFiles, equip, label, axLabels, title, outputPath, [x,0], yPsi]

    if analysisType == 'plotAOIdelta':
        axLabels = {"x": "AOI (\N{DEGREE SIGN})", "y": "$\Delta$ (\N{DEGREE SIGN})"}
        suffix   = " - delta AOI"
        vars     = [nFiles, equip, label, axLabels, title, outputPath, [x,0], yDelta]

    if analysisType == 'plotTimePsi':
        axLabels = {"x": "Time (s)", "y": "$\Psi$ (\N{DEGREE SIGN})"}
        suffix   = " - psi Time"
        vars     = [nFiles, equip, label, axLabels, title, outputPath, [x,0], yPsi]

    if analysisType == 'plotTimeDelta':
        axLabels = {"x": "Time (s)", "y": "$\Delta_{Lipids}$ - $\Delta_{Buffer}$ (\N{DEGREE SIGN})"}
        suffix   = " - delta Time"
        vars     = [nFiles, equip, label, axLabels, title, outputPath, [x,0], yDelta]

    return key, vars, suffix



if __name__ == '__main__':
    print("~Running plotEllips.py~\n")
    main()
