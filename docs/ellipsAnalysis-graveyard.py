" Plots ellipsometry data from beaglehole ellipsometer "

import sys
import csv
import re
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from statistics import mean
import warnings
import cmath as cm
import math
import scipy.signal as sig
import scipy.optimize as opt
import config
import genPlot
from genFunc import modSelection, getEllipsometryFile, getFile


def importSampleInfo(instructionsFile):

    # number of isotherms to plot
    nFiles  = len(instructionsFile)

    # create ID and count structures for later parsing
    AOI_ID  = []; nAOI  = 0
    time_ID = []; nTime = 0

    # separate input data into relevant ID lists of dicts
    for i in range(nFiles):

        if instructionsFile["refname"][i].upper() == "NULL":
            AOI_ID.append({'fname': instructionsFile["filename"][i], 'label': instructionsFile["label"][i]})
            nAOI += 1

        else:
            time_ID.append({'fname': instructionsFile["filename"][i], 'refname': instructionsFile["refname"][i], 'label': instructionsFile["label"][i]})
            nTime += 1

    return AOI_ID, nAOI, time_ID, nTime



def importAOIData(inputDIR, AOI_ID, nAOI):

    # initialise y axis variables for angle of incidence plots
    AOI       = {new_list: [] for new_list in range(nAOI)}
    psi_AOI   = {new_list: [] for new_list in range(nAOI)}
    delta_AOI = {new_list: [] for new_list in range(nAOI)}
    label_AOI = {new_list: [] for new_list in range(nAOI)}

    # extract angle of incidence data
    for i in range(nAOI):

        fileDIR = inputDIR + '/' + AOI_ID[i].get('fname') +".txt"
        AOI_df  = getEllipsometryFile(fileDIR)

        # for each element in AOI data frame, parse data
        for j in range(0,len(AOI_df)):
            AOI[i]      .append(float(AOI_df[j][5]))
            psi_AOI[i]  .append(float(AOI_df[j][0]))
            delta_AOI[i].append(float(AOI_df[j][1]))

        label_AOI[i] = AOI_ID[i].get('label')

    return AOI, psi_AOI, delta_AOI, label_AOI



def importTimeData(inputDIR, time_ID, nTime):

    # initialise variables for time plots
    psi_ref   = []
    delta_ref = []
    t         = {new_list: [] for new_list in range(nTime)}
    psi_t     = {new_list: [] for new_list in range(nTime)}
    delta_t   = {new_list: [] for new_list in range(nTime)}
    label_t   = {new_list: 0 for new_list in range(nTime)}


    for i in range(nTime):

        ## Extract reference data
        fileDIR = inputDIR + '/' + time_ID[i].get('refname') +".ds.dat"#".txt"
        #ref_df  = getEllipsometryFile(fileDIR)
        ref_df  = getFile(fileDIR,0,"\t")

        # extract raw ref data from dataframe
        #for j in range(0,len(ref_df)):
            #psi_ref  .append(float(ref_df[j][0]))
            #delta_ref.append(float(ref_df[j][1]))

        psi_ref   = ref_df[ref_df.columns.values[4]]
        delta_ref = ref_df[ref_df.columns.values[3]]

        # average reference values; each list ele is an averaged ref
        avRefPsi = mean(psi_ref)
        avRefDel = mean(delta_ref)

        ## Extract time series data
        fileDIR = inputDIR + '/' + time_ID[i].get('fname') + ".ds.dat"#".txt"
        #time_df = getEllipsometryFile(fileDIR)
        time_df  = getFile(fileDIR,0,"\t")

        # subtract ref data from these to give lipid-only data
        #for j in range(0,len(time_df)):
            #t[i]      .append(float(time_df[j][6]))
        t[i] = time_df[time_df.columns.values[11]]


        if config.ellipSubtractBufferRef == True:
                #psi_t[i]  .append(float(time_df[j][0]) - avRefPsi)
                #delta_t[i].append(float(time_df[j][1]) - avRefDel)
            psi_t[i]  = time_df[time_df.columns.values[4]] - avRefPsi
            delta_t[i] = time_df[time_df.columns.values[3]] - avRefDel

        else:
                #psi_t[i]  .append(float(time_df[j][0]) )
                #delta_t[i].append(float(time_df[j][1]) )
            psi_t[i]   = time_df[time_df.columns.values[4]]
            delta_t[i] = time_df[time_df.columns.values[3]]

        label_t[i] = time_ID[i].get('label')

    return t, psi_t, delta_t, label_t



def main(instructionsFile, title, inputDIR, outputPath):

    # filter warnings
    warnings.filterwarnings("ignore")


    # give analysis choice to user
    analysisOptions = ['plotAOIpsi','plotAOIdelta','plotTimePsi','plotTimeDelta']

    analysisRunning = True
    while analysisRunning:
        analysisType, analysisRunning = modSelection(analysisOptions)

        if analysisRunning == False:
            break


        # user input & sample instructions information
        AOI_ID, nAOI, time_ID, nTime = importSampleInfo(instructionsFile)


        # process data
        if analysisType == 'plotAOIpsi' or analysisType == 'plotAOIdelta':

            if nAOI < 1:
                print("Error: No AOI files in input.")
                sys.exit()

            # import data
            AOI, psi_AOI, delta_AOI, label_AOI = importAOIData(inputDIR, AOI_ID, nAOI)

            # converting delta variables to degrees
            for i in range(nAOI):
                for j in range(len(delta_AOI.get(i))):
                    delta_AOI[i][j] = delta_AOI.get(i)[j] * 180 / np.pi



        if analysisType == 'plotTimePsi' or analysisType == 'plotTimeDelta':

            if nTime < 1:
                print("Error: No Time files in input.")
                sys.exit()

            # import data
            t, psi_t, delta_t, label_t = importTimeData(inputDIR, time_ID, nTime)

            # converting delta variables to degrees
            #for i in range(nTime):
            #    for j in range(len(delta_t.get(i))):
            #        delta_t[i][j] = delta_t.get(i)[j] * 180 / np.pi


        # plot data
        equip = 'null'
        key   = (1,1)
        if analysisType == 'plotAOIpsi':
            axLabels = {"x": "AOI (\N{DEGREE SIGN})", "y": "$\Psi$ (\N{DEGREE SIGN})"}
            suffix   = " - psi AOI"
            vars     = (nAOI, equip, label_AOI, axLabels, title, outputPath, (AOI,0), psi_AOI)
            genPlot.main(key,vars,suffix)

        if analysisType == 'plotAOIdelta':
            axLabels = {"x": "AOI (\N{DEGREE SIGN})", "y": "$\Delta$ (\N{DEGREE SIGN})"}
            suffix   = " - delta AOI"
            vars     = (nAOI, equip, label_AOI, axLabels, title, outputPath, (AOI,0), delta_AOI)
            genPlot.main(key,vars,suffix)

        if analysisType == 'plotTimePsi':
            axLabels = {"x": "Time (s)", "y": "$\Psi$ (\N{DEGREE SIGN})"}
            suffix   = " - psi Time"
            vars     = (nTime, equip, label_t, axLabels, title, outputPath, (t,0), psi_t)
            genPlot.main(key,vars,suffix)

        if analysisType == 'plotTimeDelta':
            axLabels = {"x": "Time (s)", "y": "$\Delta_{Lipids}$ - $\Delta_{Buffer}$ (\N{DEGREE SIGN})"}
            suffix   = " - delta Time"
            vars     = (nTime, equip, label_t, axLabels, title, outputPath, (t,0), delta_t)
            genPlot.main(key,vars,suffix)


        # program executed
        print('\nAnalysis Complete!\n')
        sys.exit()
    return



if __name__ == '__main__':
    print("~Running plotEllips.py~\n")
    main()
