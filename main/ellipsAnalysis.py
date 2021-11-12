"Reads ellipsometry data"

import sys
import csv
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from statistics import mean
import warnings
import config
import genPlot


def importSampleInfo(info):

    # number of header rows
    nHeaders = 2

    # number of isotherms to plot
    nFiles  = len(info) - nHeaders

    # create ID and count structures for later parsing
    AOI_ID  = []; nAOI  = 0
    time_ID = []; nTime = 0
    ref_ID  = []; nRef  = 0

    # separate input data into relevant ID lists of dicts
    for i in range(nHeaders,nHeaders+nFiles):

        if info[i][1] in ["AOI", "AoI", "aoi"]:
            AOI_ID.append({'fname': info[i][0], 'label': info[i][2]})
            nAOI += 1

        if info[i][1] in ["T", "t", "Time", 'time']:
            time_ID.append({'fname': info[i][0], 'label': info[i][2]})
            nTime += 1

        if info[i][1] in ["Ref", "ref"]:
            ref_ID.append({'fname': info[i][0], 'label': info[i][2]})
            nRef += 1

    return nFiles, AOI_ID, nAOI, time_ID, nTime, ref_ID, nRef



def getFile(fileDIR):

    df = []
    with open(fileDIR, newline = '') as f:

        # ignore rows starting with '#'
        rdr = csv.DictReader(filter(lambda row: row[0]!='#', f))

        for row in rdr:
            df.append( row[None][0].split())

    return df



def importAOIData(inputDIR, AOI_ID, nAOI):

    # initialise y axis variables for angle of incidence plots
    AOI       = {new_list: [] for new_list in range(nAOI)}
    psi_AOI   = {new_list: [] for new_list in range(nAOI)}
    delta_AOI = {new_list: [] for new_list in range(nAOI)}
    label_AOI = {new_list: [] for new_list in range(nAOI)}

    # extract angle of incidence data
    for i in range(nAOI):

        fileDIR = inputDIR + '/' + AOI_ID[i].get('fname') +".txt"
        AOI_df  = getFile(fileDIR)

        # for each element in AOI data frame, parse data
        for j in range(0,len(AOI_df)):
            AOI[i]      .append(float(AOI_df[j][5]))
            psi_AOI[i]  .append(float(AOI_df[j][0]))
            delta_AOI[i].append(float(AOI_df[j][1]))

        label_AOI[i] = AOI_ID[i].get('label')

    return AOI, psi_AOI, delta_AOI, label_AOI



def importTimeData(inputDIR, time_ID, nTime, ref_ID, nRef):

    # initialise variables for time plots
    t         = {new_list: [] for new_list in range(nTime)}
    psi_t     = {new_list: [] for new_list in range(nTime)}
    delta_t   = {new_list: [] for new_list in range(nTime)}
    label_t   = {new_list: 0 for new_list in range(nTime)}

    psi_ref   = {new_list: [] for new_list in range(nRef)}
    delta_ref = {new_list: [] for new_list in range(nRef)}

    # averaged reference values for system i->nTime
    avRefPsi = []
    avRefDel = []


    # extract reference data
    for i in range(nRef):

        fileDIR = inputDIR + '/' + ref_ID[i].get('fname') +".txt"
        ref_df  = getFile(fileDIR)

        # extract data; subtract ref data from lipid measurements
        for j in range(0,len(ref_df)):
            t[i]        .append(float(ref_df[j][6]))
            psi_ref[i]  .append(float(ref_df[j][0]))
            delta_ref[i].append(float(ref_df[j][1]))


    # average reference values
    for i in range(nTime):
        avRefPsi.append(mean(psi_ref.get(i)))
        avRefDel.append(mean(delta_ref.get(i)))


    # extract time series data
    for i in range(nTime):

        fileDIR = inputDIR + '/' + time_ID[i].get('fname') +".txt"
        time_df = getFile(fileDIR)

        # subtract ref data from these to give lipid-only data
        for j in range(0,len(time_df)):
            t[i]      .append(float(time_df[j][6]))
            psi_t[i]  .append(float(time_df[j][0]) - avRefPsi[i])
            delta_t[i].append(float(time_df[j][1]) - avRefDel[i])

        label_t[i] = time_ID[i].get('label')

    return t, psi_t, delta_t, label_t



def main(info, title, inputDIR, outputPath):

    # filter warnings
    warnings.filterwarnings("ignore")

    # set equip to null - remove in future version?
    equip = 'null'

    # user input & sample instructions information
    nFiles, AOI_ID, nAOI, time_ID, nTime, ref_ID, nRef = importSampleInfo(info)


    # process data
    if config.plotAOI == True:

        if nAOI < 1:
            print("Error: No AOI files in input.")
            sys.exit()

        # import data
        AOI, psi_AOI, delta_AOI, label_AOI = importAOIData(inputDIR, AOI_ID, nAOI)

        # converting delta variables to degrees
        for i in range(nAOI):
            for j in range(len(delta_AOI.get(i))):
                delta_AOI[i][j] = delta_AOI.get(i)[j] * 180 / np.pi

        # plot data
        key      = (1,1)
        axLabels = {"x": "AOI (\N{DEGREE SIGN})", "y": "$\Psi$ (\N{DEGREE SIGN})"}
        suffix   = " - psi AOI"
        vars     = (nFiles, equip, label_AOI, axLabels, title, outputPath, (AOI,0), psi_AOI)
        genPlot.main(key,vars,suffix)

        axLabels = {"x": "AOI (\N{DEGREE SIGN})", "y": "$\Delta$ (\N{DEGREE SIGN})"}
        suffix   = " - delta AOI"
        vars     = (nFiles, equip, label_AOI, axLabels, title, outputPath, (AOI,0), delta_AOI)
        genPlot.main(key,vars,suffix)


    if config.plotTime == True:

        if nTime < 1:
            print("Error: No Time files in input.")
            sys.exit()

        # import data
        t, psi_t, delta_t, label_t = importTimeData(inputDIR, time_ID, nTime, ref_ID, nRef)

        # converting delta variables to degrees
        for i in range(nTime):
            for j in range(len(delta_t.get(i))):
                delta_t[i][j] = delta_t.get(i)[j] * 180 / np.pi

        # plot data
        key      = (1,1)
        axLabels = {"x": "Time (s)", "y": "$\Psi$ (\N{DEGREE SIGN})"}
        suffix   = " - psi Time"
        vars     = (nFiles, equip, label_t, axLabels, title, outputPath, (t,0), psi_t)
        genPlot.main(key,vars,suffix)

        axLabels = {"x": "Time (s)", "y": "$\Delta_{Lipids}$ - $\Delta_{Buffer}$ (\N{DEGREE SIGN})"}
        suffix   = " - delta Time"
        vars     = (nFiles, equip, label_t, axLabels, title, outputPath, (t,0), delta_t)
        genPlot.main(key,vars,suffix)


    # program executed
    print('\nAnalysis Complete! Have a nice day :)')

    return



if __name__ == '__main__':
    print("~Running plotEllips.py~\n")
    main()
