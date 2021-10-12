"Reads and plots ellipsometry data"

import csv
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from statistics import mean
import config
import genPlot


def importSampleData(info):

    # number of header rows
    nHeaders = 2

    # number of isotherms to plot
    nFiles  = len(info) - nHeaders

    # variable initialisations
    ID      = {init: 0 for init in range(nFiles)}
    tag     = {init: 0 for init in range(nFiles)}
    uniqIDs = {"AOI": [], "t": [], "ref": [],}

    # assign data
    for i in range(nHeaders,len(info)):

        # correct for indexing
        j = i - 1

        ID[j]  = info[i][0]
        tag[j] = info[i][2]


    for i in range(nFiles):

        if tag.get(i) == "AOI":
            uniqIDs["AOI"].append(int(ID.get(i)))

        elif tag.get(i) == "t":
            uniqIDs["t"].append(int(ID.get(i)))

        elif tag.get(i) == "ref":
            uniqIDs["ref"].append(int(ID.get(i)))


    # determine number of each types of files for plotting
    nAOI  = len(uniqIDs.get("AOI"))
    nTime = len(uniqIDs.get("t"))
    nRef  = len(uniqIDs.get("ref"))

    return nFiles, uniqIDs, tag, nAOI, nTime, nRef, info


def importAOI(uniqIDs, info, nAOI):

    # initialise y axis variables for angle of incidence plots
    AOI       = {new_list: [] for new_list in range(nAOI)}
    psi_AOI   = {new_list: [] for new_list in range(nAOI)}
    delta_AOI = {new_list: [] for new_list in range(nAOI)}


    # extract angle of incidence data
    for i in range(nAOI):

        # initialise data list for each file
        d_AOI = []

        # read file into memory d
        with open( "../input/00/" + info[uniqIDs.get("AOI")[i]][1] +".txt", newline = '') as f:

            # ignore rows starting with '#'
            rdr = csv.DictReader(filter(lambda row: row[0]!='#', f))

            # split rows into iterable lists
            for row in rdr:
                d_AOI.append( row[None][0].split())

        for j in range(0,len(d_AOI)):
            AOI[i].append(float(d_AOI[j][5]))
            psi_AOI[i].append(float(d_AOI[j][0]))
            delta_AOI[i].append(float(d_AOI[j][1]))

    return AOI, psi_AOI, delta_AOI


def importTime(uniqIDs, info, nTime, nRef):

    # initialise y axis variables for time plots
    t         = {new_list: [] for new_list in range(nTime)}
    psi_t     = {new_list: [] for new_list in range(nTime)}
    delta_t   = {new_list: [] for new_list in range(nTime)}
    psi_ref   = {new_list: [] for new_list in range(nRef)}
    delta_ref = {new_list: [] for new_list in range(nRef)}

    # extract reference data
    for i in range(nRef):

        # initialise data list for each file
        d_ref = []

        # read file into memory d
        with open( "../input/00/" + info[uniqIDs.get("ref")[i]][1] +".txt", newline = '') as f:

            # ignore rows starting with '#'
            rdr = csv.DictReader(filter(lambda row: row[0]!='#', f))

            # split rows into iterable lists
            for row in rdr:
                d_ref.append( row[None][0].split())

        # extract data; subtract ref data from lipid measurements
        for j in range(0,len(d_ref)):
            t[i].append(float(d_ref[j][6]))
            psi_ref[i].append(float(d_ref[j][0]))
            delta_ref[i].append(float(d_ref[j][1]))


    # stores averaged reference values for system i->nTime
    avRefPsi = []
    avRefDel = []

    # average reference values
    for i in range(nTime):
        avRefPsi.append(mean(psi_ref.get(i)))
        avRefDel.append(mean(delta_ref.get(i)))


    # extract time series data
    for i in range(nTime):

        # initialise data list for each file
        d_t = []

        # read file into memory d
        with open( "../input/00/" + info[uniqIDs.get("t")[i]][1] +".txt", newline = '') as f:

            # ignore rows starting with '#'
            rdr = csv.DictReader(filter(lambda row: row[0]!='#', f))

            # split rows into iterable lists
            for row in rdr:
                d_t.append( row[None][0].split())

        # subtract ref data from these to give lipid-only data
        for j in range(0,len(d_t)):
            t[i].append(float(d_t[j][6]))
            psi_t[i].append(float(d_t[j][0])   - avRefPsi[i])
            delta_t[i].append(float(d_t[j][1]) - avRefDel[i])

    return t, psi_t, delta_t


def main(info, title, outputPath):

    # filter warnings
    # warnings.filterwarnings("ignore")


    # user input & sample instructions information
    nFiles, uniqIDs, tag, nAOI, nTime, nRef, info = importSampleData(info)


    # process data
    if config.plotAOI == True:

        # import data
        AOI, psi_AOI, delta_AOI = importAOI(uniqIDs, info, nAOI)

        # converting delta variables to degrees
        for i in range(nAOI):
            for j in range(len(delta_AOI.get(i))):
                delta_AOI[i][j] = delta_AOI.get(i)[j] * 180 / np.pi

        # plot data
        key      = 1
        axLabels = {"x": "AOI (\N{DEGREE SIGN})", "y": "$\Psi$ (\N{DEGREE SIGN})"}
        suffix   = " - psi AOI"
        vars     = (nFiles, equip, l, axLabels, suffix, title, plotDIR, (AOI,0), psi_AOI)
        genPlot.main(key,vars)

        axLabels = {"x": "AOI (\N{DEGREE SIGN})", "y": "$\Delta$ (\N{DEGREE SIGN})"}
        suffix   = " - delta AOI"
        vars     = (nFiles, equip, l, axLabels, suffix, title, plotDIR, (AOI,0), delta_AOI)
        genPlot.main(key,vars)


    if config.plotTime == True:

        # import data
        t, psi_t, delta_t = importTime(uniqIDs, info, nTime, nRef)

        # converting delta variables to degrees
        for i in range(nTime):
            for j in range(len(delta_t.get(i))):
                delta_t[i][j] = delta_t.get(i)[j] * 180 / np.pi

        # plot data
        key      = 1
        axLabels = {"x": "Time (s)", "y": "$\Psi$ (\N{DEGREE SIGN})"}
        suffix   = " - psi Time"
        vars     = (nFiles, equip, l, axLabels, suffix, title, plotDIR, (t,0), psi_t)
        genPlot.main(key,vars)

        axLabels = {"x": "Time (s)", "y": "$\Delta_{Lipids}$ - $\Delta_{Buffer}$ (\N{DEGREE SIGN})"}
        suffix   = " - delta Time"
        vars     = (nFiles, equip, l, axLabels, suffix, title, plotDIR, (t,0), delta_t)
        genPlot.main(key,vars)


    # program executed
    print('\nAnalysis Complete! Have a nice day :)')

    return



if __name__ == '__main__':
    print("~Running plotEllips.py~\n")
    main()
