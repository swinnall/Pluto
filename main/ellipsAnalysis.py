" Plot Ellipsometry "
" Plots both AOI and t plots at N number of pressures"

import csv
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from statistics import mean
import config


def importInstructions(info):

    # order of the polynomial fit
    nPoly = 40

    # read plotInstructions
    with open("../input/"+ info +".txt", newline = '') as f:
        reader = csv.reader(f, delimiter=",")
        info = list(reader)

    # number of files to read
    nFiles = len(info)-1

    # variable initialisations
    ID      = {init: 0 for init in range(nFiles)}
    tag     = {init: 0 for init in range(nFiles)}
    uniqIDs = {"AOI": [], "t": [], "ref": [],}


    # assign data
    for i in range(1,len(info)):

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
    nAOI = len(uniqIDs.get("AOI"))
    nSys = len(uniqIDs.get("ref"))

    return nPoly, nFiles, uniqIDs, tag, nAOI, nSys, info



def importData(uniqIDs, info, nAOI, nSys):

    # initialise y axis variables for angle of incidence plots
    AOI       = {new_list: [] for new_list in range(nAOI)}
    psi_AOI   = {new_list: [] for new_list in range(nAOI)}
    delta_AOI = {new_list: [] for new_list in range(nAOI)}

    # initialise y axis variables for time plots
    t         = {new_list: [] for new_list in range(nSys)}
    psi_ref   = {new_list: [] for new_list in range(nSys)}
    delta_ref = {new_list: [] for new_list in range(nSys)}
    psi_t     = {new_list: [] for new_list in range(nSys)}
    delta_t   = {new_list: [] for new_list in range(nSys)}


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


    # extract reference data
    for i in range(nSys):

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


    # stores averaged reference values for system i->nSys
    avRefPsi = []
    avRefDel = []

    # average reference values
    for i in range(nSys):
        avRefPsi.append(mean(psi_ref.get(i)))
        avRefDel.append(mean(delta_ref.get(i)))


    # extract time series data
    for i in range(nSys):

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

    return AOI, psi_AOI, delta_AOI, t, psi_t, delta_t



def plotData(uniqIDs, info, title, nPoly, nAOI, nSys, AOI, psi_AOI, delta_AOI, t, psi_t, delta_t):

    # conveting delta variables to degrees
    for i in range(nAOI):
        for j in range(len(delta_AOI.get(i))):
            delta_AOI[i][j] = delta_AOI.get(i)[j] * 180 / np.pi

    for i in range(nSys):
        for j in range(len(delta_t.get(i))):
            delta_t[i][j] = delta_t.get(i)[j] * 180 / np.pi


    # Plot Style #
    # fontsize
    fs = 10


    # FIGURE 1 #
    # initialise figure
    fig = plt.figure()
    ax = plt.axes()


    # Psi vs AOI
    for i in range(nAOI):
        ax.plot(AOI.get(i), psi_AOI.get(i), '+', label = info[uniqIDs.get("AOI")[i]][3])

    ax.set_xlabel("AOI (\N{DEGREE SIGN})",fontsize=fs)
    ax.set_ylabel("$\Psi$ (\N{DEGREE SIGN})", color="black", fontsize=fs)

    # legend
    ax.legend()

    # tight layout function
    plt.tight_layout()

    # show plot
    plt.show()

    # save the plot as a file
    fig.savefig( '../output/00/' + title + '_psi_AOI.png',
            format='png',
            dpi=100,
            bbox_inches='tight')


    # FIGURE 2 #
    fig2 = plt.figure()
    ax = plt.axes()

    # Delta vs AoI
    for i in range(nAOI):
        ax.plot(AOI.get(i), delta_AOI.get(i), '+', label = info[uniqIDs.get("AOI")[i]][3])

    ax.set_xlabel("AOI (\N{DEGREE SIGN})",fontsize=fs)
    ax.set_ylabel("$\Delta$ (\N{DEGREE SIGN})", color="black", fontsize=fs)

    # legend
    ax.legend()

    # tight layout function
    plt.tight_layout()

    # show plot
    plt.show()

    # save the plot as a file
    fig2.savefig( '../output/00/' + title + '_delta_AOI.png',
            format='png',
            dpi=100,
            bbox_inches='tight')



    # FIGURE 3 #
    fig3 = plt.figure()
    ax = plt.axes()

    # Psi vs t
    for i in range(nSys):
        ax.plot(t.get(i)[0:len(delta_t.get(i))], psi_t.get(i), '+', label = info[uniqIDs.get("t")[i]][3])

    ax.set_xlabel("Time (s)",fontsize=fs)
    ax.set_ylabel("$\Psi$ (\N{DEGREE SIGN})", color="black", fontsize=fs)

    # legend
    ax.legend()

    # tight layout function
    plt.tight_layout()

    # show plot
    plt.show()

    # save the plot as a file
    fig3.savefig( '../output/00/' + title + '_psi_t.png',
            format='png',
            dpi=100,
            bbox_inches='tight')


    # FIGURE 4 #
    # initialise figure
    fig4 = plt.figure()
    ax = plt.axes()

    # Delta vs t
    for i in range(nSys):
        ax.plot(t.get(i)[0:len(delta_t.get(i))], delta_t.get(i), '+', label = info[uniqIDs.get("t")[i]][3])

    ax.set_xlabel("Time (s)",fontsize=fs)
    ax.set_ylabel("$\Delta_{lipids}$ - $\Delta_{PBS}$ (\N{DEGREE SIGN})", color="black", fontsize=fs)

    # legend
    ax.legend()

    # tight layout function
    plt.tight_layout()

    # show plot
    plt.show()

    # save the plot as a file
    fig4.savefig( '../output/00/' + title + '_delta_t.png',
            format='png',
            dpi=100,
            bbox_inches='tight')


    return



def main(info):

    # filter warnings
    # warnings.filterwarnings("ignore")

    # import title from global variables
    title = config.title

    # user input & sample instructions information
    nPoly, nFiles, uniqIDs, tag, nAOI, nSys, info = importInstructions(info)


    # import data
    AOI, psi_AOI, delta_AOI, t, psi_t, delta_t = importData(uniqIDs, info, nAOI, nSys)


    # pass data for plotting
    plotData(uniqIDs, info, title, nPoly, nAOI, nSys, AOI, psi_AOI, delta_AOI, t, psi_t, delta_t)


    # program executed
    print('\nFigures generated! ~ Have a nice day ~')

    return



if __name__ == '__main__':
    print("~Running plotEllips.py~\n")
    main()
