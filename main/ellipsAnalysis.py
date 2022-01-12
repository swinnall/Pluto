"Reads ellipsometry data"

import sys
import csv
import re
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from statistics import mean
import warnings
import cmath as cm
import scipy.signal as sig
import scipy.optimize as opt
import config
import genPlot


def modSelection(analysisOptions):

    # ask user to pick one of the analysisOptions
    print("~~~\nElasticity analysis Options:\n %s" %analysisOptions)

    analysisChoice = input("Which analysis would you like to do? Pick the associated number (0-%d) or 'q' to reutrn to landing page:\n  " %(len(analysisOptions)-1) )

    if analysisChoice == 'q':
        print("Returning to Pluto landing page.\n\n")
        analysisType    = 'n/a'
        analysisRunning = False

    elif analysisChoice == 'Q':
        print("Session closed.")
        sys.exit()

    elif analysisChoice in [str(i) for i in range(len(analysisOptions))]:
        analysisType = analysisOptions[int(analysisChoice)]
        print("You picked %s.py\n" %analysisType)
        analysisRunning = True

    else:
        print("Not a valid response. Returning to Pluto landing page.\n\n")
        analysisType    = 'n/a'
        analysisRunning = False

    return analysisType, analysisRunning


def importSampleInfo(info):

    # number of header rows
    nHeaders = 2

    # number of isotherms to plot
    nFiles  = len(info) - nHeaders

    # create ID and count structures for later parsing
    AOI_ID  = []; nAOI  = 0
    time_ID = []; nTime = 0

    # separate input data into relevant ID lists of dicts
    for i in range(nHeaders,nHeaders+nFiles):

        if info[i][1].upper() == "NULL":
            AOI_ID.append({'fname': info[i][0], 'label': info[i][2]})
            nAOI += 1

        else:
            time_ID.append({'fname': info[i][0], 'refname': info[i][1], 'label': info[i][2]})
            nTime += 1

    return AOI_ID, nAOI, time_ID, nTime



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
        fileDIR = inputDIR + '/' + time_ID[i].get('refname') +".txt"
        ref_df  = getFile(fileDIR)

        # extract raw ref data from dataframe
        for j in range(0,len(ref_df)):
            psi_ref  .append(float(ref_df[j][0]))
            delta_ref.append(float(ref_df[j][1]))

        # average reference values; each list ele is an averaged ref
        avRefPsi = mean(psi_ref)
        avRefDel = mean(delta_ref)


        ## Extract time series data
        fileDIR = inputDIR + '/' + time_ID[i].get('fname') +".txt"
        time_df = getFile(fileDIR)

        # subtract ref data from these to give lipid-only data
        for j in range(0,len(time_df)):
            t[i]      .append(float(time_df[j][6]))
            psi_t[i]  .append(float(time_df[j][0]) )#- avRefPsi) ## for now subtract psi ref
            delta_t[i].append(float(time_df[j][1]) )#- avRefDel) ## comment for bare interface

        label_t[i] = time_ID[i].get('label')

    return t, psi_t, delta_t, label_t



# read in experiment data for AOI sets, Psi and Delta
# looking to find parameters; par = [n1, d1, n2, d2]
# intial conditions: x0 = [n1_0, d1_0, n2_0 (water), d2_0 (water)] = [1, 1, 1.33 0];
# bounds = [(0,2), (0,2), (1.3,1.4), (-0.1,0.1)]
# test = opt.differential_evolution(residuals, bounds, init='sobol', maxiter=10000)
# solution = test.x
# print(solution)
# print(test.nit)
# lstsq = residuals(solution); to check how good solution was
# out = np.append(solution, lstsq)
# fit = np.vstack((fit,out))


def ellipsModel(data, par):

    # unpack experimental data
    AOI       = data[0]
    Delta_exp = data[1]
    Psi_exp   = data[2]


    # convert AOI to radians
    theta = []
    for i in range(len(AOI)):
        theta = AOI * np.pi / 180


    # wavelength of laser
    lmbda = 632.8E-9


    # refractive index of air
    N0 = 1


    # film 1
    k1 = 0
    n1 = par[0]
    d1 = par[1]

    Nt1 = n1 - cm.j*k1


    # film 2
    k2 = 0
    n2 = par[2]
    d2 = par[3]

    Nt1 = n2 - cm.j*k2


    # substrate
    k3 = 0
    n3 = 1.337

    Nt3 = n3 - cm.j*k3


    # unknown
    ct2 = -12.5663706143592*cm.j


    Rho = []
    for i in range(len(theta)):

        C0    = np.cos(theta[i])
        S0    = np.sin(theta[i])

        C1    = np.sqrt( 1 - (np.sin(theta[i])*(N0/Nt1))**2  )
        rp10  = (Nt1*C0 - N0*C1) / (Nt1*C0 + N0*C1)
        rn10  = (N0*C0 - Nt1*C1) / (N0*C0 + Nt1*C1)

        C2    = np.sqrt( 1 - (np.sin(theta[i])*(N0/Nt2))**2  )
        rp21  = (Nt2*C1 - Nt1*C2) / (Nt2*C1 + Nt1*C2)
        rn21  = (Nt1*C1 - Nt2*C2) / (Nt1*C1 + Nt2*C2)

        T2    = np.exp(ct2*Nt2*C2*d2/lmbda)

        C3    = np.sqrt( 1 - (np.sin(theta[i])*(N0/Nt3))**2  )
        rp32  = (Nt3*C2 - Nt2*C3) / (Nt3*C2 + Nt2*C3)
        rn32  = (Nt2*C2 - Nt3*C3) / (Nt2*C2 + Nt3*C3)

        crp21 = (rp21 + rp32*T2) / (1 + rp21*rp32*T2)
        crn21 = (rn21 + rn32*T2) / (1 + rn21*rn32*T2)

        T1    = np.exp(ct2*Nt1*C1*d1/lambda);
        crp10 = (rp10 + crp21*T1) / (1 + rp10*crp21*T1)
        crn10 = (rn10 + crn21*T1) / (1 + rn10*crn21*T1)

        Rho.append(crp10/crn10)


    Psi = []
    for i in range(len(Rho)):
        Psi.append( np.arctan(abs(Rho[i])) )


    Delta = []
    for i in range(len(Rho)):
        Delta.append( np.angle(Rho[i]) * 180/np.pi )

        if Delta[i] < 0:
            Delta[i] == Delta[i] + 360



    return Rho, Psi, Delta




# least square condition
def leastsquare(y, test):
    return np.sum((y-test)**2)


# residual function for genetic algorithm
def residuals(para):

    # where ellipsModel is the ellipsometry theory function
    test = ellipsModel(x0, *para)

    while np.size(y0) > np.size(test):
        test = np.append(test,120)

        if np.size(y0) == np.size(test):
            break

    return leastsquare(y0,test)












def main(info, title, inputDIR, outputPath):

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
        AOI_ID, nAOI, time_ID, nTime = importSampleInfo(info)


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
            for i in range(nTime):
                for j in range(len(delta_t.get(i))):
                    delta_t[i][j] = delta_t.get(i)[j] * 180 / np.pi


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

    return



if __name__ == '__main__':
    print("~Running plotEllips.py~\n")
    main()
