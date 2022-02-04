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
import math
import scipy.signal as sig
import scipy.optimize as opt
import config
import genPlot
from genFunc import modSelection, getEllipsometryFile


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
        fileDIR = inputDIR + '/' + time_ID[i].get('refname') +".txt"
        ref_df  = getEllipsometryFile(fileDIR)

        # extract raw ref data from dataframe
        for j in range(0,len(ref_df)):
            psi_ref  .append(float(ref_df[j][0]))
            delta_ref.append(float(ref_df[j][1]))

        # average reference values; each list ele is an averaged ref
        avRefPsi = mean(psi_ref)
        avRefDel = mean(delta_ref)


        ## Extract time series data
        fileDIR = inputDIR + '/' + time_ID[i].get('fname') +".txt"
        time_df = getEllipsometryFile(fileDIR)

        # subtract ref data from these to give lipid-only data
        for j in range(0,len(time_df)):
            t[i]      .append(float(time_df[j][6]))
            psi_t[i]  .append(float(time_df[j][0]) )#- avRefPsi) ## for now subtract psi ref
            delta_t[i].append(float(time_df[j][1]) )#- avRefDel) ## comment for bare interface

        label_t[i] = time_ID[i].get('label')

    return t, psi_t, delta_t, label_t




def sinModel(x,par):

    # initial guess parameters; amplitude, wavenumber, phase
    A      = par[0]
    lmbda1 = par[1]
    phi1   = par[2]
    B      = par[3]
    lmbda2 = par[4]
    phi2   = par[5]
    c      = par[6]

    # basic model, might add second sine function and/or time dependence
    sinY_Model = A*np.sin(2*np.pi*np.array(x)/lmbda1 + phi1) + B*np.sin(2*np.pi*np.array(x)/lmbda2 + phi2) + c

    return sinY_Model


# least square condition
def leastsquare(delta_exp, sinY_Model):

    diffSq = []
    for i, j in zip(delta_exp,sinY_Model):
        diffSq.append( np.power(i-j,2) )

    return np.sum(diffSq)


# residual function for genetic algorithm
def residuals(par, t, delta_exp):

    # where sinModel is the sinusoidal function
    sinY_Model = sinModel(t, par)

    return leastsquare(delta_exp, sinY_Model)


def fitSine(t, delta_t):

    # update variable names for easier processing
    t         = t.get(0)
    Delta_exp = delta_t.get(0)


    # high pass filter: remove values less than 0.01 (final smoothing of function)
    tolerance  = 0.01
    smoothT = []; smoothList = []
    for i in range(0,len(Delta_exp)-1):

        if abs(Delta_exp[i] - Delta_exp[i+1]) < tolerance and Delta_exp[i] < 180.60 and Delta_exp[i] > 179.81:
            smoothList.append(Delta_exp[i])
            smoothT.append(t[i])

    # define boundaries of k; 2pi/lambda
    lmbda1_lower = 5000
    lmbda1_upper = 15000
    lmbda1_guess = 6.429e+03
    lmbda2_lower = 60000
    lmbda2_upper = 80000
    lmbda2_guess = 7.98488484e+04

    # define input parameters; [A, lmbda1, phi1, B, lmbda2, phi2, c]
    #par = [0.4, lmbda1_guess, 0, 180]
    par = [0.4, lmbda1_guess, 0, 0.1, lmbda2_guess, 0, 180]

    # associated parameter bounds; could fix pars by defining in model function
    #bounds = [(0.2,0.5), (lmbda1_lower,lmbda1_upper), (-np.pi,np.pi), (179.5,180.5)]
    bounds = [(0.2,0.5), (lmbda1_lower,lmbda1_upper), (-np.pi,np.pi), (0.0,0.5), (lmbda2_lower,lmbda2_upper), (-np.pi,np.pi), (178,182)]

    # genetic algorithm
    geneticOutput = opt.differential_evolution(residuals, bounds, args=(smoothT, smoothList), maxiter=1000) # might need args=*par or make par global

    # parameter solution
    solution = geneticOutput.x
    print("\nParameter solution (A, lmbda1, phi1, B, lmbda2, phi2, c):\n %s" %solution)

    # associated cost
    lstsq = residuals(solution, t, Delta_exp)
    print("\nCost of chosen solution: %f" %lstsq)

    # number of iterations
    print("\nNumber of iterations: %s" %geneticOutput.nit)

    # bool of success
    print("\nOptimisation status: %s" %geneticOutput.success)

    # termination message
    print("\nTermination message: %s" %geneticOutput.message)


    # generate second sine set
    sinY_Model = sinModel(t,solution)

    plt.plot(t, Delta_exp, '.', label='Raw Data', c='#1643A2')
    plt.plot(smoothT, smoothList, '-', label='Smoothed Data', c='#3CA8AB')
    plt.plot(t, sinY_Model, label='Fit', linewidth=3, c='#E9BC5D')
    plt.legend()
    plt.xlabel('Time (s)')
    plt.ylabel('Delta (deg)')
    plt.show()


    return




def main(instructionsFile, title, inputDIR, outputPath):

    ## this is to calibrate the periodicity subtraction of the data
    ## Create simpleFit option that analyses the periodicity of sample data - could extend to comparison between multiple fits
    ## this should then generate parameters
    ## generate sinusoidal data (general function) with these parameters
    ## plot original data, fit, and subtracted data
    ## ask user if they would like to update default parameters with the new ones (in practice it's most like the program's default will be set and rarely updated as the same should be used continuously)
    ## ask user if they want to subtract the sin curve from this data

    ### Points to think about
    ## maybe each time the experimental data is read it should be analysed for periodicity [have to select region]
    ## would have to have a noise filter
    ## could callibrate against Andreas' data if he still has it


    # filter warnings
    warnings.filterwarnings("ignore")


    # give analysis choice to user
    analysisOptions = ['plotAOIpsi','plotAOIdelta','plotTimePsi','plotTimeDelta', 'parOpt','fitSin']

    analysisRunning = True
    while analysisRunning:
        analysisType, analysisRunning = modSelection(analysisOptions)

        if analysisRunning == False:
            break


        # user input & sample instructions information
        AOI_ID, nAOI, time_ID, nTime = importSampleInfo(instructionsFile)


        # process data
        if analysisType == 'plotAOIpsi' or analysisType == 'plotAOIdelta' or analysisType == 'parOpt':

            if nAOI < 1:
                print("Error: No AOI files in input.")
                sys.exit()

            # import data
            AOI, psi_AOI, delta_AOI, label_AOI = importAOIData(inputDIR, AOI_ID, nAOI)

            # converting delta variables to degrees
            for i in range(nAOI):
                for j in range(len(delta_AOI.get(i))):
                    delta_AOI[i][j] = delta_AOI.get(i)[j] * 180 / np.pi


            #if analysisType == 'parOpt':
            #    parameterOptimisation(AOI, psi_AOI, delta_AOI)




        if analysisType == 'plotTimePsi' or analysisType == 'plotTimeDelta' or analysisType == 'fitSin':

            if nTime < 1:
                print("Error: No Time files in input.")
                sys.exit()

            # import data
            t, psi_t, delta_t, label_t = importTimeData(inputDIR, time_ID, nTime)

            # converting delta variables to degrees
            for i in range(nTime):
                for j in range(len(delta_t.get(i))):
                    delta_t[i][j] = delta_t.get(i)[j] * 180 / np.pi

            if analysisType == 'fitSin':
                fitSine(t, delta_t)


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
