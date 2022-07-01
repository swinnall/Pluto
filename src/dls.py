" Process DLS data "

# import packages
import sys, os, csv
import warnings
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot

# import Pluto modules
import config
from genFunc import modSelection, getFile


def importSampleInfo(instructionsFile):

    # number of isotherms to plot
    nSamples  = len(instructionsFile)

    # separate input data into relevant ID lists of dicts
    sampleInfoList  = []
    for i in range(nSamples):
       sampleInfoList.append({'fname': instructionsFile["filename"][i], \
       'I_abs_std': instructionsFile["I_abs_std"][i], \
       'I_std': instructionsFile["I_std"][i], \
       'I_solv': instructionsFile["I_solv"][i], \
       'label': instructionsFile["label"][i]})

    return nSamples, sampleInfoList


def initDataStruct(nSamples, sampleInfoList):

    ## hardcoding paths; this is the general file directory
    inputRoot = '../../UoM-Data-Repository/input/TR-DLS/'

    ## store the nSamples amount of file paths; sample specific file directory
    inputPaths = []
    for sampleNum in range(nSamples):
        inputPaths.append(inputRoot + sampleInfoList[sampleNum].get("fname"))

    # initialise variables for time plots
    x        = {new_list: [] for new_list in range(nSamples)}
    y        = {new_list: [] for new_list in range(nSamples)}
    label    = {new_list: 0 for new_list in range(nSamples)}

    return inputPaths, x, y, label


def getTime(timeRow):
    date_time = timeRow[1].strip().strip('""')
    return sum(x * int(t) for x, t in zip([3600, 60, 1], date_time.split(":")))


def analyseCountRate(nSamples, sampleInfoList):

    inputPaths, x, y, label = initDataStruct(nSamples, sampleInfoList)

    # for every sample specific file directory, list out ASC files and extract data
    for sampleNum in range(nSamples):

        time   = []
        meanCR = []
        fileCount = 0
        for filename in os.listdir(inputPaths[sampleNum]):

            if filename.endswith(".ASC"):
                filepath = inputPaths[sampleNum] + '/' + filename

                with open(filepath) as f:
                    fileData = f.readlines()[0:30]

                    timeRow = fileData[2].strip("\n").split("\t")
                    sec_time = getTime(timeRow)

                    if fileCount == 0:
                        startTime = sec_time
                    time.append(sec_time - startTime)

                    meanCR0 = float(fileData[24].strip("\n").split("\t")[1])
                    meanCR1 = float(fileData[25].strip("\n").split("\t")[1])
                    meanCR.append(meanCR0+meanCR1)

                f.close()
                fileCount += 1

        x[sampleNum] = time
        y[sampleNum] = meanCR
        label[sampleNum] = sampleInfoList[sampleNum].get("label")

    return x, y, label


def closest_value(input_list, input_value): # find the closest value in the model array
    arr = np.asarray(input_list)
    idx = (np.abs(arr - input_value)).argmin()
    return arr[idx], idx


def objectiveMono(x,*popt):
    B, beta, Gamma, mu2 = popt # where x is the tau variable; mu stuff is for polydispersity
    return B + beta*np.exp(-2*Gamma*x)*(1+(mu2/2)*(x**2))**2

def objectiveBi(x,*popt):
    B, beta, Gamma, beta2, Gamma2 = popt
    return B + beta*np.exp(-2*Gamma*x) + beta2*np.exp(-2*Gamma2*x)

def objectiveTri(x,*popt):
    B, beta, Gamma, beta2, Gamma2, beta3, Gamma3 = popt
    return B + beta*np.exp(-2*Gamma*x) + beta2*np.exp(-2*Gamma2*x) + beta3*np.exp(-2*Gamma3*x)


def plotFit(tau, g1, g1LimVal, g1LimIdx, *popt):

    print('\ng1LimVal = %f; g1LimIdx = %d' %(g1LimVal,g1LimIdx))

    # plot input vs output
    pyplot.scatter(tau, g1)

    # define a sequence of inputs between the smallest and largest known inputs
    x_line = np.linspace(start=min(tau[0:g1LimIdx]), stop=max(tau[0:g1LimIdx]), num=g1LimIdx*10)

    # calculate the output for the range
    y_line = objectiveBi(x_line,*popt)
    pyplot.plot(x_line, y_line, '--', color='red')
    pyplot.xscale('log')
    pyplot.show()

    return


def analyseRadius(nSamples, sampleInfoList):

    inputPaths, x, y, label = initDataStruct(nSamples, sampleInfoList)

    # for every sample specific file directory, list out ASC files and extract data
    for sampleNum in range(nSamples):
        time      = []
        Rlist     = []
        fileCount = 0
        chiCount = 0
        for filename in os.listdir(inputPaths[sampleNum]):
            if filename.endswith(".ASC"):
                filepath = inputPaths[sampleNum] + '/' + filename

                # initialise lists to store values per file
                tau = []
                g1  = []
                with open(filepath) as f:
                    fileData = f.readlines()[0:205]

                    timeRow = fileData[2].strip("\n").split("\t")
                    sec_time = getTime(timeRow)

                    if fileCount == 0:
                        startTime = sec_time
                    time.append(sec_time - startTime)

                    # isolate region for correlation data only; could make more general
                    correlationData = fileData[30:197]

                    # split all the data corresponding to
                    for row in correlationData:
                        rowData = row.split("\t")

                        try:
                            tauVal, g1Val = float(rowData[0]), float(rowData[1])

                            # store values in lists
                            tau.append(tauVal)

                            if g1Val < 0:
                                g1.append(-np.sqrt(abs(g1Val)))
                            else:
                                g1.append(np.sqrt(g1Val))

                        except ValueError:
                            print("\nValueError: Could not convert string to float.\nString = %s\n" %rowData)

                f.close()
            fileCount += 1


            # function that calculates the reduced chiSq
            def calcChiSq(yModel, yExp, nPoints, nPar):
                chisq = 0
                for i in range(nPoints):
                    chisq += (yExp[i] - yModel[i])**2
                return chisq / (nPoints-nPar)


            # region extraction
            g1LimVal, g1LimIdx = closest_value(g1, 0.2*g1[0])


        #if fileCount < 40: # fitType.upper() == 'MONO' and
            initialparameters  = (0,0.4,1,0.1)
            popt, pcov = curve_fit(objectiveMono, tau[0:g1LimIdx], g1[0:g1LimIdx], p0=initialparameters)
            B, beta, Gamma, mu2 = popt

            tauModel = np.linspace(start=min(tau[0:g1LimIdx]), stop=max(tau[0:g1LimIdx]), num=g1LimIdx)
            g1Model  = objectiveMono(tauModel,*popt)
            chisq_red1 = calcChiSq(g1Model, g1, g1LimIdx, len(popt))


        #if fileCount >= 40: # fitType.upper() == 'BI' and
            initialparameters  = (0,0.4,1,0.1,3)
            popt, pcov = curve_fit(objectiveBi, tau[0:g1LimIdx], g1[0:g1LimIdx], p0=initialparameters)
            B,beta,Gamma,beta2,Gamma2 = popt

            tauModel = np.linspace(start=min(tau[0:g1LimIdx]), stop=max(tau[0:g1LimIdx]), num=g1LimIdx)
            g1Model  = objectiveBi(tauModel,*popt)
            chisq_red2 = calcChiSq(g1Model, g1, g1LimIdx, len(popt))

        #if fitType.upper() == 'TRI':
            # initialparameters  = (0,0.4,2,0.1,3,0.1,3)
            # popt, pcov = curve_fit(objectiveTri, tau[0:g1LimIdx], g1[0:g1LimIdx], p0=initialparameters)
            # B, beta, Gamma, beta2, Gamma2, beta3, Gamma3 = popt
            #
            # tauModel = np.linspace(start=min(tau[0:g1LimIdx]), stop=max(tau[0:g1LimIdx]), num=g1LimIdx)
            # g1Model  = objectiveTri(tauModel,*popt)
            chisq_red3 = 1000# calcChiSq(g1Model, g1, g1LimIdx, len(popt))

            def calcRh(Gamma):

                # calculate diffusion coefficient
                lmbda = 632.8 # [nm]
                q = 4*np.pi*np.sin(90/2)/lmbda # [1/nm]
                D = 1E-15 * Gamma / q**2  # (1/ms) / (1/nm)**2 = nm**2 / ms = 1E-18 / 1E-3 [m2/s] = 1E-15 m^2/s

                # calculate hydrodynamic radius from diffusion coefficient
                k    = 1.380649e-23 #Boltzmann constant, [J/K]
                T    = 293.15 # [K]
                visc = 1.002 * 1E-3 # [cp] = 10−3 Pa⋅s ; @20 deg; source: Lide, David R. CRC Handbook of Chemistry and Physics. Boca Raton, FL: CRC Press, 2005. http://www.hbcpnetbase.com.

                return k*T/(6*np.pi*visc*D) * (1/1E-9) # take out nm as a factor to get nice units


            if fileCount < 40: # chisq_red1 < chisq_red2 and chisq_red1 < chisq_red3: #
                #print("Fit = Mono")
                initialparameters  = (0,0.4,1,0.1)
                popt, pcov = curve_fit(objectiveMono, tau[0:g1LimIdx], g1[0:g1LimIdx], p0=initialparameters)
                B, beta, Gamma, mu2 = popt
                #print(popt)

                R = calcRh(Gamma)

            elif fileCount >= 40:  # chisq_red2 < chisq_red1 and chisq_red2 < chisq_red3: #
                initialparameters  = (0,0.4,1,0.1)
                popt, pcov = curve_fit(objectiveMono, tau[0:g1LimIdx], g1[0:g1LimIdx], p0=initialparameters)
                B, beta, gGamma, mu2 = popt

                #print("Fit = Bi")
                initialparameters  = (0,0.4,1,0.2,gGamma)
                popt, pcov = curve_fit(objectiveBi, tau[0:g1LimIdx], g1[0:g1LimIdx], p0=initialparameters, bounds=([0,0,0,0,(gGamma-0.01)],[0.1,1,5,1,(gGamma+0.01)]))
                B,beta,Gamma,beta2,Gamma2 = popt
                # print(popt)

                R1 = calcRh(Gamma) # big species
                R2 = calcRh(Gamma2) # small speices (does depend on how you setup the fits...)
                R = R1

            elif chisq_red3 < chisq_red1 and chisq_red3 < chisq_red2:
                print("Fit = Tri")
                sys.exit()

            ##### try fitting triExp for everything
            ##### take gamma from monoExp and use just use in everything
            ##### maybe try contin to try and get another idea of size
            #####

            Rlist.append(R)



            checkFit = False # set to true and choose which fit to analyse
            if checkFit == True and fileCount == 80:
                print("\nFile Number: %d" %fileCount)
                print(pcov)
                perr = np.sqrt(np.diag(pcov))
                print(perr)
                # generate model data for checking if fit is good
                tauModel = np.linspace(start=min(tau[0:g1LimIdx]), stop=max(tau[0:g1LimIdx]), num=g1LimIdx)
                g1Model  = objectiveBi(tauModel,*popt)


                chisq_red = calcChiSq(g1Model, g1, g1LimIdx, len(popt))
                print("\nChiSq: %.6f" %chisq_red)

                plotFit(tau, g1, g1LimVal, g1LimIdx, *popt)


        # store list of taus in x variables
        x[sampleNum] = time
        y[sampleNum] = Rlist
        label[sampleNum] = sampleInfoList[sampleNum].get("label")
    #print(y)
    return x, y, label




def analyseRaleighRatio(nSamples, sampleInfoList):

    # initialise data structures (main part here is for initialising y)
    inputPaths, x, y, label = initDataStruct(nSamples, sampleInfoList)

    # get all count rate data in dict structure; x = time = t, y = meanCR0+meanCR1 = I_solution at each timestep
    t, I_sol, label = analyseCountRate(nSamples, sampleInfoList)


    for sampleNum in range(nSamples):

        # intensity from solvent
        I_solv = sampleInfoList[sampleNum].get('I_solv')

        # absolute intensity from standard (RR value)
        I_std_abs = sampleInfoList[sampleNum].get('I_abs_std')

        # intensity from standard
        I_std = sampleInfoList[sampleNum].get('I_std')

        RR = []
        nPoints = len(t[sampleNum])
        for val in range(nPoints):

            RR.append((I_sol.get(sampleNum)[val] - I_solv)*(I_std_abs/I_std))

        # store list of taus in x variables
        y[sampleNum] = RR
        label[sampleNum] = sampleInfoList[sampleNum].get("label")

    return t, y, label



def main(instructionsFile, title, inputDIR, outputPath):

    # filter warnings
    warnings.filterwarnings("ignore")

    # give analysis choice to user
    analysisOptions = ['plotCount', 'plotRadius', 'plotRaleigh']

    analysisType, analysisRunning = modSelection(analysisOptions)

    # user input & sample instructions information
    nSamples, sampleInfoList = importSampleInfo(instructionsFile)

    if nSamples < 1:
        print("Fatal Error: No file input.")
        sys.exit()

    # plot initialisations
    equip = 'null'
    key   = (1,1)

    # analysis functions with associated plotting params
    if analysisType == 'plotCount':
        axLabels = {"x": "T", "y": "Mean Count Rate"}
        suffix   = " - TR DLS countRate"
        x, y, label = analyseCountRate(nSamples, sampleInfoList)

    if analysisType == 'plotRadius':
        axLabels = {"x": "Time", "y": "R (nm)"}
        suffix   = " - TR DLS radius"
        x, y, label = analyseRadius(nSamples, sampleInfoList)

    if analysisType == 'plotRaleigh':
        axLabels = {"x": "Time", "y": "RR"}
        suffix   = " - TR DLS RR"
        x, y, label = analyseRaleighRatio(nSamples, sampleInfoList)


    # reduce plot parameters to list of variables for genPlot module
    vars = [len(x), equip, label, axLabels, title, outputPath, [x,0], y]

    return key, vars, suffix



if __name__ == '__main__':
    main()
