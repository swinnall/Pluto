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
    fileInfoList  = []
    for i in range(nSamples):
       fileInfoList.append({'fname': instructionsFile["filename"][i], \
       'I_abs_std': instructionsFile["I_abs_std"][i], \
       'I_std': instructionsFile["I_std"][i], \
       'I_solv': instructionsFile["I_solv"][i], \
       'label': instructionsFile["label"][i]})

    return nSamples, fileInfoList


def initDataStruct(nSamples, fileInfoList):

    ## hardcoding paths; this is the general file directory
    inputRoot = '../../UoM-Data-Repository/input/TR-DLS/'

    ## store the nSamples amount of file paths; sample specific file directory
    inputPaths = []
    for sampleNum in range(nSamples):
        inputPaths.append(inputRoot + fileInfoList[sampleNum].get("fname"))

    # initialise variables for time plots
    x        = {new_list: [] for new_list in range(nSamples)}
    y        = {new_list: [] for new_list in range(nSamples)}
    label    = {new_list: 0 for new_list in range(nSamples)}

    return inputPaths, x, y, label


def getTime(timeRow):
    date_time = timeRow[1].strip().strip('""')
    return sum(x * int(t) for x, t in zip([3600, 60, 1], date_time.split(":")))


def analyseCountRate(nSamples, fileInfoList):

    inputPaths, x, y, label = initDataStruct(nSamples, fileInfoList)

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
        label[sampleNum] = fileInfoList[sampleNum].get("label")

    return x, y, label


def closest_value(input_list, input_value): # find the closest value in the model array
    arr = np.asarray(input_list)
    idx = (np.abs(arr - input_value)).argmin()
    return arr[idx], idx


def objectiveMono(x,B,beta,Gamma,mu2): # where x is the tau variable; mu stuff is for polydispersity
    return B + beta*np.exp(-2*Gamma*x)*(1+(mu2/2)*(x**2))**2


def objectiveBi(x,B,beta,Gamma,beta2,Gamma2): # where x is the tau variable
    return B + beta*np.exp(-2*Gamma*x) + beta2*np.exp(-2*Gamma2*x)


def plotFit(popt, tau, g1, g1LimVal, g1LimIdx):

    print('\ng1LimVal = %f; g1LimIdx = %d' %(g1LimVal,g1LimIdx))

    # unpack parameters
    B,beta,Gamma,mu2 = popt

    # plot input vs output
    pyplot.scatter(tau, g1)

    # define a sequence of inputs between the smallest and largest known inputs
    x_line = np.arange(min(tau[0:g1LimIdx]), max(tau[0:g1LimIdx]), 0.0001)

    # calculate the output for the range
    y_line = objectiveMono(x_line,B,beta,Gamma,mu2)
    pyplot.plot(x_line, y_line, '--', color='red')
    pyplot.xscale('log')
    pyplot.show()

    return


def analyseRadius(nSamples, fileInfoList):

    inputPaths, x, y, label = initDataStruct(nSamples, fileInfoList)

    # for every sample specific file directory, list out ASC files and extract data
    for sampleNum in range(nSamples):
        time      = []
        Dlist     = []
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

            initialparameters = (0,0.4,1,0.1)
            g1LimVal, g1LimIdx = closest_value(g1, 0.2*g1[0])

            popt, pcov = curve_fit(objectiveMono, tau[0:g1LimIdx], g1[0:g1LimIdx], p0 = initialparameters)

            # unpack parameters
            B, beta, Gamma, mu2 = popt
            nParMono = len(popt)

            # generate model data for checking if fit is good
            tauModel = np.linspace(start=min(tau[0:g1LimIdx]), stop=max(tau[0:g1LimIdx]), num=g1LimIdx)
            g1Model  = objectiveMono(tauModel,B,beta,Gamma,mu2)

            chisq_red = 0
            for i in range(g1LimIdx):
                chisq_red += (g1Model[i] - g1[i])**2
            chisq_red = chisq_red / (g1LimIdx-nParMono)

            checkFit = False # set to true and choose which fit to analyse
            if checkFit == True and fileCount == 20: # and chisq_red < 0.05:
                print("\nFile Number: %d" %fileCount)
                print("\nChiSq: %.6f" %chisq_red)
                plotFit(popt, tau, g1, g1LimVal, g1LimIdx)


            # relate Gamma (decay time) to diffusion coefficient; account for refractive index?
            lmbda = 632.8 # [nm]
            q = 4*np.pi*np.sin(90/2)/lmbda # [1/nm]
            D = 1E-15 * Gamma / q**2  # (1/ms) / (1/nm)**2 = nm**2 / ms = 1E-18 / 1E-3 [m2/s] = 1E-15 m^2/s
            Dlist.append(D)

            # calculate hydrodynamic radius from diffusion coefficient
            k    = 1.380649e-23 #Boltzmann constant, [J/K]
            T    = 293.15 # [K]
            visc = 1.002 * 1E-3 # [cp] = 10−3 Pa⋅s ; @20 deg; source: Lide, David R. CRC Handbook of Chemistry and Physics. Boca Raton, FL: CRC Press, 2005. http://www.hbcpnetbase.com.
            R = k*T/(6*np.pi*visc*D) * (1/1E-9) # take out nm as a factor to get nice units
            Rlist.append(R)

        # store list of taus in x variables
        x[sampleNum] = time
        y[sampleNum] = Rlist
        label[sampleNum] = fileInfoList[sampleNum].get("label")

    return x, y, label




def analyseRaleighRatio(nSamples, fileInfoList):

    # initialise data structures (main part here is for initialising y)
    inputPaths, x, y, label = initDataStruct(nSamples, fileInfoList)

    # get all count rate data in dict structure; x = time = t, y = meanCR0+meanCR1 = I_solution at each timestep
    t, I_sol, label = analyseCountRate(nSamples, fileInfoList)


    for sampleNum in range(nSamples):

        # intensity from solvent
        I_solv = fileInfoList[sampleNum].get('I_solv')

        # absolute intensity from standard (RR value)
        I_std_abs = fileInfoList[sampleNum].get('I_abs_std')

        # intensity from standard
        I_std = fileInfoList[sampleNum].get('I_std')

        RR = []
        nPoints = len(t[sampleNum])
        for val in range(nPoints):

            RR.append((I_sol.get(sampleNum)[val] - I_solv)*(I_std_abs/I_std))

        # store list of taus in x variables
        y[sampleNum] = RR
        label[sampleNum] = fileInfoList[sampleNum].get("label")

    return t, y, label



def main(instructionsFile, title, inputDIR, outputPath):

    # filter warnings
    warnings.filterwarnings("ignore")

    # give analysis choice to user
    analysisOptions = ['plotCount', 'plotRadius', 'plotRaleigh']

    analysisType, analysisRunning = modSelection(analysisOptions)

    # user input & sample instructions information
    nSamples, fileInfoList = importSampleInfo(instructionsFile)

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
        x, y, label = analyseCountRate(nSamples, fileInfoList)

    if analysisType == 'plotRadius':
        axLabels = {"x": "Time", "y": "R (nm)"}
        suffix   = " - TR DLS radius"
        x, y, label = analyseRadius(nSamples, fileInfoList)

    if analysisType == 'plotRaleigh':
        axLabels = {"x": "Time", "y": "RR"}
        suffix   = " - TR DLS RR"
        x, y, label = analyseRaleighRatio(nSamples, fileInfoList)


    # reduce plot parameters to list of variables for genPlot module
    vars = [len(x), equip, label, axLabels, title, outputPath, [x,0], y]

    return key, vars, suffix



if __name__ == '__main__':
    main()
