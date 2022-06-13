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
       fileInfoList.append({'fname': instructionsFile["filename"][i], 'standard': instructionsFile["standard"][i], 'solvent': instructionsFile["solvent"][i], 'label': instructionsFile["label"][i]})

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



def analyseCorrelation(nSamples, fileInfoList):

    inputPaths, x, y, label = initDataStruct(nSamples, fileInfoList)

    # for every sample specific file directory, list out ASC files and extract data
    for sampleNum in range(nSamples):
        tauList = []
        g1List  = []

        for filename in os.listdir(inputPaths[sampleNum]):
            if filename.endswith(".ASC"):
                filepath = inputPaths[sampleNum] + '/' + filename

                # initialise lists to store values per file
                tau = []
                g1  = []
                with open(filepath) as f:

                    # isolate region for correlation data only
                    # need to identify general method for determining region
                    correlationData = f.readlines()[30:205]

                    # split all the data corresponding to
                    for row in correlationData:
                        rowData = row.split("\t")

                        try:
                            tauVal, g1Val = float(rowData[0]), float(rowData[1])

                            # store values in lists
                            tau.append(tauVal) # np.log(tauVal)

                            if g1Val < 0:
                                g1.append(-np.sqrt(abs(g1Val)))
                            else:
                                g1.append(np.sqrt(g1Val))

                        except ValueError:
                            print("\nValueError: Could not convert string to float.\nString = %s\n" %rowData)

                f.close()

            # store tau in list of taus
            tauList.append(tau)
            g1List.append(g1)

            # where x is the tau variable
            def objective(x,B,beta,Gamma,mu2):
                return B + beta*np.exp(-2*Gamma*x)*(1+(mu2/2)*(x**2))**2


            initialparameters = (0,0.4,1,0.1)
            tauLim = 100
            popt, pcov = curve_fit(objective, tau[0:tauLim], g1[0:tauLim], p0 = initialparameters)
            print(popt)
            print(pcov)
            B,beta,Gamma,mu2 = popt
            print('y = %.5f + %.5f*np.exp(-2*%.5f*x)*(1+(%.5f/2)*(x**2))**2' %(B,beta,Gamma,mu2))
            # plot input vs output
            pyplot.scatter(tau, g1)
            # define a sequence of inputs between the smallest and largest known inputs
            x_line = np.arange(min(tau[0:tauLim]), max(tau[0:tauLim]), 0.0001)

            # calculate the output for the range
            y_line = objective(x_line,B,beta,Gamma,mu2)
            pyplot.plot(x_line, y_line, '--', color='red')
            pyplot.xscale('log')
            pyplot.show()
            sys.exit()

            ## relate Gamma (decay time) to diffusion coefficient
            ## relate diffusion coefficient to hydrodynamic radius
            ## store these values as a function of time (using date_time) and pass to x, y for plotting
            
            ## viscosity of the buffer (use water to start with at 20 deg = 1.002 cp) at 30 deg it's 0.7977 cp
            ## source: Lide, David R. CRC Handbook of Chemistry and Physics. Boca Raton, FL: CRC Press, 2005. http://www.hbcpnetbase.com.
            ## Add 30 seconds of dead time to the beginning

        # store list of taus in x variables
        x[sampleNum] = tauList[0]
        y[sampleNum] = g1List[0]
        label[sampleNum] = fileInfoList[sampleNum].get("label")

    return x, y, label




def analyseRaleighRatio(nSamples, fileInfoList):

    print("\n\nError: Function not written yet.\n\n")
    sys.exit()

    ## for each ASC file, plot col2 (count?) against col1 (time)
    ## with some calculation to get the absolute intensity

    return x, y, label



def main(instructionsFile, title, inputDIR, outputPath):

    # filter warnings
    warnings.filterwarnings("ignore")

    # give analysis choice to user
    analysisOptions = ['plotCount', 'plotCorrelation', 'plotRaleigh']

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
        suffix   = " - TR DLS"
        x, y, label = analyseCountRate(nSamples, fileInfoList)

    if analysisType == 'plotCorrelation':
        axLabels = {"x": "Tau", "y": "g1"}
        suffix   = " - TR DLS Correlation"
        x, y, label = analyseCorrelation(nSamples, fileInfoList)

    if analysisType == 'plotRaleigh':
        x, y, label = analyseRaleighRatio(nSamples, fileInfoList)


    # reduce plot parameters to list of variables for genPlot module
    vars = [len(x), equip, label, axLabels, title, outputPath, [x,0], y]

    return key, vars, suffix



if __name__ == '__main__':
    main()
