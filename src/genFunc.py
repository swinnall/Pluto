" General functions for Pluto "

# import packages
import csv, sys
import pandas as pd
import numpy as np
from numpy import linalg
from scipy.signal import savgol_filter

# import Pluto modules
import config


def getFile(path,nSkip,delim):
    return pd.read_csv(path, skiprows=nSkip, sep=delim, comment='#', na_values =' ', skip_blank_lines=True, encoding = "utf-8") # on_bad_lines='skip',


def modSelection(analysisOptions):

    # ask user to pick one of the analysisOptions
    print("\n~~~\nAnalysis Options:")
    for i,option in enumerate(analysisOptions):
        print("%d: %s" %(i+1,option))
    print("~~~\n")

    analysisChoice = input("Which analysis would you like to do? Pick the associated number (1-%d):\n  " %len(analysisOptions) )

    if analysisChoice.upper() == 'Q':
        print("Session closed.")
        sys.exit()

    elif analysisChoice in [str(i) for i in range(len(analysisOptions)+1)]:
        analysisType = analysisOptions[int(analysisChoice)-1]
        print("You picked %s.py\n" %analysisType)
        analysisRunning = True

    elif analysisChoice not in [str(i) for i in range(len(analysisOptions)+1)] and analysisOptions[0]!="surfacePressure":
        print("Not a valid response. Returning to Pluto landing page.\n\n")
        analysisType    = 'n/a'
        analysisRunning = False

    else:
        print("Not a valid response. Session closed.")
        sys.exit()

    return analysisType, analysisRunning


def reducePoints(x, y):

    nLists = len(x)

    redX = {init: [] for init in range(nLists)}
    redY = {init: [] for init in range(nLists)}
    for i in range(nLists):
        nPoint = len(x.get(i))
        for j in range(0,nPoint,config.nth):
            redX[i].append(x.get(i)[j])
            redY[i].append(y.get(i)[j])

    return redX, redY


def polySmoothData(y):

    # gets number of lists within the y dictionary
    nLists = len(y)

    # initialise list for func output
    smoothList = {init: [] for init in range(nLists)}

    for item in y.items():

        # unpacks tuple item
        key     = item[0]
        genList = item[1]

        # if length is even, truncate as savgol must have odd window length
        if len(genList) % 2 == 0:
            del genList[len(genList)-1]

        # remove noise introduced via derivative perturbations
        try:
            smoothList = savgol_filter(genList, len(genList), config.nPoly, mode='interp')  # variable, window size (length), polynomial order, mode

        # remove infinity and nan values due to Nima having repeated area values
        except linalg.LinAlgError:
            print("LinAlgError")
            for j in range(len(genList)):
                if genList[j] == float("inf") or genList[j] == float("-inf") or mt.isnan(genList[j]) == True:
                    genList[j] = 0

        except ValueError:
            print("ValueError: Polyorder must be less than window length. Attempting adjustment at 0.5*nPoly = %d." %int(config.nPoly/2))
            smoothList = savgol_filter(genList, len(genList), int(config.nPoly/2), mode='interp')  # variable, window size (length), polynomial order, mode


        # high pass filter: remove values less than 0.01 (final smoothing of function)
        deleteMe = []
        for i in range(0,len(smoothList)-1):
            if smoothList[i] < config.lowerLimit:
                deleteMe.append(i)

        # delete elements from the lists
        smoothList = np.delete(smoothList, deleteMe)

        # updates the key-value pair with the smoothed data
        y[key] = smoothList

    return y
