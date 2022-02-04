" Pluto "
" Author: @S.Winnall "

import glob, os, sys
import pandas as pd
import csv
import shutil
from shutil import copyfile
import config
import isoAnalysis
import ellipsAnalysis
import chemFormulations
import sldAnalysis
import SurfaceExcess


def modSelection(analysisOptions):

    # root input directory
    inputDir  = config.inputDir

    # root output directory
    outputDir = config.outputDir


    # ask user to pick one of the analysisOptions
    print("\n~~~\nAnalysis Options:")
    for i,option in enumerate(analysisOptions):
        print("%d: %s" %(i+1,option))
    print("~~~\n")

    analysisChoice = input("Which analysis would you like to do? Pick the associated number (1-%d) or 'q' to exit:\n  " %len(analysisOptions))

    if analysisChoice == 'q':
        print("Session closed.")
        sys.exit()

    elif analysisChoice in [str(i) for i in range(len(analysisOptions)+1)]:
        analysisType = analysisOptions[int(analysisChoice)-1]
        print("You picked %s.py\n" %analysisType)

    else:
        print("Not a valid response. Session closed.")
        sys.exit()


    # calls from config database
    instructionsName = config.pathNames.get(analysisType)[0]


    # input/output folder names to be accessed within root in/output directory
    inputOutputStem  = config.pathNames.get(analysisType)[1]


    # instructions file path
    instructionsPath = "" + inputDir + instructionsName + ".txt"


    # input data files path
    inputDataPath    = '' + inputDir + '/' + inputOutputStem + ''


    # output data files path
    outputDataPath   = '' + outputDir + '/' + inputOutputStem + ''


    return analysisType, instructionsName, instructionsPath, inputDataPath, outputDataPath



def organisePaths(analysisType, instructionsName, instructionsPath, outputDataPath):

    # get analysis title
    with open(instructionsPath, newline = '') as f:
        title = list(csv.reader(f))[0][0].split('=')[1]


    # read instructions file with pandas
    instructionsFile = pd.read_csv(instructionsPath, header=1, comment='#')


    # update output path folder dir to include title of chosen analysis
    # e.g. output/chemFormulations/analysisTitle/
    outputDataPath = '' + outputDataPath + '/' + title + ''


    try:
        # delete folder if exists and create it again
        shutil.rmtree(outputDataPath)
        os.mkdir(outputDataPath)
    except FileNotFoundError:
        os.mkdir(outputDataPath)


    try:
        # copy input instructions instructionsFilermation to outputDataPath directory
        shutil.copyfile(instructionsPath, outputDataPath + '/' + instructionsName + '.txt')

        # rename copied instructions file to include analysis title
        # this is the file which will be appended in the relevant analysis
        old_name              = outputDataPath + '/' + instructionsName + '.txt'
        outputInstructionFile = outputDataPath + '/' + title + ' - ' + instructionsName + '.txt'
        os.rename(old_name, outputInstructionFile)


    # if source and destination are same
    except shutil.SameFileError:
        print("Instructions Copying Error:\n  Source and destination represents the same file.")

    # if destination is a directory
    except IsADirectoryError:
        print("Instructions Copying Error:\n  Destination is a directory.")

    # if there is any permission issue
    except PermissionError:
        print("Instructions Copying Error:\n  Permission denied.")

    # other errors
    except:
        print("Instructions Copying Error:\n  Error occurred while copying file.")


    return instructionsFile, title, outputDataPath, outputInstructionFile




def main():

    analysisOptions = ['isoAnalysis','ellipsAnalysis', 'SurfaceExcess',        \
                        'chemFormulations', 'sldAnalysis']

    PlutoRunning = True
    while PlutoRunning:

        # selects module and creates relevant file paths
        analysisType, instructionsName, instructionsPath, inputDataPath, outputDataPath = modSelection(analysisOptions)

        # reads instructions instructionsFile, gets title and paths
        instructionsFile, title, outputDataPath, outputInstructionFile = organisePaths(analysisType, instructionsName, instructionsPath, outputDataPath)

        # calls analysis module
        if analysisType == analysisOptions[0]:
            isoAnalysis.main(instructionsFile, title, inputDataPath, outputDataPath)

        if analysisType == analysisOptions[1]:
            ellipsAnalysis.main(instructionsFile, title, inputDataPath, outputDataPath)

        if analysisType == analysisOptions[2]:
            SurfaceExcess.main(instructionsFile, title, inputDataPath, outputDataPath)

        if analysisType == analysisOptions[3]:
            chemFormulations.main(instructionsFile, outputInstructionFile)

        if analysisType == analysisOptions[4]:
            sldAnalysis.main(instructionsFile, outputInstructionFile)

    return


if __name__ == '__main__':
    print("\n\nLaunching...")
    print(" ___ _      _    ")
    print("| _ \ |_  _| |_ ___")
    print("|  _/ | || |  _/ _ \\")
    print("|_| |_|\_,_|\__\___/")
    print("\nAuthor: Samuel Winnall \nEmail:  winnall@ill.fr\n\n\n")

    main()