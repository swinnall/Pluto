" Pluto "
" Author: @S.Winnall "
# Reads config input and prepares information for subsequent analysis

import glob, os, sys
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
    inputDir  = "../UoM-Data-Repository/input/"

    # root output directory
    outputDir = "../UoM-Data-Repository/output/"


    # ask user to pick one of the analysisOptions
    print("\n~~~\nAnalysis Options:\n %s" %analysisOptions)

    analysisChoice = input("Which analysis would you like to do? Pick the associated number (0-%d) or 'q' to exit:\n  " %(len(analysisOptions)-1))

    if analysisChoice == 'q':
        print("Session closed.")
        sys.exit()

    elif analysisChoice in [str(i) for i in range(len(analysisOptions))]:
        analysisType   = analysisOptions[int(analysisChoice)]
        print("You picked %s.py\n" %analysisType)

    else:
        print("Not a valid response. Session closed.")
        sys.exit()


    # calls from config database
    instructionsName = config.pathNames.get(analysisType)[0]


    # instructions file path
    instructionsPath = "" + inputDir + instructionsName + ".txt"


    # create input path for data source
    inputDataPath = '' + inputDir + '/00 - ' + config.pathNames.get(analysisType)[1] + ''


    return analysisType, outputDir, instructionsName, instructionsPath, inputDataPath



def organisePaths(analysisType, outputDir, instructionsName, instructionsPath, inputDataPath):

    # read instructions
    with open(instructionsPath, newline = '') as f:
        reader = csv.reader(f, delimiter=",")
        info = list(reader)


    # filter out rows that start with '#'
    info = [x for x in info if not x[0].startswith('#')]


    # gets title of the analysis
    title = info[0][0].split('=')[1]


    # create output path for analysis
    outputPath = '' + outputDir + '/' + config.pathNames.get(analysisType)[1] + '/' + title + ''


    # delete folder if exists and create it again
    try:
        shutil.rmtree(outputPath)
        os.mkdir(outputPath)
    except FileNotFoundError:
        os.mkdir(outputPath)


    try:
        # copy input instructions information to outputPath directory
        shutil.copyfile(instructionsPath, outputPath + '/' + instructionsName + '.txt')

        # rename the copied instructions file to include the title of the analysis
        old_name           = outputPath + '/' + instructionsName + '.txt'
        outputDataFilePath = outputPath + '/' + title + ' - ' + instructionsName + '.txt'
        os.rename(old_name, outputDataFilePath)


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

    return info, title, outputPath, outputDataFilePath




def main():

    analysisOptions = ['isoAnalysis','ellipsAnalysis','chemFormulations',      \
                        'sldAnalysis','SurfaceExcess']

    PlutoRunning = True
    while PlutoRunning:

        # get name of instructions file
        analysisType, outputDir, instructionsName, instructionsPath, inputDataPath = modSelection(analysisOptions)

        # get file and path information
        info, title, outputPath, outputDataFilePath = organisePaths(analysisType, outputDir, instructionsName, instructionsPath)

        # calls analysis module
        if analysisType == analysisOptions[0]:
            isoAnalysis.main(info, title, inputDataPath, outputPath)

        if analysisType == analysisOptions[1]:
            ellipsAnalysis.main(info, title, inputDataPath, outputPath)

        if analysisType == analysisOptions[2]:
            chemFormulations.main(info, outputDataFilePath)

        if analysisType == analysisOptions[3]:
            sldAnalysis.main(info, outputDataFilePath)

        if analysisType == analysisOptions[4]:
            SurfaceExcess.main(info, title, inputDataPath, outputPath)

    return


if __name__ == '__main__':
    print("\n\nLaunching...")
    print(" ___ _      _    ")
    print("| _ \ |_  _| |_ ___")
    print("|  _/ | || |  _/ _ \\")
    print("|_| |_|\_,_|\__\___/")
    print("\nAuthor: Samuel Winnall \nEmail:  winnall@ill.fr\n\n\n")

    main()
