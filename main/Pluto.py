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

    # ask user to pick one of the analysisOptions
    print("Analysis Options: %s" %analysisOptions)

    analysisChoice = input("Which analysis would you like to do? Pick the associated number (0-4) or 'q' to exit:\n  ")

    if analysisChoice == 'q':
        print("Session closed.")
        sys.exit()

    elif analysisChoice in [str(i) for i in range(len(analysisOptions))]:
        analysisType   = analysisOptions[int(analysisChoice)]
        print("You picked %s.py\n" %analysisType)

    else:
        print("Not a valid response. Session closed.")
        sys.exit()


    # root input directory
    root = "../input/"

    # calls from config database
    fname = config.pathNames.get(analysisType)[0]

    # instructions file path
    source = "" + root + fname + ".txt"

    return analysisType, fname, source



def organisePaths(analysisType, fname, source):

    # read instructions
    with open(source, newline = '') as f:
        reader = csv.reader(f, delimiter=",")
        info = list(reader)


    # filter out rows that start with '#'
    info = [x for x in info if not x[0].startswith('#')]


    # gets title of the analysis
    title = info[0][0].split('=')[1]


    # create input path for data source
    inputPath = '../input/00 - ' + config.pathNames.get(analysisType)[1] + ''


    # create output path for analysis
    outputPath = '../output/' + config.pathNames.get(analysisType)[1] + '/' + title + ''


    # delete folder if exists and create it again
    try:
        shutil.rmtree(outputPath)
        os.mkdir(outputPath)
    except FileNotFoundError:
        os.mkdir(outputPath)


    # copy input analysis information to outputPath
    try:
        shutil.copyfile(source, outputPath + '/' + fname + '.txt')
        old_name = outputPath + '/' + fname + '.txt'
        new_name = outputPath + '/' + title + ' - ' + fname + '.txt'
        os.rename(old_name, new_name)


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

    return info, title, inputPath, outputPath




def main():

    analysisOptions = ['isoAnalysis','ellipsAnalysis','chemFormulations',      \
                        'sldAnalysis','SurfaceExcess']

    PlutoRunning = True
    while PlutoRunning:

        # get name of instructions file
        analysisType, fname, source = modSelection(analysisOptions)

        # get file and path information
        info, title, inputPath, outputPath = organisePaths(analysisType, fname, source)

        # calls analysis module
        if analysisType == analysisOptions[0]:
            isoAnalysis.main(info, title, inputPath, outputPath)

        if analysisType == analysisOptions[1]:
            ellipsAnalysis.main(info, title, inputPath, outputPath)

        if analysisType == analysisOptions[2]:
            chemFormulations.main(info, new_name)

        if analysisType == analysisOptions[3]:
            sldAnalysis.main(info, new_name)

        if analysisType == analysisOptions[4]:
            SurfaceExcess.main(info, title, inputPath, outputPath)

    return


if __name__ == '__main__':
    print("\n\nLaunching...")
    print(" ___ _      _    ")
    print("| _ \ |_  _| |_ ___")
    print("|  _/ | || |  _/ _ \\")
    print("|_| |_|\_,_|\__\___/")
    print("\nAuthor: Samuel Winnall \nEmail:  winnall@ill.fr\n\n\n")

    main()
