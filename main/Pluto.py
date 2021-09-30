" Pluto - S.Winnall "
# Reads config input and prepares information for subsequent analysis

import glob, os, sys
import csv
import shutil
from shutil import copyfile

import isoAnalysis
import chemFormulations
import sldAnalysis


def modSelection():

    count = 0

    if config.doIsoAnalysis == True:
        analysisType = 'Isotherm'
        count += 1

    if config.doChemFormulations == True:
        analysisType = 'chemFormulations'
        count += 1

    if config.doSLDAnalysis == True:
        analysisType = 'SLD'
        count += 1


    if count > 1:
        print("Error: More than one analysis type selected.")
        sys.exit()
    elif count == 0:
        print("Error: No analysis type selected.")
        sys.exit()


    # root input directory
    root = "../input/"

    # calls from config database
    fname = config.pathNames.get(analysisType)[0]

    # instructions file path
    source = "" + root + fname + ".txt"

    return fname, source



def main():

    # get name of instructions file
    fname, source = modSelection()


    # read instructions
    with open(source, newline = '') as f:
        reader = csv.reader(f, delimiter=",")
        info = list(reader)


    # gets title of the analysis
    title = info[0][0].split('=')[1]


    # create path for analysis
    outputFolder = '../output/' + config.pathNames.get(analysisType)[1] + '/' + title + ''


    # delete folder if exists and create it again
    try:
        shutil.rmtree(outputFolder)
        os.mkdir(outputFolder)
    except FileNotFoundError:
        os.mkdir(outputFolder)


    # copy input analysis information to outputPath
    try:
        shutil.copyfile(source, outputFolder + '/' + fname + '.txt')
        old_name = outputFolder + '/' + fname + '.txt'
        new_name = outputFolder + '/' + title + ' - ' + fname + '.txt'
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


    global outputPath
    outputPath = outputFolder + '/' + title + ' - ' + fname + '.txt'


    # calls analysis module
    if config.doIsoAnalysis == True:
        isoAnalysis.main(info, title, outputPath)

    if config.doChemFormulations == True:
        chemAnalysis.main(info, outputPath)

    if config.doSLDAnalysis == True:
        sldAnalysis.main(info)


    return


if __name__ == '__main__':
    print("\n\nLaunching Pluto...\n\n")
    main()
