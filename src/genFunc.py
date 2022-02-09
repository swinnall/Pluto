" Module for general functions within Pluto "
" Next: add in geneneral smooth data functions "

import csv
import sys
import pandas as pd


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

    elif analysisChoice not in [str(i) for i in range(len(analysisOptions)+1)] and analysisOptions[0]!="isoAnalysis":
        print("Not a valid response. Returning to Pluto landing page.\n\n")
        analysisType    = 'n/a'
        analysisRunning = False

    else:
        print("Not a valid response. Session closed.")
        sys.exit()


    return analysisType, analysisRunning



# general function for isotherm and surface excess modules
def getFile(fileDIR):

    df = pd.read_csv(fileDIR, sep='\t', na_values =' ')

    return df


# eventually combine with the more general getFile
def getEllipsometryFile(fileDIR):

    df = []
    with open(fileDIR, newline = '') as f:

        # ignore rows starting with '#'
        rdr = csv.DictReader(filter(lambda row: row[0]!='#', f))

        for row in rdr:
            df.append( row[None][0].split())

    return df
