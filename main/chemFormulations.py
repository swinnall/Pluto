" Calculate Lipid Volume Contributions "

import csv
import glob, os
import config
import chemFormulationsMain

def importSampleData(info):

    # number of header rows
    nHeader = 2

    # number of isotherms to plot
    nLipids  = len(info) - nHeader

    # variable initialisations
    lipidType   = {init: 0 for init in range(nLipids)}
    lipidAmount = {init: 0 for init in range(nLipids)}
    stockConc   = {init: 0 for init in range(nLipids)}

    # assign data
    for i in range(nHeader,len(info)):

        # correct for indexing
        j = i - nHeader

        lipidType[j]   =       info[i][0]
        lipidAmount[j] =       info[i][1]
        stockConc[j]   = float(info[i][2])

    # different info location for weight total and final concentration
    wTotal    = float(info[1][3].split('=')[1])
    finConc   = float(info[1][4].split('=')[1])
    ratioType =       info[1][5].split('=')[1]


    return nLipids, lipidType, lipidAmount, stockConc, wTotal, finConc, ratioType




def convertToRatio(nLipids, lipidAmount):

    # initalised variable for summing
    ratioTotal = 0

    # sum the ratios
    for i in range(nLipids):
        ratioTotal += int(lipidAmount[i])

    # divide by total ratio to give fraction
    for i in range(nLipids):
        lipidAmount[i] = int(lipidAmount[i]) / ratioTotal

    return lipidAmount




def calc(nLipids, lipidMW, lipidType, lipidFrac, stockConc, wTotal, finConc):

    # initialise variable for summing
    val = 0


    # divide wTotal by the sum of Mw*fraction
    for i in range(nLipids):
        val += lipidMW.get( str(lipidType[i]) ) * lipidFrac[i]


    # molar total
    molTot = wTotal / val # [mmol]


    volAdd = []
    for i in range(nLipids):

        # molar fraction = molar total * lipid fraction
        molFrac = molTot * lipidFrac[i]

        # weight fraction = molar fraction * molecular weight
        wFrac = molFrac * lipidMW.get( str(lipidType[i]) )

        # volume to add of each lipid to get the desired wTotal [mg]
        volAdd.append( (wFrac/stockConc[i] ) * 1000 ) # [mL -> uL]

        # volume total
        Vtot = ( wTotal / finConc ) * 1000 # [mL -> uL]

        # volume of CHCl3 to get to desired concentration
        V_CHCl3 = Vtot - sum(volAdd)


    return volAdd, V_CHCl3


def output(nLipids, lipidType, stockConc, wTotal, finConc, volAdd, V_CHCl3, outputPath):

    #outputPath = chemFormulationsMain.outputPath

    ## Write to file
    with open(outputPath, 'a', newline = '') as f:

        f.write('\nCalculated Quantities:\n')
        for i in range(nLipids):
            f.write("%.2f uL of %s\n" %(volAdd[i], lipidType[i]))
        f.write("%.2f uL of chloroform" %V_CHCl3)


    ## Print to terminal
    # Input information
    print("\nInput Lipids:\n~~~~~~~~~~~~~~~")

    for i in range(nLipids):
        print( "%s of stockConc %.2f mg/ml" %(lipidType[i], stockConc[i]))

    print("\n\nChosen Weight & Concentration:\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print( "Weight total:        %.2f mg\nFinal Concentration: %.2f mg/ml" %(wTotal, finConc))


    # Ouput information
    print("\n\nRequired Amounts:\n~~~~~~~~~~~~~~~~~~~~")

    for i in range(nLipids):
        print( "%.2f uL of %s" %(volAdd[i], lipidType[i]))

    print("%.2f uL of chloroform" %V_CHCl3)

    return


def main(info, outputPath):

    # import molecular weight database
    lipidMW = config.lipidMw


    # import chemical analysis info
    nLipids, lipidType, lipidAmount, stockConc, wTotal, finConc, ratioType = importSampleData(info)


    # if lipids are given in ratio (molar percentage) convert to fractional
    if ratioType == "ratio":
        lipidAmount = convertToRatio(nLipids, lipidAmount)

    # if lipids are given as fractions, ensure float type
    elif ratioType == "frac":
        for i in range(nLipids):
            lipidAmount[i] = float(lipidAmount[i])


    # calculations function
    volAdd, V_CHCl3 = calc(nLipids, lipidMW, lipidType, lipidAmount, stockConc, wTotal, finConc)


    # print output of program to terminal
    output(nLipids, lipidType, stockConc, wTotal, finConc, volAdd, V_CHCl3, outputPath)

    return



if __name__ == '__main__':
    print("~Running chemAnalysis.py~")
    main()
