" Calculate SLD profiles "

import csv
import glob, os
import re
import sys
import config

class Membrane:

    def __init__(self, lipidNames, molRatios, outputPath):
        self.lipidNames = lipidNames
        self.molRatios  = molRatios
        self.nLipids    = len(self.lipidNames)

        # import databases
        self.lipidStruct = config.lipidStruct
        self.atomSL      = config.atomSL
        self.lipidMolVol = config.lipidMolVol

        # import filepath
        self.outputPath = outputPath

        # initialise new variables
        self.headVolFrac   = 0
        self.volFrac       = {}
        self.lipidSL       = {}
        self.lipidSLD      = {}
        self.totalLipidVol = {'head': 0, 'tails': 0}
        self.avSL          = {'head': 0, 'tails': 0}
        self.avSLD         = {'head': 0, 'tails': 0}



    # calculates the total volume by summing lipid component volumes
    def totalVol(self):

        global monolayerMolVol

        for i, lipid in enumerate(self.lipidNames):

            # check lipid exists in database
            if self.lipidStruct.get(lipid) == None and lipid == "Monolayer": pass
            elif self.lipidStruct.get(lipid) == None:
                print("\nError: Lipid type not found in Lipid Molecular Formula Database.")
                print("Lipid: %s" %lipid)
                sys.exit()

            for j, struct in enumerate(['head','tails']):

                if lipid == "Monolayer":
                    self.totalLipidVol[struct] += (float(self.molRatios[i])/100) * monolayerMolVol[struct]
                else:
                    self.totalLipidVol[struct] += (float(self.molRatios[i])/100) * self.lipidMolVol.get(lipid)[j]


        # save Monolayer monolayer struct volumes on the first iteration
        if "Monolayer" not in self.lipidNames:
            monolayerMolVol['head']  = self.totalLipidVol.get('head')
            monolayerMolVol['tails'] = self.totalLipidVol.get('tails')

        if config.verbose == True and "Monolayer" not in self.lipidNames:
            print('\nTotal Volume:\n%s' %self.totalLipidVol)

        if config.verbose == True and "Monolayer" in self.lipidNames:
            print('\nTotal Volume:\n%s' %self.totalLipidVol)



    # calculates the volume fractions (replacing input molar fractions)
    def calcVolFrac(self):

        global monolayerMolVol

        for i, lipid in enumerate(self.lipidNames):
            self.volFrac[lipid] = {'head': 0, 'tails': 0}

            for j, struct in enumerate(['head','tails']):

                if lipid == "Monolayer":
                    self.volFrac[lipid][struct] = (float(self.molRatios[i])/100) * monolayerMolVol.get(struct) / self.totalLipidVol[struct]
                else:
                    self.volFrac[lipid][struct] = (float(self.molRatios[i])/100) * self.lipidMolVol.get(lipid)[j] / self.totalLipidVol[struct]

        if config.verbose == True:
            print('\nComponent Volume Fraction:\n%s' %self.volFrac)



    # calculates the scattering length of each lipid component
    def calcSL(self):

        global monolayerSL

        for i, lipid in enumerate(self.lipidNames):
            self.lipidSL[lipid] = {'head': 0, 'tails': 0}

            for j, struct in enumerate(['head','tails']):

                # if Monolayer isn't in lipidStruct (it's not)
                # store data from global variable and break loop
                if lipid == "Monolayer":
                    self.lipidSL[lipid][struct] = monolayerSL.get(struct)
                    #if struct == 'tails': break

                else:
                    # splits head/tail into list of constituent
                    splitStruct = re.split('-', self.lipidStruct.get(lipid)[j])

                    # this loop iterates across elements in a given head/tail and sums the atomic scattering lengths
                    for ele in splitStruct:

                        # multiply scattering length of atom by number of atoms
                        if hasNumbers(ele) == True:

                            # x is split into ['atom','number of given atom']
                            x = list(filter(None, re.split(r'(\d+)', ele)))

                            self.lipidSL[lipid][struct] += self.atomSL.get(x[0]) * int(x[1])


                        # add the scattering length of the single atom identified
                        elif hasNumbers(ele) == False:
                            self.lipidSL[lipid][struct] += self.atomSL.get(ele[0])


                    # multiply total lipid scattering length of a given lipid's head/tail by corresponding vol frac
                    self.avSL[struct] += self.volFrac[lipid][struct] * self.lipidSL[lipid][struct]


        # save Monolayer monolayer struct volumes on the first iteration
        if "Monolayer" not in self.lipidNames:
            monolayerSL['head']  = self.avSL.get('head')
            monolayerSL['tails'] = self.avSL.get('tails')

        if config.very_verbose == True:
            print("\nLipid scattering lengths:\n%s" %self.lipidSL)

        if config.verbose == True:
            print("\nAverage SL:\n%s" %self.avSL)



    # calculates the scattering length density of each lipid component
    def calcSLD(self):

        global monolayerSLD

        for i, lipid in enumerate(self.lipidNames):
            self.lipidSLD[lipid] = {'head': 0, 'tails': 0}

            for j, struct in enumerate(['head','tails']):

                # divide the scattering length by the molecular volume for each lipid or call from monolayer
                if lipid == "Monolayer":
                    self.lipidSLD[lipid][struct] = monolayerSLD.get(struct)
                else:
                    self.lipidSLD[lipid][struct] = 10 * self.lipidSL[lipid][struct] / self.lipidMolVol.get(lipid)[j]

                # add each one to the average lipid
                self.avSLD[struct] += self.volFrac[lipid][struct] * self.lipidSLD[lipid][struct]


        # save Monolayer monolayer struct volumes on the first iteration
        if "Monolayer" not in self.lipidNames:
            monolayerSLD['head']  = self.avSLD.get('head')
            monolayerSLD['tails'] = self.avSLD.get('tails')

        if config.very_verbose == True:
            print("\nLipid scattering length densities:\n%s" %self.lipidSLD)

        if config.verbose == True:
            print("\nAverage SLD:\n%s" %self.avSLD)



    # calculates volume fraction of the head group based on SLD
    def calcHeadVolumeFraction(self):

        # set initial membrane thickness (1: tails, 2: heads)
        if "Monolayer" not in self.lipidNames:
            d1 = 10
            d2 = 10

        # set new membrane thickness (1: tails, 2: heads)
        elif "Monolayer" in self.lipidNames:
            d1 = 15
            d2 = 10

        # select SLD (1: tails, 2: heads)
        rho1 = self.avSLD.get("tails")
        rho2 = self.avSLD.get("head")

        # select SL (1: tails, 2: heads)
        SL1 = self.avSL.get("tails")
        SL2 = self.avSL.get("head")

        # calculate volume fraction
        self.headVolFrac = (rho1 * d1 * SL2 ) / ( rho2 * d2 * SL1)

        print("\nHead volume fraction: %f" %self.headVolFrac)



    # write output to file
    def appendFile(self):

        ## Write to file
        with open(self.outputPath, 'a', newline = '') as f:

            f.write('\nCalculated Membrane: %s\n' %self.lipidNames)
            f.write("Average SLD; Head = %.2f; Tail = %.2f\n" %(self.avSLD.get("head"),self.avSLD.get("tails")))
            f.write("Head vol frac = %.2f\n" %self.headVolFrac)



def importSampleData(info):

    # number of header rows
    nHeader = 2

    # number of membranes to analyse
    nMemb  = len(info) - nHeader

    # variable initialisations
    membranes   = {init: 0 for init in range(nMemb)}
    lipidRatios = {init: 0 for init in range(nMemb)}
    labels      = {init: 0 for init in range(nMemb)}


    # assign data
    for i in range(nHeader,len(info)):

        # correct for indexing
        j = i - nHeader

        membranes[j]   = info[i][0]
        lipidRatios[j] = info[i][1]
        labels[j]      = info[i][2]


    return nMemb, membranes, lipidRatios, labels



def hasNumbers(inputString):
    # function to check whether a string contains a number
    return bool(re.search(r'\d', inputString))



def main(info, outputPath):

    # initialise global variables for averageAgain func
    global monolayerMolVol
    global monolayerSL
    global monolayerSLD
    monolayerMolVol = {'head': 0, 'tails': 0}
    monolayerSL     = {'head': 0, 'tails': 0}
    monolayerSLD    = {'head': 0, 'tails': 0}


    # import txt instructions
    nMemb, membranes, lipidRatios, labels = importSampleData(info)


    # calculate component volumes
    for i in range(nMemb):

        print("\n\n\nMembrane %d - %s" %(i+1,labels.get(i)))

        lipids = re.split(':',membranes.get(i))
        ratios = re.split(':',lipidRatios.get(i))

        # create class instance with input variables
        m = Membrane(lipids, ratios, outputPath)

        # calculate the total volume of the lipid components
        m.totalVol()

        # convert molar fraction to component volume fraction
        m.calcVolFrac()

        # calculate coherent scattering lengths
        m.calcSL()

        # divide scattering length by the molecular volume
        m.calcSLD()

        # calculate volume fraction of the headgroups
        m.calcHeadVolumeFraction()

        # update copied input instructions with output
        m.appendFile()

        # repeat for incorporating an injected lipid component into the membrane
        if config.addLipid == True:

            print("\n\n\nAdding injected component -")

            lipids = ["Monolayer", "DLin-MC3-DMA"]
            ratios = [80, 20]

            m = Membrane(lipids, ratios, outputPath)
            m.totalVol()
            m.calcVolFrac()
            m.calcSL()
            m.calcSLD()
            m.calcHeadVolumeFraction()
            m.appendFile()


    return



if __name__ == '__main__':
    print("~Running sldAnalysis.py~")
    main()
