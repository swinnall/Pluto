" Calculate SLD profiles "

import csv
import glob, os
import re
import sys
import config

class Membrane:

    def __init__(self, lipidNames, molRatios, thickness, outputFilePath):
        self.lipidNames  = lipidNames
        self.molRatios   = molRatios
        self.thickness   = thickness
        self.nLipids     = len(self.lipidNames)

        # import databases
        self.lipidStruct = config.lipidStruct
        self.atomSL      = config.atomSL
        self.lipidMolVol = config.lipidMolVol

        # import filepath
        self.outputFilePath = outputFilePath

        # initialise new variables
        self.headVolFrac   = 0
        self.twoSolv       = 0
        self.d3            = 0
        self.twoSLD_H2O    = 0
        self.twoSLD_D2O    = 0
        self.normMolRatios = []
        self.volFrac       = {}
        self.lipidSL       = {}
        self.lipidSLD      = {}
        self.totalLipidVol = {'head': 0, 'tails': 0}
        self.avLipidVol    = {'head': 0, 'tails': 0}
        self.avSL          = {'head': 0, 'tails': 0}
        self.avSLD         = {'head': 0, 'tails': 0}


    # converts molar ratio txt input to a normalised version e.g. 4:5 ==> 4/9:5/9
    def normaliseMolarRatios(self):

        totMol = 0
        for i, value in enumerate(self.molRatios):
            totMol += float(value)

        for i, value in enumerate(self.molRatios):
            self.normMolRatios.append( float(value) / totMol  )


    # calculates the total volume by summing lipid component volumes
    def calcVolFrac(self):

        global monolayerMolVol

        for i, lipid in enumerate(self.lipidNames):

            # check lipid exists in database
            if self.lipidStruct.get(lipid) == None and lipid == "Monolayer": pass
            elif self.lipidStruct.get(lipid) == None:
                print("\nError: Lipid type not found in Lipid Molecular Formula Database.")
                print("Lipid: %s" %lipid)
                sys.exit()

            # calculate total lipid volumes
            for j, struct in enumerate(['head','tails']):

                if lipid == "Monolayer":
                    self.totalLipidVol[struct] += self.normMolRatios[i] * monolayerMolVol[struct]
                else:
                    self.totalLipidVol[struct] += self.normMolRatios[i] * self.lipidMolVol.get(lipid)[j]


        # component volume calculation
        for i, lipid in enumerate(self.lipidNames):
            self.volFrac[lipid] = {'head': 0, 'tails': 0}

            for j, struct in enumerate(['head','tails']):
                if lipid == "Monolayer":
                    self.volFrac[lipid][struct] = self.normMolRatios[i] * monolayerMolVol.get(struct) / self.totalLipidVol[struct]
                else:
                    self.volFrac[lipid][struct] = self.normMolRatios[i] * self.lipidMolVol.get(lipid)[j] / self.totalLipidVol[struct]



        # save Monolayer monolayer struct volumes on the first iteration
        if "Monolayer" not in self.lipidNames:
            monolayerMolVol['head']  = self.totalLipidVol.get('head')
            monolayerMolVol['tails'] = self.totalLipidVol.get('tails')

        if config.verbose == True and "Monolayer" not in self.lipidNames:
            print('\nTotal Volume:\n%s' %self.totalLipidVol)

        if config.verbose == True and "Monolayer" in self.lipidNames:
            print('\nTotal Volume:\n%s' %self.totalLipidVol)

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
                    if config.useVolFrac == True:
                        self.avSL[struct] += self.volFrac[lipid][struct] * self.lipidSL[lipid][struct]
                    else:
                        self.avSL[struct] += self.normMolRatios[i] * self.lipidSL[lipid][struct]



        # save Monolayer monolayer struct volumes on the first iteration
        if "Monolayer" not in self.lipidNames:
            monolayerSL['head']  = self.avSL.get('head')
            monolayerSL['tails'] = self.avSL.get('tails')

        if config.very_verbose == True:
            print("\nLipid scattering lengths:\n%s" %self.lipidSL)

        if config.verbose == True:
            print("\nAverage SL:\n%s" %self.avSL)



    def calcAvLipidVol(self):

        for i, lipid in enumerate(self.lipidNames):

            # calculate total lipid volumes
            for j, struct in enumerate(['head','tails']):

                if config.useVolFrac == True:
                    if lipid == "Monolayer":
                        self.avLipidVol[struct] += self.volFrac[lipid][struct] * monolayerMolVol[struct]
                    else:
                        self.avLipidVol[struct] += self.volFrac[lipid][struct] * self.lipidMolVol.get(lipid)[j]

                else:
                    if lipid == "Monolayer":
                        self.avLipidVol[struct] += self.normMolRatios[i] * monolayerMolVol[struct]
                    else:
                        self.avLipidVol[struct] += self.normMolRatios[i] * self.lipidMolVol.get(lipid)[j]


        # if accounting for chain compaction
        if config.compactChains == True:
            chainCompactFactor = 0.85
            self.avLipidVol['tails'] = chainCompactFactor * self.avLipidVol['tails']


        if config.verbose == True:
            print("\nAverage Lipid Head/Tail Volume:\n%s" %self.avLipidVol)



    # calculates the scattering length density of each lipid component
    def calcSLD(self):

        global monolayerSLD

        # take avSL and divide by avStructVol (tails, heads)
        for i, struct in enumerate(['head','tails']):
            self.avSLD[struct] = 10 * self.avSL[struct] / self.avLipidVol[struct]


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

        # call SL
        SL1 = self.avSL.get("tails")
        SL2 = self.avSL.get("head")

        # call SLD
        rho1 = self.avSLD.get("tails")
        rho2 = self.avSLD.get("head")

        # call membrane thickness
        if "Monolayer" not in self.lipidNames:
            d1 = self.thickness.get("tails")
            d2 = self.thickness.get("head")

        # set new membrane thickness
        elif "Monolayer" in self.lipidNames:
            d1 = 15
            d2 = 10

        # calculate volume fraction
        self.headVolFrac = (rho1 * d1 * SL2 ) / ( rho2 * d2 * SL1)
        self.twoSolv     = (1-self.headVolFrac)*100

        # solvent volume in head group
        print("\nHead volume fraction: %f" %self.headVolFrac)
        print("\n2-solv = %f" %self.twoSolv)


    # calculates SLD2 for Igor Motofit; used when merging L3 with L2
    def calcSLD2(self):

        self.d3 = 20 - self.thickness.get("head")

        threeSolv     = 91.374
        vf_polyA      = (100 - threeSolv) / 100
        SLD_polyA_H2O = 3.67
        SLD_polyA_D2O = 4.46

        self.twoSLD_H2O = (self.headVolFrac*self.avSLD.get("head") + vf_polyA*SLD_polyA_H2O ) / (self.headVolFrac + vf_polyA)
        self.twoSLD_D2O = (self.headVolFrac*self.avSLD.get("head") + vf_polyA*SLD_polyA_D2O ) / (self.headVolFrac + vf_polyA)

        print("\n3-thick = %f" %self.d3)
        print("twoSLD_H2O = %f" %self.twoSLD_H2O)
        print("twoSLD_D2O = %f" %self.twoSLD_D2O)


    # write output to file
    def appendFile(self):

        ## Write to file
        with open(self.outputFilePath, 'a', newline = '') as f:

            f.write('\nCalculated Membrane: %s\n' %self.lipidNames)
            f.write("Average SL;  Head = %.4f; Tail = %.4f\n" %(self.avSL.get("head"),self.avSL.get("tails")))
            f.write("Average SLD; Head = %.4f; Tail = %.4f\n" %(self.avSLD.get("head"),self.avSLD.get("tails")))
            f.write("Thickness;   Head = %.4f; Tail = %.4f\n" %(self.thickness.get("head"),self.thickness.get("tails")))
            f.write("Head vol frac = %.4f\n" %self.headVolFrac)
            f.write("2-solv = %.4f\n" %self.twoSolv)

            if config.useL2_L3 == True:
                f.write("3-thick = %f\n" %self.d3)
                f.write("twoSLD_H2O = %f\n" %self.twoSLD_H2O)
                f.write("twoSLD_D2O = %f\n" %self.twoSLD_D2O)




def importSampleData(instructionsFile):

    # number of header rows
    nHeader = 2

    # number of membranes to analyse
    nMemb  = len(instructionsFile) - nHeader

    # variable initialisations
    membranes   = {init: 0 for init in range(nMemb)}
    lipidRatios = {init: 0 for init in range(nMemb)}
    t_thick     = {init: 0 for init in range(nMemb)}
    h_thick     = {init: 0 for init in range(nMemb)}
    labels      = {init: 0 for init in range(nMemb)}


    # assign data
    for i in range(nHeader,len(instructionsFile)):

        # correct for indexing
        j = i - nHeader

        membranes[j]   = instructionsFile[i][0]
        lipidRatios[j] = instructionsFile[i][1]
        t_thick[j]     = instructionsFile[i][2]
        h_thick[j]     = instructionsFile[i][3]
        labels[j]      = instructionsFile[i][4]



    return nMemb, membranes, lipidRatios, t_thick, h_thick, labels



def hasNumbers(inputString):
    # function to check whether a string contains a number
    return bool(re.search(r'\d', inputString))



def main(instructionsFile, outputFilePath):

    # initialise global variables for averageAgain func
    global monolayerMolVol
    global monolayerSL
    global monolayerSLD
    monolayerMolVol = {'head': 0, 'tails': 0}
    monolayerMolVol = {'head': 0, 'tails': 0}
    monolayerSL     = {'head': 0, 'tails': 0}
    monolayerSLD    = {'head': 0, 'tails': 0}


    # import txt instructions
    nMemb, membranes, lipidRatios, t_thick, h_thick, labels = importSampleData(instructionsFile)


    # calculate component volumes
    for i in range(nMemb):

        print("\n\n\nMembrane %d - %s" %(i+1,labels.get(i)))

        lipids = re.split(':',membranes.get(i))
        ratios = re.split(':',lipidRatios.get(i))

        thickness = {'head': 0, 'tails': 0}
        thickness['tails'] = float(t_thick.get(i))
        thickness['head']  = float(h_thick.get(i))


        # create class instance with input variables
        m = Membrane(lipids, ratios, thickness, outputFilePath)

        # converts 3:5 to 3/8:5/8
        m.normaliseMolarRatios()

        # convert molar fraction to component volume fraction
        if config.useVolFrac == True: m.calcVolFrac()

        # calculate coherent scattering lengths
        m.calcSL()

        # calculate average lipid structure volumes
        m.calcAvLipidVol()

        # divide scattering length by the molecular volume
        m.calcSLD()

        # calculate volume fraction of the headgroups
        m.calcHeadVolumeFraction()

        # calculates SLD2 for mixings layers 2 and 3
        if config.useL2_L3 == True: m.calcSLD2()

        # update copied input instructions with output
        m.appendFile()

        # repeat for incorporating an injected lipid component into the membrane
        if config.addLipid == True:

            print("\n\n\nAdding injected component -")

            lipids = ["Monolayer", "DLin-MC3-DMA"]
            ratios = [80, 20]

            m = Membrane(lipids, ratios, thicknesses, outputFilePath)
            m.normaliseMolarRatios()
            if config.useVolFrac == True: m.calcVolFrac()
            m.calcSL()
            m.calcAvLipidVol()
            m.calcSLD()
            m.calcHeadVolumeFraction()
            m.appendFile()


    return



if __name__ == '__main__':
    print("~Running sldAnalysis.py~")
    main()
