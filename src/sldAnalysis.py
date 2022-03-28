" Calculate SLD profiles "

import csv
import glob, os
import re
import sys
import config

class Membrane:

    def __init__(self, lipidNames, molRatios, thickness, monolayerPar, outputFilePath):

        # read system parameters
        self.lipidNames  = lipidNames
        self.molRatios   = molRatios
        self.thickness   = thickness
        self.nLipids     = len(self.lipidNames)

        # unpack monolayer parameters (0s if 1st run; each tuple ele is dict[struct])
        self.monolayerMolVol = monolayerPar[0]
        self.monolayerSL     = monolayerPar[1]
        self.monolayerSLD    = monolayerPar[2]

        # import filepath
        self.outputFilePath = outputFilePath

        # import databases
        self.lipidStruct = config.lipidStruct
        self.atomSL      = config.atomSL
        self.lipidMolVol = config.lipidMolVol

        # initialise new variables
        self.headVolFrac   = 0
        self.twoSolv       = 0
        self.d3            = 0
        self.threeSolv     = 0
        self.twoSLD_H2O    = 0
        self.twoSLD_D2O    = 0
        self.normMolRatios = []
        self.volFrac       = {}
        self.lipidSL       = {}
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

        if config.verbose == True and "Monolayer" not in self.lipidNames:
            print("\nLipid names:\n%s\n\nInput molar ratios:\n%s" %(self.lipidNames, self.molRatios))

        if config.very_verbose == True:
            print("\nNormalised molar ratios:\n%s" %self.normMolRatios)


    # calculates the total lipid volume
    def calcTotalLipidVol(self):

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
                    self.totalLipidVol[struct] += self.normMolRatios[i] * self.monolayerMolVol[struct]
                else:
                    self.totalLipidVol[struct] += self.normMolRatios[i] * self.lipidMolVol.get(lipid)[j]


        # save Monolayer monolayer struct volumes on the first iteration
        if "Monolayer" not in self.lipidNames:
            self.monolayerMolVol['head']  = self.totalLipidVol.get('head')
            self.monolayerMolVol['tails'] = self.totalLipidVol.get('tails')

        if config.very_verbose == True:
            print('\nTotal Volume:\n%s' %self.totalLipidVol)



    # calculates volume fraction
    def calcVolFrac(self):

        warningChoice = input("\nWarning: You have set volFrac = True. Do you want to continue? (y/n)\n ")

        if warningChoice.upper() == 'N':
            print("Session closed, you must change config parameter.")
            sys.exit()


        # component volume calculation
        for i, lipid in enumerate(self.lipidNames):
            self.volFrac[lipid] = {'head': 0, 'tails': 0}

            for j, struct in enumerate(['head','tails']):
                if lipid == "Monolayer":
                    self.volFrac[lipid][struct] = self.normMolRatios[i] * self.monolayerMolVol.get(struct) / self.totalLipidVol[struct]
                else:
                    self.volFrac[lipid][struct] = self.normMolRatios[i] * self.lipidMolVol.get(lipid)[j] / self.totalLipidVol[struct]

        if config.verbose == True:
            print('\nComponent Volume Fraction:\n%s' %self.volFrac)



    # calculates the scattering length of each lipid component
    def calcSL(self):

        for i, lipid in enumerate(self.lipidNames):
            self.lipidSL[lipid] = {'head': 0, 'tails': 0}

            for j, struct in enumerate(['head','tails']):

                if lipid == "Monolayer":
                    self.lipidSL[lipid][struct] = self.monolayerSL.get(struct)

                else:
                    # splits head/tail into list of constituent atoms
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
            self.monolayerSL['head']  = self.avSL.get('head')
            self.monolayerSL['tails'] = self.avSL.get('tails')

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
                        self.avLipidVol[struct] += self.volFrac[lipid][struct] * self.monolayerMolVol[struct]
                    else:
                        self.avLipidVol[struct] += self.volFrac[lipid][struct] * self.lipidMolVol.get(lipid)[j]

                else:
                    if lipid == "Monolayer":
                        self.avLipidVol[struct] += self.normMolRatios[i] * self.monolayerMolVol[struct]
                    else:
                        self.avLipidVol[struct] += self.normMolRatios[i] * self.lipidMolVol.get(lipid)[j]


        # if accounting for chain compaction
        if config.compactChains == True:
            chainCompactFactor = config.chainCompactFactor
            self.avLipidVol['tails'] = chainCompactFactor * self.avLipidVol['tails']


        if config.verbose == True:
            print("\nAverage Lipid Head/Tail Volume:\n%s" %self.avLipidVol)



    # calculates the scattering length density of each lipid component
    def calcSLD(self):

        # take avSL and divide by avStructVol (tails, heads)
        for i, struct in enumerate(['head','tails']):
            self.avSLD[struct] = 10 * self.avSL[struct] / self.avLipidVol[struct]

        # save Monolayer monolayer struct volumes on the first iteration
        if "Monolayer" not in self.lipidNames:
            self.monolayerSLD['head']  = self.avSLD.get('head')
            self.monolayerSLD['tails'] = self.avSLD.get('tails')

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

        # set new membrane thickness if intended
        elif "Monolayer" in self.lipidNames:
            if config.updateMonolayerThickness == True:
                d1 = config.new_d1
                d2 = config.new_d2
                print("\nYou have selected to update the monolayer thickness such that d1 = %.3f and d2 = %.3f." %(d1,d2))

            else:
                d1 = self.thickness.get("tails")
                d2 = self.thickness.get("head")
                print("\nMonolayer thicknesses are as in instructions file; d1 = %.3f and d2 = %.3f." %(d1,d2))

        # calculate volume fraction
        self.headVolFrac = (rho1 * d1 * SL2 ) / ( rho2 * d2 * SL1)
        self.twoSolv     = (1-self.headVolFrac)*100

        # solvent volume in head group
        if config.verbose == True:
            print("\nHead volume fraction: %f" %self.headVolFrac)
        print("\n2-solv = %f" %self.twoSolv)



    # calculates 2-SLD for Igor Motofit; used when merging L3 drug with L2 headgroup
    def calcSLDMergedL2_L3(self):

        # update thickness as part of the drug enters the second layer
        self.d3 = config.drugSize - self.thickness.get("head")

        # get solvent amount in third layer from neutron fit, set in config
        self.threeSolv = config.threeSolv

        # calculate the volume fraction of drug in the third layer
        vf_drug = (100 - self.threeSolv) / 100

        # update 2-Solv where it is assumed the vol.frac. of drug is the same in both layers
        self.twoSolv = self.twoSolv - (vf_drug*100)

        # get SLD values of drug in H2O and D2O
        SLD_drug_H2O  = config.SLD_drug_H2O
        SLD_drug_D2O  = config.SLD_drug_D2O

        # averages the existing SLDs of layer 2 (headgroup layer) with the added drug
        self.twoSLD_H2O = (self.headVolFrac*self.avSLD.get("head") + vf_drug*SLD_drug_H2O ) / (self.headVolFrac + vf_drug)
        self.twoSLD_D2O = (self.headVolFrac*self.avSLD.get("head") + vf_drug*SLD_drug_D2O ) / (self.headVolFrac + vf_drug)

        # calculate 3-SLD, average against buffer
        

        # print values to terminal
        print("\n\nYou have chosen to add the drug (%s) to the second layer:" %config.drugName)
        print("\n3-thick = %f" %self.d3)
        print("\n2-solv = %f\n3-solv = %f" %(self.twoSolv, self.threeSolv))
        print("\n2-SLD_H2O = %f\n2-SLD_D2O = %f" %(self.twoSLD_H2O, self.twoSLD_D2O))
        print("\n3-SLD_H2O = %f\n3-SLD_D2O = %f" %(SLD_drug_H2O, SLD_drug_D2O))


    # method that returns monolayer information for adding extra lipid components
    def getMonolayerPar(self):
        return (self.monolayerMolVol, self.monolayerSL, self.monolayerSLD)


    # write output to file
    def appendFile_Membrane(self):

        with open(self.outputFilePath, 'a', newline = '') as f:
            f.write('\n~~~~\nCalculated Membrane: %s\n' %self.lipidNames)
            f.write("Average SL;  Head = %.4f; Tail = %.4f\n" %(self.avSL.get("head"),self.avSL.get("tails")))
            f.write("Average SLD; Head = %.4f; Tail = %.4f\n" %(self.avSLD.get("head"),self.avSLD.get("tails")))
            f.write("Thickness;   Head = %.4f; Tail = %.4f\n" %(self.thickness.get("head"),self.thickness.get("tails")))
            f.write("Head vol frac = %.4f\n" %self.headVolFrac)
            f.write("2-solv = %.4f\n" %self.twoSolv)


    def appendFile_lipidExchangedMembrane(self):

        with open(self.outputFilePath, 'a', newline = '') as f:
            f.write('\nCalculated Membrane: %s\n' %self.lipidNames)
            f.write("Average SL;  Head = %.4f; Tail = %.4f\n" %(self.avSL.get("head"),self.avSL.get("tails")))
            f.write("Average SLD; Head = %.4f; Tail = %.4f\n" %(self.avSLD.get("head"),self.avSLD.get("tails")))
            f.write("Thickness;   Head = %.4f; Tail = %.4f\n" %(self.thickness.get("head"),self.thickness.get("tails")))
            f.write("Head vol frac = %.4f\n" %self.headVolFrac)
            f.write("2-solv = %.4f\n" %self.twoSolv)


    def appendFile_bindingToMembrane(self):

        with open(self.outputFilePath, 'a', newline = '') as f:
            f.write("\nDrug (%s) added to layer 2 of monolayer:\n" %config.drugName)
            f.write("d3 = %f\n" %self.d3)
            f.write("2-solv = %f\n" %self.twoSolv)
            f.write("3-solv = %.4f\n" %self.threeSolv)
            f.write("2-SLD_H2O = %f\n" %self.twoSLD_H2O)
            f.write("2-SLD_D2O = %f\n" %self.twoSLD_D2O)



def importSampleData(instructionsFile, sampleNum):

    membrane   = instructionsFile["membranes"][sampleNum]
    lipidRatio = instructionsFile["lipidRatios"][sampleNum]
    t_thick    = instructionsFile["d1"][sampleNum]
    h_thick    = instructionsFile["d2"][sampleNum]
    label      = instructionsFile["label"][sampleNum]

    return membrane, lipidRatio, t_thick, h_thick, label


# function to check whether a string contains a number
def hasNumbers(inputString):
    return bool(re.search(r'\d', inputString))



def main(instructionsFile, outputFilePath):

    # number of membranes to calculate
    nMemb = len(instructionsFile)

    # calculate component volumes
    for sampleNum in range(nMemb):

        # import sample data
        membrane, lipidRatio, t_thick, h_thick, label = importSampleData(instructionsFile, sampleNum)

        # print membrane label to terminal
        print("\n\n\n~~~\nMembrane %d - %s" %(sampleNum+1,label))

        # get list of lipids within membrane with ratio
        lipids = re.split(':',membrane)
        ratios = re.split(':',str(lipidRatio))

        # structure thickness information into dict
        thickness = {
            'head':  float(h_thick),
            'tails': float(t_thick)
        }

        # initialise monolayer information
        monolayerMolVol = {'head': 0, 'tails': 0}
        monolayerSL     = {'head': 0, 'tails': 0}
        monolayerSLD    = {'head': 0, 'tails': 0}
        monolayerPar    = (monolayerMolVol, monolayerSL, monolayerSLD)

        # create class instance with input variables
        m = Membrane(lipids, ratios, thickness, monolayerPar, outputFilePath)

        # converts 3:5 to 3/8:5/8
        m.normaliseMolarRatios()

        # calculates the total lipid volume of the monolayer (needed for addLipidToMonolayer)
        m.calcTotalLipidVol()

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

        # append file with initial monolayer calculation information
        m.appendFile_Membrane()

        # repeat for incorporating an injected lipid component into the membrane
        if config.addLipidToMonolayer == True:

            # get monlayer parameter information for next iteration
            monolayerPar = m.getMonolayerPar()

            if config.very_verbose == True: print("\nMonolayerPar:\n%s" %(monolayerPar,))

            # get added lipids and associated ratios from config file
            lipids = config.injectedComponents
            ratios = config.injectedRatios

            print("\n\n\nYou have chosen to add components from injected sample to the monolayer at the following ratios:")
            sumRatiosTest = 0
            for ele, lipid in enumerate(lipids):
                sumRatiosTest += ratios[ele]
                if lipid == "Monolayer": print("Averaged %s: %d%%." %(lipid, ratios[ele]))
                else: print("Added Lipid: %s: %d%%." %(lipid, ratios[ele]))

            if sumRatiosTest != 100:
                print("Fatal Error: Molar ratios defined in config do not equal 100.")
                sys.exit()

            # creates a new class instance to pass new config params to
            m = Membrane(lipids, ratios, thickness, monolayerPar, outputFilePath)

            m.normaliseMolarRatios()

            m.calcTotalLipidVol()

            if config.useVolFrac == True: m.calcVolFrac()

            m.calcSL()

            m.calcAvLipidVol()

            m.calcSLD()

            m.calcHeadVolumeFraction()

            m.appendFile_lipidExchangedMembrane()


        # calculates SLD for mixing drug in layer 3 to layer 2
        if config.addDrugToMonolayer == True:

            # uses existing class instance from either initial monolayer or updated
            m.calcSLDMergedL2_L3()

            m.appendFile_bindingToMembrane()

    sys.exit()
    return



if __name__ == '__main__':
    print("~Running sldAnalysis.py~")
    main()
