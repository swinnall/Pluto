" Module that defines global variables for Pluto"

####################
# Analysis Module #
###################

# set intended analysis type to True
doIsoAnalysis      = False
doChemFormulations = False
doSLDAnalysis      = True


#############
# Printing #
############

# print extra info
verbose      = True
very_verbose = False


################
# IAP Analysis #
################

# plot standard isotherm
plotIsotherm = False

# run elasticity analysis & plot
plotElasticity = False

# plot percentage area as a function of time
plotArea = False

# plot pressure as a function of time
plotPressure = True

# plot isotherm compressions only
plotCompressions = False

# plot isotherm expansions only
plotExpansions = False

# colour code cycles within isotherm
plotCycles = False

# W.I.P. overplots n separate files, colour coded by cycle / compression etc.
overplot = False

# shift to P_min = 0
shiftP = False

# order of the elasticity fit
nPoly = 40

# data reduction; take every 'nth' point
nth = 50


#################
# SLD Analysis #
################

# include an injected lipid into the system
addLipid = False


##############
# Databases #
#############

# molecular weight lipid database, [g/mol]
lipidMw = {
    "DPPC": 734.039,
    "d-DPPC": 796.43,
    "POPC": 760.076,
    "d-POPC": 791.267,
    "POPS": 783.99,
    "Cholesterol": 386.654,
    "d-Cholesterol": 432,
    "DLin-KC2-DMA": 642.1,
    "DLin-MC3-DMA": 642.09,
    "d-DLin-MC3-DMA": 704.5,
    "DOPE": 744.034,
    "SM": 760.223,
    "LBPA": 792.07,
    }

# atom coherent scattering lengths [fm], Coh b from https://www.ncnr.nist.gov/resources/n-lengths/
atomSL = {
    "H": -3.739,
    "D": 6.671,
    "C": 6.646,
    "N": 9.36,
    "O": 5.803,
    "P": 5.13,
    }

# chemical structures for each lipid: (struct_head, struct_tail)
# Assumption: only hydrogens and carbons in the tails
lipidStruct = {
    "POPC": ('N-O8-P','C42-H82'),
    "d-POPC": ('N-O8-P','C42-D31-H51'),
    "DOPE": ('N-O8-P','C41-H78'),
    "SM": ('N2-O6-P','C47-H93'),
    "LBPA": ('N-O10-P','C42-H82'),
    "Cholesterol": ('O-H','C27-H45'),
    "d-Cholesterol": ('O-H','C27-D45'),
    "DLin-MC3-DMA": ('N-O2','C43-H79'),
    "d-DLin-MC3-DMA": ('N-O2','C43-H17-D62'),
    "DSPC": ('N-O8-P','C44-H88'),
    "DMG-PEG-2000": ('O50','C122-H242'),
    }

# molecular volumes for each lipid (Angstroms cubed): (head, tail)
lipidMolVol = {
    "POPC": (344,937),
    "d-POPC": (344,937),
    "DOPE": (236,969),
    "SM": (274,953),
    "LBPA": (208,624),
    "Cholesterol": (5,624),
    "d-Cholesterol": (5,624),
    "DLin-MC3-DMA": (290,1000),
    "d-DLin-MC3-DMA": (290,1000),
    "DSPC": (322,1000),
    "DMG-PEG-2000": (200,470),
    }


# dict of different instruction files for Pluto
pathNames = {
    "Isotherm": ('Instructions - Isotherm','Isotherm'),
    "chemFormulations": ('Instructions - ChemFormulation','Chem Formulations'),
    "SLD": ('Instructions - SLD','SLD'),
    }


#############
# Plotting #
############

# font size
fs = 14

# line width
lw = 2

# number of ticks [n-1]
n_xticks = 7

# y axis tick interval
y_interval = 10

# xtick interval for P vs t plots (mins)
xmin_interval = 10

# colours
colourDict = {

    # dark blue, light blue, dark orange, light orange
    "0": ['#1e81b0','#abdbe3', '#e28743','#eab676'],

    # dark blue, light blue, dark orange, light orange, dark green, light green
    "1": ['#1e81b0','#abdbe3', '#e28743','#eab676', '#32BE25', '#A3e19d'],

    # dark blue, dark orange, dark green, dark purple, dark red, dark yellow
    "2": ['#1e81b0', '#e28743', '#32BE25', '#6A0F8E', '#AB2330', '#FFCC00'],

    }

c = colourDict.get("2")
