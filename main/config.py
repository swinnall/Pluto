" Module that defines global variables for Pluto"

####################
# Analysis Module #
###################

# set intended analysis type to True
doIsoAnalysis      = True
doEllipsAnalysis   = False
doChemFormulations = False
doSLDAnalysis      = False


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
plotPressure = False

# plot normalised pressure for injection kinetics
plotNormInjection = True

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


##########################
# Ellipsometry Analysis #
########################

# plot psi and delta against angle of incidence
plotAOI = False

# plot psi and delta against time
plotTime = False


#################
# SLD Analysis #
################

# multiplies average chain vol by 0.85 to model lipid chain compaction
compactChains = True

# for mixing drug uptake to second and third layers
useL2_L3 = True

# use vol frac (molar ratio = default); True is more accurate
useVolFrac = True

# include an injected lipid into the system
addLipid = False


##############
# Databases #
#############

# molecular weight lipid database, [g/mol]
lipidMw = {
    "DPPC": 734.039,
    "d-DPPC": 796.43,
    "POPC": 760.07,  # 760.076
    "d31-POPC": 791.07, # 791.267
    "POPS": 783.99,
    "Cholesterol": 386.65,  # 386.654
    "d45-Cholesterol": 432,
    "DLin-KC2-DMA": 642.1,
    "DLin-MC3-DMA": 642.09,
    "d-DLin-MC3-DMA": 704.5,
    "DOPE": 744.034,
    "SM": 760.223,
    "LBPA": 792.07,
    "PolyA": 385.31,
    }

# atom coherent scattering lengths [fm], Coh b from https://www.ncnr.nist.gov/resources/n-lengths/
atomSL = {
    "H": -3.739,
    "D": 6.671,
    "C": 6.646,
    "N": 9.36,
    "O": 5.803,
    "P": 5.13,
    "K": 3.67,
    }

# chemical structures for each lipid: (struct_head, struct_tail)
# Assumption: only hydrogens and carbons in the tails
lipidStruct = {
    "POPC": ('C10-H18-N-O8-P','C32-H64'),
    "d31-POPC": ('C10-H18-N-O8-P','C32-D31-H33'),
    "DOPE": ('N-O8-P','C41-H78'),
    "SM": ('N2-O6-P','C47-H93'),
    "LBPA": ('N-O10-P','C42-H82'),
    "Cholesterol": ('O-H','C27-H45'),
    "d45-Cholesterol": ('O-H','C27-D45'),
    "DLin-MC3-DMA": ('N-O2-C6-H9','C37-H70'),
    "d-DLin-MC3-DMA": ('N-O2-C6-H9','C37-H8-D62'),
    "DSPC": ('N-O8-P','C44-H88'),
    "DMG-PEG-2000": ('O50','C122-H242'),
    "PolyA": ('C10-H13-K-N5-O7-P','H')
    }

# molecular volumes for each lipid (Angstroms cubed): (head, tail)
lipidMolVol = {
    "POPC": (344,937),
    "d31-POPC": (344,937),
    "DOPE": (236,969),
    "SM": (274,953),
    "LBPA": (208,624),
    "Cholesterol": (5,624),
    "d45-Cholesterol": (5,624),
    "DLin-MC3-DMA": (290,1000),
    "d-DLin-MC3-DMA": (290,1000),
    "DSPC": (322,1000),
    "DMG-PEG-2000": (200,470),
    "PolyA": (1,1)
    }


# dict of input instructions names and output file names
pathNames = {
    "Isotherm": ('Instructions - Isotherm','Isotherm'),
    "Ellipsometry": ('Instructions - Ellipsometry','Ellipsometry'),
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
y_interval = 1

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
