" Module that defines global variables for Pluto"

##############
# File Paths #
##############

inputDir  = "../../UoM-Data-Repository/input/"
outputDir = "../../UoM-Data-Repository/output/"

# dict of input instructions names and output file names
pathNames = {
    "isoAnalysis":      ('Instructions - Surface Pressure', 'Isotherm'),
    "ellipsAnalysis":   ('Instructions - Ellipsometry',     'Ellipsometry'),
    "chemFormulations": ('Instructions - ChemFormulation',  'Chem Formulations'),
    "sldAnalysis":      ('Instructions - SLD',              'SLD'),
    "surfaceExcess":    ('Instructions - SurfaceExcess',    'Surface Excess'),
    }


#################
# Print Options #
#################

verbose      = False
very_verbose = False


################
# IAP Analysis #
################

## Functions
checkT0   = False  # ensures T0 = 0
shiftP    = False  # shift to P_min = 0
smoothIso = False  # reduces number of points for final plot
overplot  = False  # not yet functional

## Override Normal Keys
plotMultiPanel = False  # split input data between subplots
key = [(1,1)]  # [row0=(nFiles0,nFiles1),row1=(nFiles0,nFiles1)]

## Parameters
useCycles  = [0]  # list of cycles to be used in isotherm
nPoly      = 40  # order of the elasticity fit
nth        = 60  # data reduction; take every 'nth' point


#################
# SLD Analysis #
################

# multiplies average chain vol by 0.85 to model lipid chain compaction
compactChains = False

# add an injected lipid into the existing monolayer
addLipidToMonolayer      = True
injectedComponents       = ["Monolayer", "DLin-MC3-DMA"]
injectedRatios           = [90, 10]
updateMonolayerThickness = False
new_d1 = 1
new_d2 = 1

# add drug sample from third to second (headgroup) layer
addDrugToMonolayer = False
drugName           = "PolyA"
drugSize           = 20
threeSolv          = 91.374
SLD_drug_H2O       = 3.67 # polyA in H2O
SLD_drug_D2O       = 4.46 # polyA in D2O

# use vol frac (molar ratio = default); True is more accurate
useVolFrac = True

###################
# Surface Excess #
##################

# reduce number of points by len()/gammaNth
gammaNth = 1

# smooth function
gamma_nPoly = 10


#############
# Plotting #
############

## Save Options
saveAsPNG = True
saveAsPDF = True

# List of plot types that use the time axis
tAxisList = [" - pressure", " - area", " - normInjPressure", " - psi Time", " - delta Time", " - gammaL", " - gammaP"]

## Colours
colourDict = {

    # dark blue, light blue, dark orange, light orange
    "0": ['#1e81b0','#abdbe3', '#e28743','#eab676'],

    # dark blue, light blue, dark orange, light orange, dark green, light green
    "1": ['#1e81b0','#abdbe3', '#e28743','#eab676', '#32BE25', '#A3e19d'],

    # light blue, light orange, light green
    "2": ['#abdbe3', '#eab676', '#A3e19d'],

    # dark blue, dark orange, dark green, dark purple, dark red, dark yellow, persian pink, medium grey
    "3": ['#1e81b0', '#e28743', '#32BE25', '#6A0F8E', '#AB2330', '#FFCC00', '#F77FBE', '#71716F'],

    # 'Silver Lake' , 'Sea Serpent', 'Fuzzy Wuzzy', 'Cinnamon Satin'
    "4": ['#6083D0', '#60BBD0', '#D07560', '#D06083'],

    # 'green sheen', 'Turkish Rose'
    "5": ['#6FBBA6', '#BB6F84'],

    # 'Fuzzy Wuzzy', 'Cinnamon Satin'
    "6": ['#D07560', '#D06083'],

    # light blue, dark blue, light orange, dark orange, light green, dark green
    "7": ['#abdbe3', '#1e81b0', '#eab676', '#e28743', '#A3e19d', '#32BE25'],


    }
c = colourDict.get("3")

## Markers
markerDict = {

    # point
    "0": ['.','.','.','.'],

    # circles
    "1": ['o','o','o','o'],

    # squares
    "2": ['s','s','s','s'],

    # diamonds
    "3": ['D','D','D','D'],

    # cicrle-square repeat
    "4": ['o','s','o','s'],

    # cicrle-square-diamond repeat
    "5": ['o','s','d','o','s','d']

    }
markerType = markerDict.get("0")


##############
# Databases #
#############

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

# molecular weight lipid database, [g/mol]
lipidMw = {
    "DPPC":            734.039,
    "d-DPPC":          796.43,
    "POPC":            760.07,  # 760.076
    "d31-POPC":        791.07,  # 791.267
    "POPS":            783.99,
    "Cholesterol":     386.65,  # 386.654
    "d45-Cholesterol": 432,
    "DLin-KC2-DMA":    642.1,
    "DLin-MC3-DMA":    642.09,
    "d-DLin-MC3-DMA":  704.5,
    "DOPE":            744.034,
    "SM":              760.223,
    "LBPA":            792.07,
    "PolyA":           385.31,
    "DMG-PEG-2000":    2509.200,
    }

# chemical structures for each lipid: (struct_head, struct_tail)
lipidStruct = {
    "POPC":            ('N-O8-P-C11-H20', 'C31-H62'),
    "d31-POPC":        ('N-O8-P-C11-H20', 'C31-D31-H31'),
    "DOPE":            ('N-O8-P-C8-H14', 'C33-H64'),
    "SM":              ('N2-O5-P-C8-H19', 'O1-C33-H64'),
    "LBPA":            ('N-O4-P-C4-H11', 'O6-C38-H71'),
    "Cholesterol":     ('O-H','C27-H45'),
    "d45-Cholesterol": ('O-H','C27-D45'),
    "DLin-MC3-DMA":    ('N-O2-C7-H13', 'C36-H66'),
    "d-DLin-MC3-DMA":  ('N-O2-C7-H13', 'C36-H4-D62'),
    "DSPC":            ('N-O8-P-C11-H20','C33-H68'),
    "DMG-PEG-2000":    ('O5-C6-H7','C25-H52'), # polymer: ([O-C2-H4]_44 + O-C3-H7); total: O50-C122-H242
    "PolyA":           ('C10-H13-K-N5-O7-P','H'),
    }

# molecular volumes for each lipid (Angstroms cubed): (head, tail)
lipidMolVol = {
    "POPC":            (344,937),
    "d31-POPC":        (344,937),
    "DOPE":            (236,969),
    "SM":              (274,953),
    "LBPA":            (208,624),
    "Cholesterol":     (5,624),
    "d45-Cholesterol": (5,624),
    "DLin-MC3-DMA":    (260, 1030),
    "d-DLin-MC3-DMA":  (260, 1030),
    "DSPC":            (322,1000),
    "DMG-PEG-2000":    (256,767), # From Marianna: DMPE (head 0.25% total vol. 1023) PEG unit = 670
    "PolyA":           (1,1),
    }
