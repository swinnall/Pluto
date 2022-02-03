" Module that defines global variables for Pluto"


################
# IAP Analysis #
################

## Functions
checkT0   = False                     # ensures T0 = 0
shiftP    = False                     # shift to P_min = 0
smoothIso = False                     # reduces number of points for final plot
overplot  = False                    # not yet functional

## Override Normal Keys
plotMultiPanel = False              # split input data between subplots
key = [(1,1)]                 # [row0=(nFiles0,nFiles1),row1=(nFiles0,nFiles1)]

## Parameters
useCycles  = [0,1,2,3,4]                # list of cycles to be used in isotherm
nPoly      = 40                      # order of the elasticity fit
nth        = 60                      # data reduction; take every 'nth' point


##########################
# Ellipsometry Analysis #
########################




#################
# SLD Analysis #
################

# multiplies average chain vol by 0.85 to model lipid chain compaction
compactChains = False

# for mixing drug uptake to second and third layers
useL2_L3 = False

# use vol frac (molar ratio = default); True is more accurate
useVolFrac = True

# include an injected lipid into the system
addLipid = False


###################
# Surface Excess #
##################

# reduce number of points by len()/gammaNth
gammaNth = 1

# smooth function
gamma_nPoly = 10

# overplot L and P; to be implemented
overPlotGammaLP = False


#############
# Printing #
############

# print extra info
verbose      = True
very_verbose = False


#############
# Plotting #
############

# List of plot types that use the time axis
tAxisList = [" - pressure", " - area", " - normInjPressure", " - psi Time", " - delta Time"]

## Save Options
saveAsPNG = True
saveAsPDF = True

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
    "POPC":            ('C10-H18-N-O8-P','C32-H64'),
    "d31-POPC":        ('C10-H18-N-O8-P','C32-D31-H33'),
    "DOPE":            ('N-O8-P','C41-H78'),
    "SM":              ('N2-O6-P','C47-H93'),
    "LBPA":            ('N-O10-P','C42-H82'),
    "Cholesterol":     ('O-H','C27-H45'),
    "d45-Cholesterol": ('O-H','C27-D45'),
    "DLin-MC3-DMA":    ('N-O2-C6-H12','C37-H67'),
    "d-DLin-MC3-DMA":  ('N-O2-C6-H12','C37-H5-D62'),
    "DSPC":            ('N-O8-P','C44-H88'),
    "DMG-PEG-2000":    ('O50','C122-H242'),
    "PolyA":           ('C10-H13-K-N5-O7-P','H')
    }

## Modelling MC3 Chol Monolayer (full | split) ##

## C2H4 MC3 Model (full | split)
#  h-MC3:  ('N-O2-C7-H13', 'C36-H66')     | ('N-O2-C7-H13', 'C2-H4', 'C34-H62')
#  d-MC3:  ('N-O2-C7-H13', 'C36-H4-D62')  | ('N-O2-C7-H13', 'C2-H4', 'C34-D62')

## C3H5 MC3 Model (full | split)
#  h-MC3:  ('N-O2-C6-H12', 'C37-H67')     | ('N-O2-C6-H12', 'C3-H5', 'C34-H62')
#  d-MC3:  ('N-O2-C6-H12', 'C37-H5-D62')  | ('N-O2-C6-H12', 'C3-H5', 'C34-D62')

# C2H3 Cholesterol Model
#  h-chol:  ('O-H', 'C27-H45')   | ('O-H', 'C2-H3', 'C25-H42')
#  d-chol:  ('O-H', 'C27-D45')   | ('O-H', 'C2-D3', 'C25-D42')



# molecular volumes for each lipid (Angstroms cubed): (head, tail)
lipidMolVol = {
    "POPC":            (344,937),
    "d31-POPC":        (344,937),
    "DOPE":            (236,969),
    "SM":              (274,953),
    "LBPA":            (208,624),
    "Cholesterol":     (5,624),
    "d45-Cholesterol": (5,624),
    "DLin-MC3-DMA":    (260,1030), # (290,1000)
    "d-DLin-MC3-DMA":  (260,1030),
    "DSPC":            (322,1000),
    "DMG-PEG-2000":    (200,470),
    "PolyA":           (1,1)
    }

# dict of input instructions names and output file names
pathNames = {
    "isoAnalysis":      ('Instructions - Isotherm','Isotherm'),
    "ellipsAnalysis":   ('Instructions - Ellipsometry','Ellipsometry'),
    "chemFormulations": ('Instructions - ChemFormulation','Chem Formulations'),
    "sldAnalysis":      ('Instructions - SLD','SLD'),
    "SurfaceExcess":    ('Instructions - SurfaceExcess', 'Surface Excess'),
    }




#
