" Module that defines global variables for Pluto"


################
# IAP Analysis #
################

## Plot Isotherms
plotIsotherm     = False            # plot standard isotherm
plotCompressions = False            # plot isotherm compressions only
plotExpansions   = False            # plot isotherm expansions only
plotCycles       = False            # colour code cycles within isotherm

## Plot Y vs Time
plotPressure      = True           # plot pressure
plotNormInjection = False           # plot normalised pressure
plotArea          = False           # plot percentage area

## Elasticity Analysis
plotElasticity = False              # plots elasticity analyis

## Other Functions
checkT0   = False                     # ensures T0 = 0
shiftP    = False                     # shift to P_min = 0
smoothIso = False                     # reduces number of points for final plot
overplot  = False                    # not yet functional

## Override Normal Keys
plotMultiPanel = False              # split input data between subplots
key = [(1,1)]                 # [row0=(nFiles0,nFiles1),row1=(nFiles0,nFiles1)]

## Parameters
useCycles  = [0,1,2,3,4]                     # list of cycles to be used in isotherm
nPoly      = 40                      # order of the elasticity fit
nth        = 60                      # data reduction; take every 'nth' point


##########################
# Ellipsometry Analysis #
########################

# plot psi and/or delta against angle of incidence
plotAOIpsi   = False
plotAOIdelta = False

# plot psi and/or delta against time
plotTimePsi   = False
plotTimeDelta = True


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

# select which quantities to plot
plotGammaL = True
plotGammaP = True

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

analysisOptions = ['isoAnalysis','ellipsAnalysis','chemFormulations',          \
                    'sldAnalysis','SurfaceExcess']

analysisType = analysisOptions[0]

def plotParameters(analysisOptions, analysisType):

    if analysisType == analysisOptions[0]:

        ## Plot Types, default line plot if false
        plotWithScatter    = False
        plotLineWithMarker = False

        ## Marker Parameters
        lw            = 3
        scatterSize   = 50
        markEdgeWidth = 1

        ## Font Size
        fs                  = 14
        legend_fs_reduction = 5

        ## Axis Adjustment Parameters
        setX_AxInt  = -1
        setY_AxInt  = 0
        xAxisMinAdj = 0
        xAxisMaxAdj = 0
        yAxisMaxAdj = 0.0000001

        ## Manual Override Parameters
        overrideNoP = False
        n0 = [0, 2400, 0, 0]
        nf = [600, 17000, 16, 16]

        overrideAxisLim = True
        xmin = 0
        xmax = 9000
        ymin = 0
        ymax = 30

        overrideTickLocation = True
        n_xticks      = 3               # number of x axis ticks in time plots (s); [n-1]
        xTickInterval = 10              # x, y axis tick interval for P vs t plots; x is mins plot only
        yTickInterval = 5




    elif analysisType == analysisOptions[1]:

        ## Plot Types, default line plot if false
        plotWithScatter    = False
        plotLineWithMarker = False

        ## Marker Parameters
        lw            = 3
        scatterSize   = 50
        markEdgeWidth = 1

        ## Font Size
        fs                  = 14
        legend_fs_reduction = 5

        ## Axis Adjustment Parameters
        setX_AxInt  = -1
        setY_AxInt  = 0
        xAxisMinAdj = 0
        xAxisMaxAdj = 0
        yAxisMaxAdj = 0.0000001

        ## Manual Override Parameters
        overrideNoP = False
        n0 = [0, 2400, 0, 0]
        nf = [600, 17000, 16, 16]

        overrideAxisLim = True
        xmin = 0
        xmax = 9000
        ymin = 0
        ymax = 6

        overrideTickLocation = True
        n_xticks      = 3               # number of x axis ticks in time plots (s); [n-1]
        xTickInterval = 10              # x, y axis tick interval for P vs t plots; x is mins plot only
        yTickInterval = 1

    elif analysisType == analysisOptions[4]:

        ## Plot Types, default line plot if false
        plotWithScatter    = False
        plotLineWithMarker = False

        ## Marker Parameters
        lw            = 3
        scatterSize   = 50
        markEdgeWidth = 1

        ## Font Size
        fs                  = 14
        legend_fs_reduction = 5

        ## Axis Adjustment Parameters
        setX_AxInt  = -1
        setY_AxInt  = 0
        xAxisMinAdj = 0
        xAxisMaxAdj = 0
        yAxisMaxAdj = 0.0000001

        ## Manual Override Parameters
        overrideNoP = False
        n0 = [0, 2400, 0, 0]
        nf = [600, 17000, 16, 16]

        overrideAxisLim = True
        xmin = 0
        xmax = 9000
        ymin = 0
        ymax = 6

        overrideTickLocation = True
        n_xticks      = 3               # number of x axis ticks in time plots (s); [n-1]
        xTickInterval = 10              # x, y axis tick interval for P vs t plots; x is mins plot only
        yTickInterval = 1

    return plotWithScatter, plotLineWithMarker, lw, scatterSize, markEdgeWidth,\
            fs, legend_fs_reduction, setX_AxInt, setY_AxInt, xAxisMinAdj,      \
            xAxisMaxAdj, yAxisMaxAdj, overrideNoP, n0, nf, overrideAxisLim,    \
            xmin, xmax, ymin, ymax, overrideTickLocation, n_xticks,            \
            xTickInterval, yTickInterval


plotWithScatter, plotLineWithMarker, lw, scatterSize, markEdgeWidth,           \
 fs, legend_fs_reduction, setX_AxInt, setY_AxInt, xAxisMinAdj, xAxisMaxAdj,    \
 yAxisMaxAdj, overrideNoP, n0, nf, overrideAxisLim, xmin, xmax, ymin,          \
 ymax, overrideTickLocation, n_xticks, xTickInterval, yTickInterval = plotParameters(analysisOptions, analysisType)


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
    "DLin-MC3-DMA": ('N-O2-C8-H13','C37-H66'),
    "d-DLin-MC3-DMA": ('N-O2-C8-H13','C37-H4-D62'),
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
    "DLin-MC3-DMA": (260,1030), # (290,1000)
    "d-DLin-MC3-DMA": (260,1030),
    "DSPC": (322,1000),
    "DMG-PEG-2000": (200,470),
    "PolyA": (1,1)
    }

# dict of input instructions names and output file names
pathNames = {
    "isoAnalysis": ('Instructions - Isotherm','Isotherm'),
    "ellipsAnalysis": ('Instructions - Ellipsometry','Ellipsometry'),
    "chemFormulations": ('Instructions - ChemFormulation','Chem Formulations'),
    "sldAnalysis": ('Instructions - SLD','SLD'),
    "SurfaceExcess": ('Instructions - SurfaceExcess', 'Surface Excess'),
    }
