" Module that defines variables for Pluto"

# Keep program running after first analysis
moreAnalysis = False

##############
# File Paths #
##############

inputDir  = "../../UoM-Data-Repository/input/"
outputDir = "../../UoM-Data-Repository/output"

# dict of input instructions names and output file names
pathNames = {
    "surfacePressure": ('Instructions - Surface Pressure', 'Surface Pressure'),
    "Ellipsometry":    ('Instructions - Ellipsometry',     'Ellipsometry'),
    "surfaceExcess":   ('Instructions - SurfaceExcess',    'Surface Excess'),
    }


############################
# Surface Pressure Module #
###########################

checkT0   = False  # ensures T0 = 0
shiftP    = False  # shift to P_min = 0
useCycles = [0]  # list of cycles to be used in isotherm; e.g. [0, 1, 3, 4]

###############
# Smooth Data #
###############

smoothByPoints = False # for reducing the number of points in the dataset
nth = 10 # keep every nth point

smoothByPolyFit = True # turns data cleaning on/off
nPoly = 15 # the order of the fitted polynomial (no params produced, just gives trend)
lowerLimit = 0.01 # filters out any anomalous produced datapoints during polynomial fitting


#############
# Plotting #
############

## Ask for user input plot parameters while program is running
askUserPlotPar = False

## Show options
showFig = True

## Save Options
saveAsPNG = True
saveAsPDF = True

## Multi panel options
plotMultiPanel = False   # split input data between subplots
key = [[4],[4],[4]]  # [row0=[nFiles0,nFiles1],row1=[nFiles0,nFiles1]]

## Vertical line option
plotVerticalLine = False
x0Line = 100
y0Line = 0
y1Line = 10

## Set default plot parameters
defaultLw = 2

plotWithScatter = False
scatterSize     = 50

plotLineWithMarker = False
markerSize         = 7
markEdgeWidth      = 4

## Font size and reduction values
fs                   = 24
legend_fs_reduction  = 10
x0Axis_fs_reduction  = 0
x1Axis_fs_reduction  = 0
y0Axis_fs_reduction  = 0
tick_fs_reduction    = 0

## Axis range (False=Automatic)
overrideAxisLim = True
config_xmin     = 0
config_xmax     = 7200 # 1200 (20 min) # 7200 (2 hr) # 28800 (8hr)
config_ymin     = 1
config_ymax     = 4.0

## Force tick locations
overrideTickLocation = False
n_xticks             = 6
xTickInterval        = 10
yTickInterval        = 0.5

## Tick sizes and major/minor
majorTickSize  = 5
majorTickWidth = 2
showMinorTicks = False
minorTickSize  = 0
minorTickWidth = 0

overrideXAxisLabel = False
xLabel = "Time (min)"

overrideYAxisLabel = False
yLabel = "NA"

legendOn  = True
legendLoc = 'best' # default = 'best'; lower right''

## List of plot types that use the time axis
tAxisList = [" - pressure", " - area", " - psi Time", " - delta Time", " - gammaL", " - gammaP"]

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

    # EEM Figures: aqua, azure, blue
    "8": ['#00FFFF', '#0080FF', '#0000FF'],

    # LEM Figures: orange, red
    "9": ['#FF7F00', '#FF0000'],

    # MC3 Surface Excess Figures: light blue, purple-blue, light orange, dark green
    "10": [["#3399FF", "#3333FF", "#FF9933", "#FF3333", "#32BE25"]],

	# MC3 PBS Structural Figures: blue; light -> dark
	"11": ["#CCE6FF", "#99CCFF", "#66B3FF", "#3399FF"],

	# MC3:chol PBS Structural Figures: purple; light -> dark
	"12": ["#CCCCFF", "#9999FF", "#6666FF", "#3333FF"],

	# MC3 Citrate Structural Figures: orange; light -> dark
	"13": ["#FFE5CC", "#FFCC99", "#FFB266", "#FF9933"],

    # three panel combination of 11-13; must be list of lists where each sublist is a subplot
	"14": [["#CCE6FF", "#99CCFF", "#66B3FF", "#3399FF"], ["#CCCCFF", "#9999FF", "#6666FF", "#3333FF"], ["#FFE5CC", "#FFCC99", "#FFB266", "#FF9933"]],
    }
c = colourDict.get("10")

## Markers
markerDict = {

    # point
    "0": ['.','.','.','.'],

    # circles
    "1": [['o','o','o','o']],

    # squares
    "2": ['s','s','s','s'],

    # diamonds
    "3": ['D','D','D','D'],

    # cicrle-square repeat
    "4": ['o','s','o','s'],

    # cicrle-square-diamond repeat
    "5": ['o','s','d','o','s','d'],

    # Neutron symbols
    "6": [['s', 'o', '^', 'd'],['s', 'o', '^', 'd'],['s', 'o', '^', 'd']],

    }
markerType = markerDict.get("1")
