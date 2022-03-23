" Generalised plotting module for Pluto "

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys
import ast
import config

def plotParameters(suffix):

    # ask user to pick one of the analysisOptions
    analysisChoice = input("Would you like to override default plot parameters (y/n): ")

    if analysisChoice.upper() == 'Y':

        ## Plot Types, default line plot if false
        scatterPlot = (input("Plot as scatter plot (y/n)?: "))
        if scatterPlot.upper() == 'Y':
            plotWithScatter = True
            scatterSize     = int(input("Scatter size = "))
        else: plotWithScatter = False; scatterSize = 0

        lineMarkerPlot = (input("Plot as line with maker (y/n)?: "))
        if lineMarkerPlot.upper() == 'Y':
            plotLineWithMarker = True
            markEdgeWidth     = int(input("Marker edge width = "))
        else: plotLineWithMarker = False; markEdgeWidth = 0

        ## Marker Parameters
        lw = int(input("Line width = "))

        ## Font Size
        fs                  = int(input("Fontsize = "))
        legend_fs_reduction = int(input("Reduce fontsize by factor for legend = "))

    # default parameters
    else:
        plotWithScatter    = config.plotWithScatter
        plotLineWithMarker = config.plotLineWithMarker

        lw            = config.defaultLw
        scatterSize   = config.scatterSize
        markEdgeWidth = config.markEdgeWidth

        fs                  = config.fs
        legend_fs_reduction = config.legend_fs_reduction


    ## Axis Adjustment Parameters
    setX_AxInt  = -1
    setY_AxInt  = 0
    xAxisMinAdj = 0
    xAxisMaxAdj = 0
    yAxisMaxAdj = 1



    # ask user to pick one of the analysisOptions
    analysisChoice = input("Would you like to override default plot region (y/n): ")

    if analysisChoice.upper() == 'Y':

        # Override Parameters; config_n0, config_nf
        overrideNoP = True
        config_n0   = ast.literal_eval(input("Starting data point [n0, n0, ...]:\n  "))
        config_nf   = ast.literal_eval(input("Final data point [nf, nf, ...]:\n  "))

        # Override Parameters; config_xmin, config_xmax, config_ymin, config_ymax
        overrideAxisLim = True
        config_xmin     = ast.literal_eval(input("xmin = "))
        config_xmax     = ast.literal_eval(input("xmax = "))
        config_ymin     = ast.literal_eval(input("ymin = "))
        config_ymax     = ast.literal_eval(input("ymax = "))

        # Override Parameters; n_xticks, xTickInterval, yTickInterval
        overrideTickLocation = True
        n_xticks             = ast.literal_eval(input("Number of x-axis ticks = "))  # number of x axis ticks in time plots (s); [n-1]
        xTickInterval        = ast.literal_eval(input("xTick interval = "))  # x, y axis tick interval for P vs t plots; x is mins plot only
        yTickInterval        = ast.literal_eval(input("yTick interval = "))


    # default parameters, initialise with arbitrary numbers
    else:
        overrideNoP = config.overrideNoP
        config_n0   = config.config_n0
        config_nf   = config.config_nf

        overrideAxisLim = config.overrideAxisLim
        config_xmin     = config.config_xmin
        config_xmax     = config.config_xmax
        config_ymin     = config.config_ymin
        config_ymax     = config.config_ymax

        overrideTickLocation = config.overrideTickLocation
        n_xticks             = config.n_xticks
        xTickInterval        = config.xTickInterval
        yTickInterval        = config.yTickInterval


    return plotWithScatter, plotLineWithMarker, lw, scatterSize, markEdgeWidth,\
            fs, legend_fs_reduction, setX_AxInt, setY_AxInt, xAxisMinAdj,      \
            xAxisMaxAdj, yAxisMaxAdj, overrideNoP, config_n0, config_nf, overrideAxisLim,    \
            config_xmin, config_xmax, config_ymin, config_ymax, overrideTickLocation, n_xticks,            \
            xTickInterval, yTickInterval


def isolateFiles(count, key, suffix, row, col, X, Y, LABELS):

    # extract the number of total files to be plotted
    nFilesTotal = len(X[0])

    # separate number of files per subplot if multiplot
    if config.plotMultiPanel == True and suffix == " - isotherm":
        nFilesPerPlot = key[row][col]
    else:
        nFilesPerPlot = nFilesTotal

    # initialise output dicts
    x      = {new_list: [] for new_list in range(nFilesPerPlot)}
    y      = {new_list: [] for new_list in range(nFilesPerPlot)}
    labels = {new_list: [] for new_list in range(nFilesPerPlot)}

    # count is how far through the list of total files you are
    # iterate through all files
    if config.plotMultiPanel == True and suffix == " - isotherm":
        for i in range(count, nFilesTotal):

            if i == count:
                for j in range(nFilesPerPlot):
                    x[j]      = X[0].get(i+j)
                    y[j]      = Y.get(i+j)
                    labels[j] = LABELS.get(i+j)
                break
            break

        count += nFilesPerPlot


    elif suffix == " - elasticity" and col==0:
        x = X[0]
        y = Y; labels = LABELS
    elif suffix == " - elasticity" and col==1:
        x = X[1]
        y = Y; labels = LABELS


    else:
        x = X[0]; y = Y; labels = LABELS

    return count, nFilesPerPlot, x, y, labels



def plot(key, vars, suffix):

    plotWithScatter, plotLineWithMarker, lw, scatterSize, markEdgeWidth,       \
     fs, legend_fs_reduction, setX_AxInt, setY_AxInt, xAxisMinAdj, xAxisMaxAdj,\
     yAxisMaxAdj, overrideNoP, config_n0, config_nf, overrideAxisLim, config_xmin, config_xmax, config_ymin,      \
     config_ymax, overrideTickLocation, n_xticks, xTickInterval, yTickInterval = plotParameters(suffix)



    # unpack key into rows and columns for subplot
    if config.plotMultiPanel == True and suffix == " - isotherm":
        key = config.key
        nRow = len(key)
        nCol = len(key[0]) # assumes same num columns on both rows
    else:
        nRow, nCol = key


    # Create key x 1 subplot grid
    gs = gridspec.GridSpec(nRow, nCol)

    # initialise figure
    fig = plt.figure()
    ax  = plt.axes()




    # iterate along subplots, currently just row 0 column 1-2
    count = 0
    for row in range(nRow):
        for col in range(nCol):


            # unpack variables evertime to prevent overwriting within plot code
            N, equip, LABELS, axLabels, title, plotDIR, X, Y = vars

            # iterate through files and check number of subplots, isolate files accordingly
            count, nFilesPerPlot, x, y, labels = isolateFiles(count, key, suffix, row, col, X, Y, LABELS)
            N = nFilesPerPlot


            # initialise the subplot
            ax = plt.subplot(gs[row, col]) # row 0, col 0

            # set line width of each spine of given subplot
            for axis in ['top','bottom','left','right']:
                ax.spines[axis].set_linewidth(2)

            # initialise lists for axis range params
            min_x_vals = []; max_x_vals = []
            min_y_vals = []; max_y_vals = []


    ## Set region of interest

            # default region of interest (all values)
            n0 = [0 for i in range(N)]
            nf = []
            #print(x)
            # set upper limit to be shorter of two lists to ensure same length
            for i in range(N):
                if len(x.get(i)) == len(y.get(i)):
                    nf.append(len(x.get(i)))
                elif len(x.get(i)) > len(y.get(i)):
                    nf.append(len(y.get(i)))
                elif len(x.get(i)) < len(y.get(i)):
                    nf.append(len(x.get(i)))

            # alt. region of interest, NoP = number of points
            if overrideNoP == True:
                n0 = config_n0
                nf = config_nf


    ## Plot
            for i in range(N):

                # plots scatter plot with empty circles
                if plotWithScatter == True:
                    ax.scatter(x.get(i)[n0[i]:nf[i]], y.get(i)[n0[i]:nf[i]], label = labels.get(i), s=scatterSize, edgecolors=config.c[i], linewidth=lw, facecolors='none')

                # line plot with marker
                elif plotLineWithMarker == True:
                    ax.plot(x.get(i)[n0[i]:nf[i]], y.get(i)[n0[i]:nf[i]], label = labels.get(i), color=config.c[i], linewidth=lw, marker=config.markerType[i], markerfacecolor="None", markeredgewidth=markEdgeWidth)

                # default line plot
                else:
                    ax.plot(x.get(i)[n0[i]:nf[i]], y.get(i)[n0[i]:nf[i]], label = labels.get(i), color=config.c[i], linewidth=lw)


                # store minimum and maximum values for axis scales
                #min_x_vals.append( n0[i] ); max_x_vals.append( nf[i] )
                #min_y_vals.append( min(y.get(i)) ); max_y_vals.append( max(y.get(i)) )
                min_x_vals.append( min(x.get(i)) ); max_x_vals.append( max(x.get(i)) )
                min_y_vals.append( min(y.get(i)) ); max_y_vals.append( max(y.get(i)) )




    ## Set Axis ranges / limits

            ymin = 0
            if int(round(min(min_x_vals),-1)) + xAxisMinAdj >= 0:
                xmin = int(round(min(min_x_vals),-1)) + xAxisMinAdj
            else: xmin = 0

            xmax = int(round(max(max_x_vals),setX_AxInt)) + xAxisMaxAdj
            ymax = int(round(max(max_y_vals),setY_AxInt)) + yAxisMaxAdj


            # alt. region of interest
            if overrideAxisLim == True:
                xmin = config_xmin
                xmax = config_xmax
                ymin = config_ymin
                ymax = config_ymax

            ax.set_xlim([xmin,xmax])
            ax.set_ylim([ymin,ymax])



    ## Set tick locations plots

            # thresholds for different x axis scales
            if xmax >= 7200 and suffix in config.tAxisList:
                axLabels["x"] = "Time (hr)"

                # set axis ticks
                init_xticks = np.arange(xmin, xmax+1, step=(3600))
                ax.set_xticks(init_xticks)
                ax.set_yticks(np.arange(ymin, ymax+yAxisMaxAdj, step=yTickInterval))

                # overwrite tick numbers
                new_xticks = [i for i in range(0,int(xmax/3600)+1)]
                plt.xticks(init_xticks, new_xticks)


            elif xmax < 7200 and xmax > 600 and suffix in config.tAxisList:
                axLabels["x"] = "Time (min)"

                init_xticks = np.arange(xmin, xmax+1, step=600)
                ax.set_xticks(init_xticks)
                ax.set_yticks(np.arange(ymin, ymax+yAxisMaxAdj, step=yTickInterval))

                new_xticks = [i for i in range(0,int(round(xmax/60,-1))+xTickInterval,xTickInterval)]
                plt.xticks(init_xticks, new_xticks)


            elif xmax < 600 and suffix in config.tAxisList:
                axLabels["x"] = "Time (s)"
                ax.set_xticks(np.arange(xmin, xmax+1, step=int( round(((xmin+xmax)/n_xticks),-1) )))
                ax.set_yticks(np.arange(0, ymax+1, step=yTickInterval))


            elif overrideTickLocation == True:
                ax.set_xticks(np.arange(xmin, xmax+1, step=xTickInterval))
                ax.set_yticks(np.arange(ymin, ymax+yAxisMaxAdj, step=yTickInterval))

            else: pass

    ## Axis labels

            # axis labels; in for loop as iterates along number of subplots
            if config.plotMultiPanel == True and suffix == " - isotherm" and (nRow == 2 and nCol == 2):

                if row == 0 and col == 0:
                    plt.setp(ax.get_xticklabels(), visible=False)
                    plt.setp(ax.get_yticklabels(), visible=True)
                elif row == 0 and col == 1:
                    plt.setp(ax.get_xticklabels(), visible=False)
                    plt.setp(ax.get_yticklabels(), visible=False)
                elif row == 1 and col == 0:
                    plt.setp(ax.get_xticklabels(), visible=True)
                    plt.setp(ax.get_yticklabels(), visible=True)
                elif row == 1 and col == 1:
                    plt.setp(ax.get_xticklabels(), visible=True)
                    plt.setp(ax.get_yticklabels(), visible=False)

            elif config.plotMultiPanel == True and suffix == " - isotherm" and (nRow == 1 and nCol == 2):

                if row == 0 and col == 0:
                    plt.setp(ax.get_yticklabels(), visible=True)
                elif row == 0 and col == 1:
                    plt.setp(ax.get_yticklabels(), visible=False)

            else:
                if col == 0:
                    ax.set_xlabel(axLabels.get("x"), fontsize=fs, fontweight='bold')
                    ax.set_ylabel(axLabels.get("y"), fontsize=fs, fontweight='bold')
                elif col == 1:
                    ax.set_xlabel(axLabels.get("x1"), fontsize=fs-4, fontweight='bold')
                    plt.setp(ax.get_yticklabels(), visible=False)


            # legend; plot along with every figure unless elasticity - no property plotIsotherm...
            #if config.plotIsotherm == True and config.plotElasticity == True and col == 1:
            #    pass
                #ax.legend(prop={'size': fs-legend_fs_reduction, 'weight':'bold'}, frameon = False)
            #elif config.plotIsotherm == True and config.plotElasticity == True and col == 0:
            #    pass
            #else:
            #    ax.legend(prop={'size': fs-legend_fs_reduction, 'weight':'bold'}, frameon = False)
            ax.legend(prop={'size': fs-legend_fs_reduction, 'weight':'bold'}, frameon = False)


    ## Tick label size; legend; layout; show fig; save fig

    # set axis parameters, size etc.
            ax.tick_params(axis='x', labelsize=fs-1)
            ax.tick_params(axis='y', labelsize=fs-1)

            ax.xaxis.set_tick_params(which='major', size=5, width=2, direction='in', top='on')
            ax.yaxis.set_tick_params(which='major', size=5, width=2, direction='in', right='on')
            #ax.xaxis.set_tick_params(which='minor', size=7, width=2, direction='in', top='on')
            #ax.yaxis.set_tick_params(which='minor', size=7, width=2, direction='in', right='on')


    # merge axis of multipanel isotherm plots
    if config.plotMultiPanel == True and suffix == " - isotherm":
        fig.text(0.5, -0.03, axLabels.get("x"), ha='center', fontsize=fs, fontweight='bold')
        fig.text(-0.03, 0.5, axLabels.get("y"), va='center', rotation='vertical', fontsize=fs, fontweight='bold')


    # plot vertical line
    #plt.axvline(900, 0, 6, label='pyplot vertical line', c='r')


    # tight layout function
    plt.tight_layout()
    fig.subplots_adjust(wspace=0.05, hspace=0.05)

    # show plot
    print("\nFigure: %s%s" %(title, suffix))
    plt.show()

    # save the plot as a file
    if config.saveAsPNG == True:
        fig.savefig( plotDIR + '/' + title + suffix + '.png',
            format='png',
            dpi=400,
            bbox_inches='tight')

    if config.saveAsPDF == True:
        fig.savefig( plotDIR + '/' + title + suffix + '.pdf',
            format='pdf',
            dpi=400,
            bbox_inches='tight')

    return



def main(key, vars, suffix):

    # plot either single or dual style plot depending on input key
    # accepts dict structures only
    plot(key, vars, suffix)

    return



if __name__ == '__main__':
    print("Generating Plot...\n")
    main()
