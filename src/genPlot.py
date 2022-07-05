" Generalised plotting module for Pluto "

import sys, ast
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# import Pluto modules
import config

def plotParameters():

    if config.askUserPlotPar == True:

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


    if config.askUserPlotPar == True:

        # ask user to pick one of the analysisOptions
        analysisChoice = input("Would you like to override default plot region (y/n): ")

        if analysisChoice.upper() == 'Y':

            # Override Parameters; config_xmin, config_xmax, config_ymin, config_ymax
            overrideAxisLim = True
            config_xmin     = ast.literal_eval(input("xmin = "))
            config_xmax     = ast.literal_eval(input("xmax = "))
            config_ymin     = ast.literal_eval(input("ymin = "))
            config_ymax     = ast.literal_eval(input("ymax = "))

            # Override Parameters; n_xticks, xTickInterval, yTickInterval
            n_xticks             = ast.literal_eval(input("Number of x-axis ticks = "))  # number of x axis ticks in time plots (s); [n-1]
            xTickInterval        = ast.literal_eval(input("xTick interval = "))  # x, y axis tick interval for P vs t plots; x is mins plot only
            yTickInterval        = ast.literal_eval(input("yTick interval = "))


    # default parameters, initialise with arbitrary numbers
    else:
        overrideAxisLim = config.overrideAxisLim
        config_xmin     = config.config_xmin
        config_xmax     = config.config_xmax
        config_ymin     = config.config_ymin
        config_ymax     = config.config_ymax

        n_xticks             = config.n_xticks
        xTickInterval        = config.xTickInterval
        yTickInterval        = config.yTickInterval


    return plotWithScatter, plotLineWithMarker, lw, scatterSize, markEdgeWidth,\
            fs, legend_fs_reduction, setX_AxInt, setY_AxInt, xAxisMinAdj,\
            xAxisMaxAdj, yAxisMaxAdj, overrideAxisLim,\
            config_xmin, config_xmax, config_ymin, config_ymax, n_xticks,\
            xTickInterval, yTickInterval



def isolateFiles(count, key, suffix, row, col, X, Y, LABELS):

    # extract the number of total files to be plotted
    nFilesTotal = len(X[0])

    # separate number of files per subplot if multiplot
    if config.plotMultiPanel == True:
        nFilesPerPlot = key[row][col]
    else:
        nFilesPerPlot = nFilesTotal

    # initialise output dicts
    x      = {new_list: [] for new_list in range(nFilesPerPlot)}
    y      = {new_list: [] for new_list in range(nFilesPerPlot)}
    labels = {new_list: [] for new_list in range(nFilesPerPlot)}

    # count is how far through the list of total files you are
    # iterate through all files
    if config.plotMultiPanel == True:
        for i in range(count, nFilesTotal):

            if i == count:
                for j in range(nFilesPerPlot):
                    x[j]      = X[0].get(i+j)
                    y[j]      = Y.get(i+j)
                    labels[j] = LABELS.get(i+j)
                break
            break

        count += nFilesPerPlot

    # specific case for ellipsometry files where ax1=P and ax2=Am
    elif suffix == " - elasticity" and col==0:
        x = X[0]
        y = Y; labels = LABELS
    elif suffix == " - elasticity" and col==1:
        x = X[1]
        y = Y; labels = LABELS

    # base condition information parsed as usual
    else:
        x = X[0]; y = Y; labels = LABELS

    return count, nFilesPerPlot, x, y, labels



def plot(key, vars, suffix):

    # get plot parameters from user or config
    plotWithScatter, plotLineWithMarker, lw, scatterSize, markEdgeWidth,       \
     fs, legend_fs_reduction, setX_AxInt, setY_AxInt, xAxisMinAdj, xAxisMaxAdj,\
     yAxisMaxAdj, overrideAxisLim, config_xmin, config_xmax, config_ymin,      \
     config_ymax, n_xticks, xTickInterval, yTickInterval = plotParameters()


    # update and unpack key into rows and columns for subplot
    if config.plotMultiPanel == True:
        key = config.key
        nRow = len(key)
        nCol = len(key[0]) # assumes same num columns on both rows
    else:
        nRow, nCol = key


    # Create correpsonding subplot grid
    gs = gridspec.GridSpec(nRow, nCol)

    # initialise figure
    fig = plt.figure()
    ax  = plt.axes()

    # iterate along subplots
    count = 0
    for row in range(nRow):
        for col in range(nCol):

            # unpack variables evertime to prevent overwriting within plot code
            if config.plotAvRR == True:
                nFilesPerPlot, equip, LABELS, axLabels, title, plotDIR, X, Y, y_error = vars
            else:
                nFilesPerPlot, equip, LABELS, axLabels, title, plotDIR, X, Y = vars

            # iterate through files and check number of subplots, isolate files accordingly
            count, nFilesPerPlot, x, y, labels = isolateFiles(count, key, suffix, row, col, X, Y, LABELS)

            # initialise the subplot
            ax = plt.subplot(gs[row, col]) # row 0, col 0

            # set line width of each spine of given subplot
            for axis in ['top','bottom','left','right']:
                ax.spines[axis].set_linewidth(2)

            # initialise lists for axis range params
            min_x_vals = []; max_x_vals = []
            min_y_vals = []; max_y_vals = []


    ## Set region of interest (number of points to be plotted)

            # region of interest (all values)
            n0 = [0 for i in range(nFilesPerPlot)]
            nf = []

            # set upper limit to be shorter of two lists to ensure same length
            for i in range(nFilesPerPlot):
                if len(x.get(i)) == len(y.get(i)):
                    nf.append(len(x.get(i)))
                elif len(x.get(i)) > len(y.get(i)):
                    nf.append(len(y.get(i)))
                elif len(x.get(i)) < len(y.get(i)):
                    nf.append(len(x.get(i)))


    ## Plot
            for i in range(nFilesPerPlot):

                # plots scatter plot with empty circles
                if plotWithScatter == True and config.plotWithScatterError == False:
                    ax.scatter(x.get(i)[n0[i]:nf[i]], y.get(i)[n0[i]:nf[i]], label = labels.get(i), s=scatterSize, edgecolors=config.c[row][i], linewidth=lw, facecolors='none')

                if plotWithScatter == True and config.plotWithScatterError == True:
                    try:
                        ax.errorbar(x.get(i)[n0[i]:nf[i]], y.get(i)[n0[i]:nf[i]], yerr = y_error.get(i)[n0[i]:nf[i]], color=config.c[row][i], ms=scatterSize, mec=config.c[row][i], mfc=config.c[row][i], ecolor=config.c[row][i], label = labels.get(i), linewidth=lw) # , s=scatterSize, edgecolors=config.c[row][i], facecolors='none'
                    except ValueError:

                        try:
                            ax.errorbar(x.get(i)[n0[i]:nf[i]], y.get(i)[n0[i]:nf[i]], yerr = y_error, label = labels.get(i), s=scatterSize, edgecolors=config.c[row][i], linewidth=lw, facecolors='none')
                            print("\n\ngenPlot Error: y_error was empty, plotting without error bars.\n\n")
                        except ValueError:
                            print("\n\ngenPlot Error: Unknown ValueError. System clsoing.\n\n")
                            sys.exit()

                # line plot with marker
                elif plotLineWithMarker == True:
                    ax.plot(x.get(i)[n0[i]:nf[i]], y.get(i)[n0[i]:nf[i]], label = labels.get(i), color=config.c[row][i], linewidth=lw, marker=config.markerType[row][i], markerfacecolor="None", markeredgewidth=markEdgeWidth, markersize=config.markerSize)

                # default line plot
                else:
                    ax.plot(x.get(i)[n0[i]:nf[i]], y.get(i)[n0[i]:nf[i]], label = labels.get(i), color=config.c[row][i], linewidth=lw)


                # store minimum and maximum values of pre-processed data for axis scales
                min_x_vals.append( min(x.get(i)) ); max_x_vals.append( max(x.get(i)) )
                min_y_vals.append( min(y.get(i)) ); max_y_vals.append( max(y.get(i)) )




    ## Set Axis ranges / limits (the actual values of the points)

            # automatically set axis limits
            ymin = 0
            if int(round(min(min_x_vals),-1)) + xAxisMinAdj >= 0:
                xmin = int(round(min(min_x_vals),-1)) + xAxisMinAdj
            else: xmin = 0

            xmax = int(round(max(max_x_vals),setX_AxInt)) + xAxisMaxAdj
            ymax = int(round(max(max_y_vals),setY_AxInt)) + yAxisMaxAdj

            # override axis limits/region of interest via config/user
            if overrideAxisLim == True:
                xmin = config_xmin
                xmax = config_xmax
                ymin = config_ymin
                ymax = config_ymax

            ax.set_xlim([xmin,xmax])
            ax.set_ylim([ymin,ymax])



    ## Set tick locations plots

            # ticks correspond to values rather than number of points
            # thresholds for different x axis scales:

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

            else: pass



    ## Axis labels

            if config.overrideXAxisLabel == True:
                axLabels["x"] = config.xLabel

            if config.overrideYAxisLabel == True:
                axLabels["y"] = config.yLabel

            # axis labels; in for loop as iterates along number of subplots
            if config.plotMultiPanel == True and (nRow == 2 and nCol == 2):

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

            elif config.plotMultiPanel == True and (nRow == 1 and nCol == 2):

                if row == 0 and col == 0:
                    plt.setp(ax.get_yticklabels(), visible=True)
                elif row == 0 and col == 1:
                    plt.setp(ax.get_yticklabels(), visible=False)

            elif config.plotMultiPanel == True and (nRow == 1 and nCol == 3):

                if row == 0 and col == 0:
                    ax.set_ylabel(axLabels.get("y"), fontsize=fs-config.y0Axis_fs_reduction, fontweight='bold')
                    plt.setp(ax.get_yticklabels(), visible=True)
                elif row == 0 and col == 1:
                    ax.set_xlabel(axLabels.get("x"), fontsize=fs-config.x0Axis_fs_reduction, fontweight='bold')
                    plt.setp(ax.get_yticklabels(), visible=False)
                elif row == 0 and col == 2:
                    plt.setp(ax.get_yticklabels(), visible=False)

            elif config.plotMultiPanel == True and (nRow == 3 and nCol == 1):

                if row == 0 and col == 0:
                    plt.setp(ax.get_xticklabels(), visible=False)
                elif row == 1 and col == 0:
                    ax.set_ylabel(axLabels.get("y"), fontsize=fs-config.y0Axis_fs_reduction, fontweight='bold')
                    plt.setp(ax.get_xticklabels(), visible=False)
                elif row == 2 and col == 0:
                    ax.set_xlabel(axLabels.get("x"), fontsize=fs-config.x0Axis_fs_reduction, fontweight='bold')
                    plt.setp(ax.get_xticklabels(), visible=True)

            else:
                if col == 0:
                    ax.set_xlabel(axLabels.get("x"), fontsize=fs-config.x0Axis_fs_reduction, fontweight='bold')
                    ax.set_ylabel(axLabels.get("y"), fontsize=fs-config.y0Axis_fs_reduction, fontweight='bold')
                elif col == 1:
                    ax.set_xlabel(axLabels.get("x1"), fontsize=fs-config.x1Axis_fs_reduction, fontweight='bold')
                    plt.setp(ax.get_yticklabels(), visible=False)


    ## Legend; plot on every subplot
            if config.legendOn == True:
                ax.legend(prop={'size': fs-legend_fs_reduction, 'weight':'bold'}, frameon = False, loc=config.legendLoc)


    ## Tick label size

            # tick label size
            ax.tick_params(axis='x', labelsize=fs-config.tick_fs_reduction)
            ax.tick_params(axis='y', labelsize=fs-config.tick_fs_reduction)

            # tick line size and width
            ax.xaxis.set_tick_params(which='major', size=config.majorTickSize, width=config.majorTickWidth, direction='in', top='on')
            ax.yaxis.set_tick_params(which='major', size=config.majorTickSize, width=config.majorTickWidth, direction='in', right='on')

            if config.showMinorTicks == True:
                ax.xaxis.set_tick_params(which='minor', size=config.minorTickSize, width=config.minorTickWidth, direction='in', top='on')
                ax.yaxis.set_tick_params(which='minor', size=config.minorTickSize, width=config.minorTickWidth, direction='in', right='on')


    ## Other

    # merge axis of multipanel isotherm plots
    if config.plotMultiPanel == True and suffix == " - isotherm":
        fig.text(0.5, -0.03, axLabels.get("x"), ha='center', fontsize=fs, fontweight='bold')
        fig.text(-0.03, 0.5, axLabels.get("y"), va='center', rotation='vertical', fontsize=fs, fontweight='bold')


    # plot vertical line
    if config.plotVerticalLine == True:
        plt.axvline(config.x0Line, config.y0Line, config.y1Line, label='pyplot vertical line', c='r')


    if config.plotxLog10 == True:
        ax.set_xscale('log')
    if config.plotyLog10 == True:
        ax.set_yscale('log')

    # tight layout function
    plt.tight_layout()
    fig.subplots_adjust(wspace=0.05, hspace=0.05)

    # show plot
    print("\nFigure generated: %s%s." %(title, suffix))
    if config.showFig == True:
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
    return plot(key, vars, suffix)



if __name__ == '__main__':
    main()
