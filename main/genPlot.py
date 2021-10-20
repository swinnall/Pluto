" Generalised plotting module for Pluto "

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys
import config

def isolateFiles(count, row, col, X, Y):

    nFilesPerPlot = config.key[row][col]
    nFilesTotal   = len(X[0])

    x = {new_list: [] for new_list in range(nFilesPerPlot)}
    y = {new_list: [] for new_list in range(nFilesPerPlot)}



    # count is how far through the list of total files you are
    # iterate through all files
    for i in range(count, nFilesTotal):
        #print(i)
        if i == count:
            for j in range(nFilesPerPlot):
                #print('SAVED')
                x[j] = X[0].get(i+j)
                y[j] = Y.get(i+j)
            break
        break


    count += nFilesPerPlot
    #print(x)
    return count, nFilesPerPlot, x, y



def plot(key, vars):

    ## Add in axis label options for 2x2 case


    # unpack key into rows and columns for subplot
    nRow, nCol = key

    # Create key x 1 subplot grid
    gs = gridspec.GridSpec(nRow, nCol)

    # initialise figure
    fig = plt.figure()
    ax  = plt.axes()

    # fontsize
    fs = config.fs




    # iterate along subplots, currently just row 0 column 1-2
    count = 0
    for row in range(nRow):
        for col in range(nCol):

            # unpack variables evertime to prevent overwriting within plot code
            N, equip, labels, axLabels, suffix, title, plotDIR, X, Y = vars


            ## need: for a given subplot (row, col): iterate along the given numbr of files specified
            if config.plotSpecialIsotherm == True:
                count, nFilesPerPlot, x, y = isolateFiles(count, row, col, X, Y)
                N = nFilesPerPlot
            else:
                x = X; y = Y



            # initialise the subplot
            ax = plt.subplot(gs[row, col]) # row 0, col 0

            # initialise lists for axis range params
            min_x_vals = []; max_x_vals = []
            min_y_vals = []; max_y_vals = []


    ## Set region of interest

            # default region of interest (all values)
            n0 = [0 for i in range(N)]
            nf = []

            # set upper limit to be shorter of two lists to ensure same length
            for i in range(N):
                if len(x.get(i)) == len(y.get(i)):
                    nf.append(len(x.get(i)))
                elif len(x.get(i)) > len(y.get(i)):
                    nf.append(len(y.get(i)))
                elif len(x.get(i)) < len(y.get(i)):
                    nf.append(len(x.get(i)))

            # alt. region of interest
            if config.overrideNoP == True:
                n0 = config.n0
                nf = config.nf


    ## Plot
            for i in range(N):

                # plots scatter plot with empty circles
                if suffix in config.scatterSuffixList:
                    ax.scatter(x.get(i)[n0[i]:nf[i]], y.get(i)[n0[i]:nf[i]], label = labels.get(i), s=config.scatterSize, edgecolors=config.c[i], linewidth=config.lw, facecolors='none')

                # line plot with marker
                elif suffix not in config.scatterSuffixList and config.plotWithMarker == True:
                    ax.plot(x.get(i)[n0[i]:nf[i]], y.get(i)[n0[i]:nf[i]], label = labels.get(i), color=config.c[i], linewidth=config.lw, marker=config.markerType[i], markerfacecolor="None", markeredgewidth=config.markEdgeWidth)

                # default line plot
                else:
                    ax.plot(x.get(i)[n0[i]:nf[i]], y.get(i)[n0[i]:nf[i]], label = labels.get(i), color=config.c[i], linewidth=config.lw)


                # store minimum and maximum values for axis scales
                #min_x_vals.append( n0[i] ); max_x_vals.append( nf[i] )
                #min_y_vals.append( min(y.get(i)) ); max_y_vals.append( max(y.get(i)) )
                min_x_vals.append( min(x.get(i)) ); max_x_vals.append( max(x.get(i)) )
                min_y_vals.append( min(y.get(i)) ); max_y_vals.append( max(y.get(i)) )



    ## Set Axis ranges / limits

            # set min values
            ymin = 0
            if int(round(min(min_x_vals),-1)) + config.xAxisMinAdj >= 0:
                xmin = int(round(min(min_x_vals),-1)) + config.xAxisMinAdj
            else: xmin = 0

            # set max values
            xmax = int(round(max(max_x_vals),-1)) + config.xAxisMaxAdj
            ymax = int(round(max(max_y_vals),-1)) + config.yAxisMaxAdj

            # alt. region of interest
            if config.overrideAxisLim == True:
                xmin = config.xmin
                xmax = config.xmax
                ymin = config.ymin
                ymax = config.ymax

            ax.set_xlim([xmin,xmax])
            ax.set_ylim([ymin,ymax])



    ## Set tick locations plots

            # thresholds for different x axis scales
            if xmax >= 7200 and suffix in [" - pressure", " - area", " - normInjPressure"]:
                axLabels["x"] = "Time (hr)"

                # set axis ticks
                init_xticks = np.arange(xmin, xmax+1, step=(3600))
                ax.set_xticks(init_xticks)
                ax.set_yticks(np.arange(0, ymax+1, step=config.yTickInterval))

                # overwrite tick numbers
                new_xticks = [i for i in range(0,int(xmax/3600)+1)]
                plt.xticks(init_xticks, new_xticks)


            elif xmax < 7200 and xmax > 600 and suffix in [" - pressure", " - area", " - normInjPressure"]:
                axLabels["x"] = "Time (min)"

                init_xticks = np.arange(xmin, xmax+1, step=600)
                ax.set_xticks(init_xticks)
                ax.set_yticks(np.arange(0, ymax+1, step=config.yTickInterval))

                new_xticks = [i for i in range(0,int(round(xmax/60,-1))+config.xTickInterval,config.xTickInterval)]
                plt.xticks(init_xticks, new_xticks)


            elif xmax < 600 and suffix in [" - pressure", " - area", " - normInjPressure"]:
                axLabels["x"] = "Time (s)"
                ax.set_xticks(np.arange(xmin, xmax+1, step=int( round(((xmin+xmax)/config.n_xticks),-1) )))
                ax.set_yticks(np.arange(0, ymax+1, step=config.yTickInterval))


            elif suffix in [" - isotherm"] and config.overrideTickLocation == True:
                ax.set_xticks(np.arange(xmin, xmax+1, step=config.xTickInterval))
                ax.set_yticks(np.arange(ymin, ymax+1, step=config.yTickInterval))



    ## Axis labels

            # axis labels; in for loop as iterates along number of subplots
            if config.plotSpecialIsotherm == True:

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

            else:
                if col == 0:
                    ax.set_xlabel(axLabels.get("x"), fontsize=fs, fontweight='bold')
                    ax.set_ylabel(axLabels.get("y"), fontsize=fs, fontweight='bold')
                elif col == 1:
                    ax.set_xlabel(axLabels.get("x1"), fontsize=fs-4, fontweight='bold')
                    plt.setp(ax.get_yticklabels(), visible=False)


            # legend; must plot along with every figure
            ax.legend(prop={'size': fs-5, 'weight':'bold'}, frameon = False)

    ## Tick label size; legend; layout; show fig; save fig

    # set axis parameters, size etc.
            ax.tick_params(axis='x', labelsize=fs-1)
            ax.tick_params(axis='y', labelsize=fs-1)



    #
    if config.plotSpecialIsotherm == True:
        fig.text(0.5, -0.05, axLabels.get("x"), ha='center', fontsize=fs, fontweight='bold')
        fig.text(-0.05, 0.5, axLabels.get("y"), va='center', rotation='vertical', fontsize=fs, fontweight='bold')


    # tight layout function
    plt.tight_layout()
    fig.subplots_adjust(wspace=0, hspace=0)

    # show plot
    print("\nFigure: %s%s" %(title, suffix))
    plt.show()

    # save the plot as a file
    fig.savefig( plotDIR + '/' + title + suffix + '.png',
        format='png',
        dpi=400,
        bbox_inches='tight')

    fig.savefig( plotDIR + '/' + title + suffix + '.pdf',
        format='pdf',
        dpi=400,
        bbox_inches='tight')

    return



def main(key, vars):

    # plot either single or dual style plot depending on input key
    # accepts dict structures only
    plot(key, vars)

    return



if __name__ == '__main__':
    print("Generating Plot...\n")
    main()
