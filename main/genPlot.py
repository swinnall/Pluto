" Generalised plotting module for Pluto "

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys
import config

def plot(key, vars):

    # Create key x 1 subplot grid
    gs = gridspec.GridSpec(1, key)

    # initialise figure
    fig = plt.figure()
    ax  = plt.axes()

    # fontsize
    fs = config.fs

    # unpack variables
    N, equip, labels, axLabels, suffix, title, plotDIR, x, y = vars


    # iterate along subplots, currently just row 0 column 1-2
    for k in range(key):

        # initialise the subplot
        ax = plt.subplot(gs[0, k]) # row 0, col 0

        # initialise lists for axis range params
        min_x_vals = []; max_x_vals = []
        min_y_vals = []; max_y_vals = []


    ## Set region of interest

        # default region of interest (all values)
        n0 = [0 for i in range(N)]
        nf = []

        # set upper limit to be shorter of two lists to ensure same length
        for i in range(N):
            if len(x[k].get(i)) == len(y.get(i)):
                nf.append(len(x[k].get(i)))
            elif len(x[k].get(i)) > len(y.get(i)):
                nf.append(len(y.get(i)))
            elif len(x[k].get(i)) < len(y.get(i)):
                nf.append(len(x[k].get(i)))

        # alt. region of interest [number of points]; uncomment to override default
        n0 = [0, 0, 0, 0]
        nf = [16, 16, 16, 16]


    ## Plot

        for i in range(N):

            if suffix in suffix in [" - normInjPressure"]:
                ax.scatter(x[k].get(i)[n0[i]:nf[i]], y.get(i)[n0[i]:nf[i]], label = labels.get(i), edgecolors=config.c[i], linewidth=config.lw, facecolors='none')

            else:
                ax.plot(x[k].get(i)[n0[i]:nf[i]], y.get(i)[n0[i]:nf[i]], label = labels.get(i), color=config.c[i], linewidth=config.lw)

            # store minimum and maximum values for axis scales
            min_x_vals.append( min(x[k].get(i)) ); max_x_vals.append( max(x[k].get(i)) )   #min_x_vals.append( n0[i] ); max_x_vals.append( nf[i] )
            min_y_vals.append( min(y.get(i)) ); max_y_vals.append( max(y.get(i)) )         #min_y_vals.append( min(y.get(i)) ); max_y_vals.append( max(y.get(i)) )


    ## Set Axis ranges

        # set min values
        ymin = config.yAxisMin
        if int(round(min(min_x_vals),-1))-5 >= 0:
            xmin = int(round(min(min_x_vals),-1)) + config.xAxisMinAdj
        else: xmin = config.xAxisMin

        # set max values
        xmax = int(round(max(max_x_vals),-1)) + config.xAxisMaxAdj
        ymax = int(round(max(max_y_vals),-1)) + config.yAxisMaxAdj

        ax.set_xlim([xmin,xmax])
        ax.set_ylim([ymin,ymax])



    ## Determine axis ticks for x-axis = time plots

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



    ## Axis labels; tick label size; legend; layout; show fig; save fig

        # axis labels; in for loop as iterates along number of subplots
        if k == 0:
            ax.set_xlabel(axLabels.get("x"), fontsize=fs, fontweight='bold')
            ax.set_ylabel(axLabels.get("y"), fontsize=fs, fontweight='bold')
        elif k == 1:
            ax.set_xlabel(axLabels.get("x1"), fontsize=fs-4, fontweight='bold')
            #ax.set_ylabel(axLabels.get("y"), fontsize=fs, fontweight='bold')
            plt.setp(ax.get_yticklabels(), visible=False)

    # set axis parameters, size etc.
    ax.tick_params(axis='x', labelsize=fs-1)
    ax.tick_params(axis='y', labelsize=fs-1)

    # legend
    ax.legend(prop={'size': fs, 'weight':'bold'}, frameon = False)

    # tight layout function
    plt.tight_layout()

    # show plot
    print("\nFigure: %s%s" %(title, suffix))
    plt.show()

    # save the plot as a file
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
