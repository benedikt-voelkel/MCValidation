#!/usr/bin/env python

import sys
# Use argparse to work with command line args
import argparse
# To do some os related things like joining file paths
import os
# JSON to read config files encoded in JSON format.
import json
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

'''
Some global variables:
    to handle JSON files and corresponding plotting dictionaries:
        XDATA_NAME:
        YDATA_NAME:
'''
XDATA_NAME = "xData"
XDATA_ERROR_NAME = "xDataError"
YDATA_NAME = "yData"
YDATA_ERROR_NAME = "yDataError"
FILES = "files"
# Available extraction algorithms are 'identity' and 'mean'
EXTRACTION_ALGORITHM = "extractionAlgorithm"
XAXIS_LABEL = "xAxisLabel"
YAXIS_LABEL = "yAxisLabel"
PLOT_LABEL = "label"
INFILE_DIR = "indir"
# Available plot types are 'hist' and 'scatter'
PLOT_TYPE = "plotType"



'''
Arguments:
    plotDicts: list of dictionaries holding plot info
    outputPath (string): full output file path up to extension
    xmin,xmax,ymin,ymax (float): lower and upper values for x and y axis,
                                 respectively
    ratioIndex (int): add ratio plot taking the rationPlot'th entry in plotDicts
                      as the reference. If <1 no ratio plot is added
Description:
    Extracts a plot from each entry of the list and produces an overlay plot.
'''
def plot(plotDicts, outputPath, ratioIndex, xmax, xmin, ymax, ymin):

    # Extract the data as numpy
    #xDataArray = np.array()
    #yDataArray = np.array()
    #for pd in plotDicts:
        #xDataArray.append([pd[XDATA_NAME], axis=0])
        #yDataArray.append([pd[YDATA_NAME], axis=0])

    # prepare plots, maybe need a ratio plot
    referenceDict = None
    plotRows = 13
    plotColumns = 1
    # In case of having a ratio plot this is the number of rows the main plot
    # will get.
    rowspanMain = 10
    # Just a tuple with main and ratio axes
    axes = []

    if ratioIndex < 1:
        axMain = plt.subplot2grid((plotRows, plotColumns), (0, 0), rowspan=plotRows, colspan=plotColumns)
        axes.append(axMain)
        axMain.set_xlabel(plotDicts[0][XAXIS_LABEL])
    elif len(plotDicts) == 1:
        print "WARNING: Only one plot, ratio doesn't make sense here. Skip..."
    elif ratioIndex > len(plotDicts):
        print "WARNING: Have only " + str(len(plotDicts)) + " plots but " + str(ratioIndex) + " was requested as reference. Skip..."
    else:
        axMain = plt.subplot2grid((plotRows, plotColumns), (0, 0), rowspan=rowspanMain, colspan=plotColumns)
        axRatio = plt.subplot2grid((plotRows, plotColumns), (rowspanMain, 0), rowspan=(plotRows - rowspanMain), colspan=plotColumns, sharex=axMain)
        axes.append(axMain)
        axes.append(axRatio)
        # Make clear what the single points are normalised to
        axRatio.set_ylabel("1/" + plotDicts[ratioIndex-1][PLOT_LABEL])
        axRatio.set_xlabel(plotDicts[0][XAXIS_LABEL])
        referenceDict = plotDicts[ratioIndex-1]
    # Set to main axes
    # Set x and y labels for main plot taken from first entry in plotDicts
    # \todo Needs to be more general and given as an argument to this function.

    axes[0].set_ylabel(plotDicts[0][YAXIS_LABEL])

    colors = iter(cm.rainbow(np.linspace(0, 1, len(plotDicts))))
    markers = iter([".", "^", "v"])
    if xmax:
        axes[0].set_xlim(right=xmax)
    if xmin:
        axes[0].set_xlim(left=xmin)
    if ymax:
        axes[0].set_ylim(top=ymax)
    if ymin:
        axes[0].set_ylim(bottom=ymin)
    for pd in plotDicts:
        color = next(colors)
        marker = next(markers)
        axes[0].errorbar(pd[XDATA_NAME], pd[YDATA_NAME], xerr=pd[XDATA_ERROR_NAME], yerr=pd[YDATA_ERROR_NAME], color=color, marker=marker, linestyle="", label=pd[PLOT_LABEL])
        if referenceDict and pd != referenceDict:
            ratio = np.array(pd[YDATA_NAME]) / np.array(referenceDict[YDATA_NAME])
            ratioNomError = None
            ratioError = None

            if pd[YDATA_ERROR_NAME]:
                ratioNomError = np.square(np.array(pd[YDATA_ERROR_NAME]) / np.array(referenceDict[YDATA_NAME]))

            if referenceDict[YDATA_ERROR_NAME]:
                ratioError = np.square(np.array(pd[YDATA_NAME]) * np.array(referenceDict[YDATA_ERROR_NAME])  / np.square(np.array(referenceDict[YDATA_NAME])))


            # Sum the squares if possible
            if ratioNomError is not None and ratioError is not None:
                ratioError = ratioError + ratioNomError
            # Otherwise try only the nomError
            elif ratioNomError is not None and ratioError is None:
                ratioError = ratioNomError
            if ratioError is not None:
                ratioError = np.sqrt(ratioError)
            # So now ratioError is either None or a filled numpy array

            # Plot the stuff
            axes[1].errorbar(pd[XDATA_NAME], ratio, yerr=ratioError, color=color, marker=marker, linestyle="")

    axes[0].legend(loc="best", fontsize=10)
    plt.savefig(outputPath + ".eps")
    plt.savefig(outputPath + ".png")


'''
Arguments:
    configFile: file path to a JSON file
    indir: overwrites the input directory prefix
Description:
    Reads configFile and expecting following fields:
        required:
            'files': list of input files (their paths)
            'nEvents': list of number of events belonging to a file
                       (==> len(files) == len(nEvents))
        optional:
            'label': a label which will be printed in a legend of a plot
            'xAxisLabel': ...
            'yAxisLabel': ...

    Everything but the files are directly forwarded to a doctionary. From each file
    is the mean CPU time is extracted and also added to the dictionary as a list of
    means.
    ==> len(means) == len(nEvents)
'''
def extractFiles(configFile, indir):

    # The dictionary the plotting parameters are written to.
    configDict = {}

    # Try to extract parameters from passed JSON file.
    try:
        with open(configFile, 'r') as f:
            configDict = json.load(f)
    except Exception, e:
        print "ERROR: Could not open config JSON file " + configFile
        exit(1)

    # Check whether any parameters are specified.
    if len(configDict) == 0:
        print "ERROR: No parameters found in JSON file " + configFile
        exit(1)

    # Config seems to be sane.

    # Extract axis labels if any.
    configDict[XAXIS_LABEL] = configDict.get(XAXIS_LABEL, "xAxis")
    configDict[YAXIS_LABEL] = configDict.get(YAXIS_LABEL, "yAxis")
    # Extract a label if any.
    configDict[PLOT_LABEL] = configDict.get(PLOT_LABEL, "label")

    # Check for required x data
    if not configDict.get(XDATA_NAME, None):
        print "ERROR: No x-data found"
        exit(1)
    configDict[XDATA_ERROR_NAME] = configDict.get(XDATA_ERROR_NAME, None)
    # If number of error values is different from number of nominal values, fail
    if configDict[XDATA_ERROR_NAME] and len(configDict[XDATA_NAME]) != len(configDict[XDATA_ERROR_NAME]):
        print "ERROR: Same number of nominal values and errors required for x-data"
        exit(1)

    # Check if y data are already there and just return if so.
    configDict[YDATA_ERROR_NAME] = configDict.get(YDATA_ERROR_NAME, None)
    if configDict.get(YDATA_NAME, None):
        # Check for same length of x and y
        if len(configDict[YDATA_NAME]) != len(configDict[XDATA_NAME]):
            print "ERROR: Same number of x- and y-data (including errors) required"
            exit(1)
        if configDict[YDATA_ERROR_NAME] and len(configDict[YDATA_ERROR_NAME]) != len(configDict[YDATA_NAME]):
            print "ERROR: Same number for y-data and respective errors required"
            exit(1)
        return configDict


    # Expect input files to read data from.
    if not configDict.get(FILES, None):
        print "ERROR: Input files or y-data required"
        exit(1)

    if len(configDict[XDATA_NAME]) != len(configDict[FILES]):
        print "ERROR: Need as many x-data as files " + configFile
        exit(1)

    # Get extraction algorithm for y data.
    configDict[EXTRACTION_ALGORITHM] = configDict.get(EXTRACTION_ALGORITHM, "identity")

    # Joing file names with prefix
    prefix = configDict.get(INFILE_DIR, indir)
    configDict[FILES] = [os.path.join(prefix, inFile) for inFile in configDict[FILES]]
    configDict[YDATA_NAME] = []
    configDict[YDATA_ERROR_NAME] = []
    for file in configDict[FILES]:
        value, error = extract(file, configDict[EXTRACTION_ALGORITHM])
        configDict[YDATA_NAME].append(value)
        configDict[YDATA_ERROR_NAME].append(error)

    return configDict

'''
Arguments:
    filepath: Path to input file data points should be extracted from
    algorithm: algorithm to be used ('identity' or 'mean')
Description:
    This extracts one data point as a function of other data. The function is
    given by the argument 'algorithm'
'''
def extract(filepath, algorithm):

    value = -1.
    error = 0.
    with open(filepath, 'r') as f:
        if algorithm == "identity":
            value, error = extractIdentity(f)
        elif algorithm == "mean":
            value, error = extractMean(f)
        else:
            print "ERROR: Unknown extraction algorithm " + algorithm
            exit(1)
    return value, error

'''
Arguments:
    file: Opened file to be processed
Returns:
    1. Just one value read from a file
    2. an error of 0.
Description:
    Assumes one number in the file. If there are more just the first one is
    returned.
'''
def extractIdentity(file):
    for line in file.readlines():
        lineStripped = line.strip()
        # Check for float
        try:
            lineStripped = float(lineStripped)
            return lineStripped, 0.
        except Exception, e:
            print "ERROR: Found entry which is not a float"
            exit(1)

'''
Arguments:
    file: Opened file to be processed
Return:
    1. The mean value of values found in file
    2. standard deviation derived from values in file
Description:
    Calculates and returns the mean of data points found in the file.
'''
def extractMean(file):
    values = []
    for line in file.readlines():
        lineStripped = line.strip()
        # Check for float
        try:
            lineStripped = float(lineStripped)
        # \todo That is for now a workaround but should in general not be accepted.
        except Exception, e:
            continue
        # Increment number of lines for averaging
        values.append(lineStripped)
    # Finally the actual mean and standard deviation are returned
    return np.mean(values), np.std(values)


'''
Arguments:
    configFiles: list of file paths to JSON files
    indir: Can be used as a prefix to where potential files are found (deprecated)
    outputPath: full path to file the plot should be written to
'''
def main(configFiles, indir, outputPath, ratioIndex, xmax, xmin, ymax, ymin):

    plotDicts = []

    # Fill dictionary with plotting info
    for cf in configFiles:
        plotDicts.append(extractFiles(cf, indir))
    plot(plotDicts, outputPath, ratioIndex, xmax, xmin, ymax, ymin)

    return 0


#if __name__ == "main":
import argparse
parser = argparse.ArgumentParser()
# Require a config file with some plotting info
parser.add_argument("-d", "--indir", default="./", help="directory where to find input files specified in config file")
parser.add_argument("-c", "--config", nargs="+", required=True, help="config file")
parser.add_argument("-o", "--outprefix", required=True, help="output file name incl. path")
parser.add_argument("-r", "--ratio", default=-1, type=int, help="index of plot used as reference for ratio")
parser.add_argument("--xmax", type=float, help="upper limit on x axis")
parser.add_argument("--xmin", type=float, help="lower limit on x axis")
parser.add_argument("--ymax", type=float, help="upper limit on y axis")
parser.add_argument("--ymin", type=float, help="lower limit on y axis")

args = parser.parse_args();

# Forward config to main
sys.exit(main(args.config, args.indir, args.outprefix, args.ratio, args.xmax, args.xmin, args.ymax, args.ymin))
