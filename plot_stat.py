#!/usr/bin/env python
###################################################################
# Plot stuff from multiple stat files on a single plot
#
#Copyright (C) 2013 Jon Hill, jon.hill@imperial.ac.uk
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


import os
import sys
import math
import string
import biology_util
from mld_util import *
import argparse
from fluidity_tools import stat_parser
from pylab import *

def main():

    parser = argparse.ArgumentParser(
         prog="plot stat",
         description="""Plot a variable from the stat file for several runs."""
                     )
    parser.add_argument(
            '-v', 
            '--verbose', 
            action='store_true', 
            help="Verbose output: mainly progress reports.",
            default=False
            )
    parser.add_argument(
            '-o', 
            '--output', 
            nargs='+',
            help="Fluidity stat files",
            required=True
            )
    parser.add_argument(
            '-l', 
            '--labels', 
            help="What to call the Fluidity data",
            required=True,
            nargs='+'
            )
    parser.add_argument(
            '-n',
            '--nyears',
            type=int,
            help="Number of years to plot. Default is longest time from Fluidity data",
            required=False,
            )
    parser.add_argument(
            '-s',
            '--smooth',
            type=int,
            help="Take every n data points from the statplot. Use a large number if your stat file is large.",
            required=False,
            default = 1
            )
    parser.add_argument(
            '-y',
            '--ylabel',
            help="Supply a different label for the y-axis",
            required=False,
            default = 1
            )
    parser.add_argument(
            '-m',
            '--material',
            help="Material name. Default is Fluid",
            required=False,
            default = "Fluid"
            )

    parser.add_argument(
            'variable',
            metavar='variable',
            nargs=1,
            help="The variable to pull",
            )
    parser.add_argument(
            'statistic',
            metavar='statistic',
            nargs=1,
            help="The statistic to pull",
            )
    parser.add_argument(
            'output_file', 
            metavar='output_file',
            nargs=1,
            help='The output filename'
            )

    args = parser.parse_args()
    verbose = args.verbose
    labels = args.labels
    output_file = args.output_file[0]
    variable = args.variable[0]
    statistic = args.statistic[0]
    smooth = args.smooth
    ylabel_string = args.ylabel
    material = args.material
    if (args.nyears == None):
        nDays = None
    else:
        nDays = args.nyears*365
    if (ylabel_string == None):
        ylabel_string = statistic + " " + variable
    

    if (verbose):
        print "Plot Summary: "

    data_array = []
    time_array = []

    files = []
    longest_day = 0
    for s in args.output:
        if (verbose):
            print "  Scanning "+d
        # grab data from stat file
        stat=stat_parser( s, smooth )
        time = stat["ElapsedTime"]["value"]
        days = []
        for t in time:
            days.append(t/(24*60*60))
        time_array.append(days)
        if (days[-1] > longest_day): 
            longest_day = days[-1]
        time = []
        days = []
        data = stat[material][variable][statistic]
        data_array.append(data)

    if (not nDays):
        nDays = longest_day


    colours = [
            'r',
            'b',
            'g',
            'k',
            'DarkBlue',
            'Olive',
            'DarkGoldenRod',
            ]
    params = {
      'legend.fontsize': 30,
      'xtick.labelsize': 26,
      'ytick.labelsize': 26,
      'font.size' : 26,
      'axes.labelsize' : 30,
      'text.fontsize' : 30,
      'figure.subplot.hspace' : 0.5,
      'text.usetex' : True
        }
    rcParams.update(params)

    i = 0
    fig = figure(figsize=(25,10),dpi=90)
    ax = fig.add_subplot(111) 
    for data in data_array:
        ax.plot(time_array[i][:],data[:],colours[i],label=labels[i],lw=3)
        i = i+1

    ax.yaxis.grid(True)
    xlabel('Time (days)')
    ylabel(ylabel_string)
    ax.set_xlim(0,nDays)
    #ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1),
    #      ncol=2, fancybox=True, shadow=True)
    ax.legend(loc=0,ncol=2,fancybox=True,shadow=True)
    savefig(output_file, dpi=90,format='png')

if __name__ == "__main__":
    main()
