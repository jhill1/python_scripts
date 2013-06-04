#!/usr/bin/env python

import os
import sys
import math
import string
import biology_util
from mld_util import *
import argparse
import numpy
import glob

def main():

    parser = argparse.ArgumentParser(
         prog="plot mesh",
         description="""Plot a mesh (time vs depth) of the mesh"""
                     )
    parser.add_argument(
            '-v', 
            '--verbose', 
            action='store_true', 
            help="Verbose output: mainly progress reports.",
            default=False
            )
    parser.add_argument(
            '-m', 
            '--no_mld', 
            help="Do not plot the MLD",
            required=False,
            action="store_true",
            default=False
            )
    parser.add_argument(
            '--tke',
            help="Plot MLD based on TKE, not density",
            action='store_true',
            default=False
            )
    parser.add_argument(
            '-z', 
            type=float,
            help="Maximum depth to plot",
            default=1000,
            required=False,
            )
    parser.add_argument(
            '-n',
            '--nyears',
            type=int,
            default=None,
            help="Number of years to plot. Default is longest time from Fluidity data",
            required=False,
            )
    parser.add_argument(
            'input', 
            metavar='input',
            help="Fluidity results directory",
            )
    parser.add_argument(
            'output_file', 
            metavar='output_file',
            help='The output filename'
            )

    args = parser.parse_args()
    verbose = args.verbose
    adaptive = False
    output_file = args.output_file
    max_depth = args.z
    plot_tke = args.tke
    no_mld = args.no_mld
    if (args.nyears == None):
        nDays = None
    else:
        nDays = args.nyears*365
    # actually irrelavent as we do it by day of year...
    start = datetime.strptime("1970-01-01 00:00:00", "%Y-%m-%d %H:%M:%S")

    params = {
      'legend.fontsize': 24,
      'xtick.labelsize': 20,
      'ytick.labelsize': 20,
      'font.size' : 24,
      'axes.labelsize' : 24,
      'text.fontsize' : 24,
    }
    rcParams.update(params)
    if (verbose):
        print "Plot 2D plot of : "

    if (verbose):
        print "  Scanning "+d
    files = get_vtus(args.input)

    if (not verbose):
        # print progress bar
        total_vtus = len(files)
        percentPerVtu = 100.0/float(total_vtus)

    current_file = 1
    # set up our master variables
    mld = []
    dates = []
    fig = figure(figsize=(15,8),dpi=90)
    ax = fig.add_subplot(111)
    for vtu in files:
        if (not verbose):
            progress(50,current_file*percentPerVtu)
        # obtain surface values from each dataset
        try:
            if (verbose):
                print vtu
            os.stat(vtu)
        except:
            print "No such file: %s" % file
            sys.exit(1)

        # open vtu and derive the field indices of the edge at (x=0,y=0) ordered by depth
        u=vtktools.vtu(vtu)
        pos = u.GetLocations()
        ind = get_1d_indices(pos)

        # deal with depths
        depth = [-pos[i,2] for i in ind]
        

        # handle time and convert to calendar time
        time = u.GetScalarField('Time')
        tt = [time[i] for i in ind]
        cur_dt = start + timedelta(seconds=time[0])
        days = (cur_dt - start).days
        tt = [days for i in ind]
        dates.append( tt[0] )

        # plot out mesh
        ax.plot(tt, depth, 'ok', markersize=1.5)
        
        if plot_tke:
             d = u.GetScalarField('GLSTurbulentKineticEnergy')
             den = [d[i] for i in ind]
             mld.append(calc_mld_tke(den, depth))
        else:
            # grab density profile and calculate MLD_den
            d = u.GetScalarField('Density')
            den = [d[i] * 1000 for i in ind]
            mld.append(calc_mld_den(den, depth, den0=0.125) )

        current_file = current_file + 1

    if (not nDays):
        nDays = dates[-1]

    filename = output_file
    if (not no_mld):
        ax.plot(dates, mld, color="r", lw="2")
    ax.set_xlim(0,nDays)
    ax.set_ylim(max_depth,0)
    xlabel('Days')
    ylabel("Depth")
    savefig(filename, dpi=180,format='png')
    close(fig)

    



if __name__ == "__main__":
    main()
