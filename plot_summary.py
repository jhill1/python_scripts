#!/usr/bin/env python

import os
import sys
import math
import string
import biology_util
from mld_util import *
import argparse

def main():

    parser = argparse.ArgumentParser(
         prog="plot summary",
         description="""Plot the five part summary for biology models. Add multiple directories (up to
                     4) to plot different Fluidity outputs on the same graph.
                     
                     The station option sets which files to load and the obs_directory should contain those files."""
                     )
    parser.add_argument(
            '-v', 
            '--verbose', 
            action='store_true', 
            help="Verbose output: mainly progress reports.",
            default=False
            )
    parser.add_argument(
            '-s', 
            '--station', 
            choices=['papa', 'india', 'bermuda'],
            help="Choose a station: india, papa or bermuda",
            required=True
            )
    parser.add_argument(
            '--emulator',
            help="Plot emulator run, so MLD is precalculated. Should be a file called mld_cache.csv",
            action='store_true',
            default=False
            ) 
    parser.add_argument(
            '-o', 
            '--output', 
            nargs='+',
            help="Fluidity output directories",
            required=True
            )
    parser.add_argument(
            '-d', 
            '--directory', 
            help="Where the observation data are. Defaults to ../obs_data/",
            required=False,
            default="../obs_data/"
            )
    parser.add_argument(
            '--data',
            help="Add data to the plot. Data plotted depends on station",
            action="store_true",
            default="False"
            )
    parser.add_argument(
            '-l', 
            '--labels', 
            help="What to call the Fluidity data",
            required=True,
            nargs='+'
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
            default=500,
            required=False,
            )
    parser.add_argument(
            '--startday',
            type=int,
            help="Start day to use for plotting",
            default=0,
            required=False,
            )
    parser.add_argument(
            '-n',
            '--nyears',
            type=int,
            help="Number of years to plot. Default is longest time from Fluidity data",
            required=False,
            )
    parser.add_argument(
            'output_file', 
            metavar='output_file',
            nargs=1,
            help='The output filename'
            )

    args = parser.parse_args()
    verbose = args.verbose
    adaptive = False
    labels = args.labels
    output_file = args.output_file[0]
    max_depth = args.z
    plot_tke = args.tke
    emulator = args.emulator
    plot_data = args.data
    if (args.nyears == None):
        nDays = None
    else:
        nDays = args.nyears*365
    # actually irrelavent as we do it by day of year...
    start = datetime.strptime("1970-01-01 00:00:00", "%Y-%m-%d %H:%M:%S")


    if (verbose):
        print "Plot Summary: "

    files = []
    for d in args.output:
        if (verbose):
            print "  Scanning "+d
        files.append(get_vtus(d))

    if (not verbose):
        # print progress bar
        total_vtus = 0
        for f in files:
            total_vtus += len(f)

        percentPerVtu = 100.0/float(total_vtus)

    if (args.station == 'bermuda'):
        pp_averaged = True
    else:
        pp_averaged = False

    current_file = 1
    # set up our master variables
    mlds = []
    chlr = []
    zoo = []
    pp = []
    nit = []
    times_all = []
    mld_times_all = []
    current_sim = 0
    for sim in files:
        # these are our variables from this simulation
        nutrient = []
        primprod = []
        chlorophyll = []
        mld = []
        dates = []
        times = []
        if (emulator):
            mld, mld_times = read_cache(os.path.join(args.output[current_sim],'mld_cache.csv'))
            # mls gets appended later
            mld_times_all.append(mld_times)

        primprod, times, dates = biology_util.primary_productivity(sim, start, pp_averaged)
        chlorophyll, times, dates = biology_util.surface_chlorophyll(sim, start)
        nutrient, times, dates = biology_util.surface_nutrient(sim, start)
        dates = []
        times = []
        for vtu in sim:
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

            # handle time and convert to calendar time
            time = u.GetScalarField('Time')
            tt = [time[i] for i in ind]
            if (adaptive):
                tt.extend([time[0] for i in range(len(ind),nPoints)])
            cur_dt = start + timedelta(seconds=time[0])
            days = (cur_dt - start).days
            if (days < start_day):
                continue
            dates.append( days )

            # deal with depths
            depth = [-pos[i,2] for i in ind]
            if (adaptive):
                depth.extend([maxDepth+1 for i in range(len(ind),nPoints)])


            
            if (not emulator):
                if plot_tke:
                     d = u.GetScalarField('GLSTurbulentKineticEnergy')
                     den = [d[i] for i in ind]
                     if (adaptive):
                         den.extend([1e-10 for i in range(len(ind),nPoints)])
                     mld.append(-1 * calc_mld_tke(den, depth))
                else:
                    # grab density profile and calculate MLD_den
                    d = u.GetScalarField('Density')
                    den = [d[i] * 1000 for i in ind]
                    if (adaptive):
                        den.extend([d[-1]*1000. for i in range(len(ind),nPoints)])
                    mld.append( -1 * calc_mld_den(den, depth, den0=0.125) )
           
            
            current_file = current_file + 1
            # end of this sims data gathering
        # collate data from all simulations
        mlds.append(mld)
        chlr.append(chlorophyll)
        pp.append(primprod)
        nit.append(nutrient)
        times_all.append(dates)
        current_sim += 1

    if (not nDays):
        nDays = times_all[0][-1]

    nYears = int(nDays / 365) + 1

    # get observed values
    obs_files, obs_columns = biology_util.set_data_files(args.directory,args.station)
    observations = []
    i = 0
    for d in obs_files:
        days = None
        data = None
        if (not d == ''):
            days, data = biology_util.load_observed_data(d,obs_columns[i],years=nYears)
        i += 1
        observations.append([days, data])

    filename = output_file
    if (len(mld_times_all) == 0):
        mld_times_all = None

    biology_util.plot_summary(filename, observations, times_all, mlds, chlr, pp, nit, labels, mld_times=times_all, end_day=nDays, z_end=max_depth, pp_averaged=pp_averaged, plot_data=plot_data)
 

def read_cache(file_in):
    f = open(file_in , 'r')
    mld = []
    mld_times = []
    for line in f:
        data = line.split(",")
        mld.append(-1.0*abs(float(data[0])))
        mld_times.append(abs(float(data[1]))/86400.) #days 
    f.close
    return mld, mld_times

if __name__ == "__main__":
    main()
