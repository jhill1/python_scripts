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
import libspud

def main():

    parser = argparse.ArgumentParser(
         prog="plot 2D",
         description="""Plot a 2D (time vs depth) of a variable"""
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
            '--var', 
            choices=['dz', 'chlor', 'det', 'nit', 'zoo', 'phyto', 'pp', 'temp', 'sal', 'rho', 'tke','vertdiff'],
            help="Choose a variable to plot",
            required=True
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
            default=500,
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
            '--minvalue',
            default=None,
            type=float,
            help="Set a minimum value to use in the plot",
            required=False,
            )
    parser.add_argument(
            '--maxvalue',
            default=None,
	        type=float,
            help="Set a maximum value to use in the plot",
            required=False,
            )
    parser.add_argument(
            '--logscale',
            help="Use a log colour scale",
            action='store_true',
            default=False
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
    var = args.var
    minvalue = args.minvalue
    maxvalue = args.maxvalue
    logscale = args.logscale
    if (args.nyears == None):
        nDays = None
    else:
        nDays = args.nyears*365
    # actually irrelavent as we do it by day of year...
    start = datetime.strptime("1970-01-01 00:00:00", "%Y-%m-%d %H:%M:%S")

    variable = None
    long_name = None
    # work out which variable to pull out
    # 'dz', 'chlor', 'det', 'nit', 'zoo', 'phyto', 'pp', 'temp', 'sal', 'rho', 'tke'
    if (var == 'dz'):
        variable = None # don't need one - special case
        long_name = r'$\delta z$'
    elif (var == 'chlor'):
        variable = 'Chlorophyll'
        long_name = 'Chlorophyll-a ($mmol/m^3$)'
    elif (var == 'det'):
        variable = 'Detritus'
        long_name = 'Detritus ($mmol/m^3$)'
    elif (var == 'nit'):
        variable = 'Nutrient'
        long_name = 'Nitrate ($mmol/m^3$)'
    elif (var == 'zoo'):
        variable = 'Zooplankton'
        long_name = 'Zooplankton ($mmol/m^3$)'
    elif (var == 'phyto'):
        variable = 'Phytoplankton'
        long_name = 'Phytoplankton ($mmol/m^3$)'
    elif (var == 'pp'):
        variable = 'DailyAveragedPrimaryProduction'
        long_name = 'Primary Production ($mmol/m^3$)'
    elif (var == 'temp'):
        variable = 'Temperature'
        long_name = 'Temperature ($^/circ C$)'
    elif (var == 'sal'):
        variable = 'Salinity'
        long_name = 'Salinity ($PSU$)'
    elif (var == 'rho'):
        variable = 'Density'
        long_name = 'Density ($kg/m^3$)'
    elif (var == 'tke'):
        variable = 'GLSTurbulentKineticEnergy'
        long_name = 'Turbulent Kinetic Energy ($m^2s^{-2}$)'
    elif (var == 'vertdiff'):
        variable = 'GLSVerticalDiffusivity'
        long_name = 'Vertical Diffusivity ($m^2s^{-1}$)'
    else:
        print "Unknown variable"
        sys.exit(-1)

    if (verbose):
        print "Plot 2D plot of : "

    files = get_vtus(args.input)

    if (not verbose):
        # print progress bar
        total_vtus = len(files)
        percentPerVtu = 100.0/float(total_vtus)

    if (args.station == 'bermuda'):
        pp_averaged = True
    else:
        pp_averaged = False

    # Are we adaptive - if so, set up some variables
    # We need maxdepth, min edge length, max edge length, etc
    # I think the easiest is to parse the flml...glob flml files
    # that don't have checkpoint in their name...
    flmls = (glob.glob(os.path.join(args.input, "*.flml")))
    flmls.sort()
    
    # grab bits of the flml we want
    if flmls == None or len(flmls) == 0:
        print "No flml files available: can't check for adaptivity, this might fail..."
    elif len(flmls) > 1:
        print "Found multiple FLMLs, using the first one"
    else:
        flml = flmls[0]
        # grab min_edge_length, max_edge_length, bottom depth and starting resolution
        libspud.load_options(flml)
        if (libspud.have_option('mesh_adaptivity/')):
            adaptive = True
            minEdgeLength = libspud.get_option('/mesh_adaptivity/hr_adaptivity/tensor_field::MinimumEdgeLengths/anisotropic_symmetric/constant')[-1][-1]
            maxEdgeLength = libspud.get_option('/mesh_adaptivity/hr_adaptivity/tensor_field::MaximumEdgeLengths/anisotropic_symmetric/constant')[-1][-1]
            maxDepth = libspud.get_option('/geometry/mesh::CoordinateMesh/from_mesh/extrude/regions::WholeMesh/bottom_depth/constant')
            initial_spacing = min(minEdgeLength,libspud.get_option('/geometry/mesh::CoordinateMesh/from_mesh/extrude/regions::WholeMesh/sizing_function/constant'))
            nPoints = int(maxDepth/max(initial_spacing,minEdgeLength))+1 # maximum number of points will be at start with uniform resolution
            if (verbose):
                print minEdgeLength, maxEdgeLength, maxDepth, nPoints, initial_spacing
        else:
        # If not found, it's not adaptive...
            adaptive = False

    if (verbose):
        print "Adaptive run: ", adaptive

    current_file = 1
    # set up our master variables
    mld = []
    data = []
    dates = []
    depths = []
    for vtu in files:
        if (not verbose):
            progress(50,current_file*percentPerVtu)
        # obtain surface values from each dataset
        try:
            #if (verbose):
            #    print vtu
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
        if (adaptive):
            depth.extend([maxDepth+1 for i in range(len(ind),nPoints)])
        depths.append(depth)

        # handle time and convert to calendar time
        time = u.GetScalarField('Time')
        tt = [time[i] for i in ind]
        cur_dt = start + timedelta(seconds=time[0])
        days = (cur_dt - start).days
        tt = [days for i in ind]
        if (adaptive):
            tt.extend([days for i in range(len(ind),nPoints)])
        dates.append( tt )
        
        if plot_tke:
             d = u.GetScalarField('GLSTurbulentKineticEnergy')
             den = [d[i] for i in ind]
             if (adaptive):
                 den.extend([1e-10 for i in range(len(ind),nPoints)])
             mld.append(calc_mld_tke(den, depth))
        else:
            # grab density profile and calculate MLD_den
            d = u.GetScalarField('Density')
            den = [d[i] * 1000 for i in ind]
            if (adaptive):
                den.extend([d[-1]*1000. for i in range(len(ind),nPoints)])
            mld.append(calc_mld_den(den, depth, den0=0.125) )


        if (var == 'dz'):
            z = []
            for i in range(1,len(depth)):
                if (depth[i] - depth[i-1] < minEdgeLength): # we're in the region where we've made up the edges...
                    z.append(minEdgeLength)
                else:
                    z.append(depth[i] - depth[i-1])
            z.append(100) # to make sure shape(dz) == shape(depths)
            data.append(vtktools.arr(z))
        else:
            # primprod is depth integrated or averaged, so we need the whole field
            t = u.GetScalarField(variable)
            if (var == 'pp'):
                temp = [6.7*12*8.64e4*t[i] for i in ind]
            elif (var == 'density'):
                temp = [1024.0*t[i] for i in ind]
            else:
                temp = [t[i] for i in ind]
            if (adaptive):
                temp.extend([t[-1] for i in range(len(ind),nPoints)])
            data.append(temp)

        current_file = current_file + 1

    if (not nDays):
        nDays = dates[-1][0]

    filename = output_file
    data = numpy.array(data)

    if (var == 'dz'):
        # Add or subtract a bit more than edgelength as DZ might be a bit more or less than this
        min_data = minEdgeLength - (0.1*minEdgeLength)
        max_data = maxEdgeLength + (0.1*maxEdgeLength)
    else:
        if (minvalue == None):
            min_data = data.min()
        else:
            min_data = minvalue
        if (maxvalue == None):
            max_data = data.max()
        else:
            max_data = maxvalue
    plot_2d_data(data, depths, dates, filename, long_name, finish_day=nDays,mld_data=mld,max_depth=max_depth,minimum=min_data,maximum=max_data,logscale=logscale)

if __name__ == "__main__":
    main()
