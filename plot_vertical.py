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
         prog="plot vertical",
         description="""Plot a vertical profile of a single variable from multiple VTU files"""
                     )
    parser.add_argument(
            '-v', 
            '--verbose', 
            action='store_true', 
            help="Verbose output: mainly progress reports.",
            default=False
            )
    parser.add_argument(
            '--var', 
            choices=['dz', 'chlor', 'det', 'nit', 'zoo', 'phyto', 'pp', 'temp', 'sal', 'rho', 'tke','vertdiff'],
            help="Choose a variable to plot",
            required=True
            )
    parser.add_argument(
            '-m', 
            '--mld', 
            help="Plot a horizontal line representing the MLD",
            required=False,
            action="store_true",
            default=False
            )
    parser.add_argument(
            '--tke',
            help="Plot MLD based on TKE, not density. Only useful with --mld",
            action='store_true',
            default=False
            )
    parser.add_argument(
            '-z', 
            type=float,
            help="Maximum depth to plot",
            default=100,
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
            help="Use a log scale on the horizontal axis",
            action='store_true',
            default=False
            )   
    parser.add_argument(
            '--adaptive',
            help="One of the input vtu is from an adaptive run",
            action='store_true',
            default=False
            )    
    parser.add_argument(
            '--grey',
            help="Greyscale plt",
            action='store_true',
            default=False
            )
    parser.add_argument(
            '--summer',
            metavar='summer',
            help="Fluidity VTU files for plot 1 - summer",
            required=True,
            nargs='+'
            )
    parser.add_argument(
            '--winter',
            metavar='winter',
            help="Fluidity VTU files for plot 2 - winter",
            required=True,
            nargs='+'
            )
    parser.add_argument(
            '-l', 
            '--labels', 
            help="What to call the Fluidity data on the plot",
            required=True,
            nargs='+'
            )
    parser.add_argument(
            'output_file', 
            metavar='output_file',
            help='The output filename'
            )

    args = parser.parse_args()
    verbose = args.verbose
    output_file = args.output_file
    labels = args.labels

    max_depth = args.z
    plot_mld = args.mld
    input_files_summer = args.summer
    input_files_winter = args.winter
    plot_tke = args.tke
    var = args.var
    minvalue = args.minvalue
    maxvalue = args.maxvalue
    logscale = args.logscale
    adaptive = args.adaptive
    grey = args.grey
    
    variable = None
    long_name = None
    nPoints = 1000
    # work out which variable to pull out
    # 'dz', 'chlor', 'det', 'nit', 'zoo', 'phyto', 'pp', 'temp', 'sal', 'rho', 'tke'
    if (var == 'dz'):
        variable = None # don't need one - special case
        long_name = r'$\mathrm{\delta z}$'
    elif (var == 'chlor'):
        variable = 'Chlorophyll'
        long_name = 'Chlorophyll-a ($\mathrm{mmol/m^3}$)'
    elif (var == 'det'):
        variable = 'Detritus'
        long_name = 'Detritus ($\mathrm{mmol/m^3}$)'
    elif (var == 'nit'):
        variable = 'Nutrient'
        long_name = 'Nitrate ($\mathrm{mmol/m^3}$)'
    elif (var == 'zoo'):
        variable = 'Zooplankton'
        long_name = 'Zooplankton ($\mathrm{mmol/m^3}$)'
    elif (var == 'phyto'):
        variable = 'Phytoplankton'
        long_name = 'Phytoplankton ($\mathrm{mmol/m^3}$)'
    elif (var == 'pp'):
        variable = 'DailyAveragedPrimaryProduction'
        long_name = 'Primary Production (\mathrm{$mmol/m^3}$)'
    elif (var == 'temp'):
        variable = 'Temperature'
        long_name = 'Temperature ($\mathrm{^/circ C}$)'
    elif (var == 'sal'):
        variable = 'Salinity'
        long_name = 'Salinity ($\mathrm{PSU}$)'
    elif (var == 'rho'):
        variable = 'Density'
        long_name = 'Density ($\mathrm{kg/m^3}$)'
    elif (var == 'tke'):
        variable = 'GLSTurbulentKineticEnergy'
        long_name = 'Turbulent Kinetic Energy ($mathrm{m^2s^{-2}}$)'
    elif (var == 'vertdiff'):
        variable = 'GLSVerticalDiffusivity'
        long_name = 'Vertical Diffusivity ($\mathrm{m^2s^{-1}}$)'
    else:
        print "Unknown variable"
        sys.exit(-1)

    if (verbose):
        print "Plot vertical profile of : "

    params = {
          'legend.fontsize': 34,
          'xtick.labelsize': 34,
          'ytick.labelsize': 34,
          'font.size' : 34,
          'axes.labelsize' : 34,
          'text.fontsize' : 34,
          'figure.subplot.hspace' : 0.5,
          'text.usetex'        : True,
          'xtick.major.pad' : 10,
          'ytick.major.pad' : 10,
    }
    rcParams.update(params)
    if (grey):
        obs_colour = 'sk'
        colours = [
                'k-',
                'k--',
                'k:',
                'k.-',
                '0.75',
                '0.5',
                '0.25',
                ]
    else:
        obs_colour = 'sg'
        colours = [
                'r',
                'b',
                'g',
                'k',
                'Olive',
                'DarkGoldenRod',
                ]


    fig_summary = figure(figsize=(20,15),dpi=90)
    
    # summer plot - LHS
    current_file = 1
    # set up our master variables
    mld = []
    data = []
    dates = []
    depths = []
    last_index = []
    for vtu in input_files_summer:
        # obtain surface values from each dataset
        try:
            os.stat(vtu)
        except:
            print "No such file: %s" % vtu
            # pop the label 
            continue

        # open vtu and derive the field indices of the edge at (x=0,y=0) ordered by depth
        u=vtktools.vtu(vtu)
        pos = u.GetLocations()
        ind = get_1d_indices(pos)

        # deal with depths
        depth = [-pos[i,2] for i in ind]
        if (adaptive):
                depth.extend([depth[-1] for i in range(len(ind),nPoints)])
        depths.append(depth)
        # find last index that contains the depth of interest
        depth_ind = 0
        for d in depth:
            if (d > max_depth):
                last_index.append(depth_ind+1)
                break
            depth_ind = depth_ind+1


        if plot_tke:
             d = u.GetScalarField('GLSTurbulentKineticEnergy')
             den = [d[i] for i in ind]
             mld.append(calc_mld_tke(den, depth))
        else:
            # grab density profile and calculate MLD_den
            d = u.GetScalarField('Density')
            den = [d[i] * 1000 for i in ind]
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

    filename = output_file
    data = numpy.array(data)

    if (var == 'dz'):
        # Add or subtract a bit more than edgelength as DZ might be a bit more or less than this
        min_data = minEdgeLength - (0.1*minEdgeLength)
        max_data = maxEdgeLength + (0.1*maxEdgeLength)
    
    ax1 = fig_summary.add_subplot(131)
    i = 0
    for d in data:
        ax1.plot(d[0:last_index[i]],depths[i][0:last_index[i]],colours[i],lw=2,ms=10,label=labels[i])
        i=i+1
   
    if (plot_mld):
        ax1.axhline(y=mld[0], lw=5,alpha=0.5,color='0.5')
    ax1.set_ylim(max_depth,0)
    for label in ax1.get_yticklabels():
        label.set_verticalalignment('top')

    #ax.set_xlim(min_data,max_data)
    ax1.xaxis.tick_top()
    ax1.yaxis.grid(True)
    ylabel('Depth (m)')


    # winter plot - LHS
    current_file = 1
    # set up our master variables
    mld = []
    data = []
    dates = []
    depths = []
    last_index = []
    for vtu in input_files_winter:
        # obtain surface values from each dataset
        try:
            os.stat(vtu)
        except:
            print "No such file: %s" % vtu
            # pop the label 
            continue

        # open vtu and derive the field indices of the edge at (x=0,y=0) ordered by depth
        u=vtktools.vtu(vtu)
        pos = u.GetLocations()
        ind = get_1d_indices(pos)

        # deal with depths
        depth = [-pos[i,2] for i in ind]
        if (adaptive):
                depth.extend([depth[-1] for i in range(len(ind),nPoints)])
        depths.append(depth)
        # find last index that contains the depth of interest
        depth_ind = 0
        for d in depth:
            if (d > max_depth):
                last_index.append(depth_ind+1)
                break
            depth_ind = depth_ind+1


        if plot_tke:
             d = u.GetScalarField('GLSTurbulentKineticEnergy')
             den = [d[i] for i in ind]
             mld.append(calc_mld_tke(den, depth))
        else:
            # grab density profile and calculate MLD_den
            d = u.GetScalarField('Density')
            den = [d[i] * 1000 for i in ind]
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

    filename = output_file
    data = numpy.array(data)

    if (var == 'dz'):
        # Add or subtract a bit more than edgelength as DZ might be a bit more or less than this
        min_data = minEdgeLength - (0.1*minEdgeLength)
        max_data = maxEdgeLength + (0.1*maxEdgeLength)
    
    ax2 = fig_summary.add_subplot(132, sharey=ax1)
    i = 0
    for d in data:
        ax2.plot(d[0:last_index[i]],depths[i][0:last_index[i]],colours[i],lw=2,ms=10,label=labels[i])
        i=i+1
   
    if (plot_mld):
        ax2.axhline(y=mld[0], lw=5,alpha=0.5,color='0.5')
    ax2.set_ylim(max_depth,0)
    ax2.xaxis.tick_top()
    ax2.yaxis.grid(True)
    for label in ax2.get_yticklabels():
        label.set_verticalalignment('top')

    xlabel(long_name, x=0.5, ha='right')
    ax2.xaxis.labelpad = 20
    ax2.xaxis.set_label_position('top') 
    locator_params(axis = 'x', nbins = 4)
    legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    savefig(output_file, dpi=90,format='svg')
    
if __name__ == "__main__":
    main()
