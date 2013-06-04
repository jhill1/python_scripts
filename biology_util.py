#!/usr/bin/env python

###################################################################
# Bunch of utility functions for plotting up biology simulations
# from Fluidity
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

import matplotlib
matplotlib.use("Agg")
# Above is so we don't need and X-screen to save figures
from fluidity_tools import stat_parser
import numpy
from numpy import arange,concatenate,array,argsort,append,fromfile
from scipy.interpolate import interp1d
import os
import sys
import vtktools
import math
import string
from pylab import *
from matplotlib.ticker import MaxNLocator
import re 
from datetime import datetime, timedelta
import numpy.ma as MA

# This might also be needed for the saving of figures without X
ioff()

def set_data_files(directory, station):
    """ Function to set up the observed data files and the columns to read
        for parimary production, surface nutrients, surface chlorophyll
        and mixed layer depth.

        Returns two arrays, one contianing the output; the other containing
        the column numbers for the day of year and values from the files. The
        order is primary productivity, surface nutrients, surface chlorophyl-a,
        and mixed layer depth. The filename will be '' and columns None if those
        data are not available.
    """

    if (directory == None):
        directory = '.'
    output_files = []
    output_columns = []
    if station == 'papa':
        output_files.append(os.path.join(directory,'papa_prpr.dat'))
        output_files.append(os.path.join(directory,'stnp_b_NIT.dat'))
        output_files.append(os.path.join(directory,'stnp_chl_SURFCHLA.dat'))
        output_files.append(os.path.join(directory,'mld_stationpapa.csv'))
        output_columns.append([0,1])
        output_columns.append([2,6])
        output_columns.append([2,6])
        output_columns.append([2,6])
    elif station == 'bermuda':
        output_files.append(os.path.join(directory,'bats_pprod_ml.dat'))
        output_files.append(os.path.join(directory,'bats_no3_ml.dat'))
        output_files.append(os.path.join(directory,'bats_chl_ml.dat'))
        output_files.append(os.path.join(directory,'mld_bats.csv'))
        output_columns.append([1,3])
        output_columns.append([2,3])
        output_columns.append([1,4])
        output_columns.append([2,6])
    elif station == 'india':
        output_files.append(os.path.join(directory,'india_prpr.dat'))
        output_files.append(os.path.join(directory,'india_nn.dat'))
        output_files.append(os.path.join(directory,'india_chl.dat'))
        output_files.append('')
        output_columns.append([0,1])
        output_columns.append([0,1])
        output_columns.append([0,1])
        output_columns.append(None)
    else:
        print "Unknown station. Exiting"
        sys.exit(-1)

    return  output_files, output_columns


def load_observed_data(filename, columns, years=None):
    """ Loads observed data from text file with days in columns[0]
        and the value in columns[1]. 

        If years is not None, a repeated set of measurements is created
        for n years

        Returns: 2 arrays, one with days and one with values
    """
   
    try:
        days,values=loadtxt(filename,unpack=True,skiprows=1,usecols=(columns[0],columns[1]))
    except:
        try:
            days,values=loadtxt(filename,unpack=True,skiprows=1,usecols=(columns[0],columns[1]),delimiter=',')
        except:
            print "Error reading " + filename
            sys.exit(-1)


    if (not years == None):
        full_days = []
        full_values = []
        for i in range(0,int(years)+1):
            for d in days:
                full_days.append(i*365 + d)
            for v in values:
                full_values.append(v)
        days = array(full_days)
        values = array(full_values)
   
    return days, values


def plot_summary(filename, observations, days, mlds, chlr, zoo, pp, nit, labels, mld_times=None, start_day=0, end_day=365, z_end=500, params=None, grey=False, pp_averaged=False):
    """ Plot a summary graph for a biology run. This contains 5 plots:
        A large central plot of MLD, and beneath it; a 2x2 grid of biological
        parameters.

        Input: filename - where to save the image
               observations - list of observation data in the following order, as a tuple
                         0     [[pp_days],[pp_data]
                         1      [no_days],[no_data]
                         2      [chl_days],[chl_data]
                         3      [mld_days],[mld_data]]
               days - list of day numbers for modelled data
               mlds - list of mld simulation data as tuples (length n) of days and data, repeated for each model. n can
                      be up to 4, or we'll run out of colours.
               chlr - as above, but chlorophyll
               zoo - guess
               pp - guess this too
               nit - surface nutrients
               labels - length n array of labels to give lines
               start_day - usually zero
               end_days - 365, 720, etc
               z_end - maximum depth for MLD plot. Default 500m. Positive.
               params - want to change the font size, etc? Use this.
               grey - greyscale plotting, default False
               pp_averaged - is the PP averaged or integrated (default)
        """
    
    if (params == None):
        params = {
          'legend.fontsize': 18,
          'xtick.labelsize': 16,
          'ytick.labelsize': 16,
          'font.size' : 18,
          'axes.labelsize' : 18,
          'text.fontsize' : 18,
          'figure.subplot.hspace' : 0.5,
          'text.usetex'        : True
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
                'DarkBlue',
                'Olive',
                'DarkGoldenRod',
                ]

    # summary plot
    markersize = 3
    fig_summary = figure(figsize=(15,20),dpi=180)
    #mld
    ax = fig_summary.add_subplot(511)
    if (not observations[3][1] == None):
        ax.plot(observations[3][0],-1.0*observations[3][1],obs_colour,markersize=markersize, alpha=0.5)
    i = 0
    if (mld_times == None):
        mld_times = days
    for data in mlds:
        ax.plot(mld_times[i],data,colours[i],lw=2)
        i = i+1
    ax.set_ylim(-z_end, 0)
    ax.set_xlim(start_day,end_day)
    ax.yaxis.grid(True)
    xlabel('Day')
    ylabel('UML depth (m)')
    # chlorophyll / phyto
    ax = fig_summary.add_subplot(512)
    ax.plot(observations[2][0],observations[2][1],obs_colour,markersize=markersize, alpha=0.5)
    i = 0
    for data in chlr:
        ax.plot(days[i],data,colours[i],lw=2)
        i += 1
    ax.yaxis.grid(True)
    ax.set_xlim(start_day,end_day)
    xlabel('Day')
    ylabel(r"Surface Chlorophyll""\n"r"(mmol/m$^3$)")
    # zooplankton
    ax = fig_summary.add_subplot(513)
    i = 0
    for data in zoo:
        ax.plot(days[i],data,colours[i],lw=2)
        i += 1
    ax.yaxis.grid(True)
    ax.set_xlim(start_day,end_day)
    xlabel('Day')
    ylabel(r"Surface Zooplankton""\n""(mmol/m$^3$)")
    # primary production
    ax = fig_summary.add_subplot(514)
    ax.plot(observations[0][0],observations[0][1],obs_colour,markersize=markersize, alpha=0.5)
    i = 0
    for data in pp:
        ax.plot(days[i],data,colours[i],lw=2)
        i += 1
    ax.set_xlim(start_day,end_day)
    ax.yaxis.grid(True)
    xlabel('Day')
    if (pp_averaged):
        ylabel(r"Averaged primary production""\n"r"(mgC/m$^3$/day)")
    else:
        ylabel(r"Integrated primary production""\n"r"(mgC/m$^2$/day)")
    # nutirents
    ax = fig_summary.add_subplot(515)
    ax.plot(observations[1][0],observations[1][1],obs_colour,label="Obs.",markersize=markersize, alpha=0.5)
    i = 0
    for data in nit:
        ax.plot(days[i],data,colours[i],label=labels[i],lw=2)
        i += 1
    ax.set_xlim(start_day,end_day)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25),
          ncol=3, fancybox=True, shadow=True)
    ax.yaxis.grid(True)
    xlabel('Day')
    ylabel(r"Surface Nutrients""\n"r"(mmol/m$^3$)")
    savefig(filename, dpi=90,format='png')


def write_cache(filename, data):
    """ Write a cache file"""

    import pickle

    f = open(filename,'wb')
    


def read_cache(filename):
    """ Read a cache and return the data"""
