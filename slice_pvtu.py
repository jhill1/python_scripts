#!/usr/bin/env python

###################################################################
# Slice a PVTU file in any direction, then sotre data as XYZ or
# plot it up. Still a bit rough around the edges.
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
import vtk
import numpy as np
import argparse
from lxml import etree
xml_parser = etree.XMLParser(remove_blank_text=True)
import os
import pickle
import math
import pylab
from scipy.interpolate import griddata
import spheretools as sp
import sys

# This might also be needed for the saving of figures without X
pylab.ioff()


def main():

    def coords(s):
        try:
            x, y, z = map(int, s.split(','))
            return x, y, z
        except:
            raise argparse.ArgumentTypeError("Origin and/or normal must be x,y,z")


    parser = argparse.ArgumentParser(
         prog="plot slice",
         description="""Plot the pre-generated pkl file from gen_slice.py"""
         )
    subparsers = parser.add_subparsers(help='sub-command help')

    # Create slice data command - PVTU to pkl
    parser_cm = subparsers.add_parser('create_slice',
            help='Create pkl file containing the slice data'
            )
    parser_cm.add_argument(
            '-v', 
            '--verbose', 
            action='store_true', 
            help="Verbose output: mainly progress reports.",
            default=False
            )
    parser_cm.add_argument(
            '-o',
            '--origin',
            help="Origin tuple, e.g. [0,0,0]",
            type=coords,
            required=False
            ) 
    parser_cm.add_argument(
            '-s',
            '--sphere',
            help="Data is on the sphere. A horizontal slice at depth N (see below) will be done.",
            required=False,
            action='store_true', 
            default=False
            )
    parser_cm.add_argument(
            '--sphere_radius',
            help="Radius of the sphere. Default is 6.37101e+06 m",
            required=False,
            type=float,
            default=6.37101e+06
            )
    parser_cm.add_argument(
            '--sphere_depth',
            help="Depth to create the slice. Default is 10 m. The depth is subtracted from the sphere_radius.",
            required=False,
            type=float,
            default=10.0
            )
    parser_cm.add_argument(
            '-n', 
            '--normal', 
            help="Normal direction for the slice, e.g. [0,0,1] for a horizontal slice, or [1,0,0] for a vertical slice along Y.",
            type=coords,
            required=False
            )
    parser_cm.add_argument(
            'input_file', 
            metavar='input_file',
            nargs=1,
            help='The input PVTU'
            )
    parser_cm.add_argument(
            'output_file', 
            metavar='output_file',
            nargs=1,
            help='The output pkl file'
            )
    parser_cm.add_argument(
            'variable',
            nargs=1,
            help="""The variable to plot. 
            Note that VelocityMagnitude, BedShearStressMagnitude, MaxBedShearStressMagnitude, and AveBedShearStressMagnitude 
            is also acceptable (and will calculate the appropriate magnitude). VelocityU and VelocityV will grab the U,V components
            respectively and likewise for BedShearStressU and BedShearStressV."""
            )
    parser_cm.set_defaults(func=gen_slice)


    # Plot slice command: pkl -> Png
    parser_cm = subparsers.add_parser('plot_slice',
            help='Plot slice from pkl data'
            )    
    parser_cm.add_argument(
            '-v', 
            '--verbose', 
            action='store_true', 
            help="Verbose output: mainly progress reports.",
            default=False
            )
    parser_cm.add_argument(
            '-s', 
            '--on_sphere', 
            action='store_true', 
            help="Plot the data on the sphere. If the data are on the sphere, "+ \
                 "but this option is not enabled, the data are projected into lon/lat space",
            default=False
            )
    parser_cm.add_argument(
            '-c',
            '--colourmap',
            help="Colormap to use. Default is jet",
            choices=['jet','hot','Blues','gray','rainbow'],
            default='jet'
            )
    parser_cm.add_argument(
            '--subsample',
            help="Subsample the pkl file. If your plot take > 5 mins to generate, try setting this to 5 or semething",
            default='1',
            type=int
            )
    parser_cm.add_argument(
            '--resolution',
            help="Set the resolution of the plot. Default is 0.5 degrees",
            type=float,
            default=0.5
            )
    parser_cm.add_argument(
            '--dpi',
            help="Set the DPI. Default is 180",
            default=180
            )
    parser_cm.add_argument(
            '--format',
            help="File format. Default is based on filename",
            choices=['pdf','png'],
            )
    parser_cm.add_argument(
            '--minvalue',
            type=float,
            help="Minimum value to use in plotting",
            )
    parser_cm.add_argument(
            '--maxvalue',
            type=float,
            help="Maximum value to use in plotting",
            )
    parser_cm.add_argument(
            '--addvelocity',
            nargs=1,
            help="Plot velocity vectors on top of a scalar",
            )
    parser_cm.add_argument(
            '--velocityscale',
            type=float,
            default=1,
            help="Scale the velocity vectors. Default is 1. Bigger number means smaller arrows",
            )
    parser_cm.add_argument(
            '--mask_file', 
            metavar='mask_file',
            nargs=1,
            help='A NETCDF mask. 0 is land, 1 is ocean. Variable should be called mask'
            )
    parser_cm.add_argument(
            '-l', 
            '--hardlimit', 
            help="Contourf is white where plots are over the limits specified. This puts any value > max or < min to the max/min respectively.",
            required=False,
            action='store_true', 
            default=False
            )
    parser_cm.add_argument(
            'input_file', 
            metavar='input_file',
            nargs=1,
            help='The input pkl file'
            )
    parser_cm.add_argument(
            'output_file', 
            metavar='output_file',
            nargs=1,
            help='The output PNG/PDF file'
            )
    parser_cm.add_argument(
            'variable',
            nargs=1,
            help="The name of the variable being plotted"
            )
    # Set the function that gets called when the sub command is activated
    parser_cm.set_defaults(func=plot_slice)

    # Plot slice command: pkl -> XYZ
    parser_cm = subparsers.add_parser('pkl_to_xyz',
            help='Convert pickle file to XYZ'
            )    
    parser_cm.add_argument(
            '--subsample',
            help="Subsample the pkl file. If your plot take > 5 mins to generate, try setting this to 5 or semething",
            default='1',
            type=int
            )
    parser_cm.add_argument(
            '-g', 
            '--grid', 
            action='store_true', 
            help="Grid the data onto regular grid specified by resolution. Default is to just output a point cloud in Lat/Lon space",
            default=False
            )
    parser_cm.add_argument(
            '--resolution',
            help="Set the resolution of the XYZ. Default is 0.5 degrees",
            type=float,
            default=0.5
            )
    parser_cm.add_argument(
             '--mask_file',
             metavar='mask_file',
             nargs=1,
             default=["no_mask"],
             help='A NETCDF mask. 0 is land, 1 is ocean. Variable should be called mask'
             )
    parser_cm.add_argument(
            'input_file', 
            metavar='input_file',
            nargs=1,
            help='The input pkl file'
            )
    parser_cm.add_argument(
            'output_file', 
            metavar='output_file',
            nargs=1,
            help='The output XYZ file'
            )
    parser_cm.add_argument(
            '-v', 
            '--verbose', 
            action='store_true', 
            help="Verbose output: mainly progress reports.",
            default=False
            )
    # Set the function that gets called when the sub command is activated
    parser_cm.set_defaults(func=convert_to_xyz)
    
    # parse the arguments, and strip out the common options
    args = parser.parse_args()
    verbose = args.verbose
    # All of the rest may or may not apply to the command chosen
    # The rest of this function is effectively checking all of the arguments
    # against the command specified
    # Each command is given a function to execute via the "set_default" argument
    args.func(args)

def convert_to_xyz(args):
    """ To big to plot? Convert to the PKL file to XYZ and use GMT"""
    
    verbose = args.verbose
    subsample = args.subsample
    output_file = args.output_file[0]
    input_file =  args.input_file[0]
    resolution = args.resolution
    mask_file = args.mask_file[0]
    grid = args.grid

    if (verbose):
        print "Loading data: "+input_file
    data = pickle.load( open( input_file, "rb" ) )
    xi = data[0]
    yi = data[1]
    zi = data[2]
    axes_info = data[3]
    index_1 = axes_info[0]
    index_2 = axes_info[1]
    x_axis = axes_info[2]
    y_axis = axes_info[3]
   
    if (not mask_file == "no_mask"):
        from Scientific.IO import NetCDF
        file = NetCDF.NetCDFFile(mask_file, 'r')
        longitude = file.variables['lon'][:] 
        latitude = file.variables['lat'][:] 
        mask = file.variables['z'][:, :] 
        lon,lat = np.meshgrid(longitude, latitude)
        lon = np.reshape(lon,len(longitude)*len(latitude))
        lat = np.reshape(lat,len(longitude)*len(latitude))
        mask = np.reshape(mask,len(longitude)*len(latitude))


    # subsample data
    x = xi[0::subsample]
    y = yi[0::subsample]
    z = zi[0::subsample]
    zi = z
    yi = y
    xi = x
    z =[];x=[];y=[]

    if (grid):
        if (verbose):
            print "Creating uniform mesh"
        # Need to have different numebr of points for X and Y
        # Warning or error if n points is huge
        X, Y = np.meshgrid(np.linspace( xi.min(), xi.max(), (180./resolution)+1), np.linspace(yi.min(), yi.max(), (360./resolution)+1))
        Z = pylab.griddata(xi, yi, zi, X, Y)
        if (not mask_file == "no_mask"):
            mask_grd = pylab.griddata(lat, lon, mask, X, Y)

            # loop over mask_grd. Where land, set Z to NaN, which will mask out these areas in the plot
            for i in range(0,len(mask_grd)):
                for j in range(0,len(mask_grd[0])):
                    if (mask_grd[i][j] == 0):
                        Z[i][j] = float('NaN')

        f = open(output_file,"w")
        for i in range(0,len(X)):
            for j in range(0,len(Y[0])):
                val = Z[i][j]
                # skip NaNs - GMT fills them in for us
                if (str(val) == "--"):
                        continue
                f.write(str(Y[i][j])+" "+str(X[i][j])+" "+str(Z[i][j])+"\n")
        f.close()
    else:
        f = open(output_file,"w")
        for i in range(0,len(xi)):
            val = zi[i]
            # skip NaNs - GMT fills them in for us
            if (str(val) == "--"):
                    continue
            f.write(str(yi[i])+" "+str(xi[i])+" "+str(zi[i])+"\n")
        f.close()


def plot_slice(args):

    verbose = args.verbose
    output_file = args.output_file[0]
    input_file =  args.input_file[0]
    variable =  args.variable[0]
    cmap = args.colourmap
    dpi = args.dpi
    fmt = args.format
    minimum = args.minvalue
    maximum = args.maxvalue
    vscale = args.velocityscale
    addvelocity = None
    on_sphere = args.on_sphere
    subsample = args.subsample
    if (args.mask_file):
        mask_file = args.mask_file[0]
    else:
        mask_file = False
    hardlimit = args.hardlimit  
    resolution = args.resolution
    if (not args.addvelocity == None):
        addvelocity=args.addvelocity[0]
    if (cmap == "jet"):
        colourmap = pylab.cm.jet
    elif (cmap == "hot"):
        colourmap = pylab.cm.hot
    elif (cmap == "Blues"):
        colourmap = pylab.cm.Blues
    elif (cmap == "gray"):
        colourmap = pylab.cm.gray
    elif (cmap == "rainbow"):
        colourmap = pylab.cm.gist_rainbow_r
    else:
        colourmap = pylab.cm.jet
    
    if (fmt == None):
        fmt = output_file[-3:]
        if (not str.lower(fmt) == 'png' and 
            not str.lower(fmt) == 'pdf'):
                print "Format should be one of pdf or png"
                print "Defaulting to png"
                fmt = 'png'

    if (verbose):
        print "Loading data: "+input_file
    data = pickle.load( open( input_file, "rb" ) )
    xi = data[0]
    yi = data[1]
    zi = data[2]
    axes_info = data[3]
    index_1 = axes_info[0]
    index_2 = axes_info[1]
    x_axis = axes_info[2]
    y_axis = axes_info[3]
    plotting_velocity = False
    n_comp = len(np.shape(zi[0]))
    if n_comp < 1:
        n_comp = 1
    if (mask_file):
        from Scientific.IO import NetCDF
        file = NetCDF.NetCDFFile(mask_file, 'r')
        longitude = file.variables['lon'][:] 
        latitude = file.variables['lat'][:] 
        mask = file.variables['z'][:, :] 
        lon,lat = np.meshgrid(longitude, latitude)
        lon = np.reshape(lon,len(longitude)*len(latitude))
        lat = np.reshape(lat,len(longitude)*len(latitude))
        mask = np.reshape(mask,len(longitude)*len(latitude))

    # subsample data
    x = xi[0::subsample]
    y = yi[0::subsample]
    z = zi[0::subsample]
    zi = np.array(z)
    yi = np.array(y)
    xi = np.array(x)
    z =[];x=[];y=[]


    if (n_comp > 1):
        plotting_velocity = True

    if (verbose):
        print "Creating uniform mesh"
    # Need to have different numebr of points for X and Y
    # Warning or error if n points is huge
    X, Y = np.meshgrid(np.linspace( xi.min(), xi.max(), 360./resolution ), np.linspace(yi.min(), yi.max(), 180./resolution))

    params = {
      'legend.fontsize': 18,
      'xtick.labelsize': 16,
      'ytick.labelsize': 16,
      'font.size' : 18,
      'axes.labelsize' : 18,
      'text.fontsize' : 18,
      'figure.subplot.hspace' : 0.5
    }
    pylab.rcParams.update(params)

    fig = pylab.figure(figsize=(15,8),dpi=90)
    if (index_1 == 0 and index_2 == 1):
        ax = fig.add_subplot(111, aspect='equal')
    else:
        ax = fig.add_subplot(111)

    if (verbose):
        print "Interpolating data"
    if (plotting_velocity):
        u = zi[index_1,:]
        v = zi[index_2,:]
        U = pylab.griddata(xi, yi, u, X, Y)
        V = pylab.griddata(xi, yi, v, X, Y)
        q = pylab.quiver(X,Y,U,V,scale=vscale,color='k')
    else:
        #tri = Triangulation(xi,yi)
        #interp = tri.nn_interpolator(zi)
        #zi = interp(X,Y)
        Z = pylab.griddata(xi, yi, zi, X, Y)
        if (mask_file):
            mask_grd = pylab.griddata(lat, lon, mask, X, Y)

            # loop over mask_grd. Where land, set Z to NaN, which will mask out these areas in the plot
            for i in range(0,len(mask_grd)):
                for j in range(0,len(mask_grd[0])):
                    if (mask_grd[i][j] == 0):
                        Z[i][j] = float('NaN')

        
        # define min/max and spacing of data if not given (so we see all of the data)
        spacing = None
        if (minimum == None):
            minimum = zi.min() 
            minimum = minimum - (0.05*minimum) 
        if (maximum == None):
            maximum = zi.max()
            maximum = maximum + (0.05*maximum) 
        if (spacing == None):
            spacing = (maximum - minimum) /256.

        if (hardlimit):
            if (verbose):
                print "applying hard limit"
            for i in range(0,len(Z)):
                for j in range(0,len(Z[0])):
                    if (Z[i][j] <= minimum):
                        Z[i][j] = minimum+(0.01*minimum) 
                    if (Z[i][j] >= maximum):
                        Z[i][j] = maximum-(0.01*maximum) 

        if (verbose):
            print "Plotting"
        cs=ax.contour(Y,X,Z,np.arange(minimum,maximum,spacing),cmap=colourmap)
        #cs=pylab.scatter(yi,xi,marker='.',s=2,c=zi,cmap=colourmap,edgecolors='none')
        cs=ax.contourf(Y,X,Z,np.arange(minimum,maximum,spacing),cmap=colourmap)
        #print "triangulation"
        #triang = tri.Triangulation(xi, yi)
        #print "plotting"
        #cs=pylab.tricontour(xi, yi, zi, 15, linewidths=0.5, colors='k')
        #cs=pylab.tricontourf(xi, yi, zi, 15, cmap=plt.cm.jet)
        pp=pylab.colorbar(cs,format='%.2f')
        pp.set_label(variable)
        if (not addvelocity == None):
            data = pickle.load( open( addvelocity, "rb" ) )
            xi = data[0]
            yi = data[1]
            zi = data[2]
            axes_info = data[3]
            index_1 = axes_info[0]
            index_2 = axes_info[1]
            x_axis = axes_info[2]
            y_axis = axes_info[3]
            u = zi[:,index_1]
            v = zi[:,index_2]
            U = pylab.griddata(xi, yi, u, X, Y)
            V = pylab.griddata(xi, yi, v, X, Y)
            q = pylab.quiver(X,Y,U,V,scale=vscale,color='k',zorder=100)

    
    pylab.xlabel(x_axis+' (m)')
    pylab.ylabel(y_axis+' (m)')
    pylab.savefig(output_file, dpi=dpi,format=fmt)


def gen_slice(args):

    verbose = args.verbose
    origin = args.origin
    normal = args.normal
    output_file = args.output_file[0]
    input_file =  args.input_file[0]
    variable =  args.variable[0]
    on_sphere = args.sphere
    sphere_radius = args.sphere_radius
    depth = args.sphere_radius

    xi = []
    yi = []
    zi = []
    vtu_files = []
    path = os.path.dirname(input_file)
    
    if (not on_sphere):
        if (normal == (0,0,1)):
            index_1 = 0
            index_2 = 1
            x_axis = "X"
            y_axis = "Y"
        elif(normal == (0,1,0)):
            index_1 = 0
            index_2 = 2
            x_axis = "X"
            y_axis = "Depth"
        elif(normal == (1,0,0)):
            index_1 = 1
            index_2 = 2
            x_axis = "Y"
            y_axis = "Depth"
        else:
            print "Normal should be one of (1,0,0), (0,1,0), or (0,0,1)"
            sys.exit(-1)
        # check origin
        if (origin == None):
            print "Origin should be a set of three corrdinate, eg. (0,0,-10)"
            sys.exit(-1)
    else:
        import fluidity.spheretools as st
        index_1 = 0
        index_2 = 1
        index_3 = 2
        x_axis = "Longitude"
        y_axis = "Latitude"


    # get vtus from pvtu file. They are in lines such as:
    # <Piece Source="restratA-np64-adapt-C_99/restratA-np64-adapt-C_99_0.vtu" />
    # As PVTU is XML, parse as such and grab the Piece elements
    xml_root = etree.parse(input_file,xml_parser)
    find = etree.XPath("//Piece")
    peices = find(xml_root)
    for p in peices:
        name = p.attrib['Source']
        vtu_files.append(os.path.join(path,name))
      
    if (on_sphere):
        sphere = vtk.vtkSphere()
        sphere.SetCenter(0,0,0)
        sphere.SetRadius(sphere_radius-10)
    else:
        plane = vtk.vtkPlane()
        plane.SetOrigin(origin[0], origin[1], origin[2])
        plane.SetNormal(normal[0],normal[1],normal[2])


    getMagnitude = False
    component = -1
    # Special cases of variables
    if (variable == "VelocityMagnitude"):
        getMagnitude = True
        variable = "Velocity_projection"
    if (variable == "VelocityU"):
        getMagnitude = False
        variable = "ProjectedVelocity"
        component = 0
    if (variable == "VelocityV"):
        getMagnitude = False
        variable = "ProjectedVelocity"
        component = 1   
    if (variable == "BedShearStressMagnitude"):
        getMagnitude = True
        variable = "BedShearStress_projection"
    if (variable == "MaxBedShearStressMagnitude"):
        getMagnitude = True
        variable = "MaxBedShearStress_projection"
    if (variable == "AveBedShearStressMagnitude"):
        getMagnitude = True
        variable = "AveBedShearStress_projection"
    if (variable == "MaxVelocityMagnitude"):
        getMagnitude = True
        variable = "MaxVelocity_projection"
    if (variable == "AveVelocityMagnitude"):
        getMagnitude = True
        variable = "AveVelocity_projection"
    if (variable == "AveVelocityU"):
        getMagnitude = False
        component = 0
        variable = "AveVelocity_projection"
    if (variable == "AveVelocityV"):
        getMagnitude = False
        component = 1
        variable = "AveVelocity_projection"
    if (variable == "MaxVelocityU"):
        getMagnitude = False
        component = 0
        variable = "MaxVelocity_projection"
    if (variable == "MaxVelocityV"):
        getMagnitude = False
        component = 1
        variable = "MaxVelocity_projection"
    if (variable == "AveBedShearStressU"):
        getMagnitude = False
        component = 0
        variable = "AveBedShearStress_projection"
    if (variable == "AveBedShearStressV"):
        getMagnitude = False
        component = 1
        variable = "AveBedShearStress_projection"
    if (variable == "MaxBedShearStressU"):
        getMagnitude = False
        component = 0
        variable = "MaxBedShearStress_projection"
    if (variable == "MaxBedShearStressV"):
        getMagnitude = False
        component = 1
        variable = "MaxBedShearStress_projection"
    
    comp_1 = []
    comp_2 = []
    comp_3 = []
    # the above are for later...numpy.array[[list],[list]] is sllloooowww

    #for each vtu, grab any slice data
    i = 0
    for v in vtu_files:
        if (verbose):
            print v
        # grab the number
        num = i

        # Start by loading some data.
        reader=vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(v)
        reader.Update()
        grid=reader.GetOutput()
        grid.GetPointData().SetActiveScalars(variable)

        cut = vtk.vtkCutter()
        if (on_sphere):
            cut.SetCutFunction(sphere)
        else:
            cut.SetCutFunction(plane)
        cut.SetInputConnection(reader.GetOutputPort())
        cut.Update()

        i = i+1
        
        d = cut.GetOutput().GetPointData().GetScalars()
        d.GetNumberOfTuples()
        slice_data = np.array([d.GetTuple(i) for i in range(d.GetNumberOfTuples())])
        c = cut.GetOutput()
        coords = []
        for i in range(0,c.GetNumberOfPoints()):
            coords.append(c.GetPoint(i))
        coords = np.array(coords)

        if (len(slice_data) > 0):
            n_comps = len(slice_data[0])
            if (on_sphere):
                if (n_comps  == 3):


                    # project vectors and coords
                    for i in range(0,len(coords[:,0])):
                        #print coords[i,0],coords[i,1],coords[i,2]
                        #print slice_data[i,0],slice_data[i,1],slice_data[i,2]
                        [[radial, polar, azimuthal],[lat,lon]] = \
                            sp.vector_cartesian_2_spherical_polar(slice_data[i,0],slice_data[i,1],slice_data[i,2],
                                    coords[i,0],coords[i,1],coords[i,2])
                        xi.append(math.degrees(lon))
                        yi.append(math.degrees(lat))
                        slice_data[i,0] = radial
                        slice_data[i,1] = polar
                        slice_data[i,2] = azimuthal
                elif (n_comps == 1):
                    #project points only
                    for i in range(0,len(coords[:,0])):
                        lat,lon = st.cart2polar(coords[i,:])
                        xi.append(math.degrees(lon))
                        yi.append(math.degrees(lat))
                else:
                    print("Unknown number of componenets - not a 3D vector or scalar")
                    sys.exit(-1)
            else:
                xi.extend(coords[:,index_1])
                yi.extend(coords[:,index_2])
            
            if (getMagnitude):
                # How to calculate the magnitude nicely..? Below seems rather brute force.
                mag = []
                mag = slice_data[:,0]
                for j in range(0,len(mag)):
                        mag[j] = mag[j]*mag[j]
                for i in range(1,n_comps):
                    for j in range(0,len(mag)):
                        mag[j] = mag[j] + slice_data[j,i]*slice_data[j,i]
                for j in range(0,len(mag)):
                        mag[j] =math.sqrt(mag[j])
                zi.extend(mag)
            elif (component > -1):
                zi.extend(slice_data[:,component])
            else:
                if (n_comps == 3):
                    comp_1.extend(slice_data[:,0])
                    comp_2.extend(slice_data[:,1])
                    comp_3.extend(slice_data[:,2])
                else:
                    zi.extend(slice_data[:,0])


    # THe data contain > 1 components, but we only want one
    if (component > -1 or getMagnitude) :
        n_comps = 1
        #if (component == 1):
        #    zi = comp_1
        #elif (component == 2):
        #    zi = comp_2
        #elif (component == 3):
        #    zi = comp_3
        #else:
        #    print "Error - you seem to be generating 4D data!"
        #    sys.exit(-1)
    coords = []
    # need to tidy up these arrays - lot's duplicate points from halos, DG, or whatever
    if (n_comps == 3):
        for i in range(len(xi)):
            coords.append([xi[i],yi[i],comp_1[i],comp_2[i],comp_3[i]])
    elif (n_comps == 1):
        for i in range(len(xi)):
            coords.append([xi[i],yi[i],zi[i]])
    else:
        print "Error. I expect 1 or 3 components"
        sys.exit(-1)

    if (verbose):
        print "Removing duplicated points"
    coords_checked = unique(coords)
    
    if (verbose):
        print "Saving file"

    xi = []
    yi = []
    zi = []
    for c in coords_checked:
        xi.append(c[0])
        yi.append(c[1])
        if (n_comps == 3):
            zi.append([c[2], c[3], c[4]])
        else:
            zi.append(c[2])
        i += 1

    # Store data as a pickle
    data = [xi,yi,zi]
    # Add axes labels
    axes_info = [index_1,index_2,x_axis,y_axis]
    data.append(axes_info)
    pickle.dump( data, open( output_file, "wb" ) )

## {{{ http://code.activestate.com/recipes/52560/ (r1)
def unique(s):
    """Return a list of the elements in s, but without duplicates.

    For example, unique([1,2,3,1,2,3]) is some permutation of [1,2,3],
    unique("abcabc") some permutation of ["a", "b", "c"], and
    unique(([1, 2], [2, 3], [1, 2])) some permutation of
    [[2, 3], [1, 2]].

    For best speed, all sequence elements should be hashable.  Then
    unique() will usually work in linear time.

    If not possible, the sequence elements should enjoy a total
    ordering, and if list(s).sort() doesn't raise TypeError it's
    assumed that they do enjoy a total ordering.  Then unique() will
    usually work in O(N*log2(N)) time.

    If that's not possible either, the sequence elements must support
    equality-testing.  Then unique() will usually work in quadratic
    time.
    """

    n = len(s)
    if n == 0:
        return []

    # Try using a dict first, as that's the fastest and will usually
    # work.  If it doesn't work, it will usually fail quickly, so it
    # usually doesn't cost much to *try* it.  It requires that all the
    # sequence elements be hashable, and support equality comparison.
    u = {}
    try:
        for x in s:
            u[x] = 1
    except TypeError:
        del u  # move on to the next method
    else:
        return u.keys()

    # We can't hash all the elements.  Second fastest is to sort,
    # which brings the equal elements together; then duplicates are
    # easy to weed out in a single pass.
    # NOTE:  Python's list.sort() was designed to be efficient in the
    # presence of many duplicate elements.  This isn't true of all
    # sort functions in all languages or libraries, so this approach
    # is more effective in Python than it may be elsewhere.
    try:
        t = list(s)
        t.sort()
    except (TypeError, ValueError):
        del t  # move on to the next method
    else:
        assert n > 0
        last = t[0]
        lasti = i = 1
        while i < n:
            if t[i] != last:
                t[lasti] = last = t[i]
                lasti += 1
            i += 1
        return t[:lasti]

    # Brute force is all that's left.
    u = []
    print "Warning: brute force method needed for uniqueness. This may take a while..."
    for x in s:
        if x not in u:
            u.append(x)
    return u
## end of http://code.activestate.com/recipes/52560/ }}}


if __name__ == "__main__":
    main()
