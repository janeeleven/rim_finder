### HiRISE Rim Finder Tool ####
### by Janette Levin, 2020 ####
### Python 2.7 ####
### uses profiles generated from CTX ###
###change obs ID and center Lat at the bottom of this script


import arcpy
from arcpy import env
from arcpy.ddd import *
from arcpy.sa import *
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import errno

arcpy.CheckOutExtension("3D")
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True

#ctx directory
root_dir = "X:/User/jlevin/ctx/alcove_spectral_analysis/"

#hirise directory
work_dir = "X:/User/jlevin/hirise/groundtruthing/"

def mkdir_p(path):
    """
    Allows mkdir -p functionality : make function and skip it if it already exists
    Input:
    - path: str
    """
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def getdegrees(datapnt,root,output_dir):
    """
    Extract data from .shp file
    Input:
    - datapnt: str, name of shapefile with interpolated points
    - root: str, name of file
    - output_dir: str, path
    
    """
    X_coords = []
    Y_coords = []
    Z_coords = []
    #open new file
    f = open(output_dir + "data_" + root + ".txt","w")
    f.write("X" + " Y" + " Z" + "\n")
    #get arcpy SearchCursor
    s_cur = arcpy.da.SearchCursor(datapnt, ['SHAPE@'])
    for row in s_cur:
        polyline = row[0]
        for point in polyline:
            X_coords.append(point.X)
            Y_coords.append(point.Y)
            Z_coords.append(point.Z)
    print("Unpacked")
    if len(X_coords) == 0:
        return
    else: 
        minloc = np.argmin(np.asarray(X_coords)) #find index of minimum x point (west-most point)
        i = 0
        curloc = minloc
        xcorr = np.zeros(len(X_coords)) #empty array to store X_coords
        ycorr = np.zeros(len(X_coords))
        zcorr = np.zeros(len(X_coords))
        
        #reorganizes the data to go from west counterclockwise (-pi to pi on a circle)
        while i < len(Y_coords):
            xcorr[i] = X_coords[curloc]
            ycorr[i] = Y_coords[curloc]
            zcorr[i] = Z_coords[curloc]
            curloc += 1
            #print( str(xcorr[i]) + " " + str(ycorr[i]) + " " + str(zcorr[i]))
            
            #loop around if necessary
            if curloc == len(Y_coords):
                curloc = 0
            i += 1
        try:
             for j in range(0,len(Z_coords)):
                #print(str(X_coords[j]) + " " + str(Y_coords[j]) + " " + str(Z_coords[j]))
                f.write(str(xcorr[j]) + " " + str(ycorr[j]) + " " + str(zcorr[j]) + "\n")
             f.close()
        except:
            print "Something went wrong with " + root + ". Check if it has data. Skipping it." 
            f.close()
            pass

def polyfind(obs,lat):
    """
    Find crater rim of given crater
    Input:
    - obs: str, Crater ID
    - lat: float
    """
    obs_dir = work_dir + obs + "/"
    obs_shp = work_dir + obs + "/shp/"
    obs_dir_dat = root_dir + obs + "/"
    input_dir = obs_dir + "shp/"
    alcv_dir = obs_dir_dat + "shp/alcv/"
    output_dir = obs_dir + "data/"

    #make folders if they do not yet exist
    
    mkdir_p(obs_shp)
    mkdir_p(output_dir)
    
    #make sure the initial crater shp file is in the same projection as the DEM
    inputdem = glob.glob(obs_dir + "hirise_files/*1.tif")[0].replace('\\','/')
    shppath = obs_dir_dat + "shp/clipsize.shp"
    refproj = arcpy.Describe(inputdem).spatialReference

    arcpy.Project_management(shppath,input_dir + "bufferpoly.shp", refproj)

    # Get extent of feature

    clipout = input_dir + "bufferpoly.shp"
    desc = arcpy.Describe(clipout)
    xmin = desc.extent.XMin
    print(xmin)
    xmax = desc.extent.XMax
    print(xmax)
    ymin = desc.extent.YMin
    print(ymin)
    ymax = desc.extent.YMax
    print(ymax)
    
    #save extent coords of original guess
    k = open(output_dir + "data_coords.txt","w")
    k.write("xmin ymin xmax ymax lat" + "\n")
    k.write(str(xmin) + " " + str(ymin) + " " + str(xmax) + " " + str(ymax) + " " + str(lat))
    k.close()
    
    for i in glob.glob(alcv_dir + "alcovecov*.shp"):
        i = i.replace('\\','.')
        root = i.split(".")[1]
        outputline = alcv_dir + "line_" + root + ".shp"
        outputvar = obs_shp + "3d_" + root + ".shp"
        
        #interpolate using hirise dem                  
        arcpy.InterpolateShape_3d(inputdem, outputline, outputvar, sample_distance="10", z_factor="1", method="BILINEAR", vertices_only="DENSIFY", pyramid_level_resolution="0")
        arcpy.ClearWorkspaceCache_management()
        outputpnt = obs_shp + "3dpnt_" + root
        #convert vertices to points
        arcpy.FeatureVerticesToPoints_management(outputvar,outputpnt,"ALL")
        arcpy.ClearWorkspaceCache_management()
        datapnt = outputpnt + ".shp"

        
        getdegrees(datapnt,root,output_dir)
    
    print "Contours extracted from  " + obs
    return
obs = '20-001585'
lat = '-23.124'
polyfind(obs,lat)
