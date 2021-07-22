### Rim Finder Tool ####
### by Janette Levin, 2020 ####
### Python 2.7 ####

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

root_dir = "X:/User/jlevin/ctx/alcove_spectral_analysis/"

   
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

def pathfind(obs):
    """
    Find corresponding DEM and shapefile path for a given crater id using initial_crater_data.shp
    Input:
    - obs: str, Crater ID
    Output:
    - shppath: str, path to initial circle guess
    - pathtodem: str, path to DEM
    """
    rows = arcpy.SearchCursor("X:/User/jlevin/ctx/initial_crater_data.shp", fields="CRATER_ID;shpName;PATH")
    for row in rows:
        cratid = row.getValue("CRATER_ID")
        if cratid == obs:
            shppath = row.getValue("shpName")
            pathtodem = row.getValue("PATH")
            break
    print pathtodem
    return shppath,pathtodem

def polyfind(obs):
    """
    Find crater rim of given crater
    Input:
    - obs: str, Crater ID
    """
    obs_dir = root_dir + obs + "/"
    input_dir = obs_dir + "shp/"
    alcv_dir = input_dir + "alcv/"
    output_dir = obs_dir + "data/"

    #make folders if they do not yet exist
    mkdir_p(alcv_dir)
    mkdir_p(flow_dir)
    mkdir_p(output_dir)
    
    # find necessary shapefile and DEM 
    try:
        shppath, pathtodem = pathfind(obs)
    except NameError:
        print("Crater ID not found")
    
    inputclass = shppath
    clipout = input_dir + "clipsize.shp"

    
    # Buffer to create clipping mask: change to higher buffer value for very irregular craters
    arcpy.Buffer_analysis(inputclass, clipout , "500 Meters", "FULL", "ROUND", "NONE", "","PLANAR")
    arcpy.ClearWorkspaceCache_management()
    print("Clip mask created")

    # Get extent of feature
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
    k.write("xmin ymin xmax ymax" + "\n")
    k.write(str(xmin) + " " + str(ymin) + " " + str(xmax) + " " + str(ymax))
    k.close()
    
    #clip DEM using feature extent
    arcpy.Clip_management(pathtodem, str(xmin) + " " + str(ymin) + " " + str(xmax) + " " + str(ymax) ,obs_dir + "clipdem.tif",clipout,"-3.276700e+04","ClippingGeometry","NO_MAINTAIN_EXTENT")
    arcpy.ClearWorkspaceCache_management()

    #resample at lower resolution
    arcpy.Resample_management(obs_dir + "clipdem.tif",obs_dir + "clipdemlow.tif", cell_size="40 40", resampling_type="CUBIC")
    #mosaic lower resolution under higher resolution to fill in holes
    arcpy.MosaicToNewRaster_management(obs_dir + "clipdem.tif;" + obs_dir + "clipdemlow.tif", obs_dir, "clipdemfuse.tif", coordinate_system_for_the_raster="", pixel_type="32_BIT_FLOAT", cellsize="", number_of_bands="1", mosaic_method="FIRST", mosaic_colormap_mode="FIRST")



    
    # slope map of clipped higher res dem
    arcpy.Slope_3d(obs_dir + "clipdem.tif", obs_dir + "clipslope.tif","DEGREE","1","PLANAR","METER")
    arcpy.ClearWorkspaceCache_management()
    print("Slope Map Created")
    
   
    #isolate values between 29 and 65 degrees (top of crater rim) 
    arcpy.Reclassify_3d(obs_dir + "clipslope.tif", "VALUE", "29 65 1", obs_dir + "reclassslope.tif", "NODATA")
    arcpy.ClearWorkspaceCache_management()
    
    # polygonize raster
    arcpy.RasterToPolygon_conversion(obs_dir + "reclassslope.tif", input_dir + "polygonized","SIMPLIFY","VALUE","SINGLE_OUTER_PART", "")
    arcpy.ClearWorkspaceCache_management()
    print("polygonization complete")
    
    polygon_lyr = input_dir + "polygonized.shp"
    
    # calculate area of each polyon
    arcpy.AddGeometryAttributes_management(polygon_lyr, "AREA")
    arcpy.MakeFeatureLayer_management(input_dir + "polygonized.shp", 'polygon_lyr')
    
    #select polygons intersecting initial circle
    contained = arcpy.SelectLayerByLocation_management('polygon_lyr', "INTERSECT", inputclass)
    
    # select only polygons above a certain size
    larger = arcpy.SelectLayerByAttribute_management(contained, 'SUBSET_SELECTION', '"POLY_AREA" > 1000')
    arcpy.CopyFeatures_management(larger, input_dir + "polycontained")
    arcpy.ClearWorkspaceCache_management()
    print("Contained Features Selected")

    
    # aggregate polygons and smooth the aggregate
    arcpy.AggregatePolygons_cartography(input_dir + "polycontained.shp", input_dir + "aggregate", "100 Meters","100000 SquareMeters","10000000 SquareMeters","NON_ORTHOGONAL","", input_dir + "aggregate_Tbl")
    arcpy.ClearWorkspaceCache_management()
    arcpy.SmoothPolygon_cartography(input_dir + "aggregate.shp", input_dir + "smooth_aggregate","PAEK", "300 Meters","FIXED_ENDPOINT","NO_CHECK","")
    arcpy.ClearWorkspaceCache_management()
    
    # buffer again
    arcpy.Buffer_analysis(input_dir + "smooth_aggregate.shp", input_dir + "smooth_buff_agg" , "100 Meters", "FULL", "ROUND", "ALL", "","PLANAR")
    arcpy.ClearWorkspaceCache_management()
    
    # eliminate internal holes
    arcpy.EliminatePolygonPart_management(input_dir + "smooth_buff_agg.shp", input_dir + "fullcov_buff", "AREA", "100000000 SquareMeters","0", "CONTAINED_ONLY")
    arcpy.ClearWorkspaceCache_management()
    
    #smooth again
    arcpy.SmoothPolygon_cartography(input_dir + "fullcov_buff.shp", alcv_dir + "fullcov","PAEK", "500 Meters","FIXED_ENDPOINT","NO_CHECK","")
    arcpy.ClearWorkspaceCache_management()

    for k in range(50,1001,25):
        alcnam = "alcovecov_" + str(k)
        meterstr = str(-100 - k) + " Meters"
        arcpy.Buffer_analysis(alcv_dir + "fullcov.shp", alcv_dir + alcnam , meterstr, "FULL", "ROUND", "NONE", "","PLANAR")
        arcpy.ClearWorkspaceCache_management()

    print("Concentric Profile Polygons Created")
    
    for i in glob.glob(alcv_dir + "*.shp"):
        #use hole filled DEM
        inputdem = obs_dir + "clipdemfuse.tif"
        print(i)
        i = i.replace('\\','.')
        root = i.split(".")[1]
        inputclass = alcv_dir + root + ".shp"
        outputline = alcv_dir + root + "_line.shp"

        #convert polygon to line using its outline
        arcpy.PolygonToLine_management(inputclass, outputline,"IGNORE_NEIGHBORS")
        arcpy.ClearWorkspaceCache_management() 
        outputvar = alcv_dir + root + "_3d.shp"                    
        print(inputclass)
        
        # Add Vertices to line
        arcpy.Densify_edit(outputline,"DISTANCE","10 Meters","0.1 Meters","10")
        arcpy.ClearWorkspaceCache_management()
        #interpolate                  
        arcpy.InterpolateShape_3d(inputdem, outputline, outputvar, sample_distance="10", z_factor="1", method="BILINEAR", vertices_only="DENSIFY", pyramid_level_resolution="0")
        arcpy.ClearWorkspaceCache_management()

        #convert to points
        outputpnt = alcv_dir + root + "_3dpnt"
        arcpy.FeatureVerticesToPoints_management(outputvar,outputpnt,"ALL")
        arcpy.ClearWorkspaceCache_management()
        datapnt = outputpnt + ".shp"

        #export points to textfile
        getdegrees(datapnt,root,output_dir)
        print "Done with " + root
    
    print "Finished with " + obs

obs = raw_input("Crater ID:")
polyfind(obs)
