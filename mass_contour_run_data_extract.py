#### By Janette Levin, 2020 #####
#### Python 3.8 ########


import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import errno
import scipy.fftpack
from scipy import signal

##functions imported from###
from rms_wav_from_profile import *

root_dir = "X:/User/jlevin/ctx/alcove_spectral_analysis/"

#set up crater class
class crater(object):
    def __init__(self,craterid):
        self.craterid = craterid
        self.rockyn = 0
        self.rockys = 0
    
    def northlen(self, rockyn):
        self.rockyn = rockyn
    def southlen(self,rockys):
        self.rockys = rockys



#data_zone.txt contains observation ID, rocky zone length, and N or S
        
follist = open(root_dir + "data_zone.txt")
next(follist) #skip header

cratlist = {} #dictionary of crater objects sorted by ID
for line in follist:
    fields = line.split(" ")
    obsid = fields[0]
    cardir = fields[2]
    #if crater object, add the other side to the object
    if obsid in cratlist:
        print("This exists already! " + obsid)

        if cardir == "N":
            cratlist[obsid].northlen(int(fields[1]))
        else:
            cratlist[obsid].southlen(int(fields[1]))
    #make new crater object if other side does not exist
    else:
        cratlist[obsid] = crater(obsid)
        print("Made new one")
        if cardir == "N":
            cratlist[obsid].northlen(int(fields[1]))
        else:
            cratlist[obsid].southlen(int(fields[1]))
follist.close()



for crat in cratlist:
    output_dir = root_dir + crat  + "/data/" 
    print("Processing " + crat)
    try:
        lat = center_coords(output_dir) #get latitude of center

        ws = 600 #window size
        
        #list of profiles within the rocky zone being processed
        
        fillist = []
        nbond = cratlist[crat].rockyn
        sbond = cratlist[crat].rockys
        bigger = max(nbond,sbond,375) #max rocky zone length
        for k in range(50,bigger + 1,25):
            filnam = "data_alcovecov_" + str(k) + ".txt"
            fillist.append(filnam)


        topo = plt.figure(figsize = [16,6])
        
        zdatlist = np.empty((360,len(fillist)), dtype=object) #empty array to store all z data
        distlist = np.empty((360,len(fillist)), dtype=object) #empty array to store all dist data
        radlist = np.empty((360,len(fillist)), dtype=object) #empty array to store all rad data
        avglist = np.empty((360,len(fillist))) #empty array to store averages for filtering out crossing lines

        #spread all files into 360 bins for comparing
        for file in fillist: 
            indx = fillist.index(file)
            dat = open(output_dir + file)
            rad,z,dist = extract(dat) # extract XYZ values and convert to radial coordinates and along-profile distance
    
            #split into bin by theta
            binlocs = []
            bins = np.linspace(-np.pi,np.pi,361)
            
            #find indeces that break the theta values (rad) up by bin
            for j in bins:
                dist_2 = (rad - j)**2
                binlocs.append(np.argmin(dist_2))

            #average values in each bin)
            binned = np.empty(360)
            for k in range(0,360):
                zdatlist[k,indx] = z[binlocs[k]:binlocs[k+1]]
                distlist[k,indx] = dist[binlocs[k]:binlocs[k+1]]
                radlist[k,indx] = rad[binlocs[k]:binlocs[k+1]]
                try:
                    #average z value for each bin
                    avgval = np.mean(z[binlocs[k]:binlocs[k+1]])
                except RuntimeWarning:
                    #if there are NaN values, exclude whole bin
                    print('hole')
                    avgval = float('NaN')
                binned[k] = avgval
            avglist[:,indx] = binned
        #filter out any particulalry corrupt part (within the Region Of Interest)
        
        if sbond > 0:
            #find concentric profile that bounds the south wall
            southbound = fillist.index("data_alcovecov_" + str(sbond) + ".txt")
        else:
            southbound = 0
        if nbond > 0:
            #find the concentric profile that bounds the north wall
            northbound = fillist.index("data_alcovecov_" + str(nbond) + ".txt")
        else:
            northbound = 0

        #check average difference down the slope for each bin
        #filter out bins where profiles are too far apart (polygonization) or where the average > 0 (profiles cross)

        for k in range(360):
            meandif = np.mean(np.diff(avglist[k,:])) #average difference for each bin
            if np.isnan(meandif) == False and meandif > -25.0 and meandif < 0.0:
                continue
            else:
                for i in range(len(fillist)):
                    zdatlist[k,i] = np.array((zdatlist[k,i]))*np.nan
                    radlist[k,i] = np.array(radlist[k,i])*np.nan
                    radlist[k,i] = np.array(radlist[k,i])*np.nan



        pdftnorth = []
        pdftsouth = []

        wavnorth = []
        wavsouth = []

        #output datafiles
        rms_dat_n = open(output_dir + "rms_padfix_data_n.txt","w")
        rms_dat_s = open(output_dir + "rms_padfix_data_s.txt","w")
        wav_dat_n = open(output_dir + "wav_padfix_data_n.txt","w")
        wav_dat_s = open(output_dir + "wav_padfix_data_s.txt","w")


        for j in range(len(fillist)):
            flatz = [] #rearranging data to make it easier to keep track
            flatrad = []
            flatdist = []
            for sublist in zdatlist[:,j]:
                for item in sublist:
                    flatz.append(item)
            for sublist in radlist[:,j]:
                for item in sublist:
                    flatrad.append(item)
            for sublist in distlist[:,j]:
                for item in sublist:
                    flatdist.append(item)
            plt.figure(topo.number)
            plt.plot(flatrad,flatz)
            print('file ' + str(j))
            
            #get measurements
            maxess,rmsls,pdfts,maxesn,rmsln,pdftn = perimeter(flatrad,flatz,flatdist,ws)
            
            #filters out rms values above 2std in the north section
            filtn = np.mean(pdftn) + 2*np.std(pdftn)
            indcn =  [i for i in range(len(pdftn)) if pdftn[i] < filtn]

            pdftnf = [pdftn[i] for i in indcn]
            maxesnf = [maxesn[i] for i in indcn]

            #filters out rms values above 2std in the south section
            filts = np.mean(pdfts) + 2*np.std(pdfts)
            indcs =  [i for i in range(len(pdfts)) if pdfts[i] < filts]

            pdftsf = [pdfts[i] for i in indcs]
            maxessf = [maxess[i] for i in indcs]
            lim = len(flatrad)//6 #limit for points per third - if there is less than this in a profile, ignore that contour

            #make sure the profile is not mostly holes enough
            if len(maxesn) > lim:
                wav_dat_n.write(' '.join(str (e) for e in maxesnf) + '\n') #save data for each contour on individual line
                rms_dat_n.write(' '.join(str (e) for e in pdftnf) + '\n')
                wavnorth.append(np.mean(maxesnf))
                pdftnorth.append(np.mean(pdftnf))
            if len(maxess) > lim:
                wav_dat_s.write(' '.join(str (e) for e in maxessf) + '\n')
                rms_dat_s.write(' '.join(str (e) for e in pdftsf) + '\n')
                wavsouth.append(np.mean(maxessf))
                pdftsouth.append(np.mean(pdftsf))
            plt.figure(topo.number)
            plt.plot(rad,z,label = file.split("_")[2].split(".")[0])
            plt.title("Topography profile")
            plt.xlabel('position (m)')
            plt.ylabel('elevation')
        #figure formating
        plt.figure(topo.number)
        plt.xlabel('Position(radians)')
        plt.ylabel('elevation')
        plt.legend(loc = "lower right")
        plt.close()


        topo.savefig(output_dir + 'denoised_topo_plot.png',dpi=150)

        rms_dat_n.close()
        rms_dat_s.close()
        wav_dat_n.close()
        wav_dat_s.close()
        rmsplot = plt.figure(figsize=[6,6])
        if lat > 0:
            if northbound > 0:
                plt.plot(wavnorth[:northbound +1], pdftnorth[:northbound +1],'.',label = 'north (equator facing)')
            if southbound > 0:
                plt.plot(wavsouth[:southbound +1],pdftsouth[:southbound +1],'.',label = 'south (pole facing)')
        else:
            if northbound > 0:
                plt.plot(wavsouth[:southbound +1],pdftsouth[:southbound +1],'.',label = 'south (equator facing)')
            if southbound > 0:
                plt.plot(wavnorth[:northbound +1], pdftnorth[:northbound +1],'.',label = 'north (pole facing)')

        plt.legend()
        plt.xlabel('average wavelength (m)')
        plt.ylabel('amplitude (m)')
        plt.title('Average RMS values')
        rmsplot.savefig(output_dir + 'rms_vals.png',dpi=200)
        plt.close()
        f = open(output_dir + "filt_rms_wav.txt","w")
        f.write("LAT" + "RZ_SOUTH" + " AVG_RMS_SOUTH" + " AVG_WAV_SOUTH" + "RZ_NORTH" + " AVG_RMS_NORTH" + " AVG_WAV_NORTH" + "\n")
        f.write(str(lat) + " " + str(sbond) + " " + str(np.mean(pdftsouth[:southbound +1])) + " " + str(np.mean(wavsouth[:southbound +1])) + " "+str(nbond) + " "+ str(np.mean(pdftnorth[:northbound +1])) + " " + str(np.mean(wavnorth[:northbound +1])))
        f.close()

    #some possible errors indicating polygonization trouble
    except ValueError:
        print("Incomplete Polygonization")
    #    echostr = 'echo ' + output_dir + ' >> incomplete_poly3.txt'
        #os.system(echostr)
        #incpol.append(output_dir.split("/")[5])
    except IndexError:
        print("Incomplete Polygonization")
    #    echostr = 'echo ' + output_dir + ' >> incomplete_poly3.txt'
        #os.system(echostr)

    except FileNotFoundError:
        print("File Not Found")


