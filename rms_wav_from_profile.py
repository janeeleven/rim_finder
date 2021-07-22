#### By Janette Levin, 2020 #####
#### Python 3.8 ########
### Functions for extracting peak wavelengths and rms values ####



import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import scipy.fftpack
from scipy import signal


def center_coords(root_dir):
    """
    Load extents of the clipped DEM from data_coords.txt
    Input:
    - root_dir: str, location of observation data
    Output:
    - Latitude of center: float
    """
    global xcentglob
    global ycentglob
    #get global coords
    fc = open(root_dir + "data_coords.txt")

    next(fc)
    for line in fc:
        coords = line.split(" ")
        xmin = coords[0]
        ymin = coords[1]
        xmax = coords[2]
        ymax = coords[3]
        lat = coords[4]
        rad = (float(xmax) - float(xmin))/2
        xcentglob = float(xmax) - rad #approximate center coordinates by subtracting radius
        ycentglob = float(ymax) - rad
    fc.close()
    return float(lat)

def file_len(fname):
    """
    Length of filename
    """
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
        return i + 1

def raduniv(dato,pnts):
    """
    Convert XY coordinates to radial (theta, r) coordinates and calculate distance of each point along the line
    Input:
    - dato: np array (pnts,3) of x,y,z values
    - pnts: number of points
    Output:
    - rad: np array, theta values for each point
    - z: np array, z values for each point
    - dist: np array, distance of each point along profile
    
    """
    datu,ind = np.unique(dato,axis=0,return_index = True)#eliminate duplicate points
    
    dat = datu[np.argsort(ind)] #sort back into proper order
    
    xvals = dat[0:-2,0] - xcentglob #turn into coordinates relative to the center point
    yvals = dat[0:-2,1] - ycentglob
    
    theta = np.arctan2(yvals,xvals) #calculate theta values from relative points
    
    dist = np.zeros(np.shape(yvals))
    #calculate distance along line
    for n in range(1,len(dist)):
        dist[n] = np.sqrt((xvals[n] - xvals[n-1])**2 + (yvals[n] - yvals[n-1])**2 ) + dist[n-1]

    zvals = np.zeros(len(dat[0:-2,2]))
    zvals[:] = dat[0:-2,2]
    theta, zvals = zip(*sorted(zip(theta,zvals)))
    return theta, zvals,dist

def extract(fh):
    """
    Extract XYZ values and convert
    Input:
    - fh: filehandle
    Output:
    - rad: np array, theta values for each point
    - z: np array, z values for each point
    - dist: np array, distance of each point along profile
    """
    next(fh) #skip header
    totalpnt= file_len(fh.name)+ 1 
    print("Total points: " + str(totalpnt))
    
    dataarray = np.zeros((totalpnt, 3))
    linenum = 0
    for line in fh:
        fields = line.split(" ")
        X = fields[0]
        Y = fields[1]
        Z = fields[2]
        
        dataarray[linenum,0] = X
        dataarray[linenum,1] = Y
        dataarray[linenum,2] = Z
        linenum = linenum + 1
    fh.close()
    
    rad, z,dist = raduniv(dataarray, totalpnt)
    return rad,z,dist

## data analysis functions
def linearize(rad,z):
    """
    Detrends z by subtracting best fit line
    Input:
    - rad: np array, theta values for each point
    - z: np array, z values for each point
    Output:
    - radc: np array, linearized theta
    - z: np array, linearized z
    
    """
    radc = (rad - rad[0])
    zc = np.polyfit(radc,z,1)
    zf1 = np.poly1d(zc)
    zn = z - zf1(radc)
    return radc,zn

def butterworth(radc,z):
    """
    Use butterworth filter to filter out features with wavelength > 400
    Input:
    - radc: np array
    - z: np array
    Output:
    - zfilt: np array
    """
    
    dr = np.mean(np.diff(radc)) # take mean spacing
    fs = 1.0/dr #sampling rate
    nyq = 0.5 * fs #nyquist frequency
    cutoff = 1.0/400.
    normal_cutoff = cutoff/nyq
    order = 3 
    b,a = signal.butter(order, normal_cutoff, btype='high', analog=False)
    zfilt = signal.filtfilt(b, a, z)
    return zfilt

def fouriertrans(radc,zn):
    """
    Fourier transformation
    Input:
    - radc: np array
    - zn: np array
    Output:
    - rf: np array, r in frequency space
    - zf_p: np array, z in frequency space
    """
    N = len(radc) #number of sample pts
    dr = np.mean(np.diff(radc)) # spacing
    
    #### Fourier Transform
    zf = scipy.fftpack.fft(zn)
    
    #Take the domain we are interested in 
    rf = np.linspace(0.0,1.0/(2.0*dr),N//2)
    zf_p = np.abs(zf[:int(N//2)])
    
    return rf, zf_p

def hannwindow(radc,zn,L):
    """
    Use Hann Window to normalize window
    Input:
    - radc: np array
    - zn: np array
    - L: float, length of array
    Output:
    - hwin: np array, z values multiplied by hann window
    - win: lambda function
    """
    win = lambda x: np.sqrt(2.0/3.0)*(1.0 - np.cos(2.0*np.pi*x/(L)))
    hwin = zn*win(radc)
    return hwin,win

def padzero(radc,zn,L):
    """
    Pad window with zeros at either end
    Input:
    - radc: np array
    - zn: np array
    - L: float, last value in theta array
    Output:
    - padrad: np array, theta values, extended
    - padz: np array, z values padded
    """
    dr = np.mean(np.diff(radc)) # spacing
    padz = np.concatenate([np.zeros(L*2), zn, np.zeros(L*2)])#padded on both sides by 2*L zeros
    radcor = radc + 2*L*dr
    padrad = np.concatenate([np.linspace(0,2*L*dr,L*2),radcor,np.linspace(radcor[-1],radcor[-1]+2*L*dr,L*2)])#corresponding theta values
    return  padrad,padz

def fftwind(rad,z,windowsize,totalpoint):
    """
    Analyze each window
    Input:
    - rad: np array
    - z: np array
    - windowsize: int
    - totalpoint: int
    Output:
    - rf: np array, frequency space of theta
    - zf: np array, frequency space of z
    - rms: np aray
    - pdft_corr: np array, corrected power spectral density
    """
    radc,zc = linearize(rad,z) #detrend with linear fit
    
    rms = np.sqrt(np.mean(zc**2))
    L = radc[-1]
    
    hwin,win = hannwindow(radc,zc,L) #multiply by hann window
    padrad,padz = padzero(radc,hwin,totalpoint) #pad with zeros for better fourier analysis
    padzn = butterworth(padrad,padz) # filter out large wavelengths
    
    rf, zf_p = fouriertrans(padrad,padzn) # fourier transform
    pdft = zf_p**2/(len(zf_p)) #calculate power spectral density
    pdft_corr = pdft/np.sum(win(radc)) #correct for windowed function (as per Perron et al. 2008)
    return rf,zf_p,rms,pdft_corr

def perimeter(rad,z,dist,windowsize):
    """
    Extract wavelength and rms values for upper third ("north") and lower third ("south") of each concentric profile
    Input:
    - rad: np array
    - z: np array
    - dist: np array
    - windowsize: int
    Output:
    - maxess - np array, peak wavelength south
    - rmsls - np array, rms raw south
    - pdfts - np array, rms from pdft south
    - maxesn - np array, peak wavelength north
    - rmsln - np array, rms raw north
    - pdftn - np array, rms from pdft north
    """
    
    totalpoint = len(rad)
    spac = np.mean(np.diff(dist)) #mean spacing of dist
    windunit = int(windowsize//spac) # how many data points per window

    halfunit = int(windunit//2)
    
    radsouth = rad[totalpoint//12:5*totalpoint//12] #lower third
    distsouth = dist[totalpoint//12:5*totalpoint//12]
    zsouth = z[totalpoint//12:5*totalpoint//12]


    rmsls = []
    maxess = []
    radones = []
    pdfts = []
    for n in range(totalpoint//12,5*totalpoint//12):
        distw = dist[n-halfunit:n+halfunit] #window 
        zw = z[n-halfunit:n+halfunit]
        
        #check if there are nan values in the window. if not,
        if np.isnan(np.sum(zw)) == False:
            freq,zval,rms,pdft_corr = fftwind(distw,zw,windowsize,totalpoint)
            peakwav = 1/freq[np.argmax(pdft_corr)]
            calrms = 2*np.sqrt(np.sum(pdft_corr))
            

            #filtering peak wavelengths of less than 50 (experimentally determined)
            if peakwav > 50.0:
                pdfts.append(calrms)
                maxess.append(peakwav)
                radones.append(rad[n])
                rmsls.append(rms)
            else:
                continue
        else:
            continue

    zradnorth = rad[7*totalpoint//12:11*totalpoint//12] #upper third
    distnorth = dist[7*totalpoint//12:11*totalpoint//12]
    znorth = z[7*totalpoint//12:11*totalpoint//12]

    rmsln = []
    maxesn = []
    radonen = []
    pdftn = []
    for n in range(7*totalpoint//12,11*totalpoint//12):
        distw = dist[n-halfunit:n+halfunit]
        zw = z[n-halfunit:n+halfunit]
        if np.isnan(np.sum(zw)) == False:
            freq,zval,rms,pdft_corr = fftwind(distw,zw,windowsize,totalpoint)
            calrms = 2*np.sqrt(np.sum(pdft_corr))
            peakwav = 1/freq[np.argmax(pdft_corr)]
            
            if peakwav > 50.0:
                pdftn.append(calrms)
                maxesn.append(peakwav)
                radonen.append(rad[n])
                rmsln.append(rms)
            else:
                continue
        else:
            continue



    return maxess,rmsls,pdfts,maxesn,rmsln,pdftn

