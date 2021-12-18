import numpy as np
from scipy.signal import fftconvolve
from math import factorial
import os
import re
import sys
import shutil
import gzip

from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.io import fits

#sys.path.append('/Users/Micha/Workspaces/python/spectroscopy')
from processing import ProcessingData

class DataProcessing(object):

    roi_window = 150
    distance_from_max = 150 # dstance from row with max intensity
    dark_window = 5

    masterdark = None
    masterflat = None
    masterbias = None
    fitspath = os.path.join('..','data')
    if os.path.exists(fitspath):
        pass
    else:
        fitspath = '.'

    def __init__(self,fitspath=fitspath, masterdark=masterdark, masterflat=masterflat, masterbias=masterbias):
        DataProcessing.masterdark = masterdark

        if fitspath is not None:
            DataProcessing.fitspath = fitspath

        if masterbias is not None:
            DataProcessing.masterbias = masterbias
        else:
            sh = masterdark.shape
            DataProcessing.masterbias = np.zeros((sh[0], sh[1]))

        if masterflat is not None:
            DataProcessing.masterflat = masterflat
        else:
            sh = masterdark.shape
            DataProcessing.masterflat = np.ones((sh[0], sh[1]))
        

    @staticmethod
    def getroi(data):
        # this method selects the region of interest as area on the CCD image:
        # * completely along the wavelength axis
        # * slice with size 2*roi_window, centered at the maximum row
        xsize = data.shape[1]
        ysize = data.shape[0]
        #print (0, ysize, int(xsize/2), int(xsize/2))
        column = data[0:ysize,int(xsize/2):int(xsize/2)+1]
        # 
        # maximum, row
        iy = np.argmax(column)
        #print (iy-Data.roi_window,iy+Data.roi_window,0,xsize, 'max_row=',iy,'max=',max(data[iy,:]))
        return (iy-DataProcessing.roi_window,iy+DataProcessing.roi_window,0,xsize)