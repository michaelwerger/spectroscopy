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



class DataFile(object):

    def __init__(self, filename, fits_path=None):
        if fits_path is None:
            self.filename = filename
        else:
            self.filename = os.path.join(fits_path, filename)
        self.header = None
        self.data = None

    def set_filename(self,filename, fits_path=None):
        if fits_path is None:
            self.filename = filename
        else:
            self.filename = os.path.join(fits_path, filename)
    
    def read(self):
        with fits.open(self.filename) as hdu_list:
            self.header = hdu_list[0].header.tostring('\n') # for simplicity here only 0
            self.data = hdu_list[0].data

    def get_header(self):
        with fits.open(DataFile.filename) as hdu_list:
            self.header = hdu_list[0].header.tostring('\n') # for simplicity here only 0
        return self.header

    def get_data(self):
        with fits.open(DataFile.filename) as hdu_list:
            self.data = hdu_list[0].data
        
        return self.data


    @staticmethod
    def show_info(f):
        hdul = fits.open(f)
        hdul.info() 

    def show_header(f):
        hdu_list = fits.open(f)
        print (len(hdu_list))
        print (hdu_list[0].header.tostring('\n'))

    @staticmethod
    def update_fits_header(fits_path, icl, targets, observatory):
        os.chdir(fits_path)
        for f,hdu in zip(icl.summary['file'],icl.hdus()):
            obstime = hdu.header['DATE-OBS']

            additional_keywords = { # see https://heasarc.gsfc.nasa.gov/docs/fcg/common_dict.html
                'OBSERVER':  ('Spectrum team','-'),
                'OBJNAME':   (f.split('_')[0],'-'),
                'APERTURE':  (300,"mm"),
                'INSTRUME':  ("DADOS+ASI1600MMPro","-"),
                'TELESCOP':  ("30cm Newton","-"),
                'GRATING':   ("200lin/mm 25mu mid","-"),
                'OBSERVAT':  ("Walter-Hohmann-Obs","-"),
                'FOCALLEN':  (3000.,'Telescope focal length in mm '),
                'APTAREA':   (28260000,'Aperture area mm^2 less central obs'),
                'APTDIA':    (300., 'Aperture diameter in mm')
            }
            
            if hdu.header.get('OBJNAME') and re.match('\w',hdu.header['OBJNAME']):
                target_name = hdu.header['OBJNAME']
            else:
                target_name = targets.from_filename(f)
            try:
                target = targets.targets[target_name]
                alt_az = AltAz(obstime=obstime, location = observatory['earth_location'])
                target['altaz'] = target['radec'].transform_to(alt_az)
   
                additional_keywords['RA']      = (str(target['radec'].ra),'HH:MM:SS')
                additional_keywords['DEC']     = (str(target['radec'].dec),'dd.mm.ss')
                additional_keywords['AIRMASS'] = (float(target['altaz'].secz),'sec(z)')
            except KeyError:
                pass
                
            fname = f+'s'
            
            for k in additional_keywords.keys():
                hdu.header[k] = additional_keywords[k]
            print (fname)
            #print (hdu.header.tostring('\n'))
            if os.path.exists(fname):
                os.remove(fname)
            hdu.writeto(fname)

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

    

    


        


    

    @staticmethod
    def check():
        print ("checked")

