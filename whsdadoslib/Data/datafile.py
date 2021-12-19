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



