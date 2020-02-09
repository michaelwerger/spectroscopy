import os
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import ccdproc
import numpy as np
import math
import sys
import re
sys.path.append('/Users/Micha/Workspaces/python/spectroscopy')
from whsdadoslib.targets import Targets
from whsdadoslib.data import DataFile

class Report(object):
    
    icl = None
    observatory = None
    max_n = 9999
    
    def __init__(self,icl=icl, observatory=observatory):
        
        Report.observatory = observatory    
        
            
            
    @staticmethod
    def list_measurements(icl, targets):
        
               
        format_headers = '{:3s} {:36s} {:10s} {:23s}  {:6s}     {:8s} {:5s}'
        format_rows = '{:3d} {:36s} {:10s} {:23s} {:6.1f} {:8.5f} {:5.1f}'
        n = 0
        
        print (format_headers.format(
            'n',
            'file',
            'target',
            'DATE-OBS',
            'alt',
            'secz',
            '°C')
        )
        for fname,hdu in zip(icl.summary['file'],icl.hdus()):
            if n < Report.max_n:
                
                try:
                    if hdu.header.get('OBJNAME') and re.match('\w',hdu.header['OBJNAME']):
                        target_name = hdu.header['OBJNAME']
                    else:
                        target_name = Targets.from_filename(fname)
                    target = Targets.targets[target_name]
                    obstime = hdu.header['DATE-OBS']
                    alt_az = AltAz(obstime=obstime, location = Report.observatory['earth_location'])

                    target['altaz'] = target['radec'].transform_to(alt_az)
                    print (format_rows.format(
                            n,
                            fname, 
                            target_name,
                            hdu.header['DATE-OBS'],
                            target['altaz'].alt, 
                            target['altaz'].secz,
                            hdu.header['CCD-TEMP']
                            )
                            )
                except KeyError:
                    obstime = hdu.header['DATE-OBS']
                    print (format_rows.format(
                            n,
                            fname, 
                            'UNKNOWN',
                            hdu.header['DATE-OBS'],
                            0.0, 
                            0.0,
                            hdu.header['CCD-TEMP']
                            )
                            )
            n += 1
            #print (hdu.header)
        
    
    @staticmethod
    def list_measurement_quality(icl, center_row=1386):
        
        format_headers = '{:4s} {:40s} {:12s} {:5s} {:5s} {:6s} {:10s}'
        format_rows = '{:4d} {:30s}  {:12d} {:5d} {:5.0f} {:6d} {:10.2f}'
        
        n = 0
        
        print (format_headers.format(
            'N', 
            'FILE', 
            'MAX_INTENS', 
            'MAX_POS', 
            'DELTA', 
            'WIDTH', 
            'DEVIATION')
        )
        
        for f,hdu in zip(icl.summary['file'], icl.hdus()): # over all images in catalog
            if hdu.header['IMAGETYP'] == 'Light Frame'  and n < Report.max_n:
                #traced = np.fliplr(hdu.data).sum(axis=1)
                traced = hdu.data.sum(axis=1)
                max_intensity = max(traced)
                max_pos = np.argmax(traced)
                
                normalizedtrace = traced/max(traced)
            
                x0 = np.argmax(normalizedtrace)
                pixels = np.where(normalizedtrace > 0.5)
                width = 2*(max(pixels[0])- min(pixels[0]))
                y =[math.exp(-(x-x0)*(x-x0)/2/width)  for x in range(0,3520)]
                diff = abs(normalizedtrace - y) 
                
                print (format_rows.format(
                    n, 
                    f, 
                    int(max_intensity), 
                    max_pos, 
                    int(max_pos-center_row),
                    width, 
                    diff.sum()))

            hdu.header['MAX_POS'] = max_pos
            hdu.header['DELTA  '] = max_pos-center_row
            hdu.header['WIDTH  '] = width
            
            n += 1
        return icl
        
            
def main():
    
    observatory = {
        'name':'WHS',
        'earth_location': EarthLocation(
            lat=51.393716*u.deg, 
            lon=6.978611*u.deg, 
            height=120*u.m),
        'utc_offset':1*u.hour
    }
                   
    target = {
    'name':'ALPAQL',
    'radec': SkyCoord.from_name('alpha_aquilae')
    }

                   
    fits_path  = os.path.join('/Users','Micha','data','20190920','data')
    
    icl = ccdproc.ImageFileCollection(
    fits_path,
    glob_include=target['name']+'_*')
                   
    report = Report(icl=icl, observatory=observatory)
    report.list_measurements(target)
    print ("done")

    


if __name__ == "__main__":
    
    main()       
       
    