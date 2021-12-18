import os
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import numpy as np
import math
import sys
sys.path.append('/Users/Micha/Workspaces/python/spectroscopy/whsdadoslib')

class Report(object):
    
    icl = None
    observatory = None
    max_n = 9999
    
    def __init__(self,icl=icl, observatory=observatory):
        
        Report.observatory = observatory    
            
    @staticmethod
    def list_measurements(icl):
        
        n = 0
        print ("{:3s} {:68} {:24s} {:12s} {:6s} {:7s}".format(
            'N', 'FILE', 'UTC', 'IMAGETYP', 'EXPT', 'CCDTEMP'))
        for f,hdu in zip(icl.summary['file'], icl.hdus()): # over all images in catalog
        
            tokens = f.split('_')
            if len(tokens) > 9:
                imagetyp, _label, _abbr, exptime, binning, temp, gain, obsdate, obsclock, frame = f.split('_')
                object = _label + '_' + _abbr
            elif len(tokens) == 9:
                imagetyp, object, exptime, binning, temp, gain, obsdate, obsclock, frame = f.split('_')
            else:
                raise NotImplementedError('do not know how to handle '+f)


            print ("{:3d} {:68s} {:24s} {:12s} {:6.1f} {:7.1f}".format(
                n, 
                f, 
                hdu.header['DATE-OBS'], 
                imagetyp, 
                hdu.header['EXPTIME'], 
                hdu.header['CCD-TEMP']))
            n += 1       
            
    @staticmethod
    def list_exposuredata(icl, targets, target_name, observatory):
        n = 0
        print ("{:3s} {:68}  {:6s} {:7s}  {:12s} {:8s}".format(
            'N', 'FILE',  'EXPT', 'CCDTEMP', 'ALTITUDE', 'SEC(Z)'))
        for f,hdu in zip(icl.summary['file'], icl.hdus()): # over all images in catalog
            obstime = hdu.header['DATE-OBS']
            alt_az = AltAz(obstime=obstime, 
                location = observatory['WHS']['earth_location'])
            targets[target_name]['altaz'] = targets[target_name]['radec'].transform_to(alt_az)
            
            tokens = f.split('_')
            if len(tokens) > 9:
                imagetyp, _label, _abbr, exptime, binning, temp, gain, obsdate, obsclock, frame = f.split('_')
                object = _label + '_' + _abbr
            elif len(tokens) == 9:
                imagetyp, object, exptime, binning, temp, gain, obsdate, obsclock, frame = f.split('_')
            else:
                raise NotImplementedError('do not know how to handle '+f)
            
            print ("{:3d} {:68s} {:6.1f} {:7.1f} {:8.4f} {:8.5f}".format(
                n, 
                f, 
                hdu.header['EXPTIME'], 
                hdu.header['CCD-TEMP'],
                targets[target_name]['altaz'].alt, 
                targets[target_name]['altaz'].secz
                ))
            n += 1
    
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
    
    print ("done")

    


if __name__ == "__main__":
    
    main()       
       
    