import numpy as np
import os


from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.io import fits

#sys.path.append('/Users/Micha/Workspaces/python/spectroscopy')

class DataProcessing(object):

    roi_window = 150
    distance_from_max = 150 # dstance from row with max intensity
    dark_window = 5
    
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

    def getwidth(icl, parameters=None, ccdparameters=None):

        widths = []
        ymin = parameters.slit_positions[2]
        ymax = parameters.slit_positions[3]

        n = 0
        _l_ave = 0
        for f,hdu in zip(icl.summary['file'], icl.hdus()): # over all images in catalog
            if parameters.flip == True:
                _data = np.flip(hdu.data, axis=1)
            else:
                _data = hdu.data
            data = _data - np.ones((ccdparameters.ysize,ccdparameters.xsize))*parameters.dark
            summed_data = data.sum(axis=1)
            max_value  = np.amax(summed_data)
            normlized_data = summed_data / max_value
            _tmp = 0.1
            
            indices = np.where(normlized_data > _tmp)
            
            widths.append(len(indices[0]))

        return widths
        