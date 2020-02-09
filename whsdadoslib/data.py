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

sys.path.append('/Users/Micha/Workspaces/python/spectroscopy')
from whsdadoslib.processing import ProcessingData



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
    def getdark(data):
        #
        # select a region left and right from the spectrum to determine the average 
        # value of tbe intensity to be used as dark
        #data = np.array(data_obj[0].data)
        traced = data.sum(axis=1)
        max_column = np.argmax(traced)
        #print (max_column,max_column-distance_from_max, 
        #       max_column-distance_from_max+dark_window)
        #print (max_column,max_column+distance_from_max-dark_window,
        #       max_column+distance_from_max)
        left_data = np.array(data[max_column-DataProcessing.distance_from_max:
                                  max_column-DataProcessing.distance_from_max+DataProcessing.dark_window,:])
        right_data = np.array(data[max_column+DataProcessing.distance_from_max-DataProcessing.dark_window:
                                   max_column+DataProcessing.distance_from_max,:])
        left_dark = np.mean(left_data)
        right_dark = np.mean(right_data)
        #print (left_dark)
        #print (right_dark)
        return (left_dark + right_dark)/2.0

    @staticmethod
    def getsky(data):
        data_ncols = data.shape[1]
        data_nrows = data.shape[0]
        column = data[0:data_nrows,int(data_ncols/2):int(data_ncols/2)+1]
        #print (column)
        max_column = np.argmax(column)
        print (max_column)
        min_window = 80
        max_window = 200
        sky_rows = [i for i in range(max_column-max_window,
                                     max_column-min_window)]
        sky_rows = sky_rows + [i for i in range(max_column+min_window,
                                                max_column+max_window)]
        sky = data[1325,:].copy()
        for col in range(0,data_ncols):
            sky[col] = 0.0

        row_counter = 0
        for row in sky_rows:
            row_counter = row_counter + 1
            for col in range(0,data_ncols):
                sky[col] = sky[col] + data[row,col]

        for col in range(0,data_ncols):
            sky[col] = sky[col] / row_counter

        return sky

    @staticmethod
    def correct_dark(icl, plot=False, fliplr=True, overwrite=True):

        import matplotlib.pyplot as plt
        for f,hdu in zip(icl.summary['file'],icl.hdus()):
            if fliplr:
                flipped = np.fliplr(hdu.data)
            else:
                flipped = hdu.data
            hdu.data = (flipped - DataProcessing.masterdark)# .astype(np.uint16)

            hdu.header['DARKCORR'] = True
        
            if plot:
                figure_width = 15
                figure_height = 12
                plt.rcParams['figure.figsize'] = [figure_width, figure_height]
                fig, ax = plt.subplots()
                axis=0 # along dispersion
                #axis=1 # perpendicular to dispersion
                plt.plot(flipped.sum(axis=axis), color='blue')
                plt.plot(DataProcessing.masterdark.sum(axis=axis), color='red')
                plt.plot(hdu.data.sum(axis=axis), color='black')
                plt.show()
            if f.endswith('s'):
                fname = f 
            else:
                fname = f + 's'
            
            fitsfile = os.path.join(ProcessingData.corr_path,fname)
            if os.path.exists(fitsfile) and overwrite == True:
                os.remove(fitsfile)
            hdu.writeto(fitsfile)

            os.system('gzip ' + fitsfile)
            print (fname + "written")
        
    @staticmethod
    def correct_sky(icl):

        for hdu in icl.hdus():
            if hdu.header['DARKCORR'] == True:
                pass

    @staticmethod
    def savitzky_golay(y, window_size, order, deriv=0, rate=1):
        r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
        The Savitzky-Golay filter removes high frequency noise from data.
        It has the advantage of preserving the original shape and
        features of the signal better than other types of filtering
        approaches, such as moving averages techniques.
        Parameters
        ----------
        y : array_like, shape (N,)
            the values of the time history of the signal.
        window_size : int
            the length of the window. Must be an odd integer number.
        order : int
            the order of the polynomial used in the filtering.
            Must be less then `window_size` - 1.
        deriv: int
            the order of the derivative to compute (default = 0 means only smoothing)
        Returns
        -------
        ys : ndarray, shape (N)
            the smoothed signal (or it's n-th derivative).
        Notes
        -----
        The Savitzky-Golay is a type of low-pass filter, particularly
        suited for smoothing noisy data. The main idea behind this
        approach is to make for each point a least-square fit with a
        polynomial of high order over a odd-sized window centered at
        the point.
        Examples
        --------
        t = np.linspace(-4, 4, 500)
        y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
        ysg = savitzky_golay(y, window_size=31, order=4)
        import matplotlib.pyplot as plt
        plt.plot(t, y, label='Noisy signal')
        plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
        plt.plot(t, ysg, 'r', label='Filtered signal')
        plt.legend()
        plt.show()
        References
        ----------
        .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
           Data by Simplified Least Squares Procedures. Analytical
           Chemistry, 1964, 36 (8), pp 1627-1639.
        .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
           W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
           Cambridge University Press ISBN-13: 9780521880688
        .. [3] https://scipy-cookbook.readthedocs.io/items/SavitzkyGolay.html
        """

        

        try:
            window_size = np.abs(np.int(window_size))
            order = np.abs(np.int(order))
        except ValueError:
            raise ValueError("window_size and order have to be of type int")
        if window_size % 2 != 1 or window_size < 1:
            raise TypeError("window_size size must be a positive odd number")
        if window_size < order + 2:
            raise TypeError("window_size is too small for the polynomials order")
        order_range = range(order+1)
        half_window = (window_size -1) // 2
        # precompute coefficients
        b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
        m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
        # pad the signal at the extremes with
        # values taken from the signal itself
        firstvals = y[0] - np.abs(y[1:half_window+1][::-1] - y[0] )
        lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
        y = np.concatenate((firstvals, y, lastvals))
        return np.convolve( m[::-1], y, mode='valid')
    
    @staticmethod
    def sgolay2d( z, window_size, order, derivative=None):
        """
        References
        ----------
        .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
           Data by Simplified Least Squares Procedures. Analytical
           Chemistry, 1964, 36 (8), pp 1627-1639.
        .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
           W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
           Cambridge University Press ISBN-13: 9780521880688
        .. [3] https://scipy-cookbook.readthedocs.io/items/SavitzkyGolay.html
        """
        # number of terms in the polynomial expression
        n_terms = ( order + 1 ) * ( order + 2)  / 2.0

        if  window_size % 2 == 0:
            raise ValueError('window_size must be odd')

        if window_size**2 < n_terms:
            raise ValueError('order is too high for the window size')

        half_size = window_size // 2

        # exponents of the polynomial. 
        # p(x,y) = a0 + a1*x + a2*y + a3*x^2 + a4*y^2 + a5*x*y + ... 
        # this line gives a list of two item tuple. Each tuple contains 
        # the exponents of the k-th term. First element of tuple is for x
        # second element for y.
        # Ex. exps = [(0,0), (1,0), (0,1), (2,0), (1,1), (0,2), ...]
        exps = [ (k-n, n) for k in range(order+1) for n in range(k+1) ]

        # coordinates of points
        ind = np.arange(-half_size, half_size+1, dtype=np.float64)
        dx = np.repeat( ind, window_size )
        dy = np.tile( ind, [window_size, 1]).reshape(window_size**2, )

        # build matrix of system of equation
        A = np.empty( (window_size**2, len(exps)) )
        for i, exp in enumerate( exps ):
            A[:,i] = (dx**exp[0]) * (dy**exp[1])

        # pad input array with appropriate values at the four borders
        new_shape = z.shape[0] + 2*half_size, z.shape[1] + 2*half_size
        Z = np.zeros( (new_shape) )
        # top band
        band = z[0, :]
        Z[:half_size, half_size:-half_size] =  band -  np.abs( np.flipud( z[1:half_size+1, :] ) - band )
        # bottom band
        band = z[-1, :]
        Z[-half_size:, half_size:-half_size] = band  + np.abs( np.flipud( z[-half_size-1:-1, :] )  -band )
        # left band
        band = np.tile( z[:,0].reshape(-1,1), [1,half_size])
        Z[half_size:-half_size, :half_size] = band - np.abs( np.fliplr( z[:, 1:half_size+1] ) - band )
        # right band
        band = np.tile( z[:,-1].reshape(-1,1), [1,half_size] )
        Z[half_size:-half_size, -half_size:] =  band + np.abs( np.fliplr( z[:, -half_size-1:-1] ) - band )
        # central band
        Z[half_size:-half_size, half_size:-half_size] = z

        # top left corner
        band = z[0,0]
        Z[:half_size,:half_size] = band - np.abs( np.flipud(np.fliplr(z[1:half_size+1,1:half_size+1]) ) - band )
        # bottom right corner
        band = z[-1,-1]
        Z[-half_size:,-half_size:] = band + np.abs( np.flipud(np.fliplr(z[-half_size-1:-1,-half_size-1:-1]) ) - band )

        # top right corner
        band = Z[half_size,-half_size:]
        Z[:half_size,-half_size:] = band - np.abs( np.flipud(Z[half_size+1:2*half_size+1,-half_size:]) - band )
        # bottom left corner
        band = Z[-half_size:,half_size].reshape(-1,1)
        Z[-half_size:,:half_size] = band - np.abs( np.fliplr(Z[-half_size:, half_size+1:2*half_size+1]) - band )

        # solve system and convolve
        if derivative == None:
            m = np.linalg.pinv(A)[0].reshape((window_size, -1))
            return fftconvolve(Z, m, mode='valid')
        elif derivative == 'col':
            c = np.linalg.pinv(A)[1].reshape((window_size, -1))
            return fftconvolve(Z, -c, mode='valid')
        elif derivative == 'row':
            r = np.linalg.pinv(A)[2].reshape((window_size, -1))
            return fftconvolve(Z, -r, mode='valid')
        elif derivative == 'both':
            c = np.linalg.pinv(A)[1].reshape((window_size, -1))
            r = np.linalg.pinv(A)[2].reshape((window_size, -1))
            return fftconvolve(Z, -r, mode='valid'), fftconvolve(Z, -c, mode='valid')

    @staticmethod
    def check():
        print ("checked")

