import os
import numpy as np
from data import DataProcessing

class Dark(object):

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
    def store_dark_to_file(dark):
        with open('dark.txt','wb') as darkf:
            darkf.write(dark)

    @staticmethod
    def read_dark_from_file():
        with open('dark.txt','rb') as darkf:
            dark = darkf.read()
        return dark

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
            if os.path.exists(fitsfile + '.gz') and overwrite == True:
                os.remove(fitsfile + '.gz')
                
            hdu.writeto(fitsfile)

            os.system('gzip ' + fitsfile)
            print (fname + "written")