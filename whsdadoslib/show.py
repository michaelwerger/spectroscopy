import matplotlib.pyplot as plt
import matplotlib._color_data as mcd
import matplotlib.patches as mpatch
import numpy as np
import sys
import os
import platform
import math
from scipy.interpolate import interp1d
from astropy.table import Table
import colour
from colour_demosaicing import (
    EXAMPLES_RESOURCES_DIRECTORY,
    demosaicing_CFA_Bayer_bilinear,
    demosaicing_CFA_Bayer_Malvar2004,
    demosaicing_CFA_Bayer_Menon2007,
    mosaicing_CFA_Bayer)

from .Calibration.wavelength import WavelengthCalibration
from .parameters import FigureSize, CCDParameters
from .linetable import Linetable

from . import find_nearest_index


class Show(object):

    figure_width = 15
    figure_height = 12
    #xsize = 4656
    #ysize = 3520
    xsize = 4944
    ysize = 3284
    
    @staticmethod
    def __init__(self, figure_width=None, figure_height=None, xsize=4656, ysize=3520):  # xsize = 4944, ysize= 3284
        if figure_width:
            Show.figure_width = figure_width
        if figure_height:
            Show.figure_height = figure_height

        Show.xsize = xsize
        Show.ysize = ysize

        
    

    @staticmethod
    def along_slit(icl, slit_positions=None, dark=None, flip=None, rgb=None):  ### replaced by ShowImages.slits()
        if dark is None:
            dark = 1.0
        if flip is None:
            flip = False
        if rgb:
            for f,hdu in zip(icl.summary['file'], icl.hdus()): # over all images in catalog

                plt.rcParams['figure.figsize'] = FigureSize.THIN
                fig, ax = plt.subplots()
                if flip == True:
                    _data = np.flip(hdu.data, axis=1)
                else:
                    _data = hdu.data
                
                data = _data - np.ones((CCDParameters.ysize,CCDParameters.xsize))*dark
                greydata = demosaicing_CFA_Bayer_bilinear(data,'RGGB')
                plotcolor = ['red','green','blue']
                for channel in [0,1,2]:
                    summed_data = greydata[:,:,channel].sum(axis=1)
                    max_value  = np.amax(summed_data, axis=0)
                    normlized_data = summed_data / max_value

                    plt.plot(normlized_data, colour=plotcolor[channel])
                    for sp in slit_positions:
                        plt.plot([sp,sp],[0,0.2])

                plt.show() 
        else:
            for f,hdu in zip(icl.summary['file'], icl.hdus()): # over all images in catalog

                plt.rcParams['figure.figsize'] = FigureSize.THIN
                fig, ax = plt.subplots()
                if flip == True:
                    _data = np.flip(hdu.data, axis=1)
                else:
                    _data = hdu.data
                
                data = _data - np.ones((CCDParameters.ysize,CCDParameters.xsize))*dark
                summed_data = data.sum(axis=1)
                max_value  = np.amax(summed_data, axis=0)
                normlized_data = summed_data / max_value

                plt.plot(normlized_data)
                for sp in slit_positions:
                    plt.plot([sp,sp],[0,0.2])

                plt.show() 


    @staticmethod
    def average_along_columns(twod_array):

        plt.rcParams['figure.figsize'] = FigureSize.NARROW
        fig, ax = plt.subplots()
        traced = twod_array.sum(axis=0)
        plt.plot(traced/twod_array.shape[0])
        plt.xlabel('pixel rows (y-axis)')
        plt.ylabel('averaged counts, summed along rows')
        plt.show()

    @staticmethod
    def full_columns(icl):
        
        plt.rcParams['figure.figsize'] = FigureSize.NARROW

        fig, ax = plt.subplots()
        for hdu in icl.hdus():
            xsize = hdu.data.shape[1]
            ysize = hdu.data.shape[0]
            column = hdu.data[0:ysize,int(xsize/2):int(xsize/2)+1]
            if column.size > 0:
                plt.plot(column)
        plt.xlabel('pixel rows (perpendicular to wavelength axis)')
        plt.ylabel('measured counts')
        #plt.ylim(0,100000)
        plt.show()

    @staticmethod
    def images(icl, xlim=[0,4655], ylim=[0,3519], parameters=None, ccdparameters=None):  # replacded by ShowImages.images()


        for f,hdu in zip(icl.summary['file'], icl.hdus()): # over all images in catalog
            
            #img = np.fliplr(hdu.data)
            
            plt.rcParams['figure.figsize'] = FigureSize.THIN
            fig, ax = plt.subplots()
            if parameters.flip == True:
                _data = np.flip(hdu.data, axis=1)
            else:
                _data = hdu.data
            data = _data - np.ones((ccdparameters.ysize,ccdparameters.xsize))*parameters.dark
            img = ax.imshow(data, vmin=parameters.vmin, vmax=parameters.vmax) 
            fig.colorbar(img)
            plt.title(f)
            plt.xlabel('cols')
            plt.ylabel('rows')
            plt.ylim(ylim)
            plt.xlim(xlim)
                
            plt.show()


    @staticmethod
    def plot_table(f, xlimits=[3500, 8000], colname=None, ylim=[0,1]):
        '''
        models_dir = os.path.join('/Users','Micha','Workspaces','Kurucz','data','archive.stsci.edu','hlsps','reference-atlases','cdbs','grid','ck04models')
        model_file  = os.path.join(models_dir,'ckp00','ckp00_19000.fits')
        Show.plot_table(model_file, colname='g40', ylim=[0,1e9])
        '''

        t = Table.read(f)
        plt.rcParams['figure.figsize'] = FigureSize.NARROW
        fig, ax = plt.subplots()

        w =  t.as_array(True,'WAVELENGTH')
        _waves = [float(_w[0]) for _w in w]
        indices = np.bitwise_and(np.array(w,dtype=np.int64) > xlimits[0], np.array(w, dtype=np.int64) < xlimits[1])
        
        linetable = Linetable()
        plt.xlim(xlimits[0], xlimits[1])

        d = t.as_array(True,colname)
        selected_values = d[indices].astype(float)
        values = d.astype(float)
        max_value  = np.amax(selected_values, axis=0)
        print (max_value)
        
        plt.plot(w,[v/max_value for v in values])
        plt.ylim(ylim)

        linetable = Linetable()
        for label, wavelength in zip(linetable.get_lines(), linetable.get_wavelengths()):
            if xlimits[0] < wavelength and wavelength < xlimits[1]:
                ix = find_nearest_index(_waves, float(wavelength))
                plt.text(wavelength,0.01,label,rotation=90, color='red',size=7.0,  horizontalalignment = 'center')
                plt.plot([wavelength, wavelength], [0.1, values[ix]*0.9/max_value], color='red')

        _title = "%s[%s]" % (f.replace('.fits',''),colname)
        plt.title(_title)
        plt.show()

    @staticmethod
    def plots(icl,factors = None, xlim=None, ylim=None):

        if xlim is None:
            xlim = [0, -1]
        if ylim is None:
            ylim = [0, -1]
        if factors is None:
            factors = np.ones(len(datafiles))
        plt.rcParams['figure.figsize'] = FigureSize.NARROW
        fig, ax = plt.subplots()
        for hdu, factor in zip(icl.hdus(),factors):
            traced = hdu.data[ylim[0]:ylim[1],:].sum(axis=0)
            plt.plot(traced * factor)
            print (max(traced))
        plt.show()
            


    @staticmethod
    def show_extinctioncurve(traces, airmass1=None, airmass2=None):
        z = []
        
        with open('wavelength_solution.dat','r') as f:
            for line in f:
                z.append(float(line))
        
        positions = range(0,len(traces))
        if len(z) == 1:
            waves = [pos * z for pos in positions]
        elif len(z) == 2:
            waves = [pos * z[0] + z[1] for pos in positions]
        elif len(z) == 3:
            waves = [pos*pos* z[0] + pos * z[1] + z[2] for pos in positions]
        else:
            raise NotImplementedError
        wavelengths = []
        extcoeffs = []
        currentdir = os.getcwd()
        if 'Venus' in platform.node():
            ext_path = os.path.join('/Users', 'Micha', 'data', 'ref', 'atmosphericextinction')
        else:
            ext_path = os.path.join('d:/','Workspaces','data','ref','atmosphericextinction')
        try:
            os.chdir(ext_path)
            with open('paranal_extinction.dat','r') as f:
                for line in f:

                    tokens = line.split(' ')
                    if len(tokens) == 2:
                        wavelength, extcoeff = tokens
                        #print (wave, extcoeff)
                        wavelengths.append(float(wavelength))
                        extcoeffs.append(float(extcoeff))
        finally:
            os.chdir(currentdir)

        

        ext_std = interp1d(wavelengths, extcoeffs, fill_value="extrapolate")
        i_std = ext_std(waves)

        plt.rcParams['figure.figsize'] = FigureSize.NARROW
        fig, ax = plt.subplots()

        plt.plot(waves,i_std)

        if airmass1 is not None and airmass2 is not None:
            n_1 = [(ext_std(wave) * airmass1) for wave in waves]
        plt.show()
        
    @staticmethod
    def selected_columns(icl, column_no=3600, xlim=[2000,3500], ylim=[0,100000], max_i_limit=90000):

        plt.rcParams['figure.figsize'] = FigureSize.NARROW
        colors = mcd.CSS4_COLORS
        color_names = list(colors.keys())
        n = 0  # sequence number of measurement
        color_i = 10  # color number

        dy = (ylim[1]-ylim[0])/50 

        for hdu in icl.hdus():
            if n < 21:
                cn = mcd.CSS4_COLORS[color_names[color_i]]
                
                roi = (DataProcessing.getroi(hdu.data))
                data = np.array(hdu.data[roi[0]:roi[1],:])
                mean_dark = DataProcessing.getdark(hdu.data)
                summed = data.sum(axis=0) - mean_dark * (roi[1] - roi[0])
                max_i = summed[column_no]
                #print (n,hdu.header['DATE-OBS'],hdu.header['CCD-TEMP'], max_i)
                if max_i > max_i_limit:
                    column = hdu.data[0:Show.ysize,int(Show.xsize/2):int(Show.xsize/2)+1]

                    plt.plot(column, color=cn)
                    text = "n=%2d %s, %s" % (n, hdu.header['DATE-OBS'],str(max_i)) 
                    plt.text(xlim[0],ylim[1]-dy*color_i,text, color=cn)
                    color_i += 1
            n += 1

        plt.xlabel('pixel rows (perpendicular to wavelength axis)')
        plt.ylabel('measured counts')
        plt.ylim(ylim[0], ylim[1])
        plt.xlim(xlim[0], xlim[1])
        plt.show()

    @staticmethod
    def selected_images(icl, peaks=None, delta=None, vmin=0, vmax=1000):
        center_row = (peaks[2] + peaks[3])/2
        plt.rcParams['figure.figsize'] = FigureSize.NARROW
        for f,hdu in zip(icl.summary['file'], icl.hdus()):
            traced = np.fliplr(hdu.data).sum(axis=1)
            max_intensity = max(traced)
            max_pos = np.argmax(traced)
            if abs(max_pos-center_row) < delta:
                print (peaks[2]-100,peaks[3]+100,0,Show.xsize,f, max_intensity)
                fig, ax = plt.subplots()
                plt.imshow(hdu.data[peaks[2]-100:peaks[3]+100,:],vmin=vmin,vmax=vmax)
                plt.show()


    @staticmethod
    def selected_plots(icl, peaks=None, delta=None, vmin=0, vmax=1000):
        center_row = (peaks[2] + peaks[3])/2
        plt.rcParams['figure.figsize'] = FigureSize.NARROW
        for f,hdu in zip(icl.summary['file'], icl.hdus()):
            traced = np.fliplr(hdu.data).sum(axis=1)
            max_intensity = max(traced)
            max_pos = np.argmax(traced)
            if abs(max_pos-center_row) < delta:
                print (peaks[2]-100,peaks[3]+100,0,Show.xsize,f, max_intensity)
                fig, ax = plt.subplots()
                axis=0 # along dispersion
                #axis=1 # perpendicular to dispersion 
                plt.plot(hdu.data.sum(axis=axis), color='black')
                plt.plot(hdu.data[peaks[2]-100:peaks[3]+100,:].sum(axis=axis), color='blue')
                
                right_inset_ax = fig.add_axes([.65, .6, .2, .2])
                right_inset_ax.plot(hdu.data.sum(axis=1))
                right_inset_ax.set_title('Trace')
                right_inset_ax.set_xlim(peaks[2],peaks[3])
                right_inset_ax.set_xticks([])
                right_inset_ax.set_yticks([])
                plt.show()

    @staticmethod
    def selected_rows(icl, column_no=3600, xlim=[2000,3500], ylim=[0,100000], max_i_limit=90000):

        plt.rcParams['figure.figsize'] = FigureSize.NARROW
        colors = mcd.CSS4_COLORS
        color_names = list(colors.keys())
        n = 0  # sequence number of measurement
        color_i = 10  # color number

        dy = (ylim[1]-ylim[0])/50 

        for hdu in icl.hdus():
            if n < 21:
                cn = mcd.CSS4_COLORS[color_names[color_i]]
                
                roi = (DataProcessing.getroi(hdu.data))
                data = np.array(hdu.data[roi[0]:roi[1],:])
                mean_dark = DataProcessing.getdark(hdu.data)
                summed = data.sum(axis=1) - mean_dark * (roi[1] - roi[0])
                max_i = summed[column_no]
                #print (n,hdu.header['DATE-OBS'],hdu.header['CCD-TEMP'], max_i)
                if max_i > max_i_limit:
                    column = hdu.data[0:Show.ysize,int(Show.xsize/2):int(Show.xsize/2)+1]

                    plt.plot(column, color=cn)
                    text = "n=%2d %s, %s" % (n, hdu.header['DATE-OBS'],str(max_i)) 
                    plt.text(xlim[0],ylim[1]-dy*color_i,text, color=cn)
                    color_i += 1
            n += 1

        plt.xlabel('pixel rows (perpendicular to wavelength axis)')
        plt.ylabel('measured counts')
        plt.ylim(ylim[0], ylim[1])
        plt.xlim(xlim[0], xlim[1])
        plt.show()

    
    @staticmethod
    def show_standard_flux(ref_waves,ref_fluxes):
        
        f_std = interp1d(ref_waves,ref_fluxes)
        plt.rcParams['figure.figsize'] = FigureSize.NARROW
        fig, ax = plt.subplots()
        plt.plot(ref_waves, f_std(ref_waves))

        plt.show()

    @staticmethod
    def sky_areas(icl, xlim=[0,4655], ylim=[0,3519], sky_width=50, parameters=None, ccdparameters=None):


        
        for f,hdu in zip(icl.summary['file'], icl.hdus()): # over all images in catalog
            
            #img = np.fliplr(hdu.data)
            
            plt.rcParams['figure.figsize'] = FigureSize.THIN
            fig, ax = plt.subplots()
            if parameters.flip == True:
                _data = np.flip(hdu.data, axis=1)
            else:
                _data = hdu.data
            data = _data - np.ones((ccdparameters.ysize,ccdparameters.xsize))*parameters.dark
            img = ax.imshow(data, vmin=parameters.vmin, vmax=parameters.vmax) 
            fig.colorbar(img)
            plt.title(f)
            plt.xlabel('cols')
            plt.ylabel('rows')

            plt.plot([xlim[0],xlim[1]],[parameters.slit_positions[2],parameters.slit_positions[2]], color='red')
            plt.plot([xlim[0],xlim[1]],[parameters.slit_positions[2]+sky_width,parameters.slit_positions[2]+sky_width], color='red')
            plt.plot([xlim[0],xlim[0]],[parameters.slit_positions[2],parameters.slit_positions[2]+sky_width], color='red')
            plt.plot([xlim[1],xlim[1]],[parameters.slit_positions[2],parameters.slit_positions[2]+sky_width], color='red')

            plt.plot([xlim[0],xlim[1]],[parameters.slit_positions[3],parameters.slit_positions[3]], color='green')
            plt.plot([xlim[0],xlim[1]],[parameters.slit_positions[3]-sky_width,parameters.slit_positions[3]-sky_width], color='green')
            plt.plot([xlim[0],xlim[0]],[parameters.slit_positions[3],parameters.slit_positions[3]-sky_width], color='green')
            plt.plot([xlim[1],xlim[1]],[parameters.slit_positions[3],parameters.slit_positions[3]-sky_width], color='green')


            plt.ylim(parameters.slit_positions[2], parameters.slit_positions[3])
                
            plt.show()

    @staticmethod
    def show_traces(icl, peaks=None, column=2800, ylim=[0,1000], max_n=9999):
        n = 0
        for hdu in icl.hdus():
            plt.rcParams['figure.figsize'] = FigureSize.NARROW
            fig, ax = plt.subplots()
            plt.plot(hdu.data[:,column])
            if peaks is not None:
                    for peak in peaks:
                        plt.plot([peak,peak],ylim)
            plt.ylim(ylim)
            plt.show()
            n += 1
            if n > max_n:
                break


    @staticmethod
    def wavelength_check(trace, wavelengths_file='wavelengths.txt'):  ### replaced by class ShowWavelengthCheck.plot():

        plt.rcParams['figure.figsize'] = FigureSize.NARROW

        fig, ax = plt.subplots()
        

        z = WavelengthCalibration.get_z()

        positions = range(0,len(trace))
        if len(z) == 1:
            waves = [pos * z for pos in positions]
        elif len(z) == 2:
            waves = [pos * z[0] + z[1] for pos in positions]
        elif len(z) == 3:
            waves = [pos*pos* z[0] + pos * z[1] + z[2] for pos in positions]
        else:
            raise NotImplementedError


        plt.plot(waves,trace)

        max_i = max(trace)



        with open(wavelengths_file, 'r') as f:
            f.readline()
            min_x= 10000
            max_x= 0
            for line in f:
                if line.startswith('#'):
                    pass
                else:
                    tokens = line.split(',')
                    #print (tokens)
                    if len(tokens) == 3:
                        wavelength = float(tokens[1])
                        plt.plot([wavelength,wavelength],[0,2*max_i/3], color='red')
                        plt.text(wavelength, 2*max_i/3, tokens[2], rotation=90, color='red') #,  rotation_mode='anchor')
        plt.xlabel('wavelengths (A)')
        plt.ylabel('measured counts')
        plt.show()

class CCDImage(Show):

    def __init__(self, figure_width=None, figure_height=None):
        super().__init__(figure_width, figure_height)

    @staticmethod
    def show(icl, xlim=[0,4655], ylim=[0,3519], parameters=None, ccdparameters=None):
        Show.images(icl, xlim, ylim, parameters, ccdparameters)

    @staticmethod
    def slits(icl, slit_positions=None, dark=None, flip=None, rgb=None):
        if slit_positions:
            Show.along_slit(icl, slit_positions, dark, flip, rgb)  # replaced by ShowImages.slits()
        else:
            Show.images(icl)


class WavelengthCheck(Show):

    def __init__(self, figure_width=None, figure_height=None):
        super().__init__(figure_width, figure_height)

    @staticmethod
    def plot(trace, wavelengths_file='wavelength.txt'):
        super().wavelength_check(trace, wavelengths_file)
