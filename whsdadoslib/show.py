import matplotlib.pyplot as plt
import matplotlib._color_data as mcd
import matplotlib.patches as mpatch
import numpy as np
import sys
import os
import platform
import math
from scipy.interpolate import interp1d
from whsdadoslib.data import DataProcessing


class Show(object):

    figure_width = 15
    figure_height = 12
    xsize = 4656
    ysize = 3520
    
    def __init__(self, figure_width=None, figure_height=None):
        if figure_width:
            Show.figure_width = figure_width
        if figure_height:
            Show.figure_height = figure_height
    
    @staticmethod

    def plots(icl,factors = None, xlim=None, ylim=None):

        if xlim is None:
            xlim = [0, -1]
        if ylim is None:
            ylim = [0, -1]
        if factors is None:
            factors = np.ones(len(datafiles))
        plt.rcParams['figure.figsize'] = [Show.figure_width, Show.figure_height]
        fig, ax = plt.subplots()
        for hdu, factor in zip(icl.hdus(),factors):
            traced = hdu.data[ylim[0]:ylim[1],:].sum(axis=0)
            plt.plot(traced * factor)
            print (max(traced))
        plt.show()
            
    @staticmethod
    def images(icl, xlim=[0,4655], ylim=[0,3519], max_n=9999):
        n = 0
        plt.rcParams['figure.figsize'] = [Show.figure_width, Show.figure_height]
        for hdu in icl.hdus():
            fig, ax = plt.subplots()
            plt.imshow(hdu.data[ylim[0]:ylim[1],xlim[0]:xlim[1]],vmin=500,vmax=1000)
            plt.show()
            n += 1
            if n > max_n:
                break

    @staticmethod
    def selected_images(icl, peaks=None, delta=None, vmin=0, vmax=1000):
        center_row = (peaks[2] + peaks[3])/2
        plt.rcParams['figure.figsize'] = [Show.figure_width, Show.figure_height]
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
        plt.rcParams['figure.figsize'] = [Show.figure_width, Show.figure_height]
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
    def full_columns(icl):
        
        plt.rcParams['figure.figsize'] = [Show.figure_width, Show.figure_height]

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
    def selected_columns(icl, column_no=3600, xlim=[2000,3500], ylim=[0,100000], max_i_limit=90000):

        plt.rcParams['figure.figsize'] = [Show.figure_width, Show.figure_height]
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
    def selected_rows(icl, column_no=3600, xlim=[2000,3500], ylim=[0,100000], max_i_limit=90000):

        plt.rcParams['figure.figsize'] = [Show.figure_width, Show.figure_height]
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
    def average_along_columns(twod_array):

        plt.rcParams['figure.figsize'] = [Show.figure_width, Show.figure_height]
        fig, ax = plt.subplots()
        traced = twod_array.sum(axis=0)
        plt.plot(traced/twod_array.shape[0])
        plt.xlabel('pixel rows (y-axis)')
        plt.ylabel('averaged counts, summed along rows')
        plt.show()

    @staticmethod
    def show_traces(icl, peaks=None, column=2800, ylim=[0,1000], max_n=9999):
        n = 0
        for hdu in icl.hdus():
            plt.rcParams['figure.figsize'] = [Show.figure_width, Show.figure_height]
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

        plt.rcParams['figure.figsize'] = [Show.figure_width, Show.figure_height]
        fig, ax = plt.subplots()

        plt.plot(waves,i_std)

        if airmass1 is not None and airmass2 is not None:
            n_1 = [(ext_std(wave) * airmass1) for wave in waves]
        plt.show()
        
    @staticmethod
    def show_standard_flux(hst_waves,hst_fluxes):
        plt.rcParams['figure.figsize'] = [Show.figure_width, Show.figure_height]
        fig, ax = plt.subplots()

        waves = range(3000,10000)
        plt.plot(hst_waves, hst_fluxes, '-')
        plt.plot(hst_waves, hst_fluxes, '.', color='red')
        f_vega = interp1d(hst_waves, hst_fluxes, fill_value="extrapolate")
        plt.plot(waves, f_vega(waves))
        plt.xlim(3000,10000)
        #plt.ylim(0,max(f_vega(waves)))
        #plt.plot(waves,summed*40000)
        plt.xlabel('wavelength (A)')
        plt.ylabel('flux (mJansky)')


        fig, ax = plt.subplots()

        #plt.plot(hst_waves, hst_fluxes)
        f_vega = interp1d(hst_waves, hst_fluxes, fill_value="extrapolate")
        f_vega_mjansky = f_vega(waves)
        f_watt = f_vega_mjansky

        for i in range(0,len(f_watt)):
            f_watt[i] = f_vega_mjansky[i]/waves[i]/waves[i]

        plt.plot(waves,f_watt)
        plt.xlim(3000,10000)
        plt.ylim(0,max(f_watt))
        #plt.plot(waves,summed*40000)
        plt.xlabel('wavelength (A)')
        plt.ylabel('flux (W)')
        plt.show()

plt.show()   