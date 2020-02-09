from scipy.signal import find_peaks
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import platform
from scipy.interpolate import interp1d
sys.path.append('/Users/Micha/Workspaces/python/spectroscopy')
from whsdadoslib.data import DataProcessing
from whsdadoslib.report import Report
from whsdadoslib.show import Show

class CalibrationData(object):

    '''
    a class which determines side parameters necessary for calibration
    '''

    slit_rows = None
    lower_pixel = 0
    upper_pixel =4656
    detection_level = 0.25
    tab_path = ''
    if 'Venus' in platform.node():
        tab_path = os.path.join('/Users', 'Micha', 'data', 'ref', 'hst')
    else:
        tab_path = os.path.join('d:/','Workspaces','data','ref','hst')
    z = []

    ext_path = ''
    if 'Venus' in platform.node():
        ext_path = os.path.join('/Users', 'Micha', 'data', 'ref', 'atmosphericextinction')
    else:
        ext_path = os.path.join('d:/','Workspaces','data','ref','atmosphericextinction')

    @staticmethod
    def get_extinctiondata(extinctionfile='paranal_extinction.dat'):
        wavelengths = []
        extcoeffs = []
        currentdir = os.getcwd()
        try:
            os.chdir(CalibrationData.ext_path)
            with open(extinctionfile,'r') as f:
                for line in f:

                    tokens = line.split(' ')
                    if len(tokens) == 2:
                        wavelength, extcoeff = tokens
                        #print (wave, extcoeff)
                        wavelengths.append(float(wavelength))
                        extcoeffs.append(float(extcoeff))
        finally:
            os.chdir(currentdir)
        return (wavelengths, extcoeffs)
    
    @staticmethod
    def correct_extinction(airmass_standard, airmass_object):
        pass

    @staticmethod
    def determine_slit_rows(data, window_size=39, order=4, distance=100):
        '''
        data: nd.array with calibration lamp measurement

        return lower and upper row of slit(s)
        '''

        plt.rcParams['figure.figsize'] = [Show.figure_width, Show.figure_height]

        print ('full image')
        fig, ax = plt.subplots()
        max_i = 1000 #max(max(data))
        plt.imshow(data, vmin=0, vmax=max_i/2)
        plt.show()
        # compute a trace through the data along columns

        print ('trace along columns')

        trace = data.sum(axis=1)
        fig, ax = plt.subplots()
        plt.plot(trace)
        plt.show()

        # smooth data

        print ('smoothed trace')
        filtered = DataProcessing.savitzky_golay(trace,window_size=window_size,order=order)
        fig, ax = plt.subplots()
        plt.plot(filtered)
        plt.show()

        # compute abs(1st derivative) and determine max. positions

        print ('maximum positions of abs(1st derivative)')
        diffed = np.diff(filtered,1)
        max_i = max(abs(diffed))
        height = max_i * CalibrationData.detection_level
        peaks, _ = find_peaks(abs(diffed), height=height, distance=distance)
        plt.rcParams['figure.figsize'] = [Show.figure_width, Show.figure_height]
        fig, ax = plt.subplots()
        plt.plot(abs(diffed),color='blue')
        plt.plot(peaks, abs(diffed[peaks]), "x")

        plt.show()

        # plot slit rows over ccd image

        if len(peaks) >= 2:
            fig, ax = plt.subplots()
            plt.imshow(data, vmin=0, vmax=max_i/2)
            plt.plot([CalibrationData.lower_pixel,CalibrationData.upper_pixel],
                [peaks[0],peaks[0]], color='red')
            plt.plot([CalibrationData.lower_pixel,CalibrationData.upper_pixel],
                [peaks[1],peaks[1]], color='red')
            if len(peaks) >= 4:
                plt.plot([CalibrationData.lower_pixel,CalibrationData.upper_pixel],
                    [peaks[2],peaks[2]], color='red')
                plt.plot([CalibrationData.lower_pixel,CalibrationData.upper_pixel],
                    [peaks[3],peaks[3]], color='red')
            if len(peaks) >= 6:
                plt.plot([CalibrationData.lower_pixel,CalibrationData.upper_pixel],
                    [peaks[4],peaks[4]], color='red')
                plt.plot([CalibrationData.lower_pixel,CalibrationData.upper_pixel],
                    [peaks[5],peaks[5]], color='red')
            plt.show()
        print ('slit rows %s' % (str(peaks)))
        return peaks

    @staticmethod
    def get_standard_flux(filename):
        
        current_dir = os.getcwd()
        
        os.chdir(CalibrationData.tab_path)
        std_waves = []
        std_fluxes = []
        with open(filename, 'r') as f:
            for l in f:
                l = f.readline()
                if l.startswith('*'):
                    # ignore line
                    pass
                elif l.startswith('SET'):
                    # ignore line
                    pass
                else:
                    
                    try:
                        wave = float(l[1:12])
                        flux = float(l[13:-1])
                        if 3000.0 < wave and wave < 10000.0:
                            std_waves.append(wave)
                            #std_fluxes.append(flux*wave*wave)
                            std_fluxes.append(flux)
                    except ValueError:
                        pass
        os.chdir(current_dir)
        return (std_waves, std_fluxes)

    @staticmethod
    def check_wavelengths(trace, wavelengths_file='wavelengths.txt'):
        plt.rcParams['figure.figsize'] = [Show.figure_width, Show.figure_height]

        fig, ax = plt.subplots()
        

        z = CalibrationData.get_z()

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

    @staticmethod
    def set_z(z):
        print ('CalibrationData.z: ', z)
        CalibrationData.z = z
        with open('wavelength_solution.dat','w') as f:
            for _z in z:
                f.write(str(_z)+'\n')

    @staticmethod
    def get_z():
        z = []
        with open('wavelength_solution.dat','r') as f:
            for line in f:
                z.append(float(line))
        CalibrationData.z = z
        return z

    @staticmethod
    def prepare_instrumentfunction(traces, std_wavlengths, std_fluxes, offset=0):

        z = CalibrationData.get_z()
        print (z)
        positions = range(0,len(traces))
        if len(z) == 1:
            waves = [pos * z for pos in positions]
        elif len(z) == 2:
            waves = [pos * z[0] + z[1] for pos in positions]
        elif len(z) == 3:
            waves = [pos*pos* z[0] + pos * z[1] + z[2] for pos in positions]
        else:
            raise NotImplementedError

        f_std = interp1d(std_wavlengths, std_fluxes, fill_value="extrapolate")
        f_std_mjansky = f_std(waves)
        f_watt = f_std_mjansky
        f_std_watt = f_std_mjansky.copy()
        for i in range(0,len(f_watt)):
            f_watt[i] = f_std_mjansky[i]/waves[i]/waves[i]    

        instr = waves.copy()
        for i in range(0,len(instr)):
            instr[i] = 0
        for i in range(100,len(waves)-100):
            instr[i] = traces[i] / (f_std_watt[i+offset])
            
        max_i = max(instr[500:-500])
        plt.plot(waves,instr)
        plt.xlabel('pixel column (x-axis)')
        plt.ylabel('relative intensity')
        plt.ylim(0,max_i)
        plt.show()
        with open('waves.dat','w') as f:
                for wave in waves:
                    f.write(str(wave)+'\n')

        with open('instr.dat','w') as f:
                for inst in instr:
                    f.write(str(inst)+'\n')

    @staticmethod
    def get_instrumentfunction(instrdat='instr.dat', xpos_ypos_dat = 'xpos_ypos.dat'):
        # xpositions = []
        # with open(xposdat,'r') as f:
        #     for line in f:
        #         xpos = float(line)
        #         xpositions.append(xpos)

        # ypositions = []
        # with open(yposdat,'r') as f:
        #     for line in f:
        #         ypos = float(line)
        #         ypositions.append(ypos)

        instr = []
        xpositions = []
        ypositions = []
        with open(xpos_ypos_dat,'r') as f:
            for line in f:
                _xp, _yp = line.split(' ')
                xpositions.append(float(_xp))
                ypositions.append(float(_yp))
        with open(instrdat,'r') as f:
            for line in f:
                _instr = float(line)
                instr.append(_instr)

        fig, ax = plt.subplots()
        #plt.plot(pixelnumbers,traces)
        plt.plot(xpositions,ypositions,'o')

        z = CalibrationData.get_z()
        
        positions = range(0,len(instr))
        if len(z) == 1:
            waves = [pos * z for pos in positions]
        elif len(z) == 2:
            waves = [pos * z[0] + z[1] for pos in positions]
        elif len(z) == 3:
            waves = [pos*pos* z[0] + pos * z[1] + z[2] for pos in positions]
        else:
            raise NotImplementedError
                    
        fitted = interp1d(xpositions, ypositions, 'linear', fill_value="extrapolate")
        instr_waves = fitted(waves)
        min_i = 0.0
        max_i = max(instr[500:-500])
        plt.plot(waves,instr_waves)
        plt.plot(waves,instr)
        plt.ylim(min_i, max_i)
        plt.show()
        return (waves, instr_waves)


def main():
    
    waves, extcoeffs = CalibrationData.get_extinctiondata()
    print (extcoeffs)
    


if __name__ == "__main__":
    
    main()      