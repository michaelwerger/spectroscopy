from scipy.signal import find_peaks
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import platform
from scipy.interpolate import interp1d
from Parameters.figuresize import FigureSize

class WavelengthCalibration(object):
    @staticmethod
    def check_wavelengths(trace, wavelengths_file='wavelengths.txt'):
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

    @staticmethod
    def set_z(z):
        print ('CalibrationData.z: ', z)
        WavelengthCalibration.z = z
        with open('wavelength_solution.dat','w') as f:
            for _z in z:
                f.write(str(_z)+'\n')

    @staticmethod
    def get_z():
        z = []
        with open('wavelength_solution.dat','r') as f:
            for line in f:
                z.append(float(line))
        WavelengthCalibration.z = z
        return z

    @staticmethod
    def get_waves(start, stop, step=1):

        z = WavelengthCalibration.get_z()

        positions = range(start, stop, step)
        if len(z) == 1:
            waves = [pos * z for pos in positions]
        elif len(z) == 2:
            waves = [pos * z[0] + z[1] for pos in positions]
        elif len(z) == 3:
            waves = [pos*pos* z[0] + pos * z[1] + z[2] for pos in positions]
        else:
            raise NotImplementedError

        return waves
