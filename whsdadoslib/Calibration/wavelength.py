from scipy.signal import find_peaks
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import platform
from scipy.interpolate import interp1d
from math import sqrt

class WavelengthCalibration(object):
    

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

    @staticmethod
    def get_inverse(waves):

        z = WavelengthCalibration.get_z()


        x1 = [-z[1]/2/z[0] + sqrt(w/z[0] + z[1]*z[1]/4/z[0]/z[0] - z[2]/z[0]) for w in waves ]
        x2 = [-z[1]/2/z[0] - sqrt(w/z[0] + z[1]*z[1]/4/z[0]/z[0] - z[2]/z[0]) for w in waves ]

        return x1, x2
