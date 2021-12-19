from scipy.signal import find_peaks
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import platform
from scipy.interpolate import interp1d

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
