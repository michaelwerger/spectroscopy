from scipy.signal import find_peaks
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import platform
from figuresize import FigureSize


class Slit(object):

    @staticmethod
    def determine_slit_rows(data, window_size=39, order=4, distance=100):
        '''
        data: nd.array with calibration lamp measurement

        return lower and upper row of slit(s)
        '''

        plt.rcParams['figure.figsize'] = FigureSize.NARROW

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
        plt.rcParams['figure.figsize'] = FigureSize.NARROW
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