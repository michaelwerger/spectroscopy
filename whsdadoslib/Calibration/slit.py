from scipy.signal import find_peaks
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import platform
from ..parameters import FigureSize, Parameters
from ..Data.filter import Filter
class Slit(object):

    @staticmethod
    def determine_slit_rows(data, window_size=39, order=4, distance=100, verbose=True):
        '''
        data: nd.array with calibration lamp measurement

        return lower and upper row of slit(s)
        '''

        if verbose:
            plt.rcParams['figure.figsize'] = FigureSize.NARROW

            print ('full image')
            fig, ax = plt.subplots()
            max_i = data.max()
            plt.imshow(data, vmin=0, vmax=max_i/2)
            plt.show()

            print ('trace along columns')

        # compute a trace through the data along columns
        trace = data.sum(axis=1)

        if verbose:
            fig, ax = plt.subplots()
            plt.plot(trace)
            plt.show()

            print ('smoothed trace')
    
        filtered = Filter.savitzky_golay(trace,window_size=window_size,order=order)

        if verbose:
            fig, ax = plt.subplots()
            plt.plot(filtered)
            plt.show()

            print ('maximum positions of abs(1st derivative)')

        # compute abs(1st derivative) and determine max. positions

        diffed = np.diff(filtered,1)
        max_i = max(abs(diffed))
        height = max_i * Parameters.detection_level
        peaks, _ = find_peaks(abs(diffed), height=height, distance=distance)

        if verbose:
            plt.rcParams['figure.figsize'] = FigureSize.NARROW
            fig, ax = plt.subplots()
            plt.plot(abs(diffed),color='blue')
            plt.plot(peaks, abs(diffed[peaks]), "x")

            plt.show()

            # plot slit rows over ccd image

            fig, ax = plt.subplots()
            plt.imshow(data, vmin=0, vmax=max_i/2)
            plt.show()

            # if len(peaks) >= 2:
            #     fig, ax = plt.subplots()
            #     plt.imshow(data, vmin=0, vmax=max_i/2)
            #     plt.plot([Parameters.lower_pixel,Parameters.upper_pixel],
            #         [peaks[0],peaks[0]], color='red')
            #     plt.plot([Parameters.lower_pixel,Parameters.upper_pixel],
            #         [peaks[1],peaks[1]], color='red')
            #     if len(peaks) >= 4:
            #         plt.plot([Parameters.lower_pixel,Parameters.upper_pixel],
            #             [peaks[2],peaks[2]], color='red')
            #         plt.plot([Parameters.lower_pixel,Parameters.upper_pixel],
            #             [peaks[3],peaks[3]], color='red')
            #     if len(peaks) >= 6:
            #         plt.plot([Parameters.lower_pixel,Parameters.upper_pixel],
            #             [peaks[4],peaks[4]], color='red')
            #         plt.plot([Parameters.lower_pixel,Parameters.upper_pixel],
            #             [peaks[5],peaks[5]], color='red')
            print ('slit rows %s' % (str(peaks)))
        return peaks