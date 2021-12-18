
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import platform
from scipy.interpolate import interp1d
sys.path.append('/Users/Micha/Workspaces/python/spectroscopy/whsdadoslib')
from Calibration.wavelength import WavelengthCalibration

class InstrumentFunction(object):

    xpos = []
    ypos = []

    @staticmethod
    def on_press(event):
        print('you pressed', event.button, event.xdata, event.ydata)
        plt.plot(event.xdata,event.ydata,'o')
        InstrumentFunction.xpos.append(event.xdata)
        InstrumentFunction.ypos.append(event.ydata)
        plt.plot(InstrumentFunction.xpos,InstrumentFunction.ypos,color='red')

    @staticmethod
    def fit_instrumentfunction():

        figure_width = 15
        figure_height = 12

        

            

        waves = []
        with open('waves.dat','r') as f:
            for line in f:
                wave = float(line)
                waves.append(wave)

        instr = []
        with open('instr.dat','r') as f:
            for line in f:
                inst = float(line)
                instr.append(inst)
        max_i = max(instr)

        print ('click (xpos,ypos) pairs with mouse; press [Enter] to leave.')

        fig, ax = plt.subplots()

        plt.plot(waves,instr)

        if os.path.exists('xpos.dat') and os.path.exists('ypos.dat'):
            xpositions = []
            with open('xpos.dat','r') as f:
                for line in f:
                    InstrumentFunction.xpos = float(line)
                    xpositions.append(InstrumentFunction.xpos)

            ypositions = []
            with open('ypos.dat','r') as f:
                for line in f:
                    InstrumentFunction.ypos = float(line)
                    ypositions.append(InstrumentFunction.ypos)
            plt.plot(xpositions,ypositions, '.', color='green')



        cid = fig.canvas.mpl_connect('button_press_event', InstrumentFunction.on_press)
        plt.ylim(0,0.0003)
        plt.show()




    @staticmethod
    def prepare_instrumentfunction(traces, std_wavlengths, std_fluxes, offset=0):

        waves = WavelengthCalibration.get_waves(0,len(traces))

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
    def get_factors(icl, xlim=None, ylim=None):

        factors = []
        max_i = 0

        if xlim is None:
            xlim = [0, -1]
        if ylim is None:
            ylim = [0, -1]
        for hdu in icl.hdus():
            traced = hdu.data[ylim[0]:ylim[1],:].sum(axis=0)
            max_i = max(max_i,max(traced))

        for hdu in icl.hdus():
            traced = hdu.data[ylim[0]:ylim[1],:].sum(axis=0)
            factor = max_i/ max(traced)
            factors.append(factor)
        
        return factors




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

        z = WavelengthCalibration.get_z()
        
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