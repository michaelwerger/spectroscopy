import os
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.interpolate import interp1d
#import ipywidgets as wdg  # Using the ipython notebook widgets
import sys

def on_press(event):
    print('you pressed', event.button, event.xdata, event.ydata)
    plt.plot(event.xdata,event.ydata,'o')
    xpos.append(event.xdata)
    ypos.append(event.ydata)
    plt.plot(xpos,ypos,color='red')


figure_width = 15
figure_height = 12

xpos = []
ypos = []

    

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
            xpos = float(line)
            xpositions.append(xpos)

    ypositions = []
    with open('ypos.dat','r') as f:
        for line in f:
            ypos = float(line)
            ypositions.append(ypos)
    plt.plot(xpositions,ypositions, '.', color='green')



cid = fig.canvas.mpl_connect('button_press_event', on_press)
plt.ylim(0,0.0003)
plt.show()