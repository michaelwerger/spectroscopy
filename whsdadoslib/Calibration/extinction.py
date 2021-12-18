import os
import sys

sys.path.append('/Users/Micha/Workspaces/python/spectroscopy/whsdadoslib')


class ExtinctionCorrection(object):



    @staticmethod
    def get_extinctiondata(extinctionfile='paranal_extinction.dat'):
        wavelengths = []
        extcoeffs = []
        currentdir = os.getcwd()
        try:
            os.chdir(Paths.ext_path)
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
        raise NotImplementedError('sorry')