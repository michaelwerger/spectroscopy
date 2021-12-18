import os
import sys
sys.path.append('/Users/Micha/Workspaces/python/spectroscopy/whsdadoslib')
from Parameters.paths import Paths
class FluxCalibration(object):
    
    @staticmethod
    def get_standard_flux(filename):
        
        current_dir = os.getcwd()
        
        os.chdir(Paths.tab_path)
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
