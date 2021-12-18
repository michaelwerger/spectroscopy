import json
import os
import shutil
from astropy.coordinates import SkyCoord

class Targets(object):

    file = 'targets.json'
    targets = dict()
    d = dict()
    d = {  
        # key has the following structure - unless it is not a bright star -
        # prefix ALP, BET, GAM, DEL, EPS, ETA, THE, ZET, ....
        # suffix AQL, CAS, CEP, CMI, LYR, ORI, ...
        'ALPAQL':'alpha_aquilae',
        'GAMCAS':'gamma_cassiopeiae',
        'DELCEP':'delta_cephei',
        'ALPCMI':'alpha_canis_majoris',
        'PCYG':  'P_cygni',
        'ALPLYR':'alpha_lyrae',
        'ALPORI':'alpha_orionis',
        'BETORI':'beta_orionis',
        'GAMORI':'gamma_orionis',
        'DELORI':'delta_orionis',
        'EPSORI':'epsilon_orionis',
        'ZETORI':'zeta_orionis',
        'KAPORI':'kappa_orionis',
        'GAMUMA':'gamma_ursae_majoris',
        'M31':'M31',
        'M42':'M42',
    }

    @staticmethod
    def write():
        for k in Targets.d:

            print (k),
            skycoord = SkyCoord.from_name(str(Targets.d[k]))
            print (skycoord.ra, skycoord.dec)
            _tmp = {
                'fullname':str(Targets.d[k]),
                'ra.hms':str(skycoord.ra),
                'dec.dms':str(skycoord.dec),
                'frame':'icrs'
            }
            print (_tmp)
            Targets.targets[k] = _tmp
        if os.path.exists(Targets.file):
            _bak_file = Targets.file+'.bak'
            if os.path.exists(_bak_file):
                os.remove(_bak_file)
            shutil.copyfile(Targets.file,_bak_file)
        with open(Targets.file,"w+") as f:
            json.dump(Targets.targets,f)

    @staticmethod
    def read(source_dir='.'):

        current_dir = os.getcwd()
            
        try:
            os.chdir(source_dir)
            with open(Targets.file,"r") as f:
                Targets.targets = json.load(f)
        finally:
            os.chdir(current_dir)

        return Targets.targets

    @staticmethod
    def update():
        if Targets.targets is None or Targets.targets == {}:
            Targets.targets = Targets.read()
        for key in Targets.targets.keys():
            Targets.targets[key]['radec'] = SkyCoord(
                Targets.targets[key]['ra.hms'], 
                Targets.targets[key]['dec.dms'], 
                frame=Targets.targets[key]['frame'])
        return Targets.targets

    @staticmethod
    def __str__():
        if Targets.targets is None or Targets.targets == {}:
            print ("Not defined; try Targets.read() first")
            raise ValueError
        else:
            for target in targets:
                print ("%s %s " % (target,
                    str(targets[target])
                    ))

    @staticmethod
    def from_filename(filename):
        r = None
        if Targets.targets is None or Targets.targets == {}:
            print ("Not defined; try Targets.read() first")
            raise ValueError
        keys = Targets.targets.keys()
        
        for key in keys:
            if key in filename:
                r = key
                break
        return r
        
        
def main():

    targets = Targets()
    #targets.write()
    targets.read()
    print (targets.targets)    

if __name__ == "__main__":
    
    main()      