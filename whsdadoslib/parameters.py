import json
import os
import shutil
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation
from .paths import Paths

class CCDParameters():
    #xsize = 4656
    #ysize = 3520
    xsize = 4944
    ysize = 3284

    @staticmethod
    def __init__(xsize=4656, ysize=3520):
        CCDParameters.xsize = xsize
        CCDParameters.ysize = ysize


class FigureSize():
    THIN     = [24, 3]
    NARROW   = [24, 8]
    MEDIUM   = [24, 15]
    LARGE    = [24, 24]

class Observatory(object):

    fname = 'observatory.json'

    def __init__(self,name):


        with open(os.path.join(Paths.json_path,Observatory.fname),"r") as f:
            self.name = name
            self.records = json.load(f)        

    def get_earth_location(self):
        return EarthLocation(
            lat=self.records[self.name]['earth_location']['latitude']*u.deg,
            lon=self.records[self.name]['earth_location']['longitude']*u.deg,
            height=self.records[self.name]['earth_location']['height']*u.m
           
        )

    def get_name(self):
        return self.records[self.name]['name']

    def get_utc_offset(self):
        return 1*u.hour

    def dict(self):

        return {
            'name':self.records[self.name]['name'],
            'earth_location': EarthLocation(
                lat=self.records[self.name]['earth_location']['latitude']*u.deg, 
                lon=self.records[self.name]['earth_location']['longitude']*u.deg, 
                height=self.records[self.name]['earth_location']['height']*u.m),
            'utc_offset':1*u.hour
        }
    
    def __str__(self):
        
        if self.records[self.name]['earth_location']['latitude'] > 0:
            northsouth = 'N'
        else:
            northsouth = 'S'

        if self.records[self.name]['earth_location']['longitude'] > 0:
            eastwest = 'E'
        else:
            eastwest = 'W'
        s = '{:s} @ {:f} {:s} {:f} {:s} '.format(
            self.records[self.name]['name'],
            self.records[self.name]['earth_location']['latitude'],
            northsouth,
            self.records[self.name]['earth_location']['longitude'],
            eastwest
            )
        return s

class Parameters(object):

     slit_rows = []
     lower_pixel = 0
     upper_pixel = 0
     detection_level =0.0
     detection_windowsize = 19
     detection_order = 4
     detection_distance = 8
     flip = False
     dark = None # either float or 2d array
     processed = True # if true, instrument calibration is determined already
     vmin, vmax = 0, 100

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

    whs = Observatory()

    print (whs)
    print (whs.get_earth_location())
    print (whs.dict())

if __name__ == "__main__":
    
    main()      

