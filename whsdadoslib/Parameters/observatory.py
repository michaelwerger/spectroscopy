import json
import astropy.units as u
from astropy.coordinates import EarthLocation

class Observatory(object):

    def __init__(self,fname):
        with open(fname,"r") as f:
            self.record = json.load(f)        

    def get_earth_location(self):
        return EarthLocation(
            lat=self.record['earth_location']['latitude']*u.deg,
            lon=self.record['earth_location']['longitude']*u.deg,
            height=self.record['earth_location']['height']*u.m
           
        )

    def dict(self):

        return {
            'name':self.record['name'],
            'earth_location': EarthLocation(
                lat=self.record['earth_location']['latitude']*u.deg, 
                lon=self.record['earth_location']['longitude']*u.deg, 
                height=self.record['earth_location']['height']*u.m),
            'utc_offset':1*u.hour
        }
    
    def __str__(self):
        
        if self.record['earth_location']['latitude'] > 0:
            northsouth = 'N'
        else:
            northsouth = 'S'

        if self.record['earth_location']['longitude'] > 0:
            eastwest = 'E'
        else:
            eastwest = 'W'
        s = '{:s} @ {:f} {:s} {:f} {:s} '.format(
            self.record['name'],
            self.record['earth_location']['latitude'],
            northsouth,
            self.record['earth_location']['longitude'],
            eastwest
            )
        return s

def main():

    whs = Observatory()

    print (whs)
    print (whs.get_earth_location())
    print (whs.dict())

if __name__ == "__main__":
    
    main()      