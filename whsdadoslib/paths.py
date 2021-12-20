import os
import platform
import json

class Paths(object):

    fits_path = '' # directory with FITS data
    src_path = '' # directory with src files
    product_path = '' # directory for product files
    json_path = '/Users/Micha/Workspaces/python/spectroscopy/whsdadoslib'
    overwrite = False # ensure to set to True if paths are overwritten
    
    # directory with extinction data
    if 'Venus' in platform.node():
        ext_path = os.path.join('/Users', 'Micha', 'data', 'ref', 'atmosphericextinction')
    elif 'goe-nb-wmi' in platform.node().lower():
        ext_path = os.path.join('d:/','Workspaces','data','ref','atmosphericextinction')
    else:
        raise ValueError('ext_path not defined; extinction correction not working')
    
    # reference data files
    if 'Venus' in platform.node():
        tab_path = os.path.join('/Users', 'Micha', 'data', 'ref', 'hst')
    elif 'goe-nb-wmi' in platform.node().lower():
        tab_path = os.path.join('d:/','Workspaces','data','ref','hst')
    
    else:
        raise ValueError('tab_path not defined; flux calibration not working')
    
    @staticmethod
    def write():

        paths = {
            'src_path': Paths.src_path,
            'fits_path': Paths.fits_path,
            'corr_path': Paths.corr_path,
            'tab_path': Paths.tab_path,

            'products_path': Paths.product_path,
            'overwrite':Paths.overwrite,
        }

        if Paths.overwrite == True:
            with ('Paths.json','w') as f:
                json.dump(paths,f)

    @staticmethod
    def _read():
        with ('Paths.json','r') as f:
            paths = json.load(f)
        Paths.src_path = paths['src_path']
        Paths.fits_path = paths['fits_path']
        Paths.corr_path = paths['corr_path']
        Paths.tab_path = paths['tab_path']
        Paths.products_path = paths['products_path']
        Paths.fits_path = paths['fits_path']
