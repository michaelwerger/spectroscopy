import os
import platform
import json

class Paths(object):

    src_path = '#' # directory with src files
    fits_path = '#' # directory with FITS data
    product_path = '#' # directory for product files
    corr_path = '#' # directory for correction files
    tab_path = '#' # directory with tab files
    products_path = '#' # direcotry with products
    overwrite = False # ensure to set to True if paths are overwritten

    
    
    # directory with extinction data
    if platform.node().split('.')[0] in ['Venus','Diana']:
        json_path = '/Users/Micha_1/Workspaces/python/spectroscopy/whsdadoslib'
        ext_path = os.path.join('/Users', 'Micha_1', 'data', 'ref', 'atmosphericextinction')
        tab_path = os.path.join('/Users', 'Micha_1', 'data', 'ref', 'hst')
    elif 'goe-nb-wmi' in platform.node().lower():
        json_path = '/Users/Micha_1/Workspaces/python/spectroscopy/whsdadoslib'
        ext_path = os.path.join('d:/','Workspaces','data','ref','atmosphericextinction')
        tab_path = os.path.join('d:/','Workspaces','data','ref','hst')
    else:
        raise ValueError('some paths are not defined; '+platform.node()+' unknown; likely not working')
    
    paths_json_file = os.path.join(json_path,'Paths.json')
    
    @staticmethod
    def dump():
        print (Paths.src_path)
        print (Paths.fits_path)
        print (Paths.product_path)
        print (Paths.corr_path)
        print (Paths.tab_path)
        
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

        print (paths)
        for key in paths.keys():
            print (key, paths[key])
        if Paths.overwrite == True:
            with open(Paths.paths_json_file,'w') as f:
                json.dump(paths,f)
            print (Paths.paths_json_file+" written")

    @staticmethod
    def _read():
        
        if os.path.exists(Paths.paths_json_file):
            with open(Paths.paths_json_file,'r') as f:
                paths = json.load(f)
            print ("hh")
            print (paths)
            Paths.src_path = paths['src_path']
            Paths.fits_path = paths['fits_path']
            Paths.corr_path = paths['corr_path']
            Paths.tab_path = paths['tab_path']
            Paths.products_path = paths['products_path']
            Paths.fits_path = paths['fits_path']
        else:
            raise FileNotFoundError("Paths.json does not exist")

def main():
    
    Paths.dump()
    Paths.write()
    Paths.overwrite = True
    Paths.write()
    
    Paths._read()

    


if __name__ == "__main__":
    
    main()   
