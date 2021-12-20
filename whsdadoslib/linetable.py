import os
class Linetable():


    def __init__(self):

        self.filename = 'linetable.csv'
        self.table = dict()

        if os.path.exists(self.filename):
            with open(self.filename,'r') as linetable_f:
                for line in linetable_f:
                    try:
                        label, wave = line.split(',')
                    
                        self.table[label] = float(wave)
                    except ValueError:
                        pass # ingore lines which are not correctly formatted, e.g. missing ','
        else:
            raise ValueError('line table not initalized')

    def get_lines(self):
        return self.table.keys()
        
    def get_wavelength(self,label):
        return self.table[label]