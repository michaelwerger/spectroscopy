import os
class Linetable():


    def __init__(self, filename='linetable.csv'):

        self.filename = filename
        self.lines = []
        self.wavelengths = []

        if os.path.exists(self.filename):
            with open(self.filename,'r') as linetable_f:
                for line in linetable_f:
                   
                    try:
                        if line.startswith('#'):
                            pass
                        elif ',' in line:
                            tokens = line.split(',')
                            label, wave = tokens[0], tokens[1]
                            self.wavelengths.append(float(wave))
                            self.lines.append(label)
                        elif ';' in line:
                            tokens = line.split(';')
                            label, wave = tokens[0], tokens[1]
                            self.wavelengths.append(float(wave))
                            self.lines.append(label)
                        else:
                            pass                    
                        
                    except ValueError:
                        pass # ingore lines which are not correctly formatted, e.g. missing ','
            if len(self.lines) < 2:
                raise ValueError('No lines read')
            else:
                print ("%d lines read" % (len(self.lines)))
        else:
            raise ValueError('line table not initalized')

    def get_lines(self):
        return self.lines
        
    def get_wavelengths(self):
        return self.wavelengths