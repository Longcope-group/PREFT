import numpy as np
from scipy.io import readsav


class Tarr:
    def __init__(self, file):
        save_d = readsav(file,verbose=False,python_dict=False)
        self._recarray = save_d["tarr"]
        self.variables = self._recarray.dtype.names

    def __getattr__(self, field_name):
        return np.stack(self._recarray[field_name])
    
    def print_var(self):
        print(self.variables)
