import pandas as pd
import numpy as np

def loadDenseMatrix(filename,log=True):
    #print(filename)
    data = pd.read_csv(filename, delimiter='\t', index_col=0)
    if log == True:
        return(np.log1p(data))
    else: return data

def loadWithNorm(filename,method= "RPM",log = False):
    data = pd.read_csv(filename, delimiter='\t', index_col=0)
    if method == "RPM":
        data = (10000000 * data) / np.nansum(data)
    if log:
        return np.log1p(data)
    else:
        return data
