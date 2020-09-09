import pandas as pd
import numpy as np

def loadDenseMatrix(filename):
    #print(filename)
    data = pd.read_csv(filename, delimiter='\t', index_col=0)
    return data

def loadWithNorm(filename,method= "RPM",log = False):
    data = pd.read_csv(filename, delimiter='\t', index_col=0)
    if method == "RPM":
        data = (10000000 * data) / np.nansum(data)
    if log:
        return np.log1p(data)
    else:
        return data
