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

#from generateCmap import *
#a = loadDenseMatrix("./Rad21KD_1/observed.KR.chr20.matrix.gz")
#import matplotlib.pyplot as plt
#plt.imshow(a,clim=(0,5),aspect=1,interpolation="nearest",cmap=generate_cmap(['#FFFFFF', '#d10a3f']))
