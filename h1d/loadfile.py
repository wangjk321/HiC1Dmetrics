import pandas as pd
import numpy as np
import os
import random

def loadDenseMatrix(filename,log=False):
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

def hic2matrix(path,res,chr,gt):
    if not gt: print("rawhic require genome_table file");exit(1)
    codepath = os.path.dirname(os.path.realpath(__file__))
    makeIntra = codepath+"/extract/makeMatrixIntra.sh"
    juicer = codepath+"/jc/jctool_1.11.04.jar"
    foldername = "./MatrixTemp"+str(random.random())
    os.system("sh "+makeIntra+" "+"KR"+" "+"."+" "+path+" "+str(res)+" "+gt+" "+juicer+" "+chr+" "+foldername)
    matrixpath = foldername+"/"+str(res)+"/observed.KR."+chr+".matrix.gz"
    return(matrixpath)
