import pandas as pd
import numpy as np
import os
import random

def loadDenseMatrix(filename,log=False):
    #print(filename)
    if not os.path.exists(filename):
        print("File does not exist, please check your file path")
        exit(1)

    print("loading data...")
    try:
        data = pd.read_csv(filename, delimiter='\t', index_col=0)
    except:
        print('Error input matrix, please check https://h1d.readthedocs.io/en/latest/overview.html#input-format')
        exit(1)
    print("loading finished")

    if log == True:
        return(np.log1p(data))
    else:
        return data

def loadWithNorm(filename,method= "RPM",log = False):
    print("loading data...")
    data = pd.read_csv(filename, delimiter='\t', index_col=0)
    if method == "RPM":
        data = (10000000 * data) / np.nansum(data)
    print("loading finished")
    
    if log:
        return np.log1p(data)
    else:
        return data

def hic2matrix(path,res,chr,gt):
    if not gt:
        print("rawhic require genome_table file");exit(1)
    try:
        gtfile = pd.read_csv(gt,sep="\t",header=None)
        if not int(gtfile.iloc[0,1]):
            exit(1)
    except:
        print("Wrong genome_table file.")
        print("Please check your genome_table file. Is it tab separated ?")
        exit(1)

    print("Start dump matrix from hic file")
    codepath = os.path.dirname(os.path.realpath(__file__))
    makeIntra = codepath+"/extract/makeMatrixIntra.sh"
    juicer = codepath+"/jc/jctool_1.11.04.jar"
    foldername = "./MatrixTemp"+str(random.random())
    os.system("bash "+makeIntra+" "+"KR"+" "+"."+" "+path+" "+str(res)+" "+gt+" "+juicer+" "+chr+" "+foldername + "> info.txt")
    matrixpath = foldername+"/"+str(res)+"/observed.KR."+chr+".matrix.gz"
    print("Finish dump")
    return(matrixpath)

def cool2matrix(path,res,chr,gt):
    print("Start dump matrix from cool file")
    if not gt:
        print("cool require genome_table file");exit(1)
    codepath = os.path.dirname(os.path.realpath(__file__))
    makeIntra = codepath+"/extract/coolerdump.sh"
    foldername = "./MatrixTemp"+str(random.random())
    os.system("bash "+makeIntra+" "+path+" "+str(res)+" "+chr+" "+gt+" "+foldername + "> info.txt")
    matrixpath = foldername+"/"+str(res)+"/"+chr+".matrix.gz"
    print("Finish dump")
    return(matrixpath)
