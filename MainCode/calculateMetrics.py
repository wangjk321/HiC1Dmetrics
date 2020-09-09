from loadfile import *
import numpy as np
import pandas as pd
import argparse
from scipy import ndimage

class BasePara:
    def __init__(self,path,resolution,chromosome,out_name="noName"):
        self.path = path
        self.matrix = loadDenseMatrix(path).values
        self.matrix_shape = self.matrix.shape[0]
        self.resolution = resolution
        self.chromosome = chromosome
        self.out_name = out_name

    def makeDF(self,array,metrics_name="unknown metrics"):
        array = np.round(array,6)
        df = pd.DataFrame({metrics_name:array})
        df["chr"] = self.chromosome
        df["start"] = df.index * self.resolution
        df["end"] = df["start"] + self.resolution
        df = df.loc[:,["chr","start","end",metrics_name]]
        return df

    def makeCSV(self,df):
        print("--------Export to csv file--------- \n")
        df.to_csv(self.out_name + ".bedGraph", sep="\t", header=False, index=False)

class InsulationScore(BasePara):
    def __init__(self,path,resolution,chromosome,out_name="InsulationScore",square_size=150000):
        super().__init__(path,resolution,chromosome,out_name)
        self.square_size = square_size
        #The default size in Homer IS is 150000

    #print the instanced class
    def __str__(self):
        return "InsulationScore(matrix_path = '" + str(self.path) +"',\n" + \
                "resolution = " + str(self.resolution) + ", chromosome = '" + \
                self.chromosome + "', " + "out_name = '" + self.out_name + \
                "' \n, square_size = " + str(self.square_size) + ")"
    __repr__ = __str__

    def getIS(self):
        squareBin = round(self.square_size/self.resolution)

        array = np.zeros(self.matrix_shape)
        for i in range(self.matrix_shape):
            if(i - squareBin < 0 or i + squareBin >= self.matrix_shape): continue
            score = self.matrix[i-squareBin: i, i+1: i+squareBin+1].mean()
            if np.isnan(score): continue  #skip NaN value
            array[i] = score
        array = np.log1p(array/np.nanmean(array)) #normalization

        return super().makeDF(array,"InsulationScore")

    def getCSV(self):
        super().makeCSV(self.getIS())

class ContrastIndex(BasePara):
    def __init__(self,path,resolution,chromosome,out_name="ContrastIndex",CI_size=150000):
        super().__init__(path,resolution,chromosome,out_name)
        self.CI_size = CI_size

    #print the instanced class
    def __str__(self):
        return "ContrastIndex(matrix_path = '" + str(self.path) +"',\n" + \
                "resolution = " + str(self.resolution) + ", chromosome = '" + \
                self.chromosome + "', " + "out_name = '" + self.out_name + \
                "' \n, CI_size = " + str(self.CI_size) + ")"
    __repr__ = __str__

    def getCI(self):
        CI_Bin = round(self.CI_size/self.resolution)

        array = np.zeros(self.matrix_shape)
        for i in range(self.matrix_shape):
            if(i - CI_Bin < 0 or i + CI_Bin >= self.matrix_shape): continue
            matA = self.matrix[i-CI_Bin:i,i-CI_Bin:i]
            matB = self.matrix[i+1:i+CI_Bin+1,i+1:i+CI_Bin+1]
            A = np.triu(matA,1).sum()
            B = np.triu(matB,1).sum()
            C = self.matrix[i-CI_Bin: i, i+1: i+CI_Bin+1].sum()
            if np.isnan(A) or np.isnan(B) or np.isnan(C): continue  #skip NaN value
            array[i] = np.log1p(A+B) - np.log1p(C)

        return super().makeDF(array,"ContrastIndex")

    def getCSV(self):
        super().makeCSV(self.getCI())

class DirectionalityIndex(BasePara):
    def __init__(self,path,resolution,chromosome,out_name="noName",DI_distance=1000000):
        super().__init__(path,resolution,chromosome,out_name)
        self.DI_distance = DI_distance
        #The default distance in Homer DI is 1000000

    #print the instanced class
    def __str__(self):
        return "DirectionalityIndex(matrix_path = '" + str(self.path) +"',\n" + \
                "resolution = " + str(self.resolution) + ", chromosome = '" + \
                self.chromosome + "', " + "out_name = '" + self.out_name + \
                "' \n, DI_distance = " + str(self.DI_distance) + ")"
    __repr__ = __str__

    def getDI(self):
        distanceBin = round(self.DI_distance/self.resolution)

        array = np.zeros(self.matrix_shape)
        for i in range(self.matrix_shape):
            if(i - distanceBin < 0 or i + distanceBin >= self.matrix_shape): continue
            A = self.matrix[i,i-distanceBin:i].sum()
            B = self.matrix[i,i+1:i+distanceBin+1].sum()
            E = (A+B)/2

            if np.isnan(A) or np.isnan(B) or (B-A == 0): continue #skip NaN value
            sign = (B-A)/abs(B-A)
            array[i] = sign * np.log1p(((A-E)**2)/E + ((B-E)**2)/E) #*****?

        return super().makeDF(array,"DirectionalityIndex")

    def getCSV(self):
        super().makeCSV(self.getDI())


class DirectionalRelativeFreq(BasePara):
    #compare two metrices need normalization
    def __init__(self,path,control_path,resolution,chromosome,out_name="DRF",
                start_distance=0, end_distance=2000000):
        super().__init__(path,resolution,chromosome,out_name)
        self.matrix = loadWithNorm(self.path,log = True).values
        self.control_path = control_path
        self.control_matrix = loadWithNorm(self.control_path,log = True).values
        self.start_distance = start_distance
        self.end_distance = end_distance

    def __str__(self):
        return "DirectionalRelativeFreq(matrix_path = '" + str(self.path) +"',\n" + \
                "control_path = '" + str(self.control_path) + "',\n" + \
                "resolution = " + str(self.resolution) + ", chromosome = '" + \
                self.chromosome + "', " + "out_name = '" + self.out_name + \
                "' \n, start_distance = " + str(self.start_distance) + \
                ", end_distance = " + str(self.end_distance) + ")"
    __repr__ = __str__

    def getDRF(self):
        if self.matrix.shape[0] != self.control_matrix.shape[0]:
            print("Error: Input/control matrix have different shape")
            exit(1)

        smooth = 3
        control = ndimage.median_filter(self.control_matrix,smooth)
        logratio = ndimage.median_filter(self.matrix - control, smooth)

        array = np.zeros(self.matrix_shape)
        startDistanceBin = int(self.start_distance / self.resolution)+1
        endDistanceBin = int(self.end_distance / self.resolution)

        for i in range(endDistanceBin, self.matrix_shape - endDistanceBin):
            right = logratio[i+startDistanceBin:i+endDistanceBin+1, i].mean()
            left = logratio[i, i-endDistanceBin:i-startDistanceBin+1].mean()
            if np.isnan(right - left):continue
            array[i] = right - left


        return super().makeDF(array,"DirectionalRelativeFreq")

    def getCSV(self):
        super().makeCSV(self.getDRF())
