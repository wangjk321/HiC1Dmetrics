from loadfile import *
import numpy as np
import pandas as pd
import argparse
from scipy import ndimage
from sklearn.decomposition import PCA
from sklearn import preprocessing

class BasePara:
    def __init__(self,path,resolution,chromosome,out_name="noName"):
        self.path = path
        #self.matrixNA = loadDenseMatrix(path).values
        self.matrix = np.nan_to_num(loadDenseMatrix(path).values)
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
    def __init__(self,path,resolution,chromosome,out_name="ContrastIndex",CI_size=300000):
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

class SeparationScore(BasePara):
    def __init__(self,path,resolution,chromosome,out_name="SeparationScore",TADss_size=300000):
        super().__init__(path,resolution,chromosome,out_name)
        self.TADss_size = TADss_size

    def __str__(self):
        return "SeparationScore(matrix_path = '" + str(self.path) +"',\n" + \
                "resolution = " + str(self.resolution) + ", chromosome = '" + \
                self.chromosome + "', " + "out_name = '" + self.out_name + \
                "' \n, TADss_size = " + str(self.TADss_size) + ")"
    __repr__ = __str__

    def getTADss(self):
        TADss_Bin = round(self.TADss_size/self.resolution)

        array = np.zeros(self.matrix_shape)
        for i in range(self.matrix_shape):
            if(i - TADss_Bin < 0 or i + TADss_Bin >= self.matrix_shape): continue
            matA = self.matrix[i-TADss_Bin:i,i-TADss_Bin:i]
            matB = self.matrix[i+1:i+TADss_Bin+1,i+1:i+TADss_Bin+1]
            A = np.triu(matA,1).sum()
            B = np.triu(matB,1).sum()
            C = self.matrix[i-TADss_Bin: i, i+1: i+TADss_Bin+1].sum()
            if np.isnan(A) or np.isnan(B) or np.isnan(C) or (min(A,B)==0): continue
            array[i] = C/(min(C+A,C+B))

        return super().makeDF(array,"SeparationScore")

    def getCSV(self):
        super().makeCSV(self.getTADss())

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

class DistalToLocalRatio(BasePara):
    def __init__(self,path,resolution,chromosome,out_name="noName",sizeDLR=3000000):
        super().__init__(path,resolution,chromosome,out_name)
        self.sizeDLR = sizeDLR

    def getDLR(self):
        sizeBin = round(self.sizeDLR/self.resolution)
        array = np.zeros(self.matrix_shape)

        for i in range(self.matrix_shape):
            if(i - sizeBin < 0 or i + sizeBin >= self.matrix_shape): continue
            A = self.matrix[i,i-sizeBin:i].sum()
            B = self.matrix[i,i+1:i+sizeBin+1].sum()
            Aout = self.matrix[i,0:i-sizeBin].sum()
            Bout = self.matrix[i,i+sizeBin+1:].sum()

            if np.isnan(A+B): continue #skip NaN value
            array[i] = np.log1p(Aout+Bout)-np.log1p(A+B)
        return super().makeDF(array,"DistalToLocalRatio")

    def getCSV(self):
        super.makeCSV(self.getDLR())

class CompartmentPC1(BasePara):
    #def __init__(self,path,resolution,chromosome,out_name="noName"):
    #    super().__init__(path,resolution,chromosome,out_name)
    #    self.matrix = loadWithNorm(path,log = True).values

    def makeExpect(self,df):
        num = df.shape[0]
        avg=[]
        for i in range(num):  #间隔为i，格数为num-i
            dig=[]
            for j in range(num):
                if (j+i) < num: dig.append(df[j,j+i])
            avg.append(np.mean(dig))

        for i in range(num):
            if avg[i] < 10: #25000res 为10， 50000应该是10*（50000/25000)^2
                for flank in range(num):
                    biggerBin = avg[i-flank:i+flank+1]
                    if np.sum(biggerBin)>=10:
                        avg[i] = np.mean(biggerBin)
                        break

        expected = np.zeros((num, num))
        for i in range(num):
            for j in range(num):
                if np.isnan(df[i,j]):
                    expected[i,j] = np.NaN
                else:
                    distance = abs(i-j)
                    expected[i,j] = avg[distance]

        return(expected)

    def getPC1(self):
        rawMT = np.nan_to_num(self.matrix)
        expectMT = self.makeExpect(rawMT)
        oeMT = rawMT / expectMT
        pearsonMT = np.corrcoef(oeMT)

        pca = PCA(n_components=5)
        trained = pca.fit(np.nan_to_num(pearsonMT))
        pc1 = trained.components_
        array = pc1[0,:]
        return super().makeDF(array,"CompartmentPC1")

    def getCSV(self):
        super.makeCSV(self.getPC1())
