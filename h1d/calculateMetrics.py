from .loadfile import *
import argparse
from scipy import ndimage
from scipy import stats
from sklearn.decomposition import PCA
from sklearn import preprocessing
#from callDirectionalTAD import TADcallIS
import warnings
import os

class BasePara:
    def __init__(self,path,resolution,chromosome,out_name="noName",useNA=True,datatype="matrix"):
        if datatype == "rawhic":
            codepath = os.path.dirname(os.path.realpath(__file__))
            makeIntra = codepath+"/extract/makeMatrixIntra.sh"
            juicer = codepath+"/jc/jctool_1.11.04.jar"
            os.system("sh "+makeIntra+" "+"KR"+" "+"."+" "+path+" "+str(resolution)+" "+gt+" "+juicer+" "+chr)
            path = "./MatrixTemp/"+str(resolution)+"/observed.KR."+chromosome+".matrix.gz"
        self.path = path
        #self.matrixNA = loadDenseMatrix(path).values
        self.matrix = loadDenseMatrix(path,log=False).values
        self.matrix_shape = self.matrix.shape[0]
        self.resolution = resolution
        self.chromosome = chromosome
        self.out_name = out_name
        self.useNA = useNA
        self.allsum = np.nansum(self.matrix)
        if useNA == True:
            self.blankarray = np.zeros(self.matrix_shape) * np.NaN
        elif useNA == False:
            self.blankarray = np.zeros(self.matrix_shape)

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

class InteractionFrequency:
    def __init__(self,path,resolution,chromosome,gt=None,datatype="rawhic",normIF=True,out_name="noName"):
        if not gt:
            raise ValueError("Genometable is required for the calculation of IF")
        codepath = os.path.dirname(os.path.realpath(__file__))
        soft = codepath+"/InteractionFreq.sh"
        juicer = codepath+"/jc/jctool_1.11.04.jar"
        chrnum = chromosome.split("chr")[1]
        os.system("sh '"+soft+"' '"+juicer+"' "+path+" "+chrnum+" "+str(resolution)+" "+gt+" "+"IF_"+chromosome) #in case of space
        score = pd.read_csv("IF_"+chromosome+".bedGraph",sep="\t",header=None)
        if normIF:
            beforlog = score[3].copy()
            afterlog = np.log1p(beforlog)
            score[3] = afterlog / np.mean(afterlog[afterlog>0])
        score.index = range(score.shape[0])
        score.columns = ["chr","start","end","InteractionFreq"]
        os.system("rm "+"IF_"+chromosome+".bedGraph")
        self.score = score

    def getIF(self):
        return(self.score)

class InsulationScore(BasePara):
    def __init__(self,path,resolution,chromosome,out_name="InsulationScore",useNA=True,square_size=300000,datatype="matrix"):
        super().__init__(path,resolution,chromosome,out_name,useNA)
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
        array = self.blankarray

        for i in range(self.matrix_shape):
            if(i - squareBin < 0 or i + squareBin >= self.matrix_shape): continue
            score = self.matrix[i-squareBin: i, i+1: i+squareBin+1].mean()
            if np.isnan(score) or score == 0: continue  #skip NaN value
            array[i] = score
        array = np.log1p(array/np.nanmean(array)) #normalization

        return super().makeDF(array,"InsulationScore")

    def getCSV(self):
        super().makeCSV(self.getIS())

def TADcallIS(matrixPath,resolution,chromosome,squareSize=300000,useNA=True):
    from scipy.signal import argrelextrema

    ISbedgraph = InsulationScore(matrixPath,resolution,chromosome,square_size=squareSize,useNA=useNA).getIS()
    ISoneNA = ISbedgraph.InsulationScore
    ISone = pd.Series(np.nan_to_num(ISoneNA))

    # local minimal
    localMinPos = argrelextrema(np.array(ISone), np.less)
    localMinIS = ISone.iloc[localMinPos]

    # 0< IS <0.8
    localMinIS = localMinIS[localMinIS!=0]
    localMinIS = localMinIS[localMinIS< np.mean(ISoneNA)]

    #Around TAD boundary
    binNum = int(100000/resolution)
    localMinAround = []
    aroundZero=[]
    diffrightleft=[]

    for i in localMinIS.index:
        # local minimal around 100kb
        localMinAround.append(ISone.loc[i-binNum:i+binNum].min())

        #IS_around !=0
        aroundZero.append(ISone.loc[i-binNum*2:i+binNum*2].min())

        #strong boundary
        minusbin = ISone.loc[i-binNum] - ISone.loc[i]
        plusbin = ISone.loc[i+binNum] - ISone.loc[i]
        diffrightleft.append(minusbin+plusbin)

    bool1 = (np.array(localMinIS)-np.array(localMinAround)) <= 0
    bool2 = np.array(aroundZero)>0
    bool3 = np.array(diffrightleft)>0.05
    localMinIS = localMinIS[bool1 * bool2 * bool3]

    # build a table as output
    TADnumber = len(localMinIS)
    chrlist = [chromosome] * (TADnumber-1)
    TADstart = (np.array(localMinIS.index)[:-1])*resolution
    TADend = (np.array(localMinIS.index)[1:])*resolution
    TADout = pd.DataFrame()
    TADout["chr"] = chrlist
    TADout["TADstart"] = TADstart
    TADout["TADend"] = TADend
    #Maximum TAD size 5MB, Minimum 0.3MB
    TADout = TADout[(TADout["TADend"]-TADout["TADstart"])<=5000000]
    TADout = TADout[(TADout["TADend"]-TADout["TADstart"])>=300000]

    withNA=[]
    for i in range(TADout.shape[0]):
        s = np.array(TADout.TADstart)[i] // resolution
        e = np.array(TADout.TADend)[i] // resolution
        whetherNAIS = np.isnan(sum(ISoneNA[s:e+1]))
        withNA.append(~whetherNAIS)
    TADout = TADout[withNA]

    return(TADout)

class ContrastIndex(BasePara):
    def __init__(self,path,resolution,chromosome,out_name="ContrastIndex",useNA=True,CI_size=200000,datatype="matrix"):
        super().__init__(path,resolution,chromosome,out_name,useNA)
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
        array = self.blankarray

        for i in range(self.matrix_shape):
            if(i - CI_Bin < 0 or i + CI_Bin >= self.matrix_shape): continue
            matA = self.matrix[i-CI_Bin:i,i-CI_Bin:i]
            matB = self.matrix[i+1:i+CI_Bin+1,i+1:i+CI_Bin+1]
            A = np.sum(np.triu(matA,1))
            B = np.sum(np.triu(matB,1))
            C = np.sum(self.matrix[i-CI_Bin: i, i+1: i+CI_Bin+1])
            if np.isnan(np.sum(A+B+C)) or min(A,B,C) == 0: continue  #skip NaN value
            array[i] = np.log1p(A+B) - np.log1p(C)

        return super().makeDF(array,"ContrastIndex")

    def getCSV(self):
        super().makeCSV(self.getCI())

class SeparationScore(BasePara):
    def __init__(self,path,resolution,chromosome,out_name="SeparationScore",useNA=True,TADss_size=300000,datatype="matrix"):
        super().__init__(path,resolution,chromosome,out_name,useNA)
        self.TADss_size = TADss_size

    def __str__(self):
        return "SeparationScore(matrix_path = '" + str(self.path) +"',\n" + \
                "resolution = " + str(self.resolution) + ", chromosome = '" + \
                self.chromosome + "', " + "out_name = '" + self.out_name + \
                "' \n, TADss_size = " + str(self.TADss_size) + ")"
    __repr__ = __str__

    def getTADss(self):
        TADss_Bin = round(self.TADss_size/self.resolution)
        array = self.blankarray

        for i in range(self.matrix_shape):
            if(i - TADss_Bin < 0 or i + TADss_Bin >= self.matrix_shape): continue
            matA = self.matrix[i-TADss_Bin:i,i-TADss_Bin:i]
            matB = self.matrix[i+1:i+TADss_Bin+1,i+1:i+TADss_Bin+1]
            A = np.triu(matA,1).sum()
            B = np.triu(matB,1).sum()
            C = self.matrix[i-TADss_Bin: i, i+1: i+TADss_Bin+1].sum()
            if np.isnan(np.sum(A+B+C)) or (min(A,B,C)==0): continue
            array[i] = C/(min(C+A,C+B))

        return super().makeDF(array,"SeparationScore")

    def getCSV(self):
        super().makeCSV(self.getTADss())

class DirectionalityIndex(BasePara):
    def __init__(self,path,resolution,chromosome,out_name="noName",useNA=True,DI_distance=1000000,datatype="matrix"):
        super().__init__(path,resolution,chromosome,out_name,useNA)
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
        array = self.blankarray

        for i in range(self.matrix_shape):
            if(i - distanceBin < 0 or i + distanceBin >= self.matrix_shape): continue
            A = self.matrix[i,i-distanceBin:i].sum()
            B = self.matrix[i,i+1:i+distanceBin+1].sum()
            E = (A+B)/2

            if np.isnan(np.sum(A+B)) or (B-A == 0) or (min(A,B)==0): continue #skip NaN value
            sign = (B-A)/abs(B-A)
            array[i] = sign * np.log1p(((A-E)**2)/E + ((B-E)**2)/E) #*****?

        return super().makeDF(array,"DirectionalityIndex")

    def getCSV(self):
        super().makeCSV(self.getDI())

class DistalToLocalRatio(BasePara):
    def __init__(self,path,resolution,chromosome,out_name="noName",useNA=True,sizeDLR=3000000,datatype="matrix"):
        super().__init__(path,resolution,chromosome,out_name,useNA)
        self.sizeDLR = sizeDLR

    def getDLR(self):
        sizeBin = round(self.sizeDLR/self.resolution)
        array = self.blankarray

        for i in range(self.matrix_shape):
            if(i - sizeBin < 0 or i + sizeBin >= self.matrix_shape): continue
            A = self.matrix[i,i-sizeBin:i].sum()
            B = self.matrix[i,i+1:i+sizeBin+1].sum()
            Aout = np.nansum(self.matrix[i,0:i-sizeBin])
            Bout = np.nansum(self.matrix[i,i+sizeBin+1:])

            if np.isnan(A+B) or min(A,B,Aout,Bout)==0: continue #skip NaN value
            array[i] = np.log(Aout+Bout)-np.log(A+B)
        return super().makeDF(array,"DistalToLocalRatio")

    def getCSV(self):
        super.makeCSV(self.getDLR())

'''
class interTADscore(BasePara):
    def getInterS(self,IS_size=300000,useNA=True,TADpath=None):
        if TADpath:
            usedPath = TADpath
        else:usedPath = self.path
        tad = TADcallIS(usedPath,self.resolution,self.chromosome,squareSize=IS_size,useNA=useNA)
        leftBorder =  np.array(tad.TADstart) // self.resolution
        rightBorder = np.array(tad.TADend) // self.resolution
        array = self.blankarray

        for i in range(self.matrix_shape):
            belongTAD = (i >= leftBorder) * (i < rightBorder)
            if sum(belongTAD) == 0: continue

            elif np.median(np.nan_to_num(self.matrix[i,:])) == 0:
                continue

            startBin = int(leftBorder[belongTAD])
            endBin = int(rightBorder[belongTAD])
            A = np.nansum(self.matrix[i,0:startBin-1])
            B = np.nansum(self.matrix[i,endBin+1:])
            if np.isnan(A+B) or max(A,B) == 0: continue
            array[i] = A+B

        #array = np.log1p(array/np.nanmean(array))
        array = (array/self.allsum)*1e4
        return super().makeDF(array,"interTADscore")
'''

class CompartmentPC1(BasePara):
    #def __init__(self,path,resolution,chromosome,out_name="noName"):
    #    super().__init__(path,resolution,chromosome,out_name)
    #    self.matrix = loadWithNorm(path,log = True).values

    def makeExpect(self,df,smooth):
        num = df.shape[0]
        avg=[]
        for i in range(num):  #间隔为i，格数为num-i
            dig=[]
            for j in range(num):
                if (j+i) < num: dig.append(df[j,j+i])
            avg.append(np.mean(dig))

        if smooth == True:
            for i in range(num):
                if avg[i] < df.mean(): #25000res 为10， 50000应该是10*（50000/25000)^2
                    for flank in range(num):
                        biggerBin = avg[i-flank:i+flank+1]
                        if np.sum(biggerBin)>=df.mean():
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

    def getPC1(self, logOE = False, signCorr = None, smooth = True):
        rawMT = np.nan_to_num(self.matrix)
        expectMT = self.makeExpect(rawMT,smooth)
        warnings.filterwarnings("ignore")
        oeMT = np.nan_to_num(rawMT / expectMT)
        warnings.filterwarnings("default")

        if logOE == True:
            oeMT = np.log(oeMT)
            oeMT[np.isinf(oeMT)] = 0

        allzero = oeMT.sum(axis=1) == 0
        notallzeroindex = np.array(range(self.matrix_shape))[~allzero]
        #pandas 比numpy慢很多
        #allzeroindex = np.array(range(oeMT.shape[0]))[allzero]
        #oeMT[allzeroindex,:] = np.NaN #把0值变成NA
        #oeMT[:,allzeroindex] = np.NaN
        #pearsonMT = np.array(pd.DataFrame(oeMT).corr())

        nonzeroMT = oeMT[~allzero,:][:,~allzero]
        nonzeroPearson = np.corrcoef(nonzeroMT)
        pca = PCA(n_components=5)
        trained = pca.fit(nonzeroPearson)
        pc1 = trained.components_
        noNAarray = pc1[0,:]

        array = np.zeros(self.matrix_shape)
        array[notallzeroindex] = noNAarray
        #pearsonMT = pd.DataFrame(np.zeros((self.matrix_shape,self.matrix_shape)) *np.NaN)
        #pearsonMT.iloc[notallzeroindex,notallzeroindex] = nonzeroPearson
        #pearsonMT = np.array(pearsonMT)

        #naPos = np.isnan(pearsonMT).all(axis=1)
        #pca = PCA(n_components=5)
        #trained = pca.fit(np.nan_to_num(pearsonMT))
        #pc1 = trained.components_
        #array = pc1[0,:]

        if signCorr == None:
            pass
        else:
            geneDensity = pd.read_csv(signCorr,header=None,sep='\t')
            gd = geneDensity[geneDensity[0] == self.chromosome][3]
            cor2gd = stats.spearmanr(gd,array)[0]
            if cor2gd <0: array = -array

        if self.useNA == True: array[allzero] = np.NaN
        return super().makeDF(array,"CompartmentPC1")

    def getCSV(self):
        super.makeCSV(self.getPC1())

class intraTADscore(CompartmentPC1):
    def getIntraS(self,IS_size=300000,useNA=True,TADpath=None,useOE=True,smooth=False,normTAD=True):   #this useNA is for TAD calling
        if TADpath:
            usedPath = TADpath
        else:usedPath = self.path

        tad = TADcallIS(usedPath,self.resolution,self.chromosome,squareSize=IS_size,useNA=useNA)
        leftBorder =  np.array(tad.TADstart) // self.resolution
        rightBorder = np.array(tad.TADend) // self.resolution
        array = self.blankarray

        if useOE == True:
            rawMT = np.nan_to_num(self.matrix)
            expectMT = self.makeExpect(rawMT,smooth)
            warnings.filterwarnings("ignore")
            mt = np.nan_to_num(rawMT / expectMT)
            warnings.filterwarnings("default")
        else:
            mt = self.matrix

        for i in range(self.matrix_shape):
            belongTAD = (i >= leftBorder) * (i < rightBorder)
            if sum(belongTAD) == 0: continue
            #elif np.median(np.nan_to_num(self.matrix[i,:])) == 0:
            #    continue
            startBin = int(leftBorder[belongTAD])
            endBin = int(rightBorder[belongTAD])
            A = mt[i,startBin:i].sum()
            B = mt[i,i+1:endBin+1].sum()
            if np.isnan(A+B) or max(A,B) == 0: continue
            if normTAD == True:
                array[i] = (A+B)/(endBin-startBin)
            else: array[i] = A+B

        #array = np.log1p(array/np.nanmean(array))
        if useOE == False:
            array = (array/self.allsum)*1e4
        return super().makeDF(array,"IntraTADscore")

class interTADscore(CompartmentPC1):
    def getInterS(self,IS_size=300000,useNA=True,TADpath=None,useOE=True,
                logOE=False,smooth=False,normTAD=True,smoothScore=0.8):   #this useNA is for TAD calling
        if TADpath:
            usedPath = TADpath
        else:usedPath = self.path

        tad = TADcallIS(usedPath,self.resolution,self.chromosome,squareSize=IS_size,useNA=useNA)
        leftBorder =  np.array(tad.TADstart) // self.resolution
        rightBorder = np.array(tad.TADend) // self.resolution
        array = self.blankarray

        if useOE == True:
            rawMT = np.nan_to_num(self.matrix)
            expectMT = self.makeExpect(rawMT,smooth)
            warnings.filterwarnings("ignore")
            mt = np.nan_to_num(rawMT / expectMT)
            warnings.filterwarnings("default")
        else:
            mt = self.matrix

        if logOE == True:
            mt = np.log(mt)
            mt[np.isinf(mt)] = 0

        for i in range(self.matrix_shape):
            belongTAD = (i >= leftBorder) * (i < rightBorder)
            if sum(belongTAD) == 0: continue
            #elif np.median(np.nan_to_num(self.matrix[i,:])) == 0:
            #    continue
            startBin = int(leftBorder[belongTAD])
            endBin = int(rightBorder[belongTAD])
            A = np.nansum(mt[i,0:startBin-1])
            B = np.nansum(mt[i,endBin+1:])
            if np.isnan(A+B) or max(A,B) == 0: continue
            if normTAD == True:
                array[i] = (A+B)/(self.matrix_shape-(endBin-startBin))
            else: array[i] = A+B
        if smoothScore:
            from scipy.ndimage import gaussian_filter1d
            array = gaussian_filter1d(array,smoothScore)

        #array = np.log1p(array/np.nanmean(array))
        if useOE == False:
            array = (array/self.allsum)*1e4
        return super().makeDF(array,"InterTADscore")
