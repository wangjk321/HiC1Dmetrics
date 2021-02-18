from .calculateMetrics import *

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

        array = np.zeros(self.matrix_shape) * np.NaN
        startDistanceBin = int(self.start_distance / self.resolution)+1
        endDistanceBin = int(self.end_distance / self.resolution)

        for i in range(endDistanceBin, self.matrix_shape - endDistanceBin):
            right = logratio[i+startDistanceBin:i+endDistanceBin+1, i].mean()
            left = logratio[i, i-endDistanceBin:i-startDistanceBin+1].mean()
            if np.isnan(right - left):continue
            array[i] = right - left


        return super().makeDF(array,"DirectionalRelativeFrequency")

    def getCSV(self):
        super().makeCSV(self.getDRF())

class TADScoreChange(BasePara):
    def __init__(self,path,control_path,resolution,chromosome,out_name=""):
        super().__init__(path,resolution,chromosome,out_name)
        self.control_path = control_path

    def getChange(self,mode,parameter):
        if mode == "IS":
            treat = InsulationScore(self.path,self.resolution,self.chromosome,square_size=parameter).getIS().InsulationScore
            control = InsulationScore(self.control_path,self.resolution,self.chromosome,square_size=parameter).getIS().InsulationScore
            title = "InsulationScore Change"
        elif mode == "CI":
            treat = ContrastIndex(self.path,self.resolution,self.chromosome,CI_size=parameter).getCI().ContrastIndex
            control = ContrastIndex(self.control_path,self.resolution,self.chromosome,CI_size=parameter).getCI().ContrastIndex
            title = "ContrastIndex Change"
        elif mode == "DI":
            treat = DirectionalityIndex(self.path,self.resolution,self.chromosome,DI_distance=parameter).getDI().DirectionalityIndex
            control = DirectionalityIndex(self.control_path,self.resolution,self.chromosome,DI_distance=parameter).getDI().DirectionalityIndex
            title = "DirectionalityIndexChange"
        elif mode == "TADss":
            treat = SeparationScore(self.path,self.resolution,self.chromosome,TADss_size=parameter).getTADss().SeparationScore
            control = SeparationScore(self.control_path,self.resolution,self.chromosome,TADss_size=parameter).getTADss().SeparationScore
            title = "SeparationScore Change"

        change = np.array(treat - control)
        return super().makeDF(change,title)

    def getCSV(self):
        super().makeCSV(self.getISC())

class deltaDLR(BasePara):
    def __init__(self,path,control_path,resolution,chromosome,out_name="deltaDLR",sizeDLR=3000000):
        super().__init__(path,resolution,chromosome,out_name)
        self.control_path = control_path
        self.sizeDLR = sizeDLR
        self.DLRtreat = DistalToLocalRatio(path,resolution,chromosome,out_name="DLR",
                                            sizeDLR=self.sizeDLR).getDLR()
        self.DLRcontrol= DistalToLocalRatio(control_path,resolution,chromosome,out_name="DLR",
                                            sizeDLR=self.sizeDLR).getDLR()

    def getDeltaDLR(self):
        dDLR = np.array(self.DLRtreat.DistalToLocalRatio - self.DLRcontrol.DistalToLocalRatio)
        return super().makeDF(dDLR,"delta-DistalToLocalRatio")

    def getCSV(self):
        super().makeCSV(self.getDeltaDLR())

class PC1change(BasePara):
    def __init__(self,path,control_path,resolution,chromosome,out_name="",useNA=True,corr_file="",smoothPC=True,logPC=False):
        super().__init__(path,resolution,chromosome,out_name,useNA)
        self.control_path = control_path
        self.corr_file = corr_file
        self.PC1treat = CompartmentPC1(path,resolution,chromosome,
                                    out_name="PC1").getPC1(signCorr=self.corr_file,smooth = smoothPC, logOE=logPC)
        self.PC1control= CompartmentPC1(control_path,resolution,chromosome,
                                    out_name="PC1").getPC1(signCorr=self.corr_file,smooth = smoothPC, logOE=logPC)

    def getPC1change(self):
        treat = self.PC1treat.CompartmentPC1
        control = self.PC1control.CompartmentPC1
        array = self.blankarray
        for i in range(self.matrix_shape):
            A = treat[i]
            B = control[i]
            if np.isnan(A+B): continue
            if A*B >= 0 :  #同为A/B则视为0
                array[i] = 0
            elif min(abs(treat[i-1:i+2]))<0.01 or min(abs(control[i-1:i+2]))<0.01:
                array[i] = 0
            else:
                array[i] = A - B  #if >0 compartmentBtoA; <0 AtoB

        return super().makeDF(array,"PC1change")

    def getCSV(self):
        super().makeCSV(self.getPC1change())

class CorrelationDifference(BasePara):
    def __init__(self,path,control_path,resolution,chromosome,out_name="",useNA=True,method="pearson"):
        super().__init__(path,resolution,chromosome,out_name,useNA)
        self.control_path = control_path
        self.method = method
        self.control_matrix = loadDenseMatrix(control_path).values

    def getCorrD(self):
        array = self.blankarray
        for i in range(self.matrix_shape):
            treat = pd.Series(np.log1p(self.matrix[i,:]))
            control = pd.Series(np.log1p(self.control_matrix[i,:]))

            if np.nansum(treat) == 0 or np.nansum(control) == 0:
                continue
            elif np.median(np.nan_to_num(treat)) == 0 or np.median(np.nan_to_num(control)) == 0:
                continue
            else:
                array[i] = treat.corr(control,method=self.method)
        return super().makeDF(array,"CorrelationDifference")

    def getCSV(self):
        super().makeCSV(self.getCorrD())

class intraScoreChange(BasePara):
    def __init__(self,path,control_path,resolution,chromosome,out_name="",useNA=True,IS_size=300000):
        super().__init__(path,resolution,chromosome,out_name,useNA)
        self.control_path = control_path

        self.treat = intraTADscore(path,resolution,chromosome).getIntraS(IS_size=IS_size,TADpath=control_path).IntraTADscore
        self.control = intraTADscore(control_path,resolution,chromosome).getIntraS(IS_size=IS_size,TADpath=control_path).IntraTADscore

    def getIntraSC(self):
        array = self.blankarray
        for i in range(self.matrix_shape):
            t = np.array(self.treat)[i]
            c = np.array(self.control)[i]
            if np.isnan(t+c) or min(t,c) == 0:
                continue

            array[i] = np.log2(t/c)
        return super().makeDF(array,"IntraTADscore Change")

    def getCSV(self):
        super().makeCSV(self.getIntraSC())

class interScoreChange(BasePara):
    def __init__(self,path,control_path,resolution,chromosome,out_name="",useNA=True,IS_size=300000):
        super().__init__(path,resolution,chromosome,out_name,useNA)
        self.control_path = control_path

        self.treat = interTADscore(path,resolution,chromosome).getInterS(IS_size=IS_size,TADpath=control_path).InterTADscore
        self.control = interTADscore(control_path,resolution,chromosome).getInterS(IS_size=IS_size,TADpath=control_path).InterTADscore

    def getInterSC(self):
        array = self.blankarray
        for i in range(self.matrix_shape):
            t = np.array(self.treat)[i]
            c = np.array(self.control)[i]
            if np.isnan(t+c) or min(t,c) == 0:
                continue

            array[i] = np.log2(t/c)
        return super().makeDF(array,"InterTADscore Change")

    def getCSV(self):
        super().makeCSV(self.getInterSC())

class InteractionFrequencyChange:
    def __init__(self,path,control_path,resolution,chromosome,gt=None,datatype="rawhic",normIF=True,out_name="noName"):
        treat = InteractionFrequency(path,resolution,chromosome,gt=gt,datatype=datatype,normIF=normIF).getIF()
        print(control_path)
        control = InteractionFrequency(control_path,resolution,chromosome,gt=gt,datatype=datatype,normIF=normIF).getIF()
        scoreChange= pd.DataFrame({"chr":treat.iloc[:,0],"start":treat.iloc[:,1],
                    "end":treat.iloc[:,2],"InteractionFrequencyChange":treat.iloc[:,3]-control.iloc[:,3]})
        self.score = scoreChange

    def getIFC(self):
        return(self.score)
