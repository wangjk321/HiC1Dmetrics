from .calculateMetrics import *
from scipy.signal import argrelextrema

class call_stripe(BasePara):
    def callStripe(self,squareSize=300000,useNA=True,seg=100000, strong_thresh= 0.2):
        iasDF = intraTADscore(self.path,self.resolution,self.chromosome).getIntraS()

        ias = iasDF.iloc[:,3]
        #local maximum
        localMaxIAS = ias.iloc[argrelextrema(np.array(ias), np.greater)]

        # IAS > mean
        localMaxIAS = localMaxIAS[localMaxIAS > ias.mean()]

        resolution = self.resolution
        binNum = int(seg/resolution)
        localMaxAround = []
        aroundZero=[]
        diffrightleft=[]

        for i in localMaxIAS.index:
            #maximum around 100kb
            localMaxAround.append(ias.loc[i-binNum:i+binNum].max())

            #around != NA ??

            #Strong peaks
            minusbin =  ias.loc[i] - ias.loc[i-binNum]
            plusbin = ias.loc[i] - ias.loc[i+binNum]
            diffrightleft.append(minusbin+plusbin)

        bool1 = (np.array(localMaxIAS)-np.array(localMaxAround)) >= 0
        bool2 = np.array(diffrightleft) > strong_thresh
        localMaxIAS = localMaxIAS[bool1*bool2]

        outDF = iasDF.loc[localMaxIAS.index]
        return outDF
