from calculateMetrics import *
from calculateTwoSample import *

class multiScore:
    def __init__(self,path,res,chr,control_path=""):
        self.path = path
        self.res = res
        self.chr = chr
        self.control_path = control_path

    def obtainOneScore(self,mode,parameter):
        if mode == "IS":
            score = InsulationScore(self.path,self.res,self.chr,square_size=parameter).getIS()
        elif mode == "CI":
            score = ContrastIndex(self.path,self.res,self.chr,CI_size=parameter).getCI()
        elif mode == "DI":
            score = DirectionalityIndex(self.path,self.res,self.chr,DI_distance=parameter).getDI()
        elif mode == "TADss":
            score = SeparationScore(self.path,self.res,self.chr,TADss_size=parameter).getTADss()
        elif mode == "DLR":
            score = DistalToLocalRatio(self.path,self.res,self.chr,sizeDLR=parameter).getDLR()
        elif mode == "intraS":
            score = intraTADscore(self.path,self.res,self.chr).getIntraS(IS_size = parameter)
        elif mode == "interS":
            score = interTADscore(self.path,self.res,self.chr).getInterS(IS_size = parameter)
        elif mode == "PC1":
            score = CompartmentPC1(self.path,self.res,self.chr).getPC1(signCorr = parameter)
        return(score)

    def allOneScore(self,typelist=["IS","CI","DI","TADss","DLR","intraS","interS","PC1"],
                    parameterlist=[300000,300000,1000000,300000,3000000,300000,300000,"NotSpecified"]):
        for i,type in enumerate(typelist):
            if i == 0:
                multiType = self.obtainOneScore(mode=typelist[i],parameter=parameterlist[i])
            else:
                next = self.obtainOneScore(mode=typelist[i],parameter=parameterlist[i]).iloc[:,3]
                multiType = pd.concat([multiType,next],axis=1)

        return(multiType)

    def obtainTwoScore(self,mode,parameter):
        if mode == "ISC":
            score = TADScoreChange(self.path,self.control_path,self.res,self.chr).getChange("IS",parameter)
        elif mode == "CIC":
            score = TADScoreChange(self.path,self.control_path,self.res,self.chr).getChange("CI",parameter)
        elif mode == "DIC":
            score = TADScoreChange(self.path,self.control_path,self.res,self.chr).getChange("DI",parameter)
        elif mode == "TADssC":
            score = TADScoreChange(self.path,self.control_path,self.res,self.chr).getChange("TADss",parameter)
        elif mode == "deltaDLR":
            score = deltaDLR(self.path,self.control_path,self.res,self.chr,sizeDLR=parameter).getDeltaDLR()
        elif mode == "intraSC":
            score = 





#def allTwoSample(path,res,chr,type=["ISC","CIC","DIC","TADssC","deltaDLR","intraSC","interSC","DRF","CorrD","PC1C"],
#                parameter=[300000,300000,1000000,300000,3000000,300000,300000,[200000,5000000],"pearson"],PC1parameter="NotSpecified")
