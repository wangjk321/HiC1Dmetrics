from multiprocessing import Pool
from .callDirectionalTAD import *
from .MultiTypeScore import *
import os

# calculate each dTAD of chromosome separately
class paralfunc(object):
    def __init__(self,myfun,num_processer):
        self.myfun = myfun
        self.num_processer = num_processer

        #chrlist= list(range(1,23)); chrlist.append("X")
        #chrolist = ["chr"+str(a) for a in chrlist]
        #self.chrolist = chrolist

    def run(self,pathName,pathControl,resolution,maxchr=22):
        chrlist= list(range(1,maxchr+1)); chrlist.append("X")
        chrolist = ["chr"+str(a) for a in chrlist]
        self.chrolist = chrolist

        resultlist=[]
        p = Pool(self.num_processer)
        for r in self.chrolist:
            result = p.apply_async(self.myfun, args=(r,pathName,pathControl,resolution,))
            resultlist.append(result)
        p.close()
        p.join()

        output_left = pd.DataFrame()
        output_right = pd.DataFrame()
        for i in resultlist:
            output_single = i.get().copy()
            output_left = output_left.append(output_single["leftTAD"])
            output_right = output_right.append(output_single["rightTAD"])
        return(output_left,output_right)

def call_dTAD(chrom,pathName,pathControl,resolution):
    filename = pathName+"/observed.KR."+chrom+".matrix.gz"
    controlfilename = pathControl + "/observed.KR."+chrom+".matrix.gz"

    dTAD = DirectionalTAD(filename,controlfilename, \
                      resolution,chr=chrom,startDRF=500000,sizeDRF=2000000,\
                      sizeIS=300000)
    leftTAD,rightTAD,_ = dTAD.extractRegion()
    return({"leftTAD":leftTAD,"rightTAD":rightTAD})

def alldTAD(pathName,pathControl,resolution,outname="twoSample",maxchr=22):
    ltad,rtad = paralfunc(call_dTAD,30).run(pathName,pathControl,resolution,maxchr)
    ltad.to_csv(outname+"_leftTAD.csv",sep="\t",index=False)
    rtad.to_csv(outname+"_rightTAD.csv",sep="\t",index=False)

#alldTAD("/Users/wangjiankang/figureServer/Nov2020/Rad21KD1_HiCmatrix",
#        "/Users/wangjiankang/figureServer/Nov2020/Control1_HiCmatrix",50000,maxchr=23)

#=====================================#
class paralfuncOneSample(object):
    def __init__(self,myfun,num_processer):
        self.myfun = myfun
        self.num_processer = num_processer

    def run(self,pathName,resolution,maxchr=22,type=None,parameter=None,prefix="/observed.KR.",controlpath=None):
        chrlist= list(range(1,maxchr+1)); chrlist.append("X")
        chrolist = ["chr"+str(a) for a in chrlist]
        self.chrolist = chrolist

        resultlist=[]
        p = Pool(self.num_processer)
        for r in self.chrolist:
            result = p.apply_async(self.myfun, args=(r,pathName,resolution,type,parameter,prefix,controlpath))
            resultlist.append(result)
        p.close()
        p.join()

        output_all = pd.DataFrame()
        for i in resultlist:
            output_single = i.get().copy()
            output_all = output_all.append(output_single)
        return(output_all)

def TAD1sample(chrom,pathName,resolution,type,parameter):
    filename = pathName+"/observed.KR."+chrom+".matrix.gz"
    tad = TADcallIS(filename,resolution,chrom,squareSize=300000)
    return(tad)

def runTAD1sample(pathName,resolution,outname="allTAD",maxchr=22):
    tad = paralfuncOneSample(TAD1sample,30).run(pathName,resolution,maxchr)
    tad.to_csv(outname+"_TAD.csv",sep="\t",index=False)

def oneScoreSinglechr(chrom,pathName,resolution,type,parameter,prefix="observed.KR.",controlpath=None):
    filename = pathName+"/"+prefix+chrom+".matrix.gz"
    score = multiScore(filename,resolution,chrom).obtainOneScore(type,parameter)
    return(score)

def oneScoreAllchr(pathName,resolution,type,parameter,outname="OneScore",maxchr=22,prefix="observed.KR.",num=10):
    allscore = paralfuncOneSample(oneScoreSinglechr,num).run(pathName,resolution,maxchr,type,parameter,prefix)
    return(allscore)
    #allscore.to_csv(outname+"_"+type+".csv",sep="\t",index=False)

def multiScoreSinglechr(chrom,pathName,resolution,typelist,parameterlist):
    filename = pathName+"/observed.KR."+chrom+".matrix.gz"
    multiType = multiScore(filename,resolution,chrom).allOneScore(typelist,parameterlist)
    return(multiType)

def multiScoreAllchr(pathName,resolution,typelist,parameterlist,outname="multiScore",maxchr=22):
    allmultiscore = paralfuncOneSample(multiScoreSinglechr,30).run(pathName,resolution,maxchr,typelist,parameterlist)
    return(allmultiscore)

#oneScoreAllchr("/Users/wangjiankang/figureServer/Nov2020/Rad21KD1_HiCmatrix",50000,maxchr=5)


class paralfunJuicer(object):
    def __init__(self,myfun,num_processer):
        self.myfun = myfun
        self.num_processer = num_processer

    def run(self,data,normalize,resolution,gt,outname,maxchr=22):
        chrlist= list(range(1,maxchr+1)); chrlist.append("X")
        chrolist = ["chr"+str(a) for a in chrlist]
        self.chrolist = chrolist

        resultlist=[]
        p = Pool(self.num_processer)
        for r in self.chrolist:
            p.apply_async(self.myfun, args=(r,data,normalize,resolution,gt,outname,))
        p.close()
        p.join()

def oneJuicer(chrom,data,normalize,resolution,gt,outname):
    print("Input data: ", data)
    print("Dumping contact " +chrom+ " matrix from .hic file ......")
    codepath = os.path.dirname(os.path.realpath(__file__))
    makeIntra = codepath+"/extract/makeMatrixIntra.sh"
    juicer = codepath+"/jc/jctool_1.11.04.jar"
    foldername = outname
    os.system("bash "+makeIntra+" "+normalize+" "+"."+" "+data+" "+
            str(resolution)+" "+gt+" "+juicer+" "+chrom+" "+foldername + "> info.txt")
    try: os.system("rm info.txt")
    except: pass
    print("Dump finished, output is in ./"+outname)

def allJuicer(data,normalize,resolution,gt,outname,maxchr=22,num=20):
    paralfunJuicer(oneJuicer,num).run(data,normalize,resolution,gt,outname,maxchr=maxchr)

#========================

class paralfuncTwoSample(object):
    def __init__(self,myfun,num_processer):
        self.myfun = myfun
        self.num_processer = num_processer

    def run(self,pathName,resolution,maxchr=22,type=None,parameter=None,prefix="/observed.KR.",controlpath=None):
        chrlist= list(range(1,maxchr+1)); chrlist.append("X")
        chrolist = ["chr"+str(a) for a in chrlist]
        self.chrolist = chrolist

        resultlist=[]
        p = Pool(self.num_processer)
        for r in self.chrolist:
            result = p.apply_async(self.myfun, args=(r,pathName,resolution,type,parameter,prefix,controlpath))
            resultlist.append(result)
        p.close()
        p.join()

        output_all = pd.DataFrame()
        for i in resultlist:
            output_single = i.get()[0].copy()
            output_all = output_all.append(output_single)
        return(output_all)

def twoScoreSinglechr(chrom,pathName,resolution,type,parameter,prefix="observed.KR.",controlpath=""):
    filename = pathName+"/"+prefix+chrom+".matrix.gz"
    controlfile = controlpath+"/"+prefix+chrom+".matrix.gz"
    score = multiScore(filename,resolution,chrom,control_path=controlfile).obtainTwoScore(type,parameter)
    return(score)

def twoScoreAllchr(pathName,controlpath,resolution,type,parameter,outname="OneScore",maxchr=22,prefix="observed.KR."):
    allscore = paralfuncTwoSample(twoScoreSinglechr,30).run(pathName,resolution,maxchr,type,parameter,prefix,controlpath=controlpath)
    return(allscore)

#allJuicer("/Users/wangjiankang/Documents/localrun/Control_1/inter_30.hic","KR",50000,
#            "./gd/hg38/gd50000.txt","justtest",maxchr=5)
