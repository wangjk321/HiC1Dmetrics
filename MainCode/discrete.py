from MultiTypeScore import *

def getDiscrete(path,res,chr,mode,parameter):
    ob = multiScore(path,res,chr)
    score = ob.obtainOneScore(mode,parameter)
    return(score)
    
