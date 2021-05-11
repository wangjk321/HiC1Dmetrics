import argparse
from .MultiTypeScore import *
from .plotMetrics import *
from .plotTwoSample import *
from .MultiSampleScore import *
from .callDirectionalTAD import *
from .calldTADAllchr import *
from .discrete import *

def CLI():
    parser = argparse.ArgumentParser(description="HiC1Dmetrics is python3-based tools to \
                                                calculate, visualize, and analyze 1D metrics for Hi-C data \n \
                                                (https://github.com/wangjk321/HiC1Dmetrics) ")
    #parser.set_defaults(func=lambda args: parser.print_help())
    subparsers = parser.add_subparsers(help="Choose the mode to use sub-commands")

    #Function 1
    def func_basic(args):
        args.matrix = args.data
        if args.mode == "plot":
            if args.datatype == "rawhic":
                path = hic2matrix(args.matrix,args.resolution,args.chromosome,args.gt)
                if args.controlmatrix: controlpath = hic2matrix(args.controlmatrix,args.resolution,args.chromosome,args.gt)
            else:
                path = args.matrix
                controlpath = args.controlmatrix

            if args.plottype == "tri":
                if not controlpath:
                    PlotTri(path,args.resolution,args.chromosome,startSite=args.start,endSite=args.end).draw()
                elif controlpath:
                    DiffDraw(path,controlpath,args.resolution,args.chromosome,startSite=args.start,endSite=args.end).draw_tri()
            elif args.plottype == "square":
                if not controlpath:
                    PlotSquare(path,args.resolution,args.chromosome,startSite=args.start,endSite=args.end).draw()
                elif controlpath:
                    DiffDraw(path,controlpath,args.resolution,args.chromosome,startSite=args.start,endSite=args.end).draw_square()
            elif args.plottype == "tad":
                if not controlpath:
                    PlotTAD(path,args.resolution,args.chromosome,startSite=args.start,endSite=args.end).drawTAD(squareSize=int(args.parameter))
                elif controlpath:
                    DiffDraw(path,controlpath,args.resolution,args.chromosome,startSite=args.start,endSite=args.end).drawTAD(squareSize=int(args.parameter))
            plt.savefig(args.outname+".pdf")
        elif args.mode == "dump":
            if args.chromosome == "all":
                allJuicer(args.data,args.normalize,args.resolution,args.gt,args.outname,
                        maxchr=args.maxchr,num=args.nProcesser)
                exit(0)

            if args.datatype not in ["rawhic"] or not args.gt:
                print("Error: dump requires rawhic file and genome_table file"); exit(1)
            codepath = os.path.dirname(os.path.realpath(__file__))
            makeIntra = codepath+"/extract/makeMatrixIntra.sh"
            juicer = codepath+"/jc/jctool_1.11.04.jar"
            foldername = args.outname
            os.system("sh "+makeIntra+" "+args.normalize+" "+"."+" "+args.matrix+" "+
                    str(args.resolution)+" "+args.gt+" "+juicer+" "+args.chromosome+" "+foldername)
        elif args.mode == "gd":
            codepath = os.path.dirname(os.path.realpath(__file__))
            gdcode = codepath+"/gd/makeDensity.sh"
            os.system("sh "+gdcode+" -r "+str(args.resolution)+" -g "+args.data+" -t "+args.chromosome+" -o "+args.outname)
        else:
            print("Unsupported mode"); exit(1)
    parser_basic = subparsers.add_parser("basic",help="Provide basic functions to visualize and handle Hi-C data.")
    parser_basic.add_argument('mode', type=str, help='Type of 1D metrics,,should be one of {dump,plot}')
    parser_basic.add_argument('data', type=str, help='Path of matrix file from JuicerResult')
    parser_basic.add_argument('resolution', type=int,help="Resolution of input matrix")
    parser_basic.add_argument("chromosome",type=str,help="Chromosome number.")
    parser_basic.add_argument("-o","--outname",help="output name",type=str,default="unname")
    parser_basic.add_argument('-c','--controlmatrix', type=str, help='Path of control matrix file from JuicerResult',default=None)
    parser_basic.add_argument('--datatype',type=str,help="matrix or rawhic",default="matrix")
    parser_basic.add_argument('--gt',type=str,help="genome table",default="")
    parser_basic.add_argument('--plottype',type=str,help="Type of plot, could be one of {tri,square,tad}",default="tri")
    parser_basic.add_argument('-s','--start',type=int,help="Start sites for plotting",default=0)
    parser_basic.add_argument('-e','--end',type=int,help="End sites for plotting",default=0)
    parser_basic.add_argument("-p","--parameter",type=str,help="Parameter for indicated metrics",default=None)
    parser_basic.add_argument("--normalize",type=str,help="Normalize methods {NONE/VC/VC_SQRT/KR}",default="KR")
    parser_basic.add_argument("-n","--nProcesser",type=int,help="Number of processors",default=10)
    parser_basic.add_argument('--maxchr',type=int,help="Maximum index of chromosome (human genome is 22,i.e.)",default=None)
    parser_basic.set_defaults(func=func_basic)

    #Function 2
    #=============================================================================
    def func_one(args):
        args.matrix = args.data
        if args.chromosome == "all":
            if not args.maxchr: print("Please sepcify the maximum chromosome"); exit(1)
            if not os.path.exists(args.data):
                print("path not exist"); exit(1)
            scoreAll =oneScoreAllchr(args.data,args.resolution,args.type,args.parameter,
                                    maxchr=args.maxchr,prefix=args.prefix,num=args.nProcesser)
            scoreAll.to_csv(args.outname+"_"+args.type+"_allchr.csv",sep="\t",index=False)
            exit(0)

        if args.parameter and args.type not in ["PC1","IF"]:
            args.parameter = int(args.parameter)

        score = multiScore(args.matrix,args.resolution,args.chromosome,
                            ).obtainOneScore(args.type,parameter=args.parameter,datatype=args.datatype,gt=args.gt)
        score.to_csv(args.outname + ".bedGraph", sep="\t", header=False, index=False)

        if args.draw:
            if args.chromosome=="all": print("Error: not supported"); exit(1)
            if args.type == "IF":
                args.parameter = args.outname + ".bedGraph"
            print("==========output figure==========")
            PlotBedGraph(args.matrix,args.resolution,args.chromosome,startSite=args.start,
                        endSite=args.end,datatype=args.datatype,
                        gt=args.gt).draw(args.type,UniqueParameter=args.parameter)
            plt.savefig(args.outname+".pdf")


    parser_one = subparsers.add_parser("one",help="1D metrics designed for one Hi-C sample.",
                                            description="1D metrics designed for one Hi-C sample.")
    parser_one.add_argument('type', type=str, help='Type of 1D metrics,,should be one of {IS,CI,DI,SS,DLR,PC1,IES,IAS,IF}.')
    parser_one.add_argument('data', type=str, help='Path of matrix or rawhic file.')
    parser_one.add_argument('resolution', type=int,help="Resolution of input matrix.")
    parser_one.add_argument("chromosome",type=str,help="Chromosome number.")
    parser_one.add_argument("-p","--parameter",type=str,help="Parameter for indicated metrics.",default=None)
    parser_one.add_argument("-o","--outname",help="output name (default: 'metrics').",type=str,default="unname")
    parser_one.add_argument("-d","--draw",action='store_true',help="Plot figure for candidate region.",default=False)
    parser_one.add_argument("-s",'--start',type=int,help="Start sites for plotting.",default=0)
    parser_one.add_argument("-e",'--end',type=int,help="End sites for plotting.",default=0)
    parser_one.add_argument('--datatype',type=str,help="Type of input data: matrix(default) or rawhic.",default="matrix")
    parser_one.add_argument('--gt',type=str,help="genome_table file.",default="")
    parser_one.add_argument('--prefix',type=str,help="${prefix}chr1.matrix.gz",default="observed.KR.")
    parser_one.add_argument('--maxchr',type=int,help="Maximum index of chromosome (human genome is 22,i.e.)",default=None)
    parser_one.add_argument("-n","--nProcesser",type=int,help="Number of processors",default=10)
    parser_one.set_defaults(func=func_one)

    #Function 3
    #=============================================================================
    def func_two(args):
        args.matrix = args.data
        args.controlmatrix = args.controldata

        if args.chromosome == "all":
            if not args.maxchr: print("Please sepcify the maximum chromosome"); exit(1)
            if not os.path.exists(args.data):
                print("path not exist"); exit(1)
            scoreAll =twoScoreAllchr(args.data,args.controldata,args.resolution,args.type,args.parameter,
                                    maxchr=args.maxchr,prefix=args.prefix)
            scoreAll.to_csv(args.outname+"_"+args.type+"_allchr.csv",sep="\t",index=False)
            exit(0)

        ms = multiScore(args.matrix,args.resolution,args.chromosome,control_path=args.controlmatrix)
        score,path,control_path = ms.obtainTwoScore(args.type,parameter=args.parameter,datatype=args.datatype,gt=args.gt)
        score.to_csv(args.outname + ".bedGraph", sep="\t", header=False, index=False)

        if args.draw:
            if args.type != "IFC":
                DiffDraw(path,control_path,args.resolution,chr = args.chromosome,startSite=args.start,
                        endSite=args.end,datatype="matrix",gt=args.gt).drawMetric("custom",customfile=args.outname + ".bedGraph",name=args.type)
            elif args.type == "IFC":
                DiffDraw(path,control_path,args.resolution,chr = args.chromosome,startSite=args.start,
                        endSite=args.end,datatype="rawhic",gt=args.gt).drawMetric("custom",customfile=args.outname + ".bedGraph",name=args.type)
            plt.savefig(args.outname+".pdf")

        os.system("rm -rf MatrixTemp*")

    parser_two = subparsers.add_parser("two",help="1D metrics designed for comparison of two Hi-C samples",
                                            description="1D metrics designed for comparison of two Hi-C samples")
    parser_two.add_argument('type', type=str, help='Type of 1D metrics for two-sample comparison,should be one of {ISC,CIC,SSC,deltaDLR,CD,IESC,IASC,IFC,DRF}')
    parser_two.add_argument('data', type=str, help='Path of treated file (matrix or rawhic).')
    parser_two.add_argument('controldata', type=str, help='Path of control file (matrix or rawhic).')
    parser_two.add_argument('resolution', type=int,help="Resolution of input matrix")
    parser_two.add_argument("chromosome",type=str,help="Chromosome number ('chr21',i.e).")
    parser_two.add_argument("-p","--parameter",type=str,help="Parameter for indicated metrics",default=None)
    parser_two.add_argument("-o","--outname",help="output name (default: metricsChange)",type=str,default="metricsChange")
    parser_two.add_argument("-d","--draw",action='store_true',help="Plot figure for candidate region",default=False)
    parser_two.add_argument('-s','--start',type=int,help="Start sites for plotting",default=0)
    parser_two.add_argument('-e','--end',type=int,help="End sites for plotting",default=0)
    parser_two.add_argument('--datatype',type=str,help="matrix or rawhic",default="matrix")
    parser_two.add_argument('--gt',type=str,help="genome table file",default="")
    parser_two.add_argument('--prefix',type=str,help="${prefix}chr1.matrix.gz",default="observed.KR.")
    parser_two.add_argument('--maxchr',type=int,help="Maximum index of chromosome (human genome is 22,i.e.)",default=None)
    parser_two.add_argument("-n","--nProcesser",type=int,help="Number of processors",default=10)
    parser_two.set_defaults(func=func_two)

    #Function 4
    #=============================================================================
    def func_types(args):
        args.matrix = args.data
        typelist = args.typelist.split(",")
        parameterlist = args.parameter.split(",")
        if not args.controlmatrix:
            if not set(typelist).issubset(["IS","CI","DI","SS","DLR","PC1","IES","IAS","IF"]):
                print("Error: not supported"); exit(1)
            if "IF" in typelist and args.datatype == "matrix": print("Error: IF required rawhic datatype"); exit(1)
            ms = multiScore(args.matrix,args.resolution,args.chromosome)
            if not args.draw:
                score = ms.allOneScore(typelist,parameterlist,datatype=args.datatype,gt=args.gt)
            elif args.draw:
                score = ms.plotOneScore(typelist,parameterlist,datatype=args.datatype,gt=args.gt,start=args.start,end=args.end)
                plt.savefig(args.outname+".pdf")
            print(score.iloc[550:650,:])
            score.to_csv(args.outname + ".csv", sep="\t", header=True, index=False)
        elif args.controlmatrix:
            if not set(typelist).issubset(["ISC","CIC","SSC","deltaDLR","CD","IESC","IASC","IFC","DRF"]):
                print("Error: not supported"); exit(1)
            if "IFC" in typelist and args.datatype == "matrix": print("Error: IFC required rawhic datatype"); exit(1)
            if "DRF" in typelist:
                DRFpos = typelist.index("DRF")
                parameterlist[DRFpos] = parameterlist[DRFpos].split("-")
                print(parameterlist)
            ms = multiScore(args.matrix,args.resolution,args.chromosome,control_path=args.controlmatrix)
            if not args.draw:
                score = ms.allTwoScore(typelist,parameterlist,datatype=args.datatype,gt=args.gt)
            elif args.draw:
                score = ms.plotTwoScore(typelist,parameterlist,datatype=args.datatype,gt=args.gt,start=args.start,end=args.end)
                plt.savefig(args.outname+".pdf")
            print(score.iloc[550:650,:])
            score.to_csv(args.outname + ".csv", sep="\t", header=True, index=False)

        os.system("rm -rf MatrixTemp*")

    parser_types = subparsers.add_parser("multitypes",help="Various types of 1D metrics for the same sample",
                                            description="Various types of 1D metrics for the same sample")
    parser_types.add_argument('typelist', type=str, help='Type of 1D metrics,should be {IS,CI,DI,SS,DLR,PC1,IES,IAS,IF} or {ISC,CIC,SSC,deltaDLR,CD,IESC,IASC,IFC,DRF}')
    parser_types.add_argument('data', type=str, help='Path of matrix file or raw .hic file')
    parser_types.add_argument('resolution', type=int,help="Resolution of input matrix")
    parser_types.add_argument("chromosome",type=str,help="Chromosome number.")
    parser_types.add_argument("-p","--parameter",type=str,help="Parameter for indicated metrics",default=None,required=True)
    parser_types.add_argument('-c','--controlmatrix', type=str, help='Path of control matrix file from JuicerResult',default=None)
    parser_types.add_argument("-o","--outname",help="output name (default metrics)",type=str,default="multitypes_metrics")
    parser_types.add_argument('--datatype',type=str,help="matrix or rawhic",default="matrix")
    parser_types.add_argument('--gt',type=str,help="genome table",default="")
    parser_types.add_argument("-d","--draw",action='store_true',help="Plot figure for candidate region",default=False)
    parser_types.add_argument('-s','--start',type=int,help="Start sites for plotting",default=0)
    parser_types.add_argument('-e','--end',type=int,help="End sites for plotting",default=0)
    parser_types.set_defaults(func=func_types)

    #Function 5
    #=============================================================================
    def func_samples(args):
        if args.type == "IF" and args.datatype == "matrix": print("Error: IF required rawhic datatype"); exit(1)
        datafile = pd.read_csv(args.data,sep="\t",header=None)
        labels = list(datafile.iloc[:,0])
        samplelist = list(datafile.iloc[:,1])
        if not args.corr and not args.heat and not args.line and not args.discrete:
            score = getMultiSamplesScore(samplelist,labels,args.resolution,args.chromosome,args.type,args.parameter,
                                        datatype=args.datatype,gt=args.gt)
        elif args.corr:
            ms = repQC(samplelist,labels,args.resolution,args.chromosome,args.type,args.parameter,datatype=args.datatype,gt=args.gt)
            score = ms.score
            ms.corr_plot()
            plt.savefig(args.outname+"_corr.pdf")
        elif args.discrete:
            if args.datatype == "rawhic":
                print("not supported");exit(1)
            msd = multiSampleDiscrete(samplelist,labels,args.resolution,args.chromosome,args.type,args.parameter)
            plotpath= samplelist[0]
            msd.plotMultiDiscrete(plotpath,args.start,args.end)
            plt.savefig(args.outname+"_discrete.pdf")
            exit(0)
        elif args.line or args.heat:
            ms = repQC(samplelist,labels,args.resolution,args.chromosome,args.type,args.parameter,
                        datatype=args.datatype,gt=args.gt)
            score = ms.score
            if args.datatype == "matrix":
                plotpath= samplelist[0]
            elif args.datatype == "rawhic":
                plotpath = hic2matrix(samplelist[0],args.resolution,args.chromosome,args.gt)

            if args.heat: plottype = "heat"
            elif args.line: plottype = "line"

            ms.heatmap_tri(plotpath,args.start,args.end,clmax=args.clmax,heatmin=None,plottype=plottype)
            plt.savefig(args.outname+"_"+plottype+".pdf")

        score.to_csv(args.outname + ".csv", sep="\t", header=True, index=False)

        print(score.iloc[550:650,:])
        os.system("rm -rf MatrixTemp*")

    parser_samples = subparsers.add_parser("multisamples",help="The same metrics for muliple samples",
                                            description="The same metrics for muliple samples")
    parser_samples.add_argument('type', type=str, help='Type of 1D metrics,,should be one of {IS,CI,DI,SS,DLR,PC1,IES,IAS,IF} or {ISC,CIC,SSC,deltaDLR,CD,IESC,IASC,IFC,DRF}')
    parser_samples.add_argument('data',type=str,help="a txt contain paths for all samples")
    parser_samples.add_argument('resolution', type=int,help="Resolution of input matrix")
    parser_samples.add_argument("chromosome",type=str,help="Chromosome number.")
    parser_samples.add_argument('--datatype',type=str,help="matrix or rawhic",default="matrix")
    parser_samples.add_argument('--samplelist', type=str, help='list of file path, can be rawhic or matrix')
    parser_samples.add_argument('--labels', type=str, help='list of file name')
    parser_samples.add_argument("-p","--parameter",type=str,help="Parameter for indicated metrics",default=None)
    parser_samples.add_argument("-o","--outname",help="output name (default metrics)",type=str,default="multisamples_metrics")
    parser_samples.add_argument('--gt',type=str,help="genome table",default="")
    parser_samples.add_argument("--corr",action='store_true',help="Plot correlation for all samples",default=False)
    parser_samples.add_argument("--heat",action='store_true',help="Plot raw heatmap for all samples",default=False)
    parser_samples.add_argument("--line",action='store_true',help="Plot line chart for all samples",default=False)
    parser_samples.add_argument("--discrete",action='store_true',help="Plot discrete heatmap for all samples",default=False)
    parser_samples.add_argument('-s','--start',type=int,help="Start sites for plotting",default=0)
    parser_samples.add_argument('-e','--end',type=int,help="End sites for plotting",default=0)
    parser_samples.add_argument('--clmax',type=int,help="End sites for plotting",default=None)

    parser_samples.set_defaults(func=func_samples)

    #Function 6
    #=============================================================================
    def func_call(args):
        args.matrix = args.data
        if args.datatype == "rawhic" and args.mode != "hubs":
            path = hic2matrix(args.matrix,args.resolution,args.chromosome,args.gt)
            if args.controlmatrix: controlpath = hic2matrix(args.controlmatrix,args.resolution,args.chromosome,args.gt)
        else:
            path = args.matrix
            controlpath = args.controlmatrix

        if args.mode == "TAD":
            if not args.parameter: args.parameter = 300000
            tad = TADcallIS(path,args.resolution,args.chromosome,squareSize=int(args.parameter))
            tad.to_csv(args.outname + "_TAD.csv", sep="\t", header=True, index=False)
        elif args.mode == "dTAD":
            if args.parameter:
                parameter = args.parameter.split("-")
            else:
                parameter = [200000,5000000]
            if not args.controlmatrix: print("Error: DRF requires control data"); exit(1)
            dt = DirectionalTAD(path,controlpath,args.resolution,args.chromosome,
                                startDRF=int(parameter[0]), sizeDRF=int(parameter[1]))
            leftdTAD,rightdTAD,_ = dt.extractRegion()
            leftdTAD.to_csv(args.outname + "_leftdTAD.csv", sep="\t", header=True, index=False)
            rightdTAD.to_csv(args.outname + "_rightdTAD.csv", sep="\t", header=True, index=False)
        elif args.mode == "stripe":
            st = stripeTAD(path,args.resolution,args.chromosome)
            stripe = st.callStripe(squareSize=300000)
            stripe.to_csv(args.outname + "_stripe.csv", sep="\t", header=True, index=False)
        elif args.mode == "hubs":
            if args.datatype != "rawhic": print("Error: hubs requires rawhic datatype"); exit(1)
            IF = InteractionFrequency(path,args.resolution,args.chromosome,gt=args.gt).getIF()
            thresh = np.percentile(IF.iloc[:,3],90)
            hubregion = IF[IF.iloc[:,3]>thresh]
            hubregion.to_csv(args.outname + "_hubs_IF.csv", sep="\t", header=True, index=False)
            os.system("sed '1d' " + args.outname + "_hubs_IF.csv" + " |bedtools merge -i stdin > "+ args.outname + "_hubs.csv")
        else:
            print("unsupported model");exit(1)

    parser_call = subparsers.add_parser("call",help="Extract secondary information from metrics (dTAD, stripeTAD, et.al)",
                                            description="Extract secondary information from metrics (dTAD, stripeTAD, et.al)")
    parser_call.add_argument('mode', type=str, help='Running mode,,should be one of {dTAD,stripe,TAD,hubs}')
    parser_call.add_argument('data', type=str, help='Path of matrix file or raw .hic file')
    parser_call.add_argument('resolution', type=int,help="Resolution of input matrix")
    parser_call.add_argument("chromosome",type=str,help="Chromosome number.")
    parser_call.add_argument("-o","--outname",help="output name",type=str,default="defaultname")
    parser_call.add_argument('-c','--controlmatrix', type=str, help='Path of control matrix file from JuicerResult',default=None)
    parser_call.add_argument('--datatype',type=str,help="matrix or rawhic",default="matrix")
    parser_call.add_argument('--gt',type=str,help="genome table",default="")
    parser_call.add_argument("-p","--parameter",type=str,help="Parameter for indicated mode",default=None)
    parser_call.set_defaults(func=func_call)

    parser.add_argument("-V","--version",help="Show h1d version",action='store_true',default=False)
    args = parser.parse_args()
    if args.version:
        print("h1d version 0.1.8")
        exit(0)
    try:
        func = args.func
    except AttributeError:
        parser.error("too few arguments")
    func(args)
    os.system("rm -rf MatrixTemp*")
    os.system("rm -rf info.txt")

if __name__ == '__main__':
    CLI()
