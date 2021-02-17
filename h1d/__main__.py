import argparse
from .MultiTypeScore import *
from .plotMetrics import *
from .plotTwoSample import *

def CLI():
    parser = argparse.ArgumentParser(description="HiC1Dmetrics is a easy to use tools to \
                                                calculate, visualize, and analyze 1D metrics for Hi-C")
    #parser.set_defaults(func=lambda args: parser.print_help())
    subparsers = parser.add_subparsers(help="Choose the mode to use sub-commands")

    #Function 1
    parser_basic = subparsers.add_parser("basic",help="Provide basic functions to load, handle and visualize Hi-C data.")

    #Function 2
    #=============================================================================
    def func_one(args):
        if args.allchr:
            if datatype not in ["folder","rawhic"]: print("Error: not supported"); exit(1)

        if args.parameter and args.type not in ["PC1","IF"]:
            args.parameter = int(args.parameter)

        score = multiScore(args.matrix,args.resolution,args.chromosome,
                            ).obtainOneScore(args.type,parameter=args.parameter,datatype=args.datatype,gt=args.gt)
        score.to_csv(args.outname + ".bedGraph", sep="\t", header=False, index=False)

        if args.draw:
            if args.allchr: print("Error: not supported"); exit(1)
            if args.type == "IF":
                args.parameter = args.outname + ".bedGraph"
            print("==========output figure==========")
            PlotBedGraph(args.matrix,args.resolution,args.chromosome,startSite=args.start,
                        endSite=args.end,datatype=args.datatype,
                        gt=args.gt).draw(args.type,UniqueParameter=args.parameter)
            plt.savefig(args.outname+".pdf")


    parser_one = subparsers.add_parser("one",help="1D metrics designed for one Hi-C sample",
                                            description="1D metrics designed for one Hi-C sample")
    parser_one.add_argument('type', type=str, help='Type of 1D metrics,,should be one of {IS,CI,DI,SS,DLR,PC1,IES,IAS,IF}')
    parser_one.add_argument('matrix', type=str, help='Path of matrix file from JuicerResult')
    parser_one.add_argument('resolution', type=int,help="Resolution of input matrix")
    parser_one.add_argument("chromosome",type=str,help="Chromosome number.")
    parser_one.add_argument("-p","--parameter",type=str,help="Parameter for indicated metrics",default=None)
    parser_one.add_argument("-o","--outname",help="output name (default metrics)",type=str,default="metrics")
    parser_one.add_argument("-d","--draw",action='store_true',help="Plot figure for candidate region",default=False)
    parser_one.add_argument('--start',type=int,help="Start sites for plotting",default=0)
    parser_one.add_argument('--end',type=int,help="End sites for plotting",default=0)
    parser_one.add_argument('--datatype',type=str,help="matrix or rawhic",default="matrix")
    parser_one.add_argument('--gt',type=str,help="genome table",default="")
    parser_one.add_argument('--allchr',action='store_true',help="Calculate metrics for multiple chromosomes",default=False)
    parser_one.set_defaults(func=func_one)

    #Function 2
    #=============================================================================
    def func_two(args):
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
    parser_two.add_argument('type', type=str, help='Type of 1D metrics,,should be one of {ISC,CIC,SSC,deltaDLR,CD,IESC,IASC,IFC,DRF}')
    parser_two.add_argument('matrix', type=str, help='Path of matrix file from JuicerResult')
    parser_two.add_argument('controlmatrix', type=str, help='Path of control matrix file from JuicerResult')
    parser_two.add_argument('resolution', type=int,help="Resolution of input matrix")
    parser_two.add_argument("chromosome",type=str,help="Chromosome number.")
    parser_two.add_argument("-p","--parameter",type=str,help="Parameter for indicated metrics",default=None)
    parser_two.add_argument("-o","--outname",help="output name (default metrics)",type=str,default="metricsChange")
    parser_two.add_argument("-d","--draw",action='store_true',help="Plot figure for candidate region",default=False)
    parser_two.add_argument('-s','--start',type=int,help="Start sites for plotting",default=0)
    parser_two.add_argument('-e','--end',type=int,help="End sites for plotting",default=0)
    parser_two.add_argument('--datatype',type=str,help="matrix or rawhic",default="matrix")
    parser_two.add_argument('--gt',type=str,help="genome table",default="")
    parser_two.set_defaults(func=func_two)

    #Function 3
    #=============================================================================
    def func_types(args):
        typelist = args.type.split(",")

    parser_types = subparsers.add_parser("multitypes",help="Various types of 1D metrics for the same sample",
                                            description="Various types of 1D metrics for the same sample")
    parser_types.add_argument('type', type=str, help='Type of 1D metrics,,should be one of {ISC,CIC,SSC,deltaDLR,CD,IESC,IASC,IFC,DRF}')
    parser_types.add_argument('matrix', type=str, help='Path of matrix file from JuicerResult')
    parser_types.add_argument('controlmatrix', type=str, help='Path of control matrix file from JuicerResult')
    parser_types.add_argument('resolution', type=int,help="Resolution of input matrix")
    parser_types.add_argument("chromosome",type=str,help="Chromosome number.")
    parser_types.set_defaults(func=func_types)

    args = parser.parse_args()
    try: func = args.func
    except AttributeError:
        parser.error("too few arguments")
    func(args)

if __name__ == '__main__':
    CLI()
