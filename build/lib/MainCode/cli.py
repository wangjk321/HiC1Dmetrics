import argparse

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
        score = MultiTypeScore.multiScore(args.matrix,args.resolution,args.chromosome).obtainOneScore(args.type,parameter=args.parameter)
        score.to_csv(args.outname + ".bedGraph", sep="\t", header=False, index=False)
    parser_one = subparsers.add_parser("one",help="1D metrics designed for one Hi-C sample",
                                            description="1D metrics designed for one Hi-C sample")
    parser_one.add_argument('type', type=str, help='Type of 1D metrics,,should be one of {IS,CI,DI,SS,DLR,PC1,IES,IAS,IF}')
    parser_one.add_argument('matrix', type=str, help='Path of matrix file from JuicerResult')
    parser_one.add_argument('resolution', type=int,help="Resolution of input matrix")
    parser_one.add_argument("chromosome",type=str,help="Chromosome number.")
    parser_one.add_argument("-p","--parameter",type=int,help="Parameter for indicated metrics",default=None)
    parser.add_argument("-o","--outname",help="output name (default metrics)",type=str,default="metrics")

    parser_one.set_defaults(func=func_one)


























    #Function 3
    #=============================================================================
    parser_two = subparsers.add_parser("two",help="1D metrics designed for comparison of two Hi-C samples")
    #Function 4
    #=============================================================================
    parser_call = subparsers.add_parser("call",help="Extract secondary information from metrics, (dTAD, stripe-TAD, A/B compartment et.al.)")

    parser_types = subparsers.add_parser("multitype",help="Analyze various metrics for the same sample/comparison")

    parser_samples = subparsers.add_parser("multisample",help="Analyze the same metric for various Hi-C samples.")



    args = parser.parse_args()
    try: func = args.func
    except AttributeError:
        parser.error("too few arguments")
    func(args)

#main()
