#! /usr/bin/env python
# -*- coding: utf-8 -*-
from plotMetrics import *
from calculateMetrics import *

if __name__ == '__main__':
#def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode",help="Type of hic metrics,should be one of IS,CI,DI,DFR",type=str)
    parser.add_argument("matrix",help="Path of matrix file from JuicerResult",type=str)
    parser.add_argument("resolution",help="Resolution of input matrix",type=int)
    parser.add_argument("chromosome",help="chromosome No.",type=str)

    parser.add_argument("--square_size",help="IS parameter, Square size of Insulation Score (default 150000)",type=int,default=150000)
    parser.add_argument("--CI_size",help="CI paramater, Square size of Contrast Index (default 150000)",type=int,default=150000)
    parser.add_argument("--DI_distance",help="DI parameter, Distance used to calculate DI (default 1000000)",type=int,default=1000000)
    parser.add_argument("--control_matrix",help="DRF paramater, control matrix",type=str)
    parser.add_argument("--start_distance",help="DRF paramater, start_distance",type=int,default=0)
    parser.add_argument("--end_distance",help="DRF parameter, end_distance",type=int,default=2000000)

    parser.add_argument("-O","--out_name",help="-O or out_name (default noName)",type=str,default="noName")
    parser.add_argument("--plot_from",help="Start site of plot",type=int,default=0)
    parser.add_argument("--plot_to",help="End site of plot",type=int,default=0)

    args = parser.parse_args()

    plot = PlotBedGraph(args.matrix,args.resolution,clmin=0,clmax=20,\
                    title=args.out_name,chr=args.chromosome, \
                    startSite= args.plot_from,endSite=args.plot_to)

    if args.mode == "IS":
        print("--------Insulation Score Mode--------- \n")
        print("--------Reading Files--------- \n")
        IS = InsulationScore(args.matrix,args.resolution,args.chromosome,args.out_name,args.square_size)
        print(IS)
        IS.getCSV()
        print("---------Saving pdf----------")
        plot.makePDF("IS",args.out_name)
    elif args.mode == "CI":
        print("--------Contrast Index Mode--------- \n")
        print("--------Reading Files--------- \n")
        CI = ContrastIndex(args.matrix,args.resolution,args.chromosome,args.out_name,args.CI_size)
        print(CI)
        CI.getCSV()
        print("---------Saving pdf----------")
        plot.makePDF("CI",args.out_name)
    elif args.mode == "DI":
        print("--------Directionality Index Mode--------- \n")
        print("--------Reading Files--------- \n")
        DI = DirectionalityIndex(args.matrix,args.resolution,args.chromosome,args.out_name,args.DI_distance)
        print(DI)
        DI.getCSV()
        print("---------Saving pdf----------")
        plot.makePDF("DI",args.out_name)
    elif args.mode == "DRF":
        print("--------DirectionalRelativeFreq Mode--------- \n")
        print("--------Reading Files--------- \n")
        DRF = DirectionalRelativeFreq(args.matrix,args.control_matrix,args.resolution,args.chromosome, \
                                        args.out_name,args.start_distance,args.end_distance)
        print(DRF)
        DRF.getCSV()
    else:
        print("Please specific the metrics type")
        exit(1)

    print("\n--------Finish--------- \n")
