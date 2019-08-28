#!/bin/env python
#Intro Version 1.0 
#JaeYoon Kim 
#2017 09 11

import sys, time, math, copy, random, os
import numpy as np
from itertools import *
from scipy.stats import norm
import collections
from  IntroFunc import *
import IntroCalc
from optparse import OptionParser

parser = OptionParser(usage="python Intro.py \n Usage: %prog [options]")
parser.add_option("-m","--mode",action = 'store',type = 'string',dest = 'Mode',help = "")
parser.add_option("-v","--input",action = 'store',type = 'string',dest = 'VCF_FILE',help = "")
parser.add_option("-i","--info",action = 'store',type = 'string',dest = 'INFO_FILE',help = "")
parser.add_option("-o","--popO",action = 'store',type = 'string',dest = 'OutPop_Name',help = "")
parser.add_option("-t","--popT",action = 'store',type = 'string',dest = 'ThirdPop_Name',help = "")
parser.add_option("-p","--pop1",action = 'store',type = 'string',dest = 'TargetPop1_Name',help = "")
parser.add_option("-q","--pop2",action = 'store',type = 'string',dest = 'TargetPop2_Name',help = "")
parser.add_option("-c","--chr",action = 'store',type = 'string',dest = 'CHR',help = "")
parser.add_option("-w","--window",action = 'store',type = 'int',dest = 'WINDOW_SIZE',help = "")
parser.add_option("-s","--slide",action = 'store',type = 'int',dest = 'SLIDING_SIZE',help = "")
(opt, args) = parser.parse_args()

if opt.Mode == "Dstat":
    if opt.Mode == None or opt.VCF_FILE == None or opt.INFO_FILE == None or opt.OutPop_Name == None or opt.ThirdPop_Name == None or opt.TargetPop1_Name == None or opt.TargetPop2_Name == None or opt.CHR == None or opt.WINDOW_SIZE == None or opt.SLIDING_SIZE == None:
        print('python Intro.py -m Dstat -v [VCF] -i [INFO] -o [OUT-POP] -t [THR-POP] -p [POP1] -q [POP2] -c [CHR] -w [WINDOW] -s [SLIDE]')
        print('[-v str | --input str] : Gziped VCF file')
        print('[-i str | --info str] : Info file, it consisted sample and group name columns')
        print('[-o str | --popO str] : Out-Group Population Name, O indicates alphabet O, not zero')
        print('[-t str | --popT str] : Third In-Group Population Name')
        print('[-p str | --pop1 str] : Target Population Name1')
        print('[-q str | --pop2 str] : Target Population Name2')
        print('[-c str | --chr str] : Chromosome number, integer value')
        print('[-w int | --window int] : Window Size, integer value')
        print('[-s int | --slide int] : Sliding size, integer value')
        sys.exit()
    else:
        sys.path.append('./')
        infile, ind, popO, popT, pop1, pop2, CHR, Window, Slide = opt.VCF_FILE, opt.INFO_FILE, opt.OutPop_Name, opt.ThirdPop_Name, opt.TargetPop1_Name, opt.TargetPop2_Name, opt.CHR, opt.WINDOW_SIZE, opt.SLIDING_SIZE
        Result = IntroCalc.CalcD(infile, ind, popO, popT, pop1, pop2, CHR, Window, Slide)
        if Result != "Done": print("Calulation ERROR")
elif opt.Mode == "RND":
    if opt.Mode == None or opt.VCF_FILE == None or opt.INFO_FILE == None or opt.OutPop_Name == None or opt.TargetPop1_Name == None or opt.TargetPop2_Name == None or opt.CHR == None or opt.WINDOW_SIZE == None or opt.SLIDING_SIZE == None:
        print('python Intro.py -m RND -v [VCF] -i [INFO] -o [OUT-POP] -p [POP1] -q [POP2] -c [CHR] -w [WINDOW] -s [SLIDE]')
        print('[-v str | --input str] : Gziped VCF file')
        print('[-i str | --info str] : Info file, it consisted sample and group name columns')
        print('[-o str | --popO str] : OutGroup Population Name, O indicates alphabet O, not zero')
        print('[-p str | --pop1 str] : Target Population Name1')
        print('[-q str | --pop2 str] : Target Population Name2')
        print('[-c str | --chr str] : Chromosome number, integer value')
        print('[-w int | --window int] : Window Size, integer value')
        print('[-s int | --slide int] : Sliding size, integer value')
        sys.exit()
    else:
        sys.path.append('./')
        infile, ind, popO, pop1, pop2, CHR, Window, Slide = opt.VCF_FILE, opt.INFO_FILE, opt.OutPop_Name, opt.TargetPop1_Name, opt.TargetPop2_Name, opt.CHR, opt.WINDOW_SIZE, opt.SLIDING_SIZE
        Result = IntroCalc.CalcRND(infile, ind, popO, pop1, pop2, CHR, Window, Slide)
        if Result != "Done": print("Calulation ERROR")
elif opt.Mode == "Dxy":
    if opt.Mode == None or opt.VCF_FILE == None or opt.INFO_FILE == None or opt.TargetPop1_Name == None or opt.TargetPop2_Name == None or opt.CHR == None or opt.WINDOW_SIZE == None or opt.SLIDING_SIZE == None:
        print('python Intro.py -m Dxy -v [VCF] -i [INFO] -p [POP1] -q [POP2] -c [CHR] -w [WINDOW] -s [SLIDE]')
        print('[-v str | --input str] : Gziped VCF file')
        print('[-i str | --info str] : Info file, it consisted sample and group name columns')
        print('[-p str | --pop1 str] : Target Population Name1')
        print('[-q str | --pop2 str] : Target Population Name2')
        print('[-c str | --chr str] : Chromosome number, integer value')
        print('[-w int | --window int] : Window Size, integer value')
        print('[-s int | --slide int] : Sliding size, integer value')
        sys.exit()
    else:
        sys.path.append('./')
        infile, ind, pop1, pop2, CHR, Window, Slide = opt.VCF_FILE, opt.INFO_FILE, opt.TargetPop1_Name, opt.TargetPop2_Name, opt.CHR, opt.WINDOW_SIZE, opt.SLIDING_SIZE
        Result = IntroCalc.CalcDxy(infile, ind, pop1, pop2, CHR, Window, Slide)
        if Result != "Done": print("Calulation ERROR")
else:
    print("Usage : python Intro.py -m [ Dstat | RND | Dxy ] ... ")
    sys.exit()
