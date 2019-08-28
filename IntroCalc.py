#!/bin/env python
import sys, time, math, copy, random, os
import numpy as np
from itertools import *
from scipy.stats import norm
import collections
from IntroFunc import *

def CalcD(infile, ind, popO, popT, pop1, pop2, CHR, Window, Slide):
    Window = int(Window)
    Slide = int(Slide)
    outfile_name = "Dtest." + str(Window) + "." + str(Slide) + "." + CHR + "." + popT + "-" + pop1 + "-" + pop2 + "-" + popO + ".txt"
    FolderName = "./Dtest." + str(Window) + "." + str(Slide) + "." + popT + "-" + pop1 + "-" + pop2 + "-" + popO
    try :
        os.mkdir(FolderName)
    except OSError:
        pass
    else:
        if os.path.isdir(FolderName) == True:
            pass
        else:
            che = mkdir(FolderName)
            if che == "Conti":
                pass
            else:
                os.mkdir(FolderName)
    outfile_name = FolderName + "/" + outfile_name
    print(Time() +"    " + "Making Directory...")
    print("Directroy : ", outfile_name)
    GroupIndexDic = ReadingGroup(infile, ind)
    GroupList = [popT, pop1, pop2, popO]
    CHR = CHRCheck(CHR)
    ChrDic = collections.OrderedDict()
    for chrName,group in groupby(imap(lambda x: x.decode().rstrip('\n').rstrip('\r').strip().split('\t'), gzip.open(infile)), lambda y: y[0]):
        if chrName != CHR: continue
        print(Time()+ "    " + chrName + "    Reading VCF...")
        SlideCount = 0
        PosStart = 1
        BlockSNPDic = {}
        Begin = PosStart + (Slide * SlideCount)
        End = Window + (Slide * SlideCount)

        tmpKey = ""
        BlockDic = {}
        BlockNum = 1
        for line in group:
            Chr = line[0]
            Pos = int(line[1])
            if SlideCount == 0:
                if Pos > Window:
                    BlockFreq = Freq(BlockDic[BlockKey], GroupList, GroupIndexDic)
                    Dstat, Dse, DList, Dcount = Dstatistics(BlockFreq, GroupList)
                    outlist = [[Chr, str(BlockNum), str(Begin), str(End), str(len(BlockDic[BlockKey].keys())), str(Dcount), popT, pop1, pop2, popO, str(Dstat), str(Dse)], DList]
                    ChrDic[BlockNum] = outlist
                    tmpKey = BlockKey
                    SlideCount = SlideCount + 1
                    Begin = PosStart + (Slide * SlideCount)
                    End = Window + (Slide * SlideCount)
            else:
                if Pos > End:
                    tmpKeys = [x for x in BlockDic[tmpKey].keys() if x >= Begin]
                    for tk in tmpKeys:
                        tkline = BlockDic[tmpKey][tk]
                        BlockDic.setdefault(BlockKey, {}).setdefault(tk, []).extend(tkline)
                    del BlockDic[tmpKey]
                    BlockNum = BlockNum + 1
                    tmpKey = BlockKey
                    BlockFreq = Freq(BlockDic[BlockKey], GroupList, GroupIndexDic)
                    Dstat, Dse, DList, Dcount = Dstatistics(BlockFreq, GroupList)
                    outlist = [[Chr, str(BlockNum), str(Begin), str(End), str(len(BlockDic[BlockKey].keys())), str(Dcount), popT, pop1, pop2, popO, str(Dstat), str(Dse)], DList]
                    ChrDic[BlockNum] = outlist
                    SlideCount = SlideCount + 1
                    Begin = PosStart + (Slide * SlideCount)
                    End = Window + (Slide * SlideCount)

            BlockKey = str(Begin) + "-" + str(End)
            BlockDic.setdefault(BlockKey, {}).setdefault(Pos, []).extend(line)
        else:
            BlockNum = BlockNum + 1
            BlockFreq = Freq(BlockDic[BlockKey], GroupList, GroupIndexDic)
            Dstat, Dse, DList, Dcount = Dstatistics(BlockFreq, GroupList)
            outlist = [[Chr, str(BlockNum), str(Begin), str(End), str(len(BlockDic[BlockKey].keys())), str(Dcount), popT, pop1, pop2, popO, str(Dstat), str(Dse)], DList]
            ChrDic[BlockNum] = outlist
    else:
        print(Time() + "    Calculating CHR...")

    ChrDmean, ChrDse, ChrDMeanList= StatCHR(ChrDic)
    outfile = open(outfile_name, 'w')
    outfile_head = ["CHR", "BlockNum", "SNP", "Avg_Used_SNP", "Start", "End", "ThirdPop", "Pop1", "Pop2", "OutPop", "Dstat", "Dse", "ZDstat", "CHR_Dstat", "CHR_Dstat_se"]
    print('\t'.join(outfile_head), file=outfile)

    for BNum, Rlist in ChrDic.items():
        DZ = DZtest(Rlist[1], ChrDMeanList)
        printLine = Rlist[0] + [str(DZ), str(ChrDmean), str(ChrDse)]
        print('\t'.join(printLine), file=outfile)
    else:
        print(Time() + "    Done!")
    outfile.close()
    return "Done"

def CalcRND(infile, ind, popO, pop1, pop2, CHR, Window, Slide):
    Window = int(Window)
    Slide = int(Slide)
    outfile_name = "RND." + str(Window) + "." + str(Slide) + "." + CHR + "." + pop1 + "-" + pop2 + "-" + popO + ".txt"
    FolderName = "./RND." + str(Window) + "." + str(Slide) + "." +  pop1 + "-" + pop2 + "-" + popO

    try :
        os.mkdir(FolderName)
    except OSError:
        pass
    else:
        if os.path.isdir(FolderName) == True:
            pass
        else:
            che = mkdir(FolderName)
            if che == "Conti":
                pass
            else:
                os.mkdir(FolderName)

    outfile_name = FolderName + "/" + outfile_name
    outfile = open(outfile_name, 'w')
    print(Time() +"    " + "Making Directory...")
    print("Directroy : ", outfile_name)
    GroupIndexDic = ReadingGroup(infile, ind)
    GroupList = [pop1, pop2, popO]
    CHR = CHRCheck(CHR)
    ChrDic = collections.OrderedDict()
    for chrName,group in groupby(imap(lambda x: x.decode().rstrip('\n').rstrip('\r').strip().split('\t'), gzip.open(infile)), lambda y: y[0]):
        if chrName != CHR: continue
        print(Time()+ "    " + chrName + "    Reading VCF...")
        SlideCount = 0
        PosStart = 1
        BlockSNPDic = {}
        Begin = PosStart + (Slide * SlideCount)
        End = Window + (Slide * SlideCount)

        BlockDic = {}
        tmpKey = ""
        BlockNum = 1
        for line in group:
            Chr = line[0]
            Pos = int(line[1])
            if SlideCount == 0:
                if Pos > Window:
                    BlockFreq = Freq(BlockDic[BlockKey], GroupList, GroupIndexDic)
                    RNDmean, RNDse, RNDList, RNDcount  = RND(BlockFreq, GroupList)
                    outlist = [[Chr, str(BlockNum), str(Begin), str(End), str(len(BlockDic[BlockKey].keys())), str(RNDcount), pop1, pop2, popO, str(RNDmean), str(RNDse)], RNDList]
                    ChrDic[BlockNum] = outlist
                    tmpKey = BlockKey
                    SlideCount = SlideCount + 1
                    Begin = PosStart + (Slide * SlideCount)
                    End = Window + (Slide * SlideCount)
            else:
                if Pos > End:
                    tmpKeys = [x for x in BlockDic[tmpKey].keys() if x >= Begin]
                    for tk in tmpKeys:
                        tkline = BlockDic[tmpKey][tk]
                        BlockDic.setdefault(BlockKey, {}).setdefault(tk, []).extend(tkline)
                    del BlockDic[tmpKey]
                    BlockNum = BlockNum + 1
                    tmpKey = BlockKey
                    BlockFreq = Freq(BlockDic[BlockKey], GroupList, GroupIndexDic)
                    RNDmean, RNDse, RNDList, RNDcount = RND(BlockFreq, GroupList)
                    outlist = [[Chr, str(BlockNum), str(Begin), str(End), str(len(BlockDic[BlockKey].keys())), str(RNDcount), pop1, pop2, popO, str(RNDmean), str(RNDse)], RNDList]
                    ChrDic[BlockNum] = outlist
                    SlideCount = SlideCount + 1
                    Begin = PosStart + (Slide * SlideCount)
                    End = Window + (Slide * SlideCount)

            BlockKey = str(Begin) + "-" + str(End)
            BlockDic.setdefault(BlockKey, {}).setdefault(Pos, []).extend(line)
        else:
            BlockNum = BlockNum + 1
            BlockFreq = Freq(BlockDic[BlockKey], GroupList, GroupIndexDic)
            RNDmean, RNDse, RNDList, RNDcount = RND(BlockFreq, GroupList)
            outlist = [[Chr, str(BlockNum), str(Begin), str(End), str(len(BlockDic[BlockKey].keys())), str(RNDcount), pop1, pop2, popO, str(RNDmean), str(RNDse)], RNDList]
            ChrDic[BlockNum] = outlist

    print(Time() + "    Calculating CHR...")
    ChrRNDmean, ChrRNDse, ChrRNDMeanList= StatCHR(ChrDic)
    outfile_head = ["CHR", "BlockNum", "SNP", "Avg_Used_SNP", "Start", "End", "Pop1", "Pop2", "OutPop", "RND", "RNDse", "ZRND", "CHR_RND", "CHR_RNDse"]
    print('\t'.join(outfile_head), file=outfile)

    for BNum, Rlist in ChrDic.items():
        DZ = DZtest(Rlist[1], ChrRNDMeanList)
        printLine = Rlist[0] + [str(DZ), str(ChrRNDmean), str(ChrRNDse)]
        print('\t'.join(printLine), file=outfile)
    else:
        print(Time() + "    Done!")
    outfile.close()
    return "Done"



def CalcDxy(infile, ind, pop1, pop2, CHR, Window, Slide):
    Window = int(Window)
    Slide = int(Slide)
    outfile_name = "Dxy." + str(Window) + "." + str(Slide) + "." + CHR + "." + pop1 + "-" + pop2 + ".txt"
    FolderName = "./Dxy." + str(Window) + "." + str(Slide) + "." +  pop1 + "-" + pop2

    try :
        os.mkdir(FolderName)
    except OSError:
        pass
    else:
        if os.path.isdir(FolderName) == True:
            pass
        else:
            che = mkdir(FolderName)
            if che == "Conti":
                pass
            else:
                os.mkdir(FolderName)

    outfile_name = FolderName + "/" + outfile_name
    outfile = open(outfile_name, 'w')
    print(Time() +"    " + "Making Directory...")
    print( "Directroy : ", outfile_name)
    GroupIndexDic = ReadingGroup(infile, ind)
    GroupList = [pop1, pop2]
    CHR = CHRCheck(CHR)
    ChrDic = collections.OrderedDict()
    for chrName,group in groupby(imap(lambda x: x.decode().rstrip('\n').rstrip('\r').strip().split('\t'), gzip.open(infile)), lambda y: y[0]):
        if chrName != CHR: continue
        print(Time()+ "    " + chrName + "    Reading VCF...")
        SlideCount = 0
        PosStart = 1
        BlockSNPDic = {}
        Begin = PosStart + (Slide * SlideCount)
        End = Window + (Slide * SlideCount)

        BlockDic = {}
        tmpKey = ""
        BlockNum = 1
        for line in group:
            Chr = line[0]
            Pos = int(line[1])
            if SlideCount == 0:
                if Pos > Window:
                    BlockFreq = Freq(BlockDic[BlockKey], GroupList, GroupIndexDic)
                    Dxymean, Dxyse, DxyList, Dcount = Dxy(BlockFreq, GroupList)
                    outlist = [[Chr, str(BlockNum), str(Begin), str(End), str(len(BlockDic[BlockKey].keys())), str(Dcount), pop1, pop2, str(Dxymean), str(Dxyse)], DxyList]
                    ChrDic[BlockNum] = outlist
                    tmpKey = BlockKey
                    SlideCount = SlideCount + 1
                    Begin = PosStart + (Slide * SlideCount)
                    End = Window + (Slide * SlideCount)
            else:
                if Pos > End:
                    tmpKeys = [x for x in BlockDic[tmpKey].keys() if x >= Begin]
                    for tk in tmpKeys:
                        tkline = BlockDic[tmpKey][tk]
                        BlockDic.setdefault(BlockKey, {}).setdefault(tk, []).extend(tkline)
                    del BlockDic[tmpKey]
                    BlockNum = BlockNum + 1
                    tmpKey = BlockKey
                    BlockFreq = Freq(BlockDic[BlockKey], GroupList, GroupIndexDic)
                    Dxymean, Dxyse, DxyList, Dcount = Dxy(BlockFreq, GroupList)
                    outlist = [[Chr, str(BlockNum), str(Begin), str(End), str(len(BlockDic[BlockKey].keys())), str(Dcount), pop1, pop2, str(Dxymean), str(Dxyse)], DxyList]
                    ChrDic[BlockNum] = outlist
                    SlideCount = SlideCount + 1
                    Begin = PosStart + (Slide * SlideCount)
                    End = Window + (Slide * SlideCount)

            BlockKey = str(Begin) + "-" + str(End)
            BlockDic.setdefault(BlockKey, {}).setdefault(Pos, []).extend(line)
        else:
            BlockNum = BlockNum + 1
            BlockFreq = Freq(BlockDic[BlockKey], GroupList, GroupIndexDic)
            Dxymean, Dxyse, DxyList, Dcount = Dxy(BlockFreq, GroupList)
            outlist = [[Chr, str(BlockNum), str(Begin), str(End), str(len(BlockDic[BlockKey].keys())), str(Dcount), pop1, pop2, str(Dxymean), str(Dxyse)], DxyList]
            ChrDic[BlockNum] = outlist

    print(Time() + "    Calculating CHR...")
    ChrDxymean, ChrDxyse, ChrDxyMeanList= StatCHR(ChrDic)
    outfile_head = ["CHR", "BlockNum", "SNP", "Avg_Used_SNP", "Start", "End", "Pop1", "Pop2", "Dxy", "Dxy_se", "ZDxy", "CHR_Dxy", "CHR_Dxy_se"]
    print('\t'.join(outfile_head), file=outfile)

    for BNum, Rlist in ChrDic.items():
        DZ = DZtest(Rlist[1], ChrDxyMeanList)
        printLine = Rlist[0] + [str(DZ), str(ChrDxymean), str(ChrDxyse)]
        print('\t'.join(printLine), file=outfile)
    else:
        print(Time() + "    Done!")

    outfile.close()
    return "Done"

