#!/bin/env python
import sys, time, math, copy, random, os, gzip
import numpy as np
from itertools import *
from scipy.stats import norm
import collections
try:
    from itertools import imap
except ImportError:
    imap=map

def mkdir(FolderName):
    rn = random.randint(100, 1000)
    check_list = []
    for i in range(rn):
        check = os.path.isdir(FolderName)
        check_list.append(check)
    if True in Check_list:
        return "Conti"
    else:
        return "Make"


def DZtest(Alist, Blist):
    if Alist == "NO_SNP" : return "NO_SNP"
    if Alist == "Zero_Div" : return "Zero_Div"
    Amean = np.mean(Alist)
    Bmean = np.mean(Blist)
    Asd = np.std(Alist)
    Bsd = np.std(Blist)
    Z = (Amean - Bmean) / np.sqrt(((Asd / len(Alist)) + (Bsd / len(Blist))))
    return round(Z, 4)

def Time():
    timeStamp = time.strftime("%H:%M:%S")
    return timeStamp

def Ztest(List, mu=0):
    z = (np.mean(List) - mu) / float((np.std(List)/np.sqrt(len(List))))
    pval = norm().sf(np.abs(z))*2
    return z


def DZtest(Alist, Blist):
    if Alist == "NO_SNP" : return "NO_SNP"
    if Alist == "Zero_Div" : return "Zero_Div"
    Amean = np.mean(Alist)
    Bmean = np.mean(Blist)
    Asd = np.std(Alist)
    Bsd = np.std(Blist)
    Z = (Amean - Bmean) / np.sqrt(((Asd / len(Alist)) + (Bsd / len(Blist))))
    return round(Z, 4)


def JackknifeSE(List):
    if len(set(List)) == 1 :
        JackTotalSd ="-0"
        return JackTotalSd
    JackMeanList = []
    count = 0
    for f in range(len(List)):
        count = count + 1
        JackList = copy.deepcopy(list(List))
        JackList[f:f+1] = []
        JackMean = np.mean(JackList)
        JackMeanList.append(JackMean)
    JackTotalMean = np.mean(JackMeanList)
    JackTotalMeanRep = [JackTotalMean] * len(JackMeanList)
    JackDiff = (np.array(JackMeanList) - np.array(JackTotalMeanRep))**2
    JackDiffSum = sum(JackDiff)
    JackDiffSumDiv = JackDiffSum / count
    JackTotalSd = round(np.sqrt(JackDiffSumDiv), 6)
    return JackTotalSd


def Bootstrap(Keys, Iter=100):
    BootLists = []
    for num in range(Iter):
        Resample = np.random.choice(list(Keys), len(Keys))
        BootLists.append(list(Resample))
    return BootLists

def Dstatistics(FreqDic, GroupList):
    FreqDicKeys = FreqDic.keys()
    FreqDicKeysBootList = Bootstrap(FreqDicKeys)
    DBootList = []
    tcount = []
    for Blist in FreqDicKeysBootList:
        PreD_deno = []
        PreD_no = []
        count = 0
        for key in Blist:
            SnpDataDic = FreqDic[key]
            popTAfreq = SnpDataDic[GroupList[0]][1]
            pop1Afreq = SnpDataDic[GroupList[1]][1]
            pop2Afreq = SnpDataDic[GroupList[2]][1]
            popOAfreq = SnpDataDic[GroupList[3]][1]
            if (popTAfreq + pop1Afreq + pop2Afreq + popOAfreq) == 0 : continue
            count = count + 1
            d_deno = ((1-popTAfreq)*pop1Afreq*pop2Afreq*(1-popOAfreq) - popTAfreq*(1-pop1Afreq)*pop2Afreq*(1-popOAfreq))
            d_no = ((1-popTAfreq)*pop1Afreq*pop2Afreq*(1-popOAfreq) + popTAfreq*(1-pop1Afreq)*pop2Afreq*(1-popOAfreq))
            PreD_deno.append(d_deno)
            PreD_no.append(d_no)
        if count == 0 :
            BlockD = "NO_SNP"
            BlockSe = "NO_SNP"
            DBootList = "NO_SNP"
            return BlockD, BlockSe, DBootList, count
        sum_deno = sum(PreD_deno)
        sum_no = sum(PreD_no)
        if sum_no == 0:
            BlockD = "Zero_Div"
            BlockSe = "Zero_Div"
            DBootList = "Zero_Div"
            return BlockD, BlockSe, DBootList, count
        D = float(sum_deno) / float(sum_no)
        DBootList.append(D)
        tcount.append(count)

    BlockD = round(np.mean(DBootList), 6)
    BlockSe = JackknifeSE(DBootList)
    tcountm = str(round(np.mean(tcount), 0))
    return BlockD, BlockSe, DBootList, str(tcountm)


def Dxy(FreqDic, GroupList):
    FreqDicKeys = FreqDic.keys()
    FreqDicKeysBootList = Bootstrap(FreqDicKeys)
    BootDxyList = []
    tcount = []
    for Blist in FreqDicKeysBootList:
        dxyiList = []
        count = 0
        for key in Blist:
            SnpDataDic = FreqDic[key]
            pop1Afreq = SnpDataDic[GroupList[0]][1]
            pop2Afreq = SnpDataDic[GroupList[1]][1]
            if (pop1Afreq + pop2Afreq) == 0 : continue
            count = count + 1
            dxyi = (1-pop1Afreq)*(pop2Afreq) + (pop1Afreq)*(1-pop2Afreq)

            dxyiList.append(dxyi)
        if count == 0 :
            Dxymean = "NO_SNP"
            Dxyse = "NO_SNP"
            BootDxyList = "NO_SNP"
            return Dxymean, Dxyse, BootDxyList, count
        dxy = sum(dxyiList) / float(len(dxyiList))
        BootDxyList.append(dxy)
        tcount.append(count)
    Dxymean = round(np.mean(BootDxyList), 6)
    Dxyse = JackknifeSE(BootDxyList)
    tcountm = str(round(np.mean(tcount), 0))
    return Dxymean, Dxyse, BootDxyList, tcountm


def RND(FreqDic, GroupList):
    FreqDicKeys = FreqDic.keys()
    FreqDicKeysBootList = Bootstrap(FreqDicKeys)
    BootRNDList = []
    tcount = []
    for Blist in FreqDicKeysBootList:
        dxyiList = []
        dxoiList = []
        dyoiList = []
        count = 0
        for key in Blist:
            SnpDataDic = FreqDic[key]
            pop1Afreq = SnpDataDic[GroupList[0]][1]
            pop2Afreq = SnpDataDic[GroupList[1]][1]
            popOAfreq = SnpDataDic[GroupList[2]][1]
            if (pop1Afreq + pop2Afreq + popOAfreq) == 0 : continue
            count = count + 1
            dxyi = (1-pop1Afreq)*(pop2Afreq) + (pop1Afreq)*(1-pop2Afreq)
            dxoi = (1-pop1Afreq)*(popOAfreq) + (pop1Afreq)*(1-popOAfreq)
            dyoi = (1-pop2Afreq)*(popOAfreq) + (pop2Afreq)*(1-popOAfreq)
            dxyiList.append(dxyi)
            dxoiList.append(dxoi)
            dyoiList.append(dyoi)

        if count == 0 :
            RNDmean = "NO_SNP"
            RNDse = "NO_SNP"
            BootRNDList = "NO_SNP"
            return RNDmean, RNDse, BootRNDList, count

        dxy = sum(dxyiList) / float(len(dxyiList))
        dxo = sum(dxoiList) / float(len(dxoiList))
        dyo = sum(dyoiList) / float(len(dyoiList))
        if (dxo + dyo) == 0:
            RNDmean = "NO_SNP"
            RNDse = "NO_SNP"
            BootRNDList = "NO_SNP"
            return RNDmean, RNDse, BootRNDList, count
        dout = (dxo + dyo) / 2
        RND = round(dxy / float(dout), 6)
        BootRNDList.append(RND)
        tcount.append(count)
    RNDmean = round(np.mean(BootRNDList), 6)
    RNDse = JackknifeSE(BootRNDList)
    tcountm = str(round(np.mean(tcount), 0))
    return RNDmean, RNDse, BootRNDList, tcountm



def Freq(Dic, GroupList, GroupIndexDic):
    TargetGroupIndexDic = {}
    for Glist in GroupList:
        TargetGroupIndexDic[Glist] = GroupIndexDic[Glist]
    BlockFreqDic = {}
    for SNP, Wline in Dic.items():
        for Gname, Gidx in TargetGroupIndexDic.items():
            AlleleList = [ Wline[x].split('|') for x in Gidx]
            AlleleCount = len(AlleleList) * 2
            AltAlleleCount = sum([ int(x)+int(y) for x, y, in AlleleList])
            RefAlleleCount = AlleleCount - AltAlleleCount
            AltAlleleFreq = round(AltAlleleCount / float(AlleleCount),8)
            RefAlleleFreq = round(RefAlleleCount / float(AlleleCount), 8)
            BlockFreqDic.setdefault(SNP, {}).setdefault(Gname, []).extend([RefAlleleFreq, AltAlleleFreq])
    return BlockFreqDic

def ReadingGroup(infile, ind):
    open_infile = gzip.open(infile)
    header = open_infile.readline().decode().rstrip('\n').rstrip('\r').split('\t')
    while header[0].startswith('#CHR'): break
    else: header = open_infile.readline().decode().rstrip('\n').rstrip('\r').split('\t')
    GroupDic = {}
    for line in imap(lambda x: x.rstrip('\n').rstrip('\r').strip().split('\t'), open(ind)):
        Group = line[1]
        Id = line[0]
        IdIdx = header.index(Id)
        GroupDic.setdefault(Group, []).append(IdIdx)
    return GroupDic

def CHRCheck(CHR):
    if CHR == "ALL":
        return "ALL"
    else:
        return CHR

def StatCHR(Dic):
    RNDlist = []
    for value in Dic.values():
          if value[0][-2] == "NO_SNP" : continue
          if value[0][-2] == "Zero_Div": continue
          RNDmean = float(value[0][-2])
          RNDlist.append(RNDmean)
    CHRRND = round(np.mean(RNDlist), 6)
    CHRRNDse = round(JackknifeSE(RNDlist), 6)
    return CHRRND, CHRRNDse, RNDlist


