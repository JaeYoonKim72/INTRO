#!/bin/env R
library(ggplot2)
suppressPackageStartupMessages(library("argparse"))
parser = ArgumentParser()
parser$add_argument("METHOD", nargs=1, help="Method [Dstat | RND | Dxy]")
parser$add_argument("RESULT_FILE", nargs=1, help="Result text file for the Method")
parser$add_argument("CHR", nargs=1, help="Chromosome Number, integer")

args = parser$parse_args()
mode = args$METHOD
file = args$RESULT_FILE
chr = args$CHR

##============================================================================

if ( mode == "Dstat" ){
    Index = as.integer(11)
    Lab = "D-statistic Value"
    Ymin = as.double(-1.05)
    Ymax = as.double(1.05)}
if ( mode == "RND"){
    Index = as.integer(10)
    Lab = "RND Value"
    Ymin = as.double(0)
    Ymax = as.double(2.05)}
if ( mode == "Dxy" ){
    Index = as.integer(9)
    Lab = "dxy Value"
    Ymin = as.double(0)
    Ymax = as.double(1.05)}

print("Loading file...")
data=read.table(file, header=T)
pdata = data[which((data[,Index] != "NO_SNP" & data[,Index]!= "Zero_Div")),]
pos = round((pdata[,3] + pdata[,4])/2)
d = as.double(as.vector(pdata[,Index]))
pdata = data.frame(pos, d)
pdata[,1] = pdata[,1] / 1000000
xmax=max(as.vector(pdata[,1]))
xpo = max(as.vector(pdata[,1]))
lab = paste("Chr", chr, sep="")

print("Plot ..")
outfile_name = paste(strsplit(file, ".txt")[[1]], ".jpeg", sep="")
jpeg(file=outfile_name, width=1200, height=600, quality=100)
bp = ggplot(pdata, aes(x=pdata[,1], y=pdata[,2])) + geom_line(size=0.4)  + xlab("Distance (Mb)") + ylab(Lab) + theme(plot.background = element_blank()) + theme_bw(base_size=26) + ylim(c(Ymin, Ymax))
bp +  annotate("text", x=xpo, y=Ymax, label=lab, size=6, family="Airal") +  scale_x_continuous(breaks=seq(0, xmax, 5))
dev.off()
print("Done!")
