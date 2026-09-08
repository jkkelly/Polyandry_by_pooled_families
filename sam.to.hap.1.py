# -*- coding: utf-8 -*-
import re
import sys
import os

MinMQ = 40

name = sys.argv[1]
chrom= sys.argv[2]
dirx=chrom+"files/"

snplist = open(chrom+".thin2.vcf.close.snps.50", 'r')
sam = open(dirx+name+".sam", 'r') # Sp25-1-56_S46_chr2L.sam
outfile = open(dirx+"hap1."+sys.argv[1]+".txt", 'w')

###########################################################

### First fill in a dictionary of snps of interest.
### that has structure snpDict[position] = [RefBase, AltBase]

snpDict = {}
for vline in snplist:
    vcols = vline.replace('\n', '').split('\t')
# chr2L   25708   T       C       4       25727   25729   25741   25753
# chr2L   25727   T       G       3       25729   25741   25753
    if vcols[0]==chrom:
        snpTig = vcols[0]
        snpPos = vcols[1]
        snpRef = vcols[2]
        snpAlt = vcols[3]
        snpDict[snpPos] = [snpRef, snpAlt]

#------------------------------------------------------------------------------------------

cigar_regex = re.compile(r'(\d+[MIDSHP])') # regex search for operations within cigar
cigarpart_regex = re.compile(r'(\d+)([MIDSHP])') # regex search to split up individual cigar ops


for line in sam:

# LH00663:73:23HKV2LT4:4:2145:45398:14116 83      chr2L   5831961 0       151M    =       5831800 -312    TCACTCCAAGCATAGATCTCCATTATATTGTCAATTGATCCTTTTAGTCTTTGGATTAATTCACTAAGCAGGTGAGCTGCACACAGCTCGAGTTTGGGAATTGTCTTCCTATTTTTTATAGGGTTGACTCTACTTTTGCTAGCTATTATAT       IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII       MD:Z:151        PG:Z:MarkDuplicates     RG:Z:4  NM:i:0  AS:i:151        XS:i:151        MQ:i:0  MC:Z:151M       ms:i:6040

        cols = line.replace('\n', '').split('\t')
        if cols[0][0]=="@":
                pass
        elif cols[2]==chrom: # do not read headers

            mq=float(cols[4])
            if mq>=MinMQ:
                readname = cols[0]
                contig = cols[2]
                start = int(cols[3])
    
                cigar = cols[5]
                SEQ = cols[9]
                QUAL = cols[10]
    
                # use cigar string to create a final read sequence that can be related to the reference coordinates.
    
                cigarOps = cigar_regex.findall(cigar) # divides cigar string into components
                readsite = 0 # keeps track of where we are on the read (SEQ string)
                readseq = '' # string of base calls for read
                readqual = '' # string of quality scores for read
    
                for opp in cigarOps:
                    opLength = int(cigarpart_regex.search(opp).group(1)) # the number part
                    opCode = cigarpart_regex.search(opp).group(2) # the letter part
    
                    if opCode == 'H': # if the operation is a hard clip, don't record
                        pass
    
                    elif opCode == 'I' or opCode == 'S': # if the operation is an insertion or soft clip
                        for i in range(opLength):
                            readsite += 1 # don't record information to readseq or readqual.
    
                    if opCode == 'M': # if the operation is a matching alignment
                        for i in range(opLength):
                            readseq = readseq + SEQ[readsite]
                            readqual = readqual + QUAL[readsite]
                            readsite += 1
    
                    elif opCode == 'D': # if the operation is a deletion
                        for i in range(opLength):
                            readseq = readseq + '.' # record blanks to readseq and readqual at this location
                            readqual = readqual + '.'
    
                ##### Now we look for the SNPs / indels
                LocusList = []
                AlleleList = []
    
                for j in range(len(readseq)):
                    currSite = start + j
    
                    try:
                        refBase = snpDict[str(currSite)][0]
                        altBase = snpDict[str(currSite)][1]
    
                        if readseq[j] == refBase and (ord(readqual[j])-33) >= 20:
                            LocusList.append(currSite)
                            AlleleList.append('R')
    
                        elif readseq[j] == altBase and (ord(readqual[j])-33) >= 20:
                            LocusList.append(currSite)
                            AlleleList.append('A')
    
                    except KeyError:
                        pass
    
                if len(LocusList) >= 1:
                    outfile.write(readname+'\t'+contig)
                    for k in range(len(LocusList)):
                        outfile.write('\t'+str(LocusList[k])+'\t'+str(AlleleList[k]))
                    outfile.write('\n')

outfile.close()


