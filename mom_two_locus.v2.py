# -------------------------------------------------------------------------------
# Name:
# v2:     align genotype to haplotype
# -------------------------------------------------------------------------------

import sys
import numpy as np

chrom =sys.argv[1]
mom= sys.argv[2] # Sp25_3_1_S238_L002


src = open(mom+"."+chrom+".twolocus", "r") # 
src2= open(chrom+"files/hap1."+mom+"_"+chrom+".txt", "r") #  chr2Lfiles/hap1.Sp25_4_36P_S289_L003_chr2L.txt
out1 = open(mom+"."+chrom+".twolocus.aligned", "w")

snppair={}
for line_idx, line in enumerate(src):
    cols = line.replace('\n', '').split('\t')

# Sp25_3_1_S238_L002      chr2L   25708   25727   0/1_15,5        0/1_14,6
# Sp25_3_1_S238_L002      chr2L   25727   25729   0/1_14,6        0/1_14,6
    snppair[cols[2]+"_"+cols[3]]={"s1":cols[4],"s2":cols[5],"RR":0,"RA":0,"AR":0,"AA":0}
src.close()


for line_idx, line in enumerate(src2):
    cols = line.replace('\n', '').split('\t') 
# LH00663:73:23HKV2LT4:2:2222:40285:14144 chr2L   25708   R
# LH00663:73:23HKV2LT4:2:2453:1283:6143   chr2L   25708   R
# LH00663:73:23HKV2LT4:2:1145:46692:25213 chr2L   25708   R       25727   R       25729   R
    snps = int( (len(cols)-2)/2 )
    if snps>1:
        for snp1 in range(snps-1):
            allele1 = cols[snp1*2+3]
            snp2=snp1+1
            allele2 = cols[snp2*2+3]
            try:
                uk=snppair[ cols[snp1*2+2]+"_"+cols[snp2*2+2] ]
                happy=allele1+allele2
                snppair[ cols[snp1*2+2]+"_"+cols[snp2*2+2] ][happy]+=1
            except KeyError:
                pass
src2.close()

for par in snppair:
    out1.write(mom+'\t'+chrom+'\t'+par+'\t'+snppair[par]["s1"]+'\t'+snppair[par]["s2"]+'\t'+str(snppair[par]["RR"])+'\t'+str(snppair[par]["RA"])+'\t'+str(snppair[par]["AR"])+'\t'+str(snppair[par]["AA"])+'\n')

out1.close()
