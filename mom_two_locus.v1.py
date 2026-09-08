# -------------------------------------------------------------------------------
# Name:
# Purpose:     Make list of intermediate frequency snps by position
# -------------------------------------------------------------------------------

import sys
import numpy as np

chrom =sys.argv[1]
mom= sys.argv[2] # Sp25_3_1_S238_L002


src = open(chrom+".thin2.vcf", "r") # chr2L.thin2.vcf
out1 = open(mom+"."+chrom+".twolocus", "w")

in1 = open("SP25.vcfkey", "r")
for line_idx, line in enumerate(in1):
    cols = line.replace('\n', '').split('\t')
# vcf column    name    type    collection      Family  Family_name
# 515   Sp25_3_1_S238_L002      Mother  3       1       3_1
    if cols[1]==mom:
        Families=int(cols[0])
in1.close()


snppair=[]
pull={}
in1 = open(chrom+".thin2.vcf.close.snps.50", "r")
for line_idx, line in enumerate(in1):
    cols = line.replace('\n', '').split('\t')
# chr2L   25708   T       C       4       25727   25729   25741   25753
# chr2L   25727   T       G       3       25729   25741   25753
    snppair.append([int(cols[1]),int(cols[5])])
    pull[int(cols[1])]=1
    pull[int(cols[5])]=1
in1.close()

for line_idx, line in enumerate(src):
    cols = line.replace('\n', '').split('\t')

# chr2L   7088    .       A       T       77706.7 .       DP=35863;VDB=1.80735e-10;SGB=7519.4;RPBZ=1.41804;MQBZ=0.00498717;MQSBZ=-0.00765249;BQBZ=-1.76174;SCBZ=-0.521751;MQ0F=2.78839e-05;AC=729;AN=4000;DP4=9393,15326,2350,4022;MQ=59        GT:PL:AD        ./.:0,0,0:0,0
    try:
        uk=pull[int(cols[1])]
        j=Families # mom
        vv = cols[j].split(":")  # 0/0:0,15,151:5,0
        pull[int(cols[1])]=vv[0]+"_"+vv[2]

    except KeyError:
        pass


for j in range(len(snppair)):
    pos1=snppair[j][0]    
    pos2=snppair[j][1]
    out1.write(mom+'\t'+chrom+'\t'+str(pos1)+'\t'+str(pos2)+'\t'+str(pull[pos1])+'\t'+str(pull[pos2])+'\n')
out1.close()
