# -------------------------------------------------------------------------------
# Name:
# Purpose:     Make list of intermediate frequency snps by position
# -------------------------------------------------------------------------------

import sys
import numpy as np

infile =sys.argv[1] 
src = open(infile, "r") # should be ordered vcf
out1 = open(infile+".snplist", "w")

MinMAF = 0.1
MinMoms= 300
# minmd=6 # minimum depth to call mom genotype

Families={}
Male_bulks={}
in1 = open("SP25.vcfkey", "r")
for line_idx, line in enumerate(in1):
    cols = line.replace('\n', '').split('\t')
# vcf column    name    type    collection      Family  Family_name
# 515   Sp25_3_1_S238_L002      Mother  3       1       3_1
# 514   Sp25_3_1P_merged        clutch  3       1       3_1
# 2046  Sp25_M3a_S755_L007_part_0531    Male_pool       M3a     Male    M3a_Male
# 2047  Sp25_M3b_S756_L007      Male_pool       M3b     Male    M3b_Male
    if line_idx > 0:
        if cols[2]=="Male_pool":
            try:
                Male_bulks[cols[3]].append(int(cols[0]))
            except KeyError:
                Male_bulks[cols[3]]=[int(cols[0])]
        else:
            try:
                uk=Families[cols[5]]
            except KeyError:
                Families[cols[5]]=[-9,-9]
            if cols[2]=="Mother":
                Families[cols[5]][0]=int(cols[0])
            elif cols[2]=="clutch":
                Families[cols[5]][1]=int(cols[0])
in1.close()



for line_idx, line in enumerate(src):
    cols = line.replace('\n', '').split('\t')

# chr2L   7088    .       A       T       77706.7 .       DP=35863;VDB=1.80735e-10;SGB=7519.4;RPBZ=1.41804;MQBZ=0.00498717;MQSBZ=-0.00765249;BQBZ=-1.76174;SCBZ=-0.521751;MQ0F=2.78839e-05;AC=729;AN=4000;DP4=9393,15326,2350,4022;MQ=59        GT:PL:AD        ./.:0,0,0:0,0 

    dataf=[0,0,0]
    refa=cols[3]
    alta=cols[4]
    for fam in Families:

        j=Families[fam][0] # mom
        if j>0:
            vv = cols[j].split(":")  # 0/0:0,15,151:5,0
            if len(vv) == 3:
                if vv[0] != './.':
                    depth_mom = int(vv[2].split(",")[0]) + int(vv[2].split(",")[1])
                    if vv[0]=="0/0":
                        dataf[0]+=1
                    elif vv[0]=="0/1":
                        dataf[1]+=1
                    elif vv[0]=="1/1":
                        dataf[2]+=1
    n=sum(dataf)
    if n>=MinMoms:
        Q = float(dataf[0]+0.5*dataf[1])/float(n)
        if min(Q,1-Q)>=MinMAF:
            out1.write( cols[0]+'\t'+cols[1]+'\t'+refa+'\t'+alta+'\t'+str(Q) )
            for j in range(3):
                out1.write('\t'+str(dataf[j]))
            out1.write('\n')

out1.close()

