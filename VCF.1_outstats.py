# -------------------------------------------------------------------------------
# Name:
# Purpose:     READ vcf: drop variants if MQ is too small or minor allele too rare
# -------------------------------------------------------------------------------

import sys
import numpy as np
from scipy.stats import binom

chrom = sys.argv[1]
src = open(chrom+".thin2.vcf", "r")
out1 = open(chrom+"thin2summary.txt", "w")
out2 = open(chrom+"Fam_summary.txt", "w")

PostProb_threshold = 0.99 # float(sys.argv[2]) # 0.99
Offspring_minor_threshold = 0.01 # float(sys.argv[3]) # 0.01
MinMAF = 0.05
minmd=6 # minimum depth to call mom genotype

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

momstats={}
for fam in Families:
    momstats[fam]=[0,0,0,0,0,0]


out1.write('Chrom\tpos')
for blk in Male_bulks:
    out1.write('\t'+blk)
out1.write('\tRR\tRA\tAA\tNN\tclutch_R\tclutch_A\tMomOff_test\tMomOff_fail\tp(mom)\tp(clutch)\tHW_chi2\n')

for line_idx, line in enumerate(src):
    cols = line.replace('\n', '').split('\t')

    if cols[0] == "#CHROM":
        pass  # out1.write(line)
    else:

        Qref = {"c": [], "m": [],   "Male_bulks":[]}
        Rdepth = {"c": [], "m": [], "Male_bulks":[]}
        datamb={}
        dataf={}
        for blk in Male_bulks:
            datamb[blk]=[0,0]
            for j in Male_bulks[blk]:
                vv = cols[j].split(":")  # 0/0:0,15,151:5,0
                if len(vv) == 3:
                    if vv[0] != './.':
                        datamb[blk][0]+= int(vv[2].split(",")[0])
                        datamb[blk][1]+= int(vv[2].split(",")[1])

        dataf=[0,0,0,0,0,0,0,0]
        for fam in Families:

            j=Families[fam][0] # mom
            MomG="NA"
            if j>0:
                vv = cols[j].split(":")  # 0/0:0,15,151:5,0
                if len(vv) == 3:
                    if vv[0] != './.':
                        depth_mom = int(vv[2].split(",")[0]) + int(vv[2].split(",")[1])
                        if depth_mom >= minmd and vv[0]=="0/0":
                            MomG="RR"
                            momstats[fam][0]+=1
                            dataf[0]+=1
                        elif depth_mom >= minmd and vv[0]=="0/1":
                            momstats[fam][1]+=1
                            dataf[1]+=1
                        elif depth_mom >= minmd and vv[0]=="1/1":
                            MomG="AA"
                            momstats[fam][2]+=1
                            dataf[2]+=1
                        else:
                            momstats[fam][3]+=1
                            dataf[3]+=1
                    else:
                        momstats[fam][3]+=1
                        dataf[3]+=1
            j=Families[fam][1] # clutch
            if j>0:
                vv = cols[j].split(":")  # 0/0:0,15,151:5,0
                if len(vv) == 3:
                    if vv[0] != './.':
                        R = int(vv[2].split(",")[0])
                        A = int(vv[2].split(",")[1])
                        dataf[4]+=R
                        dataf[5]+=A
                        if MomG=="RR":
                            momstats[fam][4]+=1
                            dataf[6]+=1
                            minor_allele_prob = binom.cdf(R,R+A,0.5)
                            if minor_allele_prob<=Offspring_minor_threshold:
                                momstats[fam][5]+=1
                                dataf[7]+=1
                        if MomG=="AA":
                            momstats[fam][4]+=1
                            dataf[6]+=1
                            minor_allele_prob = binom.cdf(A,R+A,0.5)
                            if minor_allele_prob<=Offspring_minor_threshold:
                                momstats[fam][5]+=1
                                dataf[7]+=1

        nmom=dataf[0]+dataf[1]+dataf[2]
        if nmom>0:
            px= float(dataf[0]+0.5*dataf[1])/float(nmom)
            pc= float(dataf[4])/float(dataf[4]+dataf[5])
            if min(px,1-px)>=MinMAF and min(pc,1-pc)>=MinMAF:
                # moms in HW?
                excc=[px*px*float(nmom),2*px*(1-px)*float(nmom),(1-px)*(1-px)*float(nmom)]
                chi2=0.0
                for u8 in range(3):
                    chi2+= (excc[u8]-dataf[u8])**2.0 / excc[u8]
                out1.write(cols[0]+'\t'+cols[1])
                for blk in Male_bulks:
                    out1.write('\t'+str(datamb[blk][0])+','+str(datamb[blk][1]))
                for j in range(8):
                    out1.write('\t'+str(dataf[j]))
                out1.write('\t'+str(px)+'\t'+str(pc)+'\t'+str(chi2)+'\n')

out1.close()

for fam in momstats:
    out2.write(fam)
    for j in range(6):
        out2.write('\t'+str(momstats[fam][j]))
    out2.write('\n')
