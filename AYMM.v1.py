#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# v2: Include maternal genotype in classification

import sys
import numpy as np
from scipy.stats import binom

chrom = sys.argv[1] # chr2L
focal_fam = sys.argv[2] # 
PostProb_threshold = 0.99 # float(sys.argv[2]) # 0.99
Offspring_minor_threshold = 0.01 # float(sys.argv[3]) # 0.01

invcf=chrom+".genotype_data.txt" # 
out1=open(chrom+"."+focal_fam+".AYMM.summary.txt","w")


################ MAIN
dist={}
for j in range(1,22):
    dist[j]={}
    for k in range(j+1):
        dist[j][k]={ "unclear":0,"RR":0,"RA":0,"AA":0  }

matfams={}
famstats={}

dubs={}
in1 = open("dubious_fams.txt", "r")
for line_idx, line in enumerate(in1):
    cols = line.replace('\n', '').split('\t')
    dubs[cols[0]]=1
in1.close()



famid={}
inx=open(invcf, "r") #
for line_idx, line in enumerate(inx):
    cols = line.replace('\n', '').split('\t')
# chrom   pos     Ref     Alt     M1a     M1b     M2a     M2b     M3a     M3b     M4a     M4b     1_100   1_101   1_113 ... 4_99    spinsters
    if line_idx==0:
        for j in range(12,len(cols)-1):
            try:
                uz = dubs[cols[j]]
                if cols[j]==focal_fam:
                    Focal_j = j
                famid[j]=cols[j]
                famstats[cols[j]]=[0,0]
            except KeyError:
                pass
    elif cols[0]=="chrom":
        pass            
    else:
    # chr2L   7156    T       C       272,53  898,135 435,58  436,46  673,102 357,68  359,62  896,111   0/0:0,21,230:7,0_19,6   0/0:0,9,100:3,0_36,0    0/0:0,21,229:7,0_115,0  0/0:0,6,75:2,0_53,0     ./.:0,0,0:0,0_72,0    0/1:255,0,255:20,19_46,14       0/0:0,33,255:11,0_40,0  0/0:0,39,255:13,0_36,0  0/1:93,0,158:4,3_31,25     
        snpstats=[0,0]
        matr=[] # obtain initial estimate of allele frequency across all maternal flies
        for j in famid:
            ma,clutch=cols[j].split("_")
            GT=ma.split(":")[0]
            if GT=="0/0":
                matr.append(1.0)
            elif GT=="0/1":
                matr.append(0.5)
            elif GT=="1/1":
                matr.append(0.0)
        pRmoms=np.average(matr)


        ma,clutch=cols[Focal_j].split("_") 
        MomG="unclear"
        GT=ma.split(":")[0] # 1/1:101,14,0:1,12 but also really weird shit:: 1/1:0,9,0:26,25    0/0:0,5,5:3,1
        pl1=ma.split(":")[1]
        ad1=ma.split(":")[2]
        if len(pl1.split(","))==3:
            fred = [float(pl1.split(",")[0]),float(pl1.split(",")[1]),float(pl1.split(",")[2])]
            gl = [ 10.0**(-0.1*fred[0]),10.0**(-0.1*fred[1]),10.0**(-0.1*fred[2]) ]
            rx = pRmoms*pRmoms*gl[0]+2*(1-pRmoms)*pRmoms*gl[1]+(1-pRmoms)*(1-pRmoms)*gl[2]
            pp = [ pRmoms*pRmoms*gl[0]/rx, 2*(1-pRmoms)*pRmoms*gl[1]/rx, (1-pRmoms)*(1-pRmoms)*gl[2]/rx ]
            if   pp[0]>=PostProb_threshold:
                MomG="RR"
            elif pp[1]>=PostProb_threshold:
                MomG="RA"
            elif pp[2]>=PostProb_threshold:
                    MomG="AA"
                    
        if MomG=="RR" or  MomG=="AA":               
            for j in famid:
                ma,clutch=cols[j].split("_")    
                R=int(clutch.split(",")[0])
                A=int(clutch.split(",")[1])
                dp = R+A
    
                if MomG=="RR":
                    famstats[famid[j]][0]+=1
                    minor_allele_prob = binom.cdf(R,dp,0.5)
                    if minor_allele_prob<=Offspring_minor_threshold:
                        famstats[famid[j]][1]+=1
                elif MomG=="AA":
                    famstats[famid[j]][0]+=1
                    minor_allele_prob = binom.cdf(A,dp,0.5)
                    if minor_allele_prob<=Offspring_minor_threshold:
                        famstats[famid[j]][1]+=1

inx.close()


for fam in famstats:
    out1.write( str(focal_fam)+'\t'+str(fam)+'\t'+str(famstats[fam][0])+'\t'+str(famstats[fam][1])+'\n'  )

out1.close()

