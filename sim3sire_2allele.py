# -*- coding: utf-8 -*-
"""
Created on Sun May 17 10:05:23 2026

@author: jackj
"""
from numpy.random import binomial,multinomial
from random import random
import sys

Q = 0.5 # float(sys.argv[5]) # ref freq at all loci
sire1= float(sys.argv[1])
sire2= float(sys.argv[2])
Famsize= int(sys.argv[3])
m = 4*Famsize # int(sys.argv[4]) # read depth for families


pdfracs=[sire1,sire2,1.0-sire1-sire2]
Sires_per = len(pdfracs)
NoFams = 1000
NoSNPs = 1000
simrep = "sim2."+sys.argv[1]+"_"+sys.argv[2]+"_"+sys.argv[3]
out2=open(simrep+".2data.txt","w")

################ MAIN # simulate data


#genotype_matrix={}
for j in range(NoFams):
    #genotype_matrix[j]=[]
    out2.write(str(j)+"\t"+str(Famsize))

    for k in range(NoSNPs):

        h7 =random()
        if h7< Q**2.0:
            mg=0
        elif h7<(Q**2.0+2*Q*(1-Q)):
            mg=1
        else:
            mg=2
        offgc=[0,0,0]
        if mg==0:
            for pgx in range(Sires_per):
                fdsize = int(pdfracs[pgx]*Famsize)
                h7 =random()
                if h7< Q**2.0:
                    offgc[0]+=fdsize
                elif h7<(Q**2.0+2*Q*(1-Q)):
                    het_cc=binomial(fdsize, 0.5)
                    offgc[0]+=(fdsize-het_cc)
                    offgc[1]+=het_cc
                else:
                    offgc[1]+=fdsize

        elif mg==1:
            for pgx in range(Sires_per):
                fdsize = int(pdfracs[pgx]*Famsize)
                h7 =random()
                if h7< Q**2.0: # dad is R/R
                    het_cc=binomial(fdsize, 0.5)
                    offgc[0]+=(fdsize-het_cc)
                    offgc[1]+=het_cc

                elif h7<(Q**2.0+2*Q*(1-Q)):
                    x11,x12,x22=multinomial(fdsize, [0.25,0.5,0.25])
                    offgc[0]+=x11
                    offgc[2]+=x22
                    offgc[1]+=x12
                else:      # dad is A/A
                    het_cc=binomial(fdsize, 0.5)
                    offgc[2]+=(fdsize-het_cc)
                    offgc[1]+=het_cc


        elif mg==2:
            for pgx in range(Sires_per):
                fdsize = int(pdfracs[pgx]*Famsize)
                h7 =random()
                if h7< Q**2.0:
                    offgc[1]+=fdsize
                elif h7<(Q**2.0+2*Q*(1-Q)):
                    het_cc=binomial(fdsize, 0.5)
                    offgc[2]+=(fdsize-het_cc)
                    offgc[1]+=het_cc
                else:
                    offgc[2]+=fdsize

        #print(mg,offgc)
        q=float(offgc[0]+0.5*offgc[1])/float(sum(offgc))
        Rcc=binomial(m, q)
        if mg==0:
            out2.write("\tRR:"+str(Rcc)+","+str(m-Rcc))
        elif mg==1:
            out2.write("\tRA:"+str(Rcc)+","+str(m-Rcc))    
        if mg==2:
            out2.write("\tAA:"+str(Rcc)+","+str(m-Rcc))            
    out2.write("\n")

