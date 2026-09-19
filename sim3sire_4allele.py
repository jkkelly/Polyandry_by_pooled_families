# -*- coding: utf-8 -*-
"""
simulate 4 allele data
"""
from numpy.random import binomial,multinomial
from random import random
import sys

sire1= float(sys.argv[1])
sire2= float(sys.argv[2])
Famsize= int(sys.argv[3])
m = 4*Famsize # int(sys.argv[4]) # read depth for families

Q = [0.25,0.25,0.25,0.25] # float(sys.argv[5]) # ref freq at all loci

pdfracs=[sire1,sire2,1.0-sire1-sire2]
Sires_per = len(pdfracs)
NoFams = 1000
NoSNPs = 200

simrep = "sim."+sys.argv[1]+"_"+sys.argv[2]+"_"+sys.argv[3]
out2=open(simrep+".4data.txt","w")

################ MAIN # simulate data

glist=[]
PopGenoFreqs={}
for j in range(4):
    for k in range(j,4):
        gx=str(j)+str(k)
        glist.append(gx)
        if j==k:
            PopGenoFreqs[gx]=Q[j]**2.0
        else:
            PopGenoFreqs[gx]=Q[j]*Q[k]*2.0            
        print(gx,PopGenoFreqs[gx])


af=[0,0,0,0]
#genotype_matrix={}
for j in range(NoFams):
    #genotype_matrix[j]=[]
    out2.write(str(j)+"\t"+str(Famsize))

    for k in range(NoSNPs):

        offgenos={}
        for z in PopGenoFreqs:
            offgenos[z]=0

        mg='na' # choose mom geno
        h7 =random()
        cdf=PopGenoFreqs["00"]
        for x in range(len(glist)):
            if h7<cdf:
                mg=glist[x]
                break
            else:
                cdf+= PopGenoFreqs[glist[x+1]]
                
        for pgx in range(Sires_per):

            fdsize = int(pdfracs[pgx]*Famsize)
            dg='na' # choose dad geno for sire pgx
            h7 =random()
            cdf=PopGenoFreqs["00"]
            for x in range(len(glist)):
                if h7<cdf:
                    dg=glist[x]
                    break
                else:
                    cdf+= PopGenoFreqs[glist[x+1]]            

            for ox in range(fdsize): #make babies
                if random()<0.5:
                    mallele=int(mg[0])
                else:
                    mallele=int(mg[1])
                if random()<0.5:
                    dallele=int(dg[0])
                else:
                    dallele=int(dg[1])                    
                    
                if mallele<dallele:
                    og=str(mallele)+str(dallele)
                else:
                    og=str(dallele)+str(mallele)                    

                offgenos[og]+=1
        
        af[0]=float( offgenos["00"]+0.5*(offgenos["01"]+offgenos["02"]+offgenos["03"]) )/float(Famsize)
        af[1]=float( offgenos["11"]+0.5*(offgenos["01"]+offgenos["12"]+offgenos["13"]) )/float(Famsize)
        af[2]=float( offgenos["22"]+0.5*(offgenos["02"]+offgenos["12"]+offgenos["23"]) )/float(Famsize)
        af[3]=float( offgenos["33"]+0.5*(offgenos["03"]+offgenos["13"]+offgenos["23"]) )/float(Famsize)

        nx0,nx1,nx2,nx3=multinomial(m, af)
        # print(af,nx0,nx1,nx2,nx3)
        out2.write("\t"+mg+":"+str(nx0)+","+str(nx1)+","+str(nx2)+","+str(nx3))
    out2.write("\n")

out2.close()

