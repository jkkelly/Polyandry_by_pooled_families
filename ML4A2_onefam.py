# -*- coding: utf-8 -*-
"""
Created on Sun May 17 10:10:04 2026
# ML2n :: normal approx, account for sampling of het dads

# Assumes ClusterBlocks = 10 for 3 sires
# Assumes ClusterBlocks = 10 for 4 sires

# autosomal inheritance

# method2 :: no discrete blocks
# method3 :: updated mean/variance calcs
# method4 :: include maternal hets

# calculate locus specific normalization factors

"""
from scipy.stats import norm
from math import log,ceil,sqrt
import numpy as np
import sys

OutputDir = "SimData/"

FILEPREFIX= sys.argv[1]
focal_fam= sys.argv[2]
out1=open(OutputDir+FILEPREFIX+"."+focal_fam+".A2.out","w")
# out2=open(FILEPREFIX+"."+focal_fam+".A2.bySNP.txt","w")
ClutchBlocks=5 # Max No Sires
granularity = 20 # determine fraction intervals
MinSingleSnpProb = 10.0**(-6.0)
x2threshold = 3.84

# src  =open(FILEPREFIX+".snpdata.txt", "r") # autosomal.snpdata.txt
src  =open(FILEPREFIX+".2data.txt", "r") # simulated data

def find_combinations(target_sum, num_parts): # this generates the possible fractions contributed by each male sire to the offspring pool
    results = []

    def backtrack(remaining_sum, current_combination, last_value):
        # Base Case: If we have selected all 4 numbers
        if len(current_combination) == num_parts:
            if remaining_sum == 0:
                results.append(list(current_combination))
            return
        # Optimization: If the remaining sum is too small to even match the decreasing order constraint, we can stop early.
        if remaining_sum < (num_parts - len(current_combination)):
            return
        # To enforce a >= b >= c >= d, we can generate them in reverse order (d, c, b, a)
        # where each subsequent number must be greater than or equal to the last.
        # This naturally avoids duplicates and enforces the order when reversed at the end.
        for next_value in range(last_value, remaining_sum + 1):
            # Remaining elements left to choose
            elements_left = num_parts - len(current_combination) - 1
            # Pruning: If the next_value is too large to leave enough room for the rest
            if next_value * elements_left > (remaining_sum - next_value):
                break
            current_combination.append(next_value)
            backtrack(remaining_sum - next_value, current_combination, next_value)
            current_combination.pop() # Backtrack

    # Start with a minimum value of 1 since they must be positive integers
    backtrack(target_sum, [], 1)
    # Reverse each combination to match the a >= b >= c >= d format
    return [comb[::-1] for comb in results]

def LL_fam_final_detailed(delt):
    llx=0.0
    for snpID in range(NoSNPs):
        L1=fam1(snpID,delt)
#        out2.write(focal_fam+'\t'+str(snpID)+'\t'+str(Q[snpID])+'\t'+str(genotype_matrix[famID][snpID][0])+'\t'+str(genotype_matrix[famID][snpID][1])+'\t'+str(genotype_matrix[famID][snpID][2])+'\t'+str(L1)+'\n')
        if L1>0.0:
            llx+=log(L1)
        else:
            print("ZeroProb Error",snpID)
    return llx


def LL_fam(delt):
    llx=0.0
    for snpID in range(NoSNPs):
        L1=fam1(snpID,delt)
        if L1>0.0:
            llx+=log(L1)
        else:
            print("ZeroProb Error",snpID)
    return llx

def fam1(snpID,delt):

    q = Q[snpID]
    Sires = len(delt)
    PriorDad={}
    PriorDad[0]= q*q
    PriorDad[1]= 2*q*(1-q)
    PriorDad[2]= 1.0-PriorDad[0]-PriorDad[1]

    MG=genotype_matrix[famID][snpID][0]
    Rcc = genotype_matrix[famID][snpID][1]
    Mreads=genotype_matrix[famID][snpID][1]+genotype_matrix[famID][snpID][2]
    if Mreads==0:
        return 1.0
    if MG=="RR":
        Qmom = 1.0
    elif MG=="AA":
        Qmom = 0.0
    elif MG=="RA":
        Qmom = 0.5
    else: # no data for mom
        return 1.0

    ttot_prob=0.0
    for mglist in posg[Sires]:
        parcc=[0,0,0]
        mgprob=1.0
        for x in range(len(mglist)):
            pgx = mglist[x] # 0, 1, 2
            mgprob*= PriorDad[pgx]
            parcc[pgx]+=delt[x]

        if abs(sum(parcc)-1.0) > 0.0001:
            print("fail 47",sum(parcc))
        QDad = parcc[0]+0.5*parcc[1]

        EQoff=0.5*QDad+0.5*Qmom
        if MG=="RA":
            VarQoff=(0.25 + QDad*(1-QDad))/(4.0*Famsize[famID])
        else:
            VarQoff=(QDad*(1-QDad))/(4.0*Famsize[famID])

        prob2=0.0
        if EQoff<0.0001 and Rcc==0:
            prob2=1.0
        elif EQoff>0.9999 and Rcc==Mreads:
            prob2=1.0
        elif EQoff>0.0001 and EQoff<0.9999:
            mu=Mreads*EQoff
            var = Mreads*EQoff*(1.0-EQoff)
            var+= float(Mreads*(Mreads-1))*VarQoff # add second term if MG = RA
            # print(parcc,Mreads*off_pa*(1.0-off_pa),float(Mreads*(Mreads-1))*freq01/float(16*Famsize))
            z1=(Rcc+0.5-mu)/sqrt(var)
            z0=(Rcc-0.5-mu)/sqrt(var)
            prob2 = norm.cdf(z1)-norm.cdf(z0)
        ttot_prob+= (prob2 * mgprob)

# NORMALIZE here
    ttot_prob=ttot_prob/maxprob[snpID]
    if ttot_prob<MinSingleSnpProb:
        ttot_prob=MinSingleSnpProb
#    print(snpID,delt,ttot_prob)
    return ttot_prob





######################################################
# Globals
Q={}
delta={}
Famsize={}
genotype_matrix={}
maxprob={}
for line_idx, line in enumerate(src):
    cols = line.replace('\n', '').split('\t')
    # 0 100     2:56,144        2:0,200 1:110,90        1:61,139        0:142,58        0:134,66
    famid =cols[0]
    Famsize[famid]= int(cols[1])

    NoSNPs = len(cols)-2
    genotype_matrix[famid]=[]
    for j in range(NoSNPs):

        vv=cols[j+2].split(":")
        ref=int(vv[1].split(",")[0])
        alt=int(vv[1].split(",")[1])
        if famid==focal_fam:
            genotype_matrix[famid].append( [vv[0],ref,alt] )
            maxprob[j]=1.0
            if ref>0 and alt>0:
                Mreads=ref+alt
                px=float(ref)/Mreads
                if px>0 and px<1.0:
                    mu=Mreads*px
                    var = Mreads*px*(1.0-px)
    
                    z1=(ref+0.5-mu)/sqrt(var)
                    z0=(ref-0.5-mu)/sqrt(var)
                    maxprob[j] = norm.cdf(z1)-norm.cdf(z0)

        if alt+ref>0:
            pr=float(ref)/float(alt+ref)
            try:
                Q[j].append(pr)
            except KeyError:
                Q[j]=[pr]
src.close()

NoFamilies = len(Famsize.keys())

for j in range(NoSNPs):
    rx = np.average(Q[j])
    Q[j]=rx # changed to float

print("no snps",NoSNPs)

posg={} # all possible sire genotype given the number of sires
posg[1]=[[0],[1],[2]]
for j in range(2,ClutchBlocks+1):
    posg[j]=[]
    for x in range(len(posg[j-1])):
        new0=[]
        new1=[]
        new2=[]
        for y in range(j-1):
            new0.append(posg[j-1][x][y])
            new1.append(posg[j-1][x][y])
            new2.append(posg[j-1][x][y])
        new0.append(0)
        new1.append(1)
        new2.append(2)
        posg[j].append(new0)
        posg[j].append(new1)
        posg[j].append(new2)



sol3 = find_combinations(granularity, 3)
sets3=[]
for triplet in sol3:
    x=[]
    for j in range(3):
        x.append(float(triplet[j])/float(granularity))
    sets3.append(x)

sol4 = find_combinations(granularity, 4)
sets4=[]
for quad in sol4:
    x=[]
    for j in range(4):
        x.append(float(quad[j])/float(granularity))
    sets4.append(x)



famID = focal_fam # famID is now locked to focal_fam
allbest={}
best={}
# One sire
# change here:: delta is now a list with proportions to sire 1, sire 2, sire 3, ...
delta[famID]=[1.0]
ll1=LL_fam(delta[famID])
best[1]=[ll1,delta[famID]]
print(famID,"one sire",best[1][0],best[1][1])
allbest=best[1]

# two sires
best[2]=[-999999999,"NA"]
for a1 in [0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95]:
    delta[famID]=[a1,1-a1]
    llx=LL_fam(delta[famID])
    if llx>best[2][0]:
        ll2=llx
        best[2]=[llx,delta[famID]]
print(famID,"2sire",best[2][0],best[2][1],ll2-ll1)

if ll2-ll1>x2threshold:
    allbest=best[2]
    best[3]=[-999999999,"NA"]
    for a1 in sets3:
        delta[famID]=a1
        llx=LL_fam(delta[famID])

        if llx>best[3][0]:
            ll3=llx
            best[3]=[llx,delta[famID]]
    print(famID,"3sire",best[3][0],best[3][1],ll3-ll2)

    if ll3-ll2>x2threshold:
        allbest=best[3]
        best[4]=[-999999999,"NA"]
        for a1 in sets4:
            delta[famID]=a1
            llx=LL_fam(delta[famID])

            if llx>best[4][0]:
                ll4=llx
                best[4]=[llx,delta[famID]]
        print(famID,"4sire",best[4][0],best[4][1],ll4-ll3)
        if ll4-ll3>x2threshold:
            allbest=best[4]

out1.write(str(famID)+'\t'+str(Famsize[famID])+'\t'+str(allbest[0])+'\t'+str(len(allbest[1]))+'\t'+str(allbest[1])+'\n')
print(str(famID)+'\tchosen\t'+str(Famsize[famID])+'\t'+str(allbest[0])+'\t'+str(allbest[1]))

delta[famID]=allbest[1]
# llz = LL_fam_final_detailed(delta[famID])

out1.close()
# out2.close()
