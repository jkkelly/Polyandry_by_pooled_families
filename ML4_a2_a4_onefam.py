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
# method4.vs :: search male fractions

"""

from scipy.stats import binom,norm,multivariate_normal
from math import log,ceil,sqrt
import numpy as np
import sys

OutputDir = "RealData/"

FILEPREFIX= sys.argv[1] # "chosen.snps1.txt"
FILEPREFIX1= sys.argv[2] # "set0" # sys.argv[1]
focal_fam= sys.argv[3] # "3_33" #

out1=open(OutputDir+FILEPREFIX+"."+focal_fam+".A2A4.out","w")
# out2=open(FILEPREFIX+"."+focal_fam+".A2A4.bySNP.txt","w")

ClutchBlocks=5 # Max No Sires
x2threshold = 3.84
MinSingleSnpProb = 10.0**(-6.0)

sm1 = 0.0001

granularity = 20 # determine fraction intervals

src  =open(FILEPREFIX+".2data.txt", "r") # autosomal.2data.txt
src1  =open(FILEPREFIX1+".4data.txt", "r") # set1.4data.txt


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
    for snpID in range(NoSNPs1):
        L1=fam4(snpID,delt)
#        out2.write(focal_fam+'\t4-allele\t'+str(snpID)+'\t'+str(Q1[snpID])+'\t'+str(genotype_matrix1[famID][snpID][0])+'\t'+str(genotype_matrix1[famID][snpID][1])+'\t'+str(genotype_matrix1[famID][snpID][2])+'\t'+str(genotype_matrix1[famID][snpID][3])+'\t'+str(genotype_matrix1[famID][snpID][4])+'\t'+str(L1)+'\n')
        llx+=log(L1)
    return llx


def LL_fam(delt):
    llx=0.0
    for snpID in range(NoSNPs):
        L1=fam1(snpID,delt)
        if L1>0.0:
            llx+=log(L1)
        else:
            print("zeroprob snp",snpID)

    for snpID in range(NoSNPs1):
        L4=fam4(snpID,delt)
        if L4>0.0:
            llx+=log(L4)
        else:
            print("zeroprob 4-allele",snpID)

    return llx

def fam4(snpID,delt): # family j probability for one 4-allele

    PriorDad={}
    for i in range(4):
        for j in range(i,4):
            if i==j:
                PriorDad[ str(i)+str(j) ]= Q1[snpID][i]**2.0
            else:
                PriorDad[ str(i)+str(j) ]= Q1[snpID][i]*Q1[snpID][j]*2.0

    Sires = len(delt)
    cv1=[ [0,0,0],[0,0,0],[0,0,0] ]
    CovQoff=[ [0,0,0],[0,0,0],[0,0,0] ]
    MG=genotype_matrix1[famID][snpID][0]
    OC = [genotype_matrix1[famID][snpID][1],genotype_matrix1[famID][snpID][2],genotype_matrix1[famID][snpID][3],genotype_matrix1[famID][snpID][4]]
    Mreads=sum(OC)
    if Mreads==0:
        return 1.0
    if MG=="NA": # no data for mom
        return 1.0
    QMOM=[0.0 for j in range(4)]
    for j in range(4):
        if MG==str(j)+str(j):
            QMOM[j]=1
        elif MG[0]==str(j) or MG[1]==str(j):
            QMOM[j]=0.5

    ttot_prob=0.0
    for mglist in posg1[Sires]: # loop over all possible collection of male genotypes
        parcc={} # count of male genotype z
        for z in glist:
            parcc[z]=0
        mgprob=1.0
        for x in range(len(mglist)): # mglist is something like ['01', '00', '00']
            pgx = mglist[x] # '00' or '01' ....
            mgprob*= PriorDad[pgx]
            parcc[pgx]+=delt[x]

        QDAD = [0.0 for j in range(4)]
        for z in glist:
            if parcc[z]>0.0:
                if z[0]==z[1]:
                    QDAD[int(z[0])]+=parcc[z]
                else:
                    QDAD[int(z[0])]+=0.5*parcc[z]
                    QDAD[int(z[1])]+=0.5*parcc[z]
        #print("da y",sum(QDAD))

        EQoff=[0.0 for j in range(4)]
        for j in range(4):
            EQoff[j]=0.5*QDAD[j]+0.5*QMOM[j]

        prob2=0.0
        done=0
        for j in range(4):
            if OC[j]>0 and EQoff[j]<sm1: # impossible outcome
                done=1
            elif OC[j]<Mreads and EQoff[j]>1-sm1: # impossible outcome
               done=1
            elif OC[j]==Mreads and EQoff[j]>1-sm1: # certain outcome
                prob2=1.0
                done=1

        if done==0:

            mu=[ Mreads*EQoff[0],Mreads*EQoff[1],Mreads*EQoff[2] ]

            for j in range(3):
                for k in range(3):
                    if j==k:
                        CovQoff[j][k]=(QMOM[j]*(1-QMOM[j]) + QDAD[j]*(1-QDAD[j]))/(4.0*Famsize[famID])
                        cv1[j][k]=Mreads*EQoff[j]*(1-EQoff[j]) + Mreads*(Mreads-1)*CovQoff[j][k]
                    else:
                        CovQoff[j][k]=(-QMOM[j]*QMOM[k] - QDAD[j]*QDAD[k])/(4.0*Famsize[famID])
                        cv1[j][k]= -Mreads*EQoff[j]*EQoff[k] + Mreads*(Mreads-1)*CovQoff[j][k]

            upper_bounds=[OC[0]+0.5,OC[1]+0.5,OC[2]+0.5]
            lower_bounds=[OC[0]-0.5,OC[1]-0.5,OC[2]-0.5]
            # print(cv1)
            prob2 = multivariate_normal.cdf(upper_bounds, mean=mu, cov=cv1, allow_singular=True, lower_limit=lower_bounds)

        ttot_prob+= (prob2 * mgprob)
    ttot_prob=ttot_prob/maxprob1[snpID]        
    if ttot_prob<MinSingleSnpProb:
        ttot_prob=MinSingleSnpProb

    return ttot_prob


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
    ttot_prob=ttot_prob/maxprob[snpID]            
    if ttot_prob<MinSingleSnpProb:
        ttot_prob=MinSingleSnpProb
    return ttot_prob



######################################################

Q={}
Q1={}
delta={}
Famsize={}
genotype_matrix={}
genotype_matrix1={}
maxprob={}
maxprob1={}

for line_idx, line in enumerate(src1):
    cols = line.replace('\n', '').split('\t')
#    0  100     22:0,0,100,0    02:28,22,50,0   23:25,0,23,52   13:0,31,41,28   11:0,59,0,41    23:0,21,22,57

    famid =cols[0]
    Famsize[famid]= int(cols[1])

    NoSNPs1 = len(cols)-2
    if line_idx==0:
        for j in range(NoSNPs1):
            Q1[j]=[[],[],[],[]]

    genotype_matrix1[famid]=[]
    for j in range(NoSNPs1):

        vv=cols[j+2].split(":")
        m=int(vv[1].split(",")[0])+int(vv[1].split(",")[1])+int(vv[1].split(",")[2])+int(vv[1].split(",")[3])
        if famid==focal_fam:
            maxprob1[j]=1.0
            genotype_matrix1[famid].append( [ vv[0],int(vv[1].split(",")[0]),int(vv[1].split(",")[1]),int(vv[1].split(",")[2]),int(vv[1].split(",")[3]) ] )
            cv1=[ [0,0,0],[0,0,0],[0,0,0] ]
            OC = [ int(vv[1].split(",")[0]),int(vv[1].split(",")[1]),int(vv[1].split(",")[2]),int(vv[1].split(",")[3]) ]
            Mreads=float(sum(OC))
            if Mreads>0:
                mu=[ OC[0],OC[1],OC[2] ]
                px=[ float(OC[0])/float(Mreads),float(OC[1])/float(Mreads),float(OC[2])/float(Mreads) ]
    
                for x1 in range(3):
                    for y1 in range(3):
                        if x1==y1:
                            cv1[x1][y1]=Mreads*px[x1]*(1-px[x1]) 
                        else:
                            cv1[x1][y1]=-Mreads*px[x1]*px[y1]
    
                upper_bounds=[OC[0]+0.5,OC[1]+0.5,OC[2]+0.5]
                lower_bounds=[OC[0]-0.5,OC[1]-0.5,OC[2]-0.5]
                maxprob1[j]= multivariate_normal.cdf(upper_bounds, mean=mu, cov=cv1, allow_singular=True, lower_limit=lower_bounds)            
        
        if m>0:
            for x0 in range(4):
                pr=float( int(vv[1].split(",")[x0]) )/float(m)
                Q1[j][x0].append(pr)

src1.close()

for j in range(NoSNPs1):
    for x0 in range(4):
        rx = np.average(Q1[j][x0])
        Q1[j][x0]=rx # changed to float

print("no 40-allele markers",NoSNPs1)

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

        genotype_matrix[famid].append( [vv[0],ref,alt] )
        maxprob[j]=1.0
        if alt+ref>0:
            pr=float(ref)/float(alt+ref)
            if pr<1 and pr>0:
                Mreads=ref+alt
                mu=Mreads*pr
                var = Mreads*pr*(1.0-pr)
                z1=(ref+0.5-mu)/sqrt(var)
                z0=(ref-0.5-mu)/sqrt(var)
                maxprob[j] = norm.cdf(z1)-norm.cdf(z0)
            
            try:
                Q[j].append(pr)
            except KeyError:
                Q[j]=[pr]
src.close()


for j in range(NoSNPs):
    rx = np.average(Q[j])
    Q[j]=rx # changed to float

print("no snps",NoSNPs)

NoFamilies = len(Famsize.keys())


# fill required vectors ######################

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

glist=[]
for i in range(4):
    for j in range(i,4):
        glist.append(str(i)+str(j))


posg1={} # all possible sire genotype given the number of sires
posg1[1]=[]
for z in glist:
    posg1[1].append([z])

for j in range(2,ClutchBlocks+1):
    posg1[j]=[]
    for z in glist:
        for x in range(len(posg1[j-1])):
            new0=[]
            for y in range(j-1):
                new0.append(posg1[j-1][x][y])

            new0.append(z)
            posg1[j].append(new0)


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

# fit models ######################
famID = focal_fam

best={}
allbest={}
# One sire
delta[famID]=[1.0]
ll1=LL_fam(delta[famID])
best[1]=[ll1,delta[famID]]
print(famID,"one sire",best[1][0],best[1][1])
allbest=best[1]


for a1 in [0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95]:
# for a1 in [0.5]:
    delta[famID]=[a1,1-a1]
    llx=LL_fam(delta[famID])
    if a1==0.5:
        ll2=llx
        best[2]=[llx,delta[famID]]
    elif llx>best[2][0]:
        ll2=llx
        best[2]=[llx,delta[famID]]
    print("\t\t\ttesting",delta[famID],llx)

print(famID,"2sire",best[2][0],best[2][1],ll2-ll1)


if ll2-ll1>x2threshold:
    allbest=best[2]
    best[3]=[-999999999,"NA"]
    for a1 in sets3:
        delta[famID]=a1
        llx=LL_fam(delta[famID])
        if best[3][1]=="NA":
            ll3=llx
            best[3]=[llx,delta[famID]]
        elif llx>best[3][0]:
            ll3=llx
            best[3]=[llx,delta[famID]]
        print("\t\t\ttesting",delta[famID],llx)
    print(famID,"3sire",best[3][0],best[3][1],ll3-ll2)


    if ll3-ll2>x2threshold:
        allbest=best[3]
        best[4]=[-999999999,"NA"]
        for a1 in sets4:
            delta[famID]=a1
            llx=LL_fam(delta[famID])
            if best[4][1]=="NA":
                ll4=llx
                best[4]=[llx,delta[famID]]
            elif llx>best[4][0]:
                ll4=llx
                best[4]=[llx,delta[famID]]
            print("\t\t\ttesting",delta[famID],llx)
        print(famID,"4sire",best[4][0],best[4][1],ll4-ll3)
        if ll4-ll3>x2threshold:
            allbest=best[4]

out1.write(str(famID)+'\t'+str(Famsize[famID])+'\t'+str(allbest[0])+'\t'+str(len(allbest[1]))+'\t'+str(allbest[1])+'\n')
# llz = LL_fam_final_detailed(allbest[1])


out1.close()
# out2.close()



