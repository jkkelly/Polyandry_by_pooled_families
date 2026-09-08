# -*- coding: utf-8 -*-

import sys
#import numpy as np

MinMinorHap = 10
Max3plus = 0
MaxMissing = 20
MaxFracWrong =0.05
OE = [0.9,1.1]
nm = sys.argv[1] 


out1 = open(nm+".4data.txt", "w")

fams={}
src = open("auto.famlist.txt", "r")
for line_idx, line in enumerate(src):
    cols = line.replace('\n', '').split('\t')
    fams[cols[0]]=cols[0]
# 4_106

src.close()

src = open("SP25_famsizes.txt", "r")
for line_idx, line in enumerate(src):
    cols = line.replace('\n', '').split('\t')
    if line_idx>0:
        try:
            uz=fams[cols[0]]
            fams[cols[0]]+= "\t"+cols[1]
    # 1_40      72
        except KeyError:
            pass
src.close()

Families={}
in1 = open("SP25.vcfkey", "r")
for line_idx, line in enumerate(in1):
    cols = line.replace('\n', '').split('\t')
# vcf column    name    type    collection      Family  Family_name
# 515   Sp25_3_1_S238_L002      Mother  3       1       3_1
# 514   Sp25_3_1P_merged        clutch  3       1       3_1
# 2046  Sp25_M3a_S755_L007_part_0531    Male_pool       M3a     Male    M3a_Male
# 2047  Sp25_M3b_S756_L007      Male_pool       M3b     Male    M3b_Male

    if line_idx > 0:
        try:
            amirelevant = fams[cols[5]]
            try:
                uk=Families[cols[5]]
            except KeyError:
                Families[cols[5]]=[-9,-9]
            if cols[2]=="Mother":
                Families[cols[5]][0]=cols[1]
            elif cols[2]=="clutch":
                Families[cols[5]][1]=cols[1]
        except KeyError:
            pass
in1.close()


keep={}
for chrom in ["chr2L","chr2R","chr3L","chr3R"]:
    print("reading",chrom)
    src = open(chrom+".twoloc.evaluation", "r")
    for line_idx, line in enumerate(src):
        cols = line.replace('\n', '').split('\t')
        snppair=cols[0]+"_"+cols[1]
    # chr2L   25708_25727     340     0       0       216     1       68      210     3       260     18
    
        nA=[int(cols[2]),int(cols[3]),int(cols[4]),int(cols[5])]
        apf=[int(cols[6]),int(cols[7]),int(cols[8]),int(cols[9])] # alleles per mom  0      1       2       3+
    
    
        if min(nA)>= MinMinorHap and apf[0]<=MaxMissing and apf[3]<=Max3plus:
            mm=[int(cols[10]),int(cols[11])]
            fracwrong = float(mm[1])/sum(mm)
            px=[]
            nx=float(sum(nA))
            ehomo=0.0
            for j in range(4):
                px.append(float(nA[j])/nx)
                ehomo+= px[j]**2.0
            ehet=1.0-ehomo
            ohet=float(apf[2])/float(apf[1]+apf[2])
            if ohet/ehet >= OE[0] and ohet/ehet <= OE[1] and fracwrong<=MaxFracWrong:
                keep[snppair]=1
    src.close()

print("NoSNPs",len(keep.keys()))


cc=0
for fm in Families:
    
    mom_name = Families[fm][0]
    datums={}
    for chrom in ["chr2L","chr2R","chr3L","chr3R"]:
    
        src = open(chrom+"files/"+mom_name+"."+chrom+".twolocus.aligned", "r") # Sp25_4_18_S263_L002.chr2L.twolocus.aligned
        for line_idx, line in enumerate(src):
            cols = line.replace('\n', '').split('\t')
            # Sp25_4_18_S263_L002     chr2L   25708_25727     0/1_11,9        0/1_12,10       11      0       0       9
            # Sp25_4_18_S263_L002     chr2L   25727_25729     0/1_12,10       0/1_12,10       12      0       0       10
            try:
                ku88 = keep[cols[1]+"_"+cols[2]]
                snppair= cols[1]+"_"+cols[2]
                nonz=[0 for j in range(4)]
                for j in range(4):
                    if int(cols[5+j])>0:
                        nonz[j]=1
                if sum(nonz)==0:
                    MG="NA"
                elif sum(nonz)==1:
                    for j in range(4):
                        if nonz[j]==1:
                            cu88=j
                    if cu88==0:
                        MG="00"
                    elif cu88==1:
                        MG="11"
                    elif cu88==2:
                        MG="22"
                    else:
                        MG="33"
                elif sum(nonz)==2:
                    cu88=''
                    for j in range(4):
                        if nonz[j]==1:
                            cu88+=str(j)
                    MG=cu88
                else:
                    MG="NA"
                datums[snppair]=MG+":"

            except KeyError:
                pass
        src.close()
    
        c_name = Families[fm][1]
        src = open(chrom+"files/"+c_name+"."+chrom+".twolocus.aligned", "r") # Sp25_4_18_S263_L002.chr2L.twolocus.aligned
        for line_idx, line in enumerate(src):
            cols = line.replace('\n', '').split('\t')
            # Sp25_4_18_S263_L002     chr2L   25708_25727     0/1_11,9        0/1_12,10       11      0       0       9
            # Sp25_4_18_S263_L002     chr2L   25727_25729     0/1_12,10       0/1_12,10       12      0       0       10
            try:
                ku88= keep[cols[1]+"_"+cols[2]]
                snppair= cols[1]+"_"+cols[2]
                cdats=cols[5]
                for j in range(1,4):
                    cdats+=','+cols[5+j]
                datums[snppair]+=cdats                     
            except KeyError:
                pass
        src.close()
    cc+=1
    print("done",fm,cc)

    out1.write(fams[fm])
    for snppair in datums:    
        out1.write("\t"+datums[snppair])
    out1.write('\n')

out1.close()
