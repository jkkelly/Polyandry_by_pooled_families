# -*- coding: utf-8 -*-

import sys

chrom = sys.argv[1] # "chr2L"

out1 = open(chrom+".twoloc.evaluation", "w")
out2 = open(chrom+".twoloc.informative_by_fam", "w")
MinHomo = 8
MinHet = 2

fams={}
src = open("auto.famlist.txt", "r")
for line_idx, line in enumerate(src):
    cols = line.replace('\n', '').split('\t')
    fams[cols[0]]=0 # 4_106
src.close()

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

stats={}
cc=0
for fm in Families:
    mom_name = Families[fm][0]
    src = open(chrom+"files/"+mom_name+"."+chrom+".twolocus.aligned", "r") # Sp25_4_18_S263_L002.chr2L.twolocus.aligned
    for line_idx, line in enumerate(src):
        cols = line.replace('\n', '').split('\t')
        # Sp25_4_18_S263_L002     chr2L   25708_25727     0/1_11,9        0/1_12,10       11      0       0       9
        # Sp25_4_18_S263_L002     chr2L   25727_25729     0/1_12,10       0/1_12,10       12      0       0       10
        snppair = cols[2]
        if cc==0:
            stats[snppair]=[0 for j in range(10)]

        g1=cols[3].split("_")[0]
        g2=cols[4].split("_")[0]
        if g1 != "./." and g2 != "./.":
            gg = g1+","+g2

        nonz=[0 for j in range(4)]
        acounts=[0 for j in range(4)]
        for j in range(4):
            if int(cols[5+j])>0:
                nonz[j]=1
                acounts[j]=int(cols[5+j])
        if sum(nonz)==0:
            stats[snppair][4]+=1
        elif sum(nonz)==1:
            if max(acounts)>=MinHomo:
                fams[fm]+=1    

            stats[snppair][5]+=1
            for j in range(4):
                if nonz[j]==1:
                    cu88=j
                    stats[snppair][j]+=2
            if (cu88==0 and gg=="0/0,0/0") or (cu88==2 and gg=="1/1,0/0") or (cu88==1 and gg=="0/0,1/1") or (cu88==3 and gg=="1/1,1/1"):
                stats[snppair][8]+=1
            else:
                stats[snppair][9]+=1
        elif sum(nonz)==2:
            stats[snppair][6]+=1
            ct1=0
            for j in range(4):
                if acounts[j]>=MinHet:
                    ct1+=1                    
            if ct1==2:
                fams[fm]+=1        
                
            cu88=''
            for j in range(4):
                if nonz[j]==1:
                    cu88+=str(j)
                    stats[snppair][j]+=1
            if (cu88=="01" and gg=="0/0,0/1") or (cu88=="02" and gg=="0/1,0/0") or (cu88=="03" and gg=="0/1,0/1") or (cu88=="12" and gg=="0/1,0/1") or (cu88=="13" and gg=="0/1,1/1") or (cu88=="23" and gg=="1/1,0/1") :
                stats[snppair][8]+=1
            else:
                stats[snppair][9]+=1
        else:
            stats[snppair][7]+=1

    src.close()
    cc+=1
    print("done",cc)

for fm in fams:
    out2.write(chrom+'\t'+fm+'\t'+str(fams[fm])+'\n')
    
for snppair in stats:
    out1.write(chrom+'\t'+snppair)
    for j in range(10):
        out1.write('\t'+str(stats[snppair][j]))
    out1.write('\n')
