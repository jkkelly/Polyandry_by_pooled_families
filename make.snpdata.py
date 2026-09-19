import sys
import numpy as np

snpfile = sys.argv[1]

minmd=6 # minimum depth to call genotype

out1 = open(snpfile+".snpdata.txt", "w")
out2 = open(snpfile+".snpdata.famgenos", "w")

fams={}
staxt={}
src = open("auto.famlist.txt", "r")
for line_idx, line in enumerate(src):
    cols = line.replace('\n', '').split('\t')
    fams[cols[0]]=cols[0]
    staxt[cols[0]]=[0,0,0]
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

snps={}
src = open(snpfile,"r")
for line_idx, line in enumerate(src):
    cols = line.replace('\n', '').split('\t')
    # chr2L   40556
    snps[ cols[0]+"_"+cols[1] ]=1
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
                Families[cols[5]][0]=int(cols[0])
            elif cols[2]=="clutch":
                Families[cols[5]][1]=int(cols[0])

        except KeyError:
            pass

in1.close()

dous=0

for infx in ["chr2L.thin2.vcf","chr2R.thin2.vcf","chr3L.thin2.vcf","chr3R.thin2.vcf"]:
    src = open(infx, "r")
    for line_idx, line in enumerate(src):
        cols = line.replace('\n', '').split('\t')
    # chr2L   7088    .       A       T       77706.7 .       DP=35863;VDB=1.80735e-10;SGB=7519.4;RPBZ=1.41804;MQBZ=0.00498717;MQSBZ=-0.00765249;BQBZ=-1.76174;SCBZ=-0.521751;MQ0F=2.78839e-05;AC=729;AN=4000;DP4=9393,15326,2350,4022;MQ=59        GT:PL:AD        ./.:0,0,0:0,0
        try:
            uk=snps[cols[0]+"_"+cols[1]]
            dous+=1
            for fam in fams:
                j=Families[fam][0] # mom
                MomG="NA"
                if j>0:
                    vv = cols[j].split(":")  # 0/0:0,15,151:5,0
                    if len(vv) == 3:
                        if vv[0] != './.':
                            depth_mom = int(vv[2].split(",")[0]) + int(vv[2].split(",")[1])
                            if depth_mom >= minmd and vv[0]=="0/0":
                                MomG="RR"
                                staxt[fam][0]+=1
                            elif depth_mom >= minmd and vv[0]=="1/1":
                                MomG="AA"
                                staxt[fam][2]+=1
                            elif vv[0]=="0/1":
                                MomG="RA"
                                staxt[fam][1]+=1
                j=Families[fam][1] # clutch
                RA="0,0"
                if j>0:
                    vv = cols[j].split(":")  # 0/0:0,15,151:5,0
                    if len(vv) == 3:
                        if vv[0] != './.':
                            RA = vv[2]

                fams[fam]+= "\t"+MomG+":"+RA

        except KeyError:
            pass
    src.close()


print("found snps",dous)

for fam in fams:
    out1.write(fams[fam]+'\n')
    out2.write(fam+'\t'+str(staxt[fam][0])+'\t'+str(staxt[fam][1])+'\t'+str(staxt[fam][2])+'\n')
