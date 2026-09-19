# -------------------------------------------------------------------------------
# Name:
# Purpose:     Identify best SNP for inference within each genomic window
# -------------------------------------------------------------------------------

import sys
import numpy as np

chrom= sys.argv[1]

windowsize = 10000
HWcut = 3.84
MinMAF=0.25

Fac1=5.0
Fac2=20.0

src = open(chrom+"thin2summary.txt" , "r")
out1 = open(chrom+".bestper_window", "w")
best={}
for line_idx, line in enumerate(src):
    cols = line.replace('\n', '').split('\t')

# Chrom   pos     M1a     M1b     M2a     M2b     M3a     M3b     M4a     M4b     RR      RA      AA      NN      clutch_R        clutch_A        MomOff_test     MomOff_fail     p(mom)  p(clutch)       HW_chi2
# chr2L   7088    205,48  652,147 285,91  300,91  476,133 212,48  248,62  658,139 164     76      14      159     18814   4884    149     2       0.7952755905511811      0.7939066587897713      1.6711040047703891

    if cols[0] == "Chrom":
        pass  # out1.write(line)
    else:
        pmom=float(cols[18])
        pclutch=float(cols[19])
        HW_chi2=float(cols[20])
        if min(pmom,1-pmom)>=MinMAF and min(pclutch,1-pclutch)>=MinMAF and HW_chi2<HWcut:
            win = int(int(cols[1])/windowsize)
            try:
                uk=best[win]
            except KeyError:
                best[win]=[100,"NA",0]
            best[win][2]+=1
            fracmissing = float(cols[13])/( float(cols[10])+float(cols[11])+float(cols[12])+float(cols[13]) )
            fracbad = float(cols[17])/( float(cols[16])+float(cols[17]) )
            mscore = Fac2*abs(pmom-pclutch) + Fac1*fracbad + fracmissing
            if mscore < best[win][0]:
                best[win][0]=mscore
                best[win][1]=line
src.close()

for win in best:
    out1.write(chrom+'\t'+str(win*windowsize)+'\t'+str(best[win][2])+'\t'+str(best[win][0])+'\t'+str(best[win][1]))
