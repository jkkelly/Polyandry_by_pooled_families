# -------------------------------------------------------------------------------
# Name:
# Purpose: Identify SNPs from *snplist that are within MaxDist
# -------------------------------------------------------------------------------
import sys
infile = sys.argv[1] # chr2L.thin2.vcf
chrom=infile.split(".")[0]
MaxDist = 50

src = open(infile+".snplist", "r") # should be ordered vcf
out1 = open(infile+".close.snps."+str(MaxDist), "w")


biglist=[] # one chrom 
bases=[]
for line_idx, line in enumerate(src):
    cols = line.replace('\n', '').split('\t')
# chr2L   7088    A       T       0.8090659340659341      244     101     19
# chr2L   7156    A       G       0.8787465940054496      282     81      4
    biglist.append(int(cols[1]))
    bases.append([cols[2],cols[3]])

src.close()

outcomes={}
for j in range(len(biglist)-1):
    pos=biglist[j]
    distforward=0
    k=j+1
    near=[]
    while distforward<=MaxDist:
        distforward=biglist[k]-pos
        if distforward<=MaxDist:
            near.append(biglist[k])
        k+=1
        if k>=len(biglist):
            break

    nx=len(near)
    try:
        outcomes[nx]+=1
    except KeyError:
        outcomes[nx]=1
    if nx>1:
        out1.write( chrom+'\t'+str(pos)+'\t'+bases[j][0]+'\t'+bases[j][1]+'\t'+str(nx) )
        for j in range(nx):
            out1.write('\t'+str(near[j]))
        out1.write('\n')

out1.close()

for nx in outcomes:
    print(nx,outcomes[nx])
    
