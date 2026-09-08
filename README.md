# Polyandry_by_pooled_families


These programs apply likelihood methods to estimate polyandry from genomic data.




Set2: Make mini-haplotype markers
The program "make.snp.list.py" identifies intermediate frequency SNPs with sufficient coverage and quality from the overall vcf file.  It produces [prefix]+".snplist" for use by subsequent programs.
The program "Determine.close.snps.py" identifies SNPs separated by less than a specified physical distance (e.g. MaxDist = 50) to be evaluated as potential components of mini-haplotype markers.




Set 3: Identify mis-matched families
The program AYMM.v1.py tests a clutch against a collection of alternative mothers to determine if one the alternatives is the true genetic mother of that clutch.
The "dubious_fams.txt" input is determine the application of xxx to the genomic data with the original mother-clutch assignments.


The program run.aymm.py creates shell scripts to run AYMM.v1.py on each suspect mother-clutch assignment identified in "dubious_fams.txt"


Set 4: Estimate polyandry given genotype files
The program ML4_a2_a4_onefam.py estimates polyandry in specified maternal family given two genotype input files (FILEPREFIX+".2data.txt", FILEPREFIX1+".4data.txt").
The program ML4A2_onefam.py performs the same operation but only uses biallelic snps (FILEPREFIX+".2data.txt").
The program ML4A4_onefam.py performs the same operation but only uses 4-allele markers (FILEPREFIX+".4data.txt").


