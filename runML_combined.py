#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

commandsFile = open("ML4.2_4_all.sh", 'w')

spx="chosen.snps1.txt"
splist = "set0"

src=open("Fams_to_estimate.txt","r")
for line_idx, line in enumerate(src):
    cols = line.replace('\n', '').split('\t')
    sampleName = cols[0]
# 1_101
    shellScriptName = 'a%s.sh' % (sampleName)
    shellScript = open(shellScriptName, 'w' )

    commandsFile.write('sbatch %s \n' % (shellScriptName))

    shellScript.write("#!/bin/bash\n" )
    shellScript.write("#SBATCH --job-name=%s\n" % (sampleName))

    shellScript.write("#SBATCH --partition=sjmac       \n" )
    shellScript.write("#SBATCH --time=0-25:59:00       \n" )
    shellScript.write("#SBATCH --mem-per-cpu=5gb      \n" )
    shellScript.write("#SBATCH --mail-type=NONE          \n" )
    shellScript.write("#SBATCH --mail-user=jkk@ku.edu    \n" )
    shellScript.write("#SBATCH --ntasks=1                   \n" )
    shellScript.write("#SBATCH --cpus-per-task=1            \n" )
    shellScript.write("#SBATCH --output=LOGS/rx_%s.log\n\n\n" % sampleName )
    shellScript.write("module load conda\n")
    shellScript.write("conda activate python3sk\n")
    shellScript.write("python ML4_a2_a4_onefam.py %s %s %s\n" % (spx,splist,sampleName))
