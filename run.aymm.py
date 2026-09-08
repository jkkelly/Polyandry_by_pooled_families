#!/usr/bin/env python3
# -*- coding: utf-8 -*-

commandsFile = open("a_all.sh", 'w')

src=open("dubious_fams.txt","r")
for line_idx, line in enumerate(src):
    cols = line.replace('\n', '').split('\t')
    sampleName = cols[0]

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
    shellScript.write("#SBATCH --output=LOGS/ra_%s.log\n\n\n" % sampleName )
    shellScript.write("module load conda\n")
    shellScript.write("conda activate python3sk\n")

    shellScript.write("python AYMM.v1.py auto %s\n" % (sampleName))
