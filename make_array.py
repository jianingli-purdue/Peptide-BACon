import random
from math import exp
from scipy.constants import Boltzmann
import box
import cgSim
import dev
import sys
from schrodinger import structure

##------------------------------------------------------------------------------------------------------------
# script to create models for use in 100ns runs
path = "/home/marlo/HDD/desmond_jobs/jobs/probacon/"

seqs = ["ECG", "KFG", "GG", "FF"]
# seqs = ["FE"]
num = 8
salt_conc = 10/1000

for seq in seqs:
    name = seq + "_" + str(int(round(num)))+ "_" + str(int(salt_conc*1000)) + "mM_AACG"
    height = num
    total = num**3
    print(name, height, total)
    b = box.box(name=name, spacing=15, height=height, total=total, rand=False, aaCG=True, path=path)


    b.addStruct(seq)
    b.build()
    b.sysBuild()
    # b.run()
    # b.analyze()
