import sys
import subprocess
import os.path
import time
from os import mkdir, chdir, getcwd
from shutil import copyfile, rmtree

# creates subprocess that checks if sysbuilder is done and starts conversion and minimization when it is
##------------------------------------------------------------------------------------------------------------
# set variables
try:
    location, jobName, trialName, path = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
except:
    location, jobName, trialName, path = sys.argv[1], sys.argv[2], None, sys.argv[3]

bashCommand = "$SCHRODINGER/utilities/multisim -JOBNAME " + jobName + "_min" + " -HOST localhost -maxjob 1 -cpu 1 " \
                    "-m mr_minimization.msj -description " + "'aacg min job' " + jobName + "_setup-out.cms"\
                    ' -set stage[1].set_family.md.jlaunch_opt=["-gpu"] -o ' + jobName + "_min-out.cms"

if trialName == None:
    dir = path + jobName +"/"
    cmsLocation = path + "/" + jobName + "/" + jobName + "_setup-out.cms"
else:
    dir = path + trialName + "/" +jobName
    cmsLocation = path + trialName + "/" + jobName + "/" + jobName + "_setup-out.cms"

# wait for sys builder to finish
print(cmsLocation)
while not os.path.exists(cmsLocation):
    time.sleep(30)

# start conversion/minimization job
copyfile(path + "/mr_minimization.msj", dir + "/mr_minimization.msj")
subprocess.Popen(bashCommand, shell=True, cwd=path+"/" + jobName)
