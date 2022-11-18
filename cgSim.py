import sys
import subprocess
import os.path
import time
import fileinput
from random import randint
from shutil import copyfile, rmtree

##------------------------------------------------------------------------------------------------------------
## class to run md simulation on built system
class cgJob:
    # constructor
    def __init__(self, jobName, trialName, path):
        self._path = path
        self._jobName = jobName
        self._trialName = trialName
        self._bashCommand = "$SCHRODINGER/utilities/multisim -JOBNAME " + self._jobName + "_aacg" + " -HOST localhost -maxjob 1 -cpu 1 " \
                            "-m mr_md.msj -description " + "'aacg job' " + self._jobName + "_min-out.cms"\
                            ' -set stage[1].set_family.md.jlaunch_opt=["-gpu"] -o ' + self._jobName + "_aacg-out.cms"

        if trialName == None:
            self._dir = self._path + self._jobName
            self._cmsLocation = self._path + self._jobName + "/" + self._jobName + "_min-out.cms"
        else:
            self._dir = self._path + self._trialName + "/" +self._jobName
            self._cmsLocation = self._path + self._trialName + "/" + self._jobName + "/" + self._jobName + "_min-out.cms"

##------------------------------------------------------------------------------------------------------------
# method to move coarse grain msj to job folder
    def moveMSJ(self):
        # wait for system builder to finish
        print(self._cmsLocation)
        while not os.path.exists(self._cmsLocation):
            time.sleep(10)
        #time.sleep(10)
        # copyfile("/home/marlo/HDD/Maestro/jobs/montecarlo/mr_minimization.msj", self._dir + "mr_minimization.msj")
        # copyfile("/home/marlo/HDD/Maestro/jobs/montecarlo/mr_md.msj", self._dir + "mr_md.msj")

        copyfile(self._path + "mr_md.msj", self._dir + "/mr_md.msj")

##------------------------------------------------------------------------------------------------------------
# method to randomize velocity seed
    def randomizeSeed(self):
        filename = self._dir + "/mr_md.msj"
        text_to_search = "2007"
        replacement_text = str(randint(0000, 9999)).rjust(4, "0")
        #print(replacement_text)
        with fileinput.FileInput(filename, inplace=True, backup='.bak') as file:
            for line in file:
                print(line.replace(text_to_search, replacement_text), end='')

##------------------------------------------------------------------------------------------------------------
# method to run job
    def run(self):
        # run command
        subprocess.Popen(self._bashCommand, shell=True, cwd=self._path + self._trialName + "/" +self._jobName)

# job = cgJob("ACR_step22_attempt1", "trip_trial_1", "/home/marlo/HDD/Maestro/jobs/montecarlo/")
# job.randomizeSeed()
