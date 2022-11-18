import sys
import subprocess
from os import mkdir, chdir, getcwd
from shutil import copyfile, rmtree

## sys builder class
# allows user to run sysbuilder on generated model
##------------------------------------------------------------------------------------------------------------ v
class sysBuilder:
    def __init__(self, location, jobName, trialName, path):
        self._path = path
        self._location = location
        self._trialName = trialName
        self._jobName = jobName
        self._bashCommand = '$SCHRODINGER/utilities/multisim -JOBNAME ' + self._jobName + "_setup" + ' -m ' \
                            + self._jobName + ".msj " + self._jobName + '.mae -o ' + self._jobName + "_setup" + \
                            "-out.cms -HOST localhost"
        # print(self._bashCommand)

        # used to convert to AACG
        self._minCommand = "$SCHRODINGER/run sysMin.py " + self._location + " " +  self._jobName + " " + self._trialName + " " + self._path
        #print(self._minCommand)

        # generate new folder
        if trialName == None:
            self._dst = self._path + "/" + self._jobName + "/" + self._jobName + ".mae"
            try:
                mkdir(self._path + "/" + self._jobName + "/")
            except:
                rmtree(self._path + "/" + self._jobName + "/")
                mkdir(self._path + "/" + self._jobName + "/")

        else:
            self._dst = self._path + trialName + "/" + self._jobName + '/' + self._jobName + ".mae"
            try:
                mkdir(self._path + trialName + "/" + self._jobName + "/")
            except:
                rmtree(self._path + trialName + "/" + self._jobName + "/")
                mkdir(self._path + trialName + "/" + self._jobName + "/")

##------------------------------------------------------------------------------------------------------------
# method to move mae to job folder
    def moveMAE(self):
        # print(self._location, self._dst)
        copyfile(self._location, self._dst)

##------------------------------------------------------------------------------------------------------------
# method to write .msj
    def writeMSJ(self):
        # write msj to job folder
        msg = '''task {
	task = "desmond:auto"
}

build_geometry {
	add_counterion = {
		 ion = Na
		 number = neutralize_system
	}
	box = {
		 shape = cubic
		 size = 7.5
		 size_type = buffer
	}
	override_forcefield = OPLS3e
        minimize_volume = true
	rezero_system = false
     salt = {
        concentration = 0.01
        negative_ion = Cl
        positive_ion = Na
     }
    solvent = SPC
}

assign_forcefield {
	forcefield = OPLS3e
}'''
        msj = open(self._path + self._trialName + "/" + self._jobName + "/" + self._jobName + ".msj", "w")
        msj.write(msg)
        msj.close()
##------------------------------------------------------------------------------------------------------------
# method to run system builder
    def run(self):
        # os.system("some_command &")
        # cmmd_file = open(self._path + self._trialName + "/" + self._jobName + "/setupCmd.sh", "w")
        # cmmd_file.write(self._bashCommand)
        # cmmd_file.close()
        #
        # subprocess.run("chmod +x setupCmd.sh", shell=True, cwd=self._path + self._trialName + "/" + self._jobName)

        subprocess.run(self._bashCommand, shell=True, cwd=self._path + self._trialName + "/" + self._jobName)
        subprocess.Popen(self._minCommand, shell=True)
