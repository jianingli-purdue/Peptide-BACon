import sys
import subprocess
import os.path
import time
import random
from os import mkdir, chdir, getcwd
from shutil import copyfile, rmtree

##------------------------------------------------------------------------------------------------------------
## class to run md simulation on built system
class mdJob:
    # constructor
    def __init__(self, jobName, trialName, path):
        self._path = path
        self._jobName = jobName
        self._trialName = trialName
        self._bashCommand = "$SCHRODINGER2/utilities/multisim -VIEWNAME desmond_molecular_dynamics_gui.MDApp -JOBNAME " + \
                            self._jobName + " -HOST localhost -maxjob 1 -cpu 1 -m " + self._jobName + ".msj -c " + self._jobName + \
                            '.cfg -description "Molecular Dynamics" ' + self._jobName + '.cms -mode umbrella -set stage[1].set_family.md.jlaunch_opt=[\"-gpu\"] ' + \
                            "-PROJ /home/marlo/.schrodinger/tmp/tproj62725a4923 -DISP append -o " + self._jobName + "-out.cms -ATTACHED"

        # generate new folder
        if trialName == None:
            self._dir = self._path +  self._jobName + "/"
            self._location = self._path + self._jobName + "/" + self._jobName + "_setup-out.cms"
        else:
            self._dir = self._path + self._trialName + "/" + self._jobName + "/"
            self._location = self._path + self._trialName + "/" + self._jobName + "/" + self._jobName + "_setup-out.cms"
        self._dst = self._dir + self._jobName + ".cms"
        # try:
        #     mkdir(self._dir + "/")
        # except:
        #     rmtree(self._dir + "/")
        #     mkdir(self._dir + "/")
##------------------------------------------------------------------------------------------------------------
# method to move out.cms to job folder
    def moveCMS(self):
        # print(self._location)
        while not os.path.exists(self._location):
            time.sleep(10)
        copyfile(self._location, self._dst)
##------------------------------------------------------------------------------------------------------------
# method to build msj
    def writeMSJ(self):
        msg = '''
# Desmond standard NPT relaxation protocol
# All times are in the unit of ps.
# Energy is in the unit of kcal/mol.
task {
   task = "desmond:auto"
   set_family = {
      desmond = {
         checkpt.write_last_step = no
      }
   }
}

simulate {
   title       = "Brownian Dynamics NVT, T = 10 K, small timesteps, and restraints on solute heavy atoms, 100ps"
   annealing   = off
   time        = 100
   timestep    = [0.001 0.001 0.003 ]
   temperature = 10.0
   ensemble = {
      class = "NVT"
      method = "Brownie"
      brownie = {
         delta_max = 0.1
      }
   }
   restrain = {
      atom = "solute_heavy_atom"
      force_constant = 50.0
   }
}

simulate {
   effect_if   = [["==" "-gpu" "@*.*.jlaunch_opt[-1]"] 'ensemble.method = Langevin']
   title       = "NVT, T = 10 K, small timesteps, and restraints on solute heavy atoms, 12ps"
   annealing   = off
   time        = 12
   timestep    = [0.001 0.001 0.003]
   temperature = 10.0
   restrain    = { atom = solute_heavy_atom force_constant = 50.0 }
   ensemble    = {
      class  = NVT
      method = Berendsen
      thermostat.tau = 0.1
   }

   randomize_velocity.interval = 1.0
   eneseq.interval             = 0.3
   trajectory.center           = []
}

simulate {
   title       = "NPT, T = 10 K, and restraints on solute heavy atoms, 12ps"
   effect_if   = [["==" "-gpu" "@*.*.jlaunch_opt[-1]"] 'ensemble.method = Langevin']
   annealing   = off
   time        = 12
   temperature = 10.0
   restrain    = retain
   ensemble    = {
      class  = NPT
      method = Berendsen
      thermostat.tau = 0.1
      barostat  .tau = 50.0
   }

   randomize_velocity.interval = 1.0
   eneseq.interval             = 0.3
   trajectory.center           = []
}

solvate_pocket {
   should_skip = true
   ligand_file = ?
}

simulate {
   title       = "NPT and restraints on solute heavy atoms, 12ps"
   effect_if   = [["@*.*.annealing"] 'annealing = off temperature = "@*.*.temperature[0][0]"'
                  ["==" "-gpu" "@*.*.jlaunch_opt[-1]"] 'ensemble.method = Langevin']
   time        = 12
   restrain    = retain
   ensemble    = {
      class  = NPT
      method = Berendsen
      thermostat.tau = 0.1
      barostat  .tau = 50.0
   }

   randomize_velocity.interval = 1.0
   eneseq.interval             = 0.3
   trajectory.center           = []
}

simulate {
   title       = "NPT and no restraints, 24ps"
   effect_if   = [["@*.*.annealing"] 'annealing = off temperature = "@*.*.temperature[0][0]"'
                  ["==" "-gpu" "@*.*.jlaunch_opt[-1]"] 'ensemble.method = Langevin']
   time        = 24
   ensemble    = {
      class  = NPT
      method = Berendsen
      thermostat.tau = 0.1
      barostat  .tau = 2.0
   }

   eneseq.interval   = 0.3
   trajectory.center = solute
}

simulate {
   cfg_file = "'''
        msg += self._jobName
        msg += '''.cfg"
   jobname  = "$MASTERJOBNAME"
   dir      = "."
   compress = ""
}

# Job launching command:
# $SCHRODINGER/utilities/multisim -VIEWNAME desmond_molecular_dynamics_gui.MDApp -JOBNAME g_rec_prism_1 -HOST localhost -maxjob 1 -cpu 1 -m g_rec_prism_1.msj -c g_rec_prism_1.cfg -description "Molecular Dynamics" g_rec_prism_1.cms -mode umbrella -set stage[1].set_family.md.jlaunch_opt=[\"-gpu\"] -PROJ /home/marlo/.schrodinger/tmp/tproj62725a4923 -DISP append -o g_rec_prism_1-out.cms -ATTACHED
        '''

        msj = open(self._path + self._trialName + "/" + self._jobName + "/" + self._jobName + ".msj", "w")
        msj.write(msg)
        msj.close()

##------------------------------------------------------------------------------------------------------------
# method to build cfg
# accepts dictionary of settings
    def writeCFG(self, job_settings):
        seed = job_settings["seed"]
        time = job_settings["time"]
        interval = job_settings["interval"]

        msg = '''annealing = false
backend = {
}
bigger_rclone = false
checkpt = {
   first = 0.0
   interval = 240.06
   name = "$JOBNAME.cpt"
   write_last_step = true
}
cpu = 1
cutoff_radius = 9.0
elapsed_time = 0.0
energy_group = false
eneseq = {
   first = 0.0
   interval = 1.2
   name = "$JOBNAME$[_replica$REPLICA$].ene"
}
ensemble = {
   barostat = {
      tau = 2.0
   }
   class = NPT
   method = MTK
   thermostat = {
      tau = 1.0
   }
}
glue = solute
maeff_output = {
   first = 0.0
   interval = 120.0
   name = "$JOBNAME$[_replica$REPLICA$]-out.cms"
   periodicfix = true
   trjdir = "$JOBNAME$[_replica$REPLICA$]_trj"
}
meta = false
meta_file = ?
pressure = [1.01325 isotropic ]
randomize_velocity = {
   first = 0.0
   interval = inf
   seed = '''
        msg += str(random.randrange(1000, 9999)) + '\n'
        msg += '''   temperature = "@*.temperature"
}
restrain = none
simbox = {
   first = 0.0
   interval = 1.2
   name = "$JOBNAME$[_replica$REPLICA$]_simbox.dat"
}
surface_tension = 0.0
taper = false
temperature = [
   [300.0 0 ]
]
time = '''
        msg += time + '\n'
        msg += '''timestep = [0.002 0.002 0.006 ]
trajectory = {
   center = []
   first = 0.0
   format = dtr
   frames_per_file = 250
   interval = '''
        msg += interval + '\n'
        msg += '''   name = "$JOBNAME$[_replica$REPLICA$]_trj"
   periodicfix = true
   write_velocity = false
}
        '''
        cfg = open(self._path + self._trialName + "/" + self._jobName + "/" + self._jobName + ".cfg", "w")
        cfg.write(msg)
        cfg.close()

##------------------------------------------------------------------------------------------------------------
# method to run job
    def run(self):
        # print("/home/marlo/HDD/Maestro/jobs/" + self._jobName)
        subprocess.run(self._bashCommand, shell=True, cwd=self._path + self._trialName + "/" + self._jobName + "/")
