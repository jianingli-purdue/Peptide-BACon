import random
import sys, numpy
import dev
import time
import pickle
import os.path
from sysBuilder import sysBuilder
from mdSim import mdJob
from cgSim import cgJob
from math import tan
from grid import grid
from shutil import rmtree
from schrodinger.structutils import transform, analyze, minimize
from schrodinger.application.desmond.packages import analysis, traj, topo, traj_util, msys
from schrodinger import structure
from schrodinger.structutils.build import add_hydrogens, connect, delete_hydrogens

## Box class accepts user input and builds populated grid

class box:
    def __init__(self, name=None, spacing=None, height=None, width=None, length=None, total=None, rand=False, concentrations=None, blocked=False, aaCG=False, trialName="", path=None):
        self._name = name
        self._spacing = spacing
        self._height = height
        self._width = width
        self._length = length
        self._total = total
        self._rand = rand
        self._blocked = blocked
        self._structures = []
        if concentrations != None:
            self._concentrations = concentrations
        else:
            self._concentrations = []
        self._aaCG = aaCG
        self._trialName = trialName
        self._path = path
##------------------------------------------------------------------------------------------------------------
## method to add beta structure without file reader
    def addStruct(self, toAdd, shape=None, list=False):
        if list == True:
            for struct in toAdd:
                self._structures.append(dev.build_helix(struct))
        else:
            try:
                if shape == None:
                    st = dev.build_sheet(toAdd)
                else:
                    st = dev.build_helix(toAdd)
                self._structures.append(st)
            except:
                self._structures.append(next(structure.StructureReader("input/" + toAdd)))

##------------------------------------------------------------------------------------------------------------
## method to read input file and extract settings/sequence
    def fileReader(self, inputfile):
        inputfile = 'input/' + inputfile
        fh = open(inputfile)
        flag = False
        for line in fh:
            temp = line.split()
            if temp[0] == "":
                continue
            if temp[0] == "#":
                if temp[1] == "!Name":
                    try:
                        self._name = temp[2]
                    except:
                        print("You didn't include a name!")
                if temp[1] == "!source":
                    box._source = temp[2]
                if temp[1] == "!conc":
                    for i in temp[2:]:
                        self._concentrations.append(i)
            if temp[0] == "*":
                for param in temp:
                    temp2 = param.split("=")
                    if "total" in temp2:
                        self._total = int(temp2[1])
                    if "height" in temp2:
                        self._height = int(temp2[1])
                    if "width" in temp2:
                        self._width = int(temp2[1])
                    if "length" in temp2:
                        self._length = int(temp2[1])
                    if "spacing" in temp2:
                        self._spacing = int(temp2[1])
                    if "rand" in temp2:
                        if temp2[1] == "True":
                            self._rand = True
                        else:
                            self._rand = False
                    if "blocked" in temp2:
                        if temp2[1] == "True":
                            self._blocked = True
                        else:
                            self._blocked = False
            if temp[0] == ">":
                flag = True
                continue
            if temp[0] == "<":
                flag = False
            if flag:
                # build structure to append
                try:
                    #print(temp[0][0])
                    #print(temp[0])
                    if temp[1] == "alpha":
                        st = dev.build_helix(temp[0])
                        self._structures.append(st)
                    elif temp[1] == "beta":
                        st = dev.build_sheet(temp[0])
                        self._structures.append(st)
                # append structure from file
                except:
                    self._structures.append(next(structure.StructureReader("input/" + str(temp[0]))))
                # add concentrations
        fh.close()

##------------------------------------------------------------------------------------------------------------
# method to build box
    def build(self):
    #    print(self._structures)
        g = grid(total=self._total, spacing=self._spacing, height=self._height, width=self._width, length=self._length, rand=self._rand, blocked=self._blocked)
        if self._rand:
            st = g.randPop(self._structures, self._concentrations)
        elif self._blocked:
            st = g.blockPop(*self._structures)
        else:
            st = g.seqPop(*self._structures)
        self._box = "out/" + self._name + ".mae"
        # delete_hydrogens(st)
        # add_hydrogens(st)
        st.write("out/" + self._name + ".mae")
        st.write("out/" + self._name + ".pdb")
        return st

##------------------------------------------------------------------------------------------------------------
# method to run sys buiilder
    def sysBuild(self):
        # sys builder
        setup = sysBuilder(self._box, self._name, self._trialName, self._path)
        ## move mae and write msj
        setup.moveMAE()
        setup.writeMSJ()
        ## run sys builder
        setup.run()
        # exit()

##------------------------------------------------------------------------------------------------------------
# method to run MD simulation
    def run(self):
        # run without coarse graining
        if not self._aaCG:
            # md job
            md = mdJob(self._name, self._trialName, self._path)
            ## move setup cms and write MSJ/CFG
            md.moveCMS()
            md.writeMSJ()
            SEED = "100"
            TIME = "100000"
            INTERVAL = "5"
            md.writeCFG({"seed":SEED, "time":TIME, "interval":INTERVAL})
            # run job
            md.run()
        # run with coarse graining
        elif self._aaCG:
            cg = cgJob(self._name, self._trialName, self._path)
            cg.moveMSJ()
            cg.randomizeSeed()
            cg.run()

##------------------------------------------------------------------------------------------------------------
# method to analyze trajectories
    def analyze(self):
        if not self._aaCG:
            cmsfile_out = self._path + self._trialName + "/" + self._name + "/" + self._name + "-out.cms"
            cmsfile_in = self._path + self._trialName + "/" + self._name + "/" + self._name + "-in.cms"
        elif self._aaCG:
            cmsfile_out = self._path + self._trialName + "/" + self._name + "/" + self._name + "_aacg-out.cms"
            cmsfile_in = self._path + self._trialName + "/" + self._name + "/" + self._name + "_aacg-in.cms"

        print(cmsfile_out)
        # print(cmsfile_in)
        # while job not done, wait
        while not os.path.exists(cmsfile_out):
            time.sleep(100)
        time.sleep(10)
        print("Calculating Relative Contacts")

        sys_out = msys.Load(cmsfile_out)
        protein_out = sys_out.selectIds("protein")
        solv_out = sys_out.selectIds("water")
        results_out = sys_out.findContactIds(3, ids=protein_out, other=solv_out)

        # calculate contacts for initial frame
        msys_model, cms_model, traj = traj_util.read_cms_and_traj(cmsfile_out)
        new_pos = traj[1].pos()
        sys_in = msys.Load(cmsfile_in)
        sys_in.setPositions(new_pos)
        protein_in = sys_in.selectIds("protein")
        solv_in = sys_in.selectIds("water")
        results_in = sys_in.findContactIds(3, ids=protein_in, other=solv_in, pos=new_pos)
        #
        # pickle.dump(results_out, open("contact_validation/" + self._name + "_out.p", "wb"))
        # pickle.dump(results_in, open("contact_validation/" + self._name + "_in.p", "wb"))

        total_out = len(results_out)
        total_in = len(results_in)
        total_rel = total_out/total_in
        return total_rel
        # # import trajectories
        # msys_model, cms_model, traj = traj_util.read_cms_and_traj(cmsfile)
        #
        # # select only protein
        # protein_aids = cms_model.select_atom('protein')
        # protein_st = cms_model.extract(protein_aids)
        # protein_gids = topo.aids2gids(cms_model, protein_aids, include_pseudoatoms=False)
        #
        # # wrap final frame
        # to_wrap = [traj[-1]]
        # wrapped = topo.center(msys_model, protein_gids, to_wrap)
        # wrapped.insert(0, traj[0]) # add first frame too
        #
        # # find SASA
        # POI_Intra = analysis.SolventAccessibleSurfaceArea(msys_model, cms_model, "protein")
        # SASA = analysis.analyze(wrapped, POI_Intra)
        # #print(SASA)
        # return SASA[-1]/SASA[0]

##------------------------------------------------------------------------------------------------------------
# method to delete maestro job folder
    def deleteJob(self):
        folder = self._path + self._trialName + "/" +self._name + "/"
        rmtree(folder)
