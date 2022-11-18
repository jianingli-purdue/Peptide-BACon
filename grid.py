import numpy as np
import math
import dev
import random
from schrodinger.structutils import transform, analyze, minimize
from schrodinger import structure
from schrodinger.structutils.build import add_hydrogens, connect, delete_hydrogens
#------------------------------------------------------------------------------------------------------------
## grid class - represents grid on box, class methods populate grid with molecules
# takes user input dimensions and builds grid (angstroms)
class grid:
    def __init__(self, total=None, spacing=None, height=None, width=None, length=None, rand=False, blocked=False, source=None):
        # initialize everything to none
        self._spacing = spacing
        self._total = total
        self._height = height
        self._width = width
        self._length = length
        self._rand = rand
        self._blocked = blocked
        self._source = source

##------------------------------------------------------------------------------------------------------------
    # private method to find dimensions of box
    def _setCoords(self, total=None, spacing=None, height=None, width=None, length=None):
        # if user gives total points & height & spacing
        if (total is not None) and (spacing is not None) and (height is not None):
            #print(self._total, self._height)
            pointsPerWidth = int(math.sqrt(self._total/self._height))
            #print(pointsPerWidth)
            self._coords = np.zeros(shape=(pointsPerWidth, pointsPerWidth, self._height))
        if (spacing is not None) and (height is not None) and (width is not None) and (length is not None):
            pointsPerHeight = height//spacing
            pointsPerWidth = width//spacing
            pointsPerLength = length//spacing
            self._coords = np.zeros(shape=(pointsPerLength, pointsPerWidth, pointsPerHeight))

##------------------------------------------------------------------------------------------------------------
    # method to test whether input is filename or structure
    def _testInput(self, input):
        if isinstance(input, str):
            return next(structure.StructureReader(self._source + input))
        elif input is None:
            return None
        else:
            return input.copy()

##------------------------------------------------------------------------------------------------------------
    # function to create list of structures to populate from
    # accepts list of molecules (structure  s/filenames) and list of concentrations
    def _createStList(self, structures, concentrations=None):
        # blocked list
        if self._blocked:
            blockedSts = []
            numBlocks = len(structures)
            # iterate  through structures and add them to list
            for st in structures:
                for i in range(self._total//numBlocks):
                    blockedSts.append(st)
            return blockedSts

        # if >1 structures
        if len(structures) > 1:
            # check conc
            total = 0
            for prob in concentrations:
                total += float(prob)
            if total > 1:
                print("The concentrations do not add up to one, please enter new concentrations now")
                # to do: user input
            # build list
            concentratedSts = []
            for i, prob in enumerate(concentrations):
                reps = float(prob) * self._total
                for j in range(int(reps)):
                    concentratedSts.append(structures[i])
            random.shuffle(concentratedSts)
            #print(concentratedSts)
            return concentratedSts

        # if only one specified structure
        elif len(structures) == 1:
            concentratedSts = []
            # add complementary probability and blank structure
            concentrations.append(1-float(concentrations[0]))
            structures.append(None)
            # build list with blanks
            for i, prob in enumerate(concentrations):
                reps = float(prob) * self._total
                #print(reps)
                for j in range(int(round(reps))):
                    concentratedSts.append(structures[i])
            random.shuffle(concentratedSts)
            # check to make sure first index is a structure
            while concentratedSts[-1] is None:
                random.shuffle(concentratedSts)
        #    print(concentratedSts)
            return concentratedSts

##------------------------------------------------------------------------------------------------------------
    # method to populates the grid sequentially with given molecule(s)
    # accepts molecule(s) as either string filenames or structures
    def seqPop(self, *argv):
        # get shape of box
        self._setCoords(total=self._total, spacing=self._spacing, height=self._height, width=self._width, length=self._length)
        x,y,z, = self._coords.shape
        res = 0


        # populate sequentially for only one molecule
        if len(argv) == 1:
            molecule = argv[0]
        #    print(molecule)
            l = 0
            # iterate through array
            for i in range(x):
                for j in range(y):
                    for k in range(z):
                        # base case
                        if (i == 0 and j == 0 and k == 0):
                            #print(i, j, k)
                            baseSt = self._testInput(molecule)
                            # for chain in baseSt.chain:
                            #     chainNum = ord(chain.name) + 1
                            #     #print(chain.name)
                            transform.translate_structure(baseSt, 0, 0, 0)
                            transform.rotate_structure(baseSt, y_angle=0.5*math.pi, z_angle=0.25*math.pi)

                            # make into zwitterion
                            # finalSt = dev.zwitter_peptide(baseSt)
                            finalSt = baseSt
                            #increment chain
                            # for chain in baseSt.chain:
                            #     chain.name = chr(chainNum)
                            #     chainNum += 1
                            # increment resnum
                            for residue in baseSt.residue:
                                residue.resnum = res
                                res += 1

                        # copy base st and move it
                        else:
                            tmp_st = baseSt.copy()
                            # measure length of molecule and add half of it to spacing
                            atoms = []
                            for atom in tmp_st.atom:
                                atoms.append(atom.xyz)
                            #print(atoms[0], atoms[-1])
                            molLength = [atoms[0][0]-atoms[-1][0], atoms[0][1]-atoms[-1][1], atoms[0][2]-atoms[-1][2]]
                            toAdd = math.sqrt(molLength[0]**2 + molLength[1]**2 + molLength[2]**2)
                        #    print(molLength)
                            transform.translate_structure(tmp_st, i*self._spacing, j*self._spacing, k*self._spacing) # PLAY WITH ME PLEASE
                            #increment chain
                            # for chain in tmp_st.chain:
                            #     print("before")
                            #     print(chain.name)
                            #     chain.name = chr(chainNum)
                            #     print("after")
                            #     print(chain.name)
                            # chainNum += 1
                            # increment resnum
                            for residue in tmp_st.residue:
                                residue.resnum = res
                                res += 1
                            # merge into final st
                            finalSt = finalSt.merge(tmp_st)
                # increment l
                        l += 1
                    l += 1
                l += 1

        # populate sequentially for > 1 molecule
        if len(argv) > 1:
            # initialize molecule iterators
            numMolecules = len(argv)
            molIterator = 0
            l = 0
            # iterate through array
            for i in range(x):
                for j in range(y):
                    for k in range(z):
                        # base case
                        if (i == 0 and j == 0 and k == 0):
                            baseSt = self._testInput(argv[molIterator])
                            transform.translate_structure(baseSt, 0, 0, 0)
                            #increment chain
                            # for chain in baseSt.chain:
                            #     chain.name = chr(res)
                            # increment resnum
                            for residue in baseSt.residue:
                                residue.resnum = res
                                res += 1
                            finalSt = baseSt
                            molIterator += 1
                        # copy base st and move it
                        else:
                            # use new molecule
                            if molIterator == numMolecules:
                                molIterator = 0
                            tmp_st = self._testInput(argv[molIterator])

                            transform.translate_structure(tmp_st, i*self._spacing, j*self._spacing, k*self._spacing)
                            #increment chain
                            # for chain in tmp_st.chain:
                            #     chain.name = chr(chainNum)
                            #     chainNum += 1
                            # increment resnum
                            for residue in tmp_st.residue:
                                residue.resnum = res
                                res += 1
                            # merge into final st
                            finalSt = finalSt.merge(tmp_st)

                            # increment to next molecule
                            molIterator += 1
                # increment l
                        l += 1
                    l += 1
                l += 1


        return finalSt

##------------------------------------------------------------------------------------------------------------
    # method to populate a grid with blocks of molecules
    # accepts molecules as either string filenames or structures
    def blockPop(self, *argv):
        # get shape of box
        self._setCoords(total=self._total, spacing=self._spacing, height=self._height, width=self._width, length=self._length)
        x,y,z, = self._coords.shape
        res = 0
        chainNum = 32

        # populate blocks
        # initialize molecule iterators
        sts = self._createStList(argv)
        l = 0
        # iterate through array
        for i in range(x):
            for j in range(y):
                for k in range(z):
                    molecule = sts.pop()
                    # base case
                    if (i == 0 and j == 0 and k == 0):
                        baseSt = self._testInput(molecule)
                        finalSt = dev.zwitter_peptide(baseSt)
                        transform.translate_structure(baseSt, 0, 0, 0)
                        #increment chain
                        # for chain in baseSt.chain:
                        #     chain.name = chr(chainNum)
                        #     chainNum += 1
                        # increment resnum
                        for residue in baseSt.residue:
                            residue.resnum = res
                            res += 1
                    # copy base st and move it
                    else:
                        tmp_st = self._testInput(molecule)

                        transform.translate_structure(tmp_st, i*self._spacing, j*self._spacing, k*self._spacing)
                        #increment chain
                        # for chain in tmp_st.chain:
                        #     chain.name = chr(chainNum)
                        #     chainNum += 1
                        # increment resnum
                        for residue in tmp_st.residue:
                            residue.resnum = res
                            res += 1
                        # merge into final st
                        finalSt = finalSt.merge(tmp_st)

        return finalSt

##------------------------------------------------------------------------------------------------------------
    # method to randomly populate the grid
    # accepts list of structures and list of concentrations
    def randPop(self, structures, concentrations):
        # get shape of box
        self._setCoords(total=self._total, spacing=self._spacing, height=self._height, width=self._width, length=self._length)
        x,y,z, = self._coords.shape
        # sequentially populate more than one molecule
        if len(structures) > 1:
            # generate list to pick from
            sts = self._createStList(structures, concentrations)
            l = 0
            res = 0
            chainNum = 32

            # iterate through array
            for i in range(x):
                for j in range(y):
                    for k in range(z):
                        # add random st from sts list
                        molecule = sts.pop()
                        # base case
                        if (i == 0 and j == 0 and k == 0):
                            # test whether input molecule is string or st
                            baseSt = self._testInput(molecule)
                            transform.translate_structure(baseSt, 0, 0, 0)
                            finalSt = dev.zwitter_peptide(baseSt)
                            #increment chain
                            # for chain in baseSt.chain:
                            #     #print(chr(chainNum))
                            #     chain.name = chr(chainNum)
                            # chainNum += 1
                            #increment resnum
                            for residue in baseSt.residue:
                                residue.resnum = res
                                res += 1
                        # copy base st and move it
                        else:
                            tmp_st = self._testInput(molecule)
                            if tmp_st == None:
                                pass
                            else:
                                transform.translate_structure(tmp_st, i*self._spacing, j*self._spacing, k*self._spacing)
                                tmp_st = dev.zwitter_peptide(tmp_st)

                            #increment chain
                            # for chain in tmp_st.chain:
                            #     #print(chr(chainNum))
                            #
                            #     chain.name = chr(chainNum)
                                # chainNum += 1
                            # increment resnum
                            for residue in tmp_st.residue:
                                residue.resnum = res
                                res += 1

                            # merge into final st
                            finalSt = finalSt.merge(tmp_st)

        # populate sequentially for >1 molecule
        if len(structures) == 1:
            # generate list to pick from
            sts = self._createStList(structures, concentrations)
            l = 0
            res = 0
            chainNum = 68
            # iterate through array
            for i in range(x):
                for j in range(y):
                    for k in range(z):
                        # add random st from sts list
                        molecule = sts.pop()
                        # base case
                        if (i == 0 and j == 0 and k == 0):
                            # test whether input molecule is string or st
                            baseSt = self._testInput(molecule)
                            transform.translate_structure(baseSt, 0, 0, 0)
                            finalSt = dev.zwitter_peptide(baseSt)
                            #increment chain
                            # for chain in finalSt.chain:
                            #     #print(chr(chainNum))
                            #     chain.name = chr(chainNum)
                            # chainNum += 1
                            # increment resnum
                            for residue in baseSt.residue:
                                residue.resnum = res
                                res += 1

                        # copy base st and move it
                        else:
                            tmp_st = self._testInput(molecule)
                            if tmp_st is None:
                                pass
                            else:
                                tmp_st = dev.zwitter_peptide(tmp_st)
                                transform.translate_structure(tmp_st, i*self._spacing, j*self._spacing, k*self._spacing)
                                #increment chain
                                # for chain in tmp_st.chain:
                                #     #print(chr(chainNum))
                                #     chain.name = chr(chainNum)
                                #     #print(chain.name)
                                # try:
                                #     finalSt.write("zwitter.pdb")
                                # except:
                                #     pass
                                chainNum += 1
                                # increment resnum
                                for residue in tmp_st.residue:
                                    residue.resnum = res
                                    res += 1

                                # merge into final st
                                finalSt = finalSt.merge(tmp_st)
        return finalSt
